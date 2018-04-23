#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2018-03-21 15:43:35
# Contact:    Nolan.Dickson@Canada.ca

import os
import ftplib
import netCDF4
import argparse
import datetime
import pyloader
import numpy as np
import pyhdf.SD as SD
import MISRutils as mu

'''
Given a coordinate and a date range, find the relevant VFM files and then
search them for the requested feature.
If feature is found, determine the geographic locations, altitudes, etc
of feature.
output the results as nc file.

The download and extract functions can also be called from another python script
'''

# TODO type, subtype arg to list?
# TODO make some examples and instructions on working with output files


ftype_list = {
    0: 'invalid',
    1: 'clear_air',
    2: 'cloud',
    3: 'aerosol',
    4: 'stratospheric',
    5: 'surface',
    6: 'subsurface',
    7: 'no_signal'
}
subtypes = {
    3: {
        0: 'not_determined',
        1: 'clean_marine',
        2: 'dust',
        3: 'polluted_continental',
        4: 'clean_continental',
        5: 'polluted_dust',
        6: 'smoke',
        7: 'other'
    },
    2: {
        0: 'low_overcast_t',
        1: 'low_overcast_o',
        2: 'transition_stratocumulus',
        3: 'broken_cumulus',
        4: 'altocumulus_t',
        5: 'altostratus_o',
        6: 'cirrus_t',
        7: 'deep_convective_o'
    },
    4: {
        0: 'not_determined',
        1: 'non_depolarizing_PSC',
        2: 'depolarizing_PSC',
        3: 'non_depolarizing_aerosol',
        4: 'depolarizing_aerosol',
        7: 'other'
    }
}


def download_cal(srch_times, savedir):
    '''
    given a list of datetimes, figure out which file will contain each time and
    download it
    '''
    ftp_url = 'ftp.icare.univ-lille1.fr'
    icare = ftplib.FTP(ftp_url, user=user, passwd=password)
    icare.cwd('SPACEBORNE/CALIOP/VFM.v4.10/')

    if debug_flag:
        _debug('Successfully accessed icare ftp site.\n')

    local_files = []

    for time in srch_times:

        ftp_dir = '{date:%Y}/{date:%Y_%m_%d}/'.format(date=time)

        try:
            icare.cwd(ftp_dir)
        except ftplib.error_perm:
            if debug_flag:
                _debug('no VFM files on {d:%Y}/{d:%Y_%m_%d}/'.format(d=time))
            continue

        file_list = icare.nlst()

        ctempl = '%Y-%m-%dT%H-%M-%S'
        vind, vfile = min(
            enumerate(file_list),
            key=lambda f:
                (time - datetime.datetime.strptime(f[1][30:49], ctempl)
                 if datetime.datetime.strptime(f[1][30:49], ctempl) < time
                 else datetime.timedelta(999))
        )

        if debug_flag:
            _debug('searching at time: {0}\nnearest file: {1}\n'.format(
                   time, vfile))

        delta_t = datetime.datetime.strptime(vfile[30:49], ctempl) - time

        if -300 < delta_t.seconds < 0:
            # file start is within 5 minutes less than time, grab earlier one
            inter_files = [file_list[vind - 1], vfile]

            if debug_flag:
                _debug('\tnear edge -> add: {}\n'.format(file_list[vind - 1]))

        else:
            inter_files = [vfile]

        for ftp_file in inter_files:
            local_path = '{dir}/{file}'.format(dir=savedir, file=ftp_file)
            local_files.append(local_path)

            if not os.path.exists(local_path):
                with open(local_path, 'wb') as f:
                    icare.retrbinary('RETR ' + ftp_file, f.write)

        icare.cwd('../../')

    icare.close()
    return local_files


def get_feature(filename, coords, ftype, subtype=None):

    vfile = SD.SD(filename, SD.SDC.READ)

    # ---------------------------------------------------------------------
    # Get geographic coordinates and find bounding indices (of ss res)
    # ---------------------------------------------------------------------

    lon_arr = vfile.select('ssLongitude')[:, 0]
    lat_arr = vfile.select('ssLatitude')[:, 0]

    range_inds = np.where(
        ((coords[0][0] < lon_arr) & (lon_arr < coords[3][0])) &
        ((coords[0][1] < lat_arr) & (lat_arr < coords[3][1]))
    )[0]

    # ---------------------------------------------------------------------
    # Get feature flags and split the 3 resolution layers into sperate
    # arrays, for reshaping
    # ---------------------------------------------------------------------

    fcf_arr = vfile.select('Feature_Classification_Flags').get()

    layer2 = fcf_arr[:, 0:165]
    layer1 = fcf_arr[:, 165:1165]
    layer0 = fcf_arr[:, 1165:5515]

    # ---------------------------------------------------------------------
    # Reshape all layers based on their resolution, and flip to ensure
    # the bottom altitude level is at the first index of each
    # ---------------------------------------------------------------------

    # TODO flip is not 2/3 compatible

    sslayer2 = np.flip(layer2.reshape(layer2.shape[0] * 3, 55), axis=1)
    sslayer1 = np.flip(layer1.reshape(layer1.shape[0] * 5, 200), axis=1)
    sslayer0 = np.flip(layer0.reshape(layer0.shape[0] * 15, 290), axis=1)

    ss_fcf_arr = np.empty((lon_arr.shape[0], 545), dtype=np.uint16)

    # ---------------------------------------------------------------------
    # Construct the complete (ss) fcf array, by stacking all layers.
    # Higher layers are repeated to fit, in line with their resolution
    # ---------------------------------------------------------------------

    ss_fcf_arr[:, 490:545] = sslayer2.repeat(5, axis=0)
    ss_fcf_arr[:, 290:490] = sslayer1.repeat(3, axis=0)
    ss_fcf_arr[:, :290] = sslayer0

    # ---------------------------------------------------------------------
    # Cut out the section of the (ss) fcf array within the bounding box
    # ---------------------------------------------------------------------

    ss_fcf_arr = ss_fcf_arr[range_inds, :]

    # ---------------------------------------------------------------------
    # determine the feature type based on the first bits. Determine the
    # indices of all the feature occurences
    # ---------------------------------------------------------------------

    type_arr = ss_fcf_arr & 7
    feat_inds = np.where((type_arr == ftype))

    # ---------------------------------------------------------------------
    # If looking specifically for a subtype, search the feature indices at
    # the 10th bit for the feature subtype
    # ---------------------------------------------------------------------

    if subtype is not None:
        subtype_arr = np.zeros_like(ss_fcf_arr)
        subtype_arr[feat_inds] = ss_fcf_arr[feat_inds]
        subtype_arr = (subtype_arr & 3584) >> 9
        feat_inds = np.where(subtype_arr == subtype)

    inds_1D = np.unique(feat_inds[0])

    # ---------------------------------------------------------------------
    # Create a boolean masked array around the feature indices
    # and change the fcf array to a binary array only containing the
    # columns with the feature in them
    # ---------------------------------------------------------------------

    mask = np.ones_like(ss_fcf_arr)
    mask[feat_inds] = 0

    mask = mask[inds_1D, :]
    fcf_ones = np.ones_like(ss_fcf_arr[inds_1D, :])

    feat_mask = np.ma.masked_array(fcf_ones, mask=mask, dtype=np.float32)

    # ---------------------------------------------------------------------
    # Create an array of altitude matching the shape of the feature mask
    # ---------------------------------------------------------------------

    alt2 = [20.2 + i * 0.18 for i in range(55)]
    alt1 = [8.2 + i * 0.06 for i in range(200)]
    alt0 = [-0.5 + i * 0.03 for i in range(290)]

    ss_alt_layer2 = np.ones((fcf_ones.shape[0], 55), dtype=np.float32) * alt2
    ss_alt_layer1 = np.ones((fcf_ones.shape[0], 200), dtype=np.float32) * alt1
    ss_alt_layer0 = np.ones((fcf_ones.shape[0], 290), dtype=np.float32) * alt0

    ss_alt = np.empty_like(fcf_ones, dtype=np.float32)

    ss_alt[:, 490:545] = ss_alt_layer2
    ss_alt[:, 290:490] = ss_alt_layer1
    ss_alt[:, :290] = ss_alt_layer0

    # ---------------------------------------------------------------------
    # Map the altitude array onto the feature mask in metres
    # ---------------------------------------------------------------------

    alt_mask = feat_mask * ss_alt * 1000.

    # ---------------------------------------------------------------------
    # Grab the correct range of lat, lon and time arrays and the maximum
    # altitude of feature in each column, for output
    # ---------------------------------------------------------------------

    lon_arr = lon_arr[range_inds][inds_1D]
    lat_arr = lat_arr[range_inds][inds_1D]
    time_arr = vfile.select('ssProfile_UTC_Time').get()[range_inds, 0][inds_1D]

    max_alt = alt_mask.max(axis=1)

    # ---------------------------------------------------------------------
    # Write debug and return
    # ---------------------------------------------------------------------

    if debug_flag:
        _debug(
            'feat_inds:\n{}\nlon_arr:\n{}\nlat_arr:\n{}\ntime_arr:\n{}\n'
            'max_alt:\n{}\nfeat_mask:\n{}\nalt_mask:\n{}\n'.format(
                feat_inds, lon_arr, lat_arr, time_arr,
                max_alt, feat_mask, alt_mask,
            )
        )

    return (lon_arr, lat_arr, time_arr, max_alt, feat_mask, alt_mask)


def _write_output(outfile, res, args):
    '''
    take in the path to an output file and a dictionary containing results for
    each day (and the starting args to get some creation metadata), and write
    the results to the file
    '''
    # TODO attributes to variables (units, descriptions, etc)
    # TODO add confidences and qa's and stuff

    out = netCDF4.Dataset(outfile, 'w')

    for date in res:
        data = res[date]  # list of tuples containing results for each file
        dt_grp = out.createGroup(date)

        # create dimensions

        numcols = sum([len(dataset[0]) for dataset in data])
        dt_grp.createDimension('numcols', numcols)
        dt_grp.createDimension('numlayers', 545)

        # create and assign variables

        flat_names = ['Longitude', 'Latitude',
                      'UTC_Profile_time', 'Max_Altitude']

        mult_names = [
            '{0}{1}'.format(
                ftype_list[args.type],
                '-' + subtypes[args.type][args.subtype] if args.subtype is not None else ''
            ),
            'Altitude'
        ]

        for fi, flat in enumerate(flat_names):
            dtype = data[0][fi].dtype
            dims = ('numcols',)
            var = dt_grp.createVariable(flat, dtype, dims, zlib=True)

            var_holder = data[0][fi]
            for dataset in data[1:]:
                var_holder = np.hstack((var_holder, dataset[fi]))

            var[:] = var_holder

        for fi, mult in enumerate(mult_names, 4):
            dtype = data[0][fi].dtype
            dims = ('numcols', 'numlayers')
            var = dt_grp.createVariable(mult, dtype, dims, zlib=True)

            var_holder = data[0][fi]
            for dataset in data[1:]:
                var_holder = np.vstack((var_holder, dataset[fi]))

            var[:] = var_holder

    # create and assign metadata

    out.longitude_range = args.coords[:2]
    out.latitude_range = args.coords[2:]

    out.date_range = [int(d.strftime('%Y%m%d')) for d in args.daterange]
    out.feature = '{0} - {1}'.format(args.type, ftype_list[args.type])
    out.subtype = '{0} - {1}'.format(
        args.subtype, subtypes[args.type][args.subtype]
    ) if args.subtype is not None else 'all'

    out.close()


if __name__ == "__main__":

# ---------------------------------------------------------------------
# (0) Parse arguments
# ---------------------------------------------------------------------

    parser = argparse.ArgumentParser(
        description="Find a specific feature type, within a coord/date range, "
                    "in the CALIOP VFM product."
    )
    parser.add_argument(
        '-d', '--daterange', required=True, nargs='*',
        help="Date, must have format YYYY.MM.DD. Can also "
             "search over range of dates, "
             "simply by providing a second, later date. Please only "
             "do a date range inside the same year."
    )
    parser.add_argument(
        '-L', '--coords', type=float, nargs='*', required=True,
        help="Longitude, latitude range to search for hdf files within: "
             "minLON maxLON minLAT maxLAT\n Please avoid crossing the "
             "poles and antemeridian with the range."
    )
    parser.add_argument(
        '-u', '--user', required=False,
        help='Your ICARE username.'
    )
    parser.add_argument(
        '-p', '--password', required=False,
        help='Your ICARE password'
    )
    parser.add_argument(
        '-t', '--type', type=int, default=3,
        help='Which feature type to extract. See Calipso Data Manual for list'
    )
    parser.add_argument(
        '-s', '--subtype', type=int, required=False, default=None,
        help='Which feature subtype to extract. See Calipso Data Manual '
        'for list. If not provided, will return entirety of feature type.'
    )
    parser.add_argument(
        '--caldir', default=None,
        help='Directory where the CALIOP VFM files will be stored.'
    )
    parser.add_argument(
        '--debug', action='store_true',
        help='Create a log file with information for use in debugging.'
    )
    args = parser.parse_args()

# ---------------------------------------------------------------------
# (1) Check all arguments
# ---------------------------------------------------------------------

    debug_flag = args.debug

    dates = [datetime.datetime.strptime(d, '%Y.%m.%d') for d in args.daterange]
    args.daterange = dates

    caldir = args.caldir if args.caldir else os.getcwd()

    if len(dates) > 2:
        parser.error("Invalid daterange. Please provide a range of two dates"
                     " in the form: YYYY.MM.DD YYYY.MM.DD")

    if len(args.coords) == 4:
        minLON, maxLON, minLAT, maxLAT = args.coords
        coords = [(minLON, minLAT), (minLON, maxLAT),
                  (maxLON, minLAT), (maxLON, maxLAT)]
    else:
        parser.error("Invalid coordinates. Please provide a range of the "
                     "form: minLON maxLON minLAT maxLAT")

    user = args.user if args.user else 'DicksonN'
    password = args.password if args.password else 'SSHSBench2015!'

    ftype = args.type
    subtype = args.subtype if args.subtype else None

# ---------------------------------------------------------------------
# (1.5) Initiate debugging log file, if requested
# ---------------------------------------------------------------------

    if debug_flag:
        debug_dir = '{}/debug'.format(os.getcwd())
        logfile = '{dir}/log_s{s:%y%m%d}_e{e:%y%m%d}_c{c}.log'.format(
            dir=debug_dir if os.path.exists(debug_dir) else os.getcwd(),
            s=dates[0],
            e=dates[-1],
            c=datetime.datetime.now().strftime('%y%m%d%H%M%S')
        )

        def _debug(mssg):
            with open(logfile, 'a') as f:
                f.write(mssg)

        mssg = '{d} CALIPSO VFM EXTRACTION {d}\n\n'.format(d='-' * 10)
        for n, a in vars(args).items():
            mssg += '{0}: {1}\n'.format(n, a)
        mssg += 'Feature type: {}\n'.format(ftype_list[ftype])
        mssg += 'Feature subtype: {}\n\n'.format(subtypes[ftype][subtype]
                                                 if subtype is not None else None)
        _debug(mssg)

# ---------------------------------------------------------------------
# (2) Determine Calipso overpasses of the ranges
# ---------------------------------------------------------------------

    load = pyloader.Load_message('> Finding overpasses')

    op_times = mu.find_overpass(coord=coords, daterange=dates,
                                swath_width=None, sat='Calipso')

    if debug_flag:
        _debug('overpass times:\n{}\n\n'.format(op_times))

    if not op_times:
        mssg = 'No Calipso overpasses match these date and coordinate ranges.'
        raise RuntimeError(mssg)

# ---------------------------------------------------------------------
# (3) Download VFM data for the relevant overpasses
# ---------------------------------------------------------------------

    load.update('> Downloading VFM files')

    file_list = download_cal(op_times, os.getcwd())

    if debug_flag:
        _debug('all VFM files downloaded:\n{}\n'.format(
               [f + '\n' for f in file_list]))

    if not file_list:
        mssg = 'No VFM files were found matching these overpass times.'
        raise RuntimeError(mssg)

# ---------------------------------------------------------------------
# (4) Gather relavant data from each VFM to store in ouput file
# ---------------------------------------------------------------------

    res = {}
    for find, filename in enumerate(file_list):

        if debug_flag:
            _debug('working on: {}\n'.format(filename))

        # -----------------------------------------------------------------
        # Find and extract feature information
        # -----------------------------------------------------------------

        ext_data = get_feature(filename, coords, ftype, subtype)

        if ext_data[0].size:

            date = os.path.basename(filename)[30:40]
            if date not in res:
                res[date] = []

            res[date].append(ext_data)

# ---------------------------------------------------------------------
# (5) Write the results, along with some metadata, to an output file
# ---------------------------------------------------------------------

    load.update('> Writing data to file')

    outfile = '{d}/CAL_{f}{sf}_s{s:%y%m%d}_e{e:%y%m%d}_c{c}.nc'.format(
        f=ftype_list[ftype],
        sf='' if not subtype else '_' + subtypes[ftype][subtype],
        d=os.getcwd(),
        s=dates[0],
        e=dates[-1],
        c=datetime.datetime.now().strftime('%y%m%d%H%M%S')
    )

    _write_output(outfile, res, args)

    load.stop()

    if debug_flag:
        _debug('wrote results to {0}\ntotal time taken: {1}\n'.format(
               outfile, load.delta_t))
