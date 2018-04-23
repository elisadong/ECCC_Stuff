#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2018-02-14 10:18:16
# Contact:    Nolan.Dickson@Canada.ca

# doubt I'll even try to document this, tis a mess
# if you're looking at this after I've left, only physical notes still exist, somewhere

# TODO alot of this swath drawing can probably be replaced
# with actual misr block checking, using the new SOM conversion
# functions

# That is an *obscene* amount of imports for such a short program
import os
import re
import sys
import math
import time
import pickle
import ftplib
import netCDF4
import argparse
import datetime
import requests
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pyorbital.orbital as porb

import Cpoly
import pyloader
import postMINX as pm
import MISRutils as mu


# problems to address:
# - the TLE drawn MISR Swath never lines up exactly with the MISR imagery exactly
# - The fire pixels are grabbed based on a time and just checked all over the swath, which leads to alot of completely irrelevant pixels
# - The non-variable drawn swath size is obviously only a perfect match at certain times in the orbit


def distance(long1, lat1, long2, lat2):
    '''
    Calculate distance in metres between
    two coordinates using the Haversine formula
    '''

    R = 6371000  # Radius of Earth, in metres
    phi_c = lat2 * math.pi / 180.
    phi_p = lat1 * math.pi / 180.
    del_phi = (lat1 - lat2) * math.pi / 180.
    del_lam = (long1 - long2) * math.pi / 180.

    a = math.sin(del_phi / 2.)**2 + math.cos(phi_p) \
        * math.cos(phi_c) * math.sin(del_lam / 2.)**2
    dist = R * 2 * math.atan2(a**0.5, (1 - a)**0.5)

    return dist


def draw_swath(path, swath_width):
    '''Given an orbit path, return a polygon of a determined width around it'''

    poly_left, poly_right = [], []
    for pind, point in enumerate(path[:-1]):
        x1, y1, time1 = point
        x2, y2, time2 = path[pind + 1]

        dx = x2 - x1
        dy = y2 - y1
        d = math.sqrt((dx * dx) + (dy * dy)) / 2.

        r = swath_width  # TODO Should be calculating based on function(index), but this is wayy easier and probably fine

        xo = x1 + math.sqrt((r * r) + (d * d))
        theta = math.atan(r / d)
        phi = math.atan(dy / dx)

        alpha_left = phi + theta if phi > 0. else phi - theta
        alpha_right = phi - theta if phi > 0. else phi + theta

        poly_time = time1 + (time2 - time1) / 2.

        poly_left.append((
            (((xo - x1) * math.cos(alpha_left)) - ((y1 - y1) * math.sin(alpha_left)) + x1),
            (((xo - x1) * math.sin(alpha_left)) + ((y1 - y1) * math.cos(alpha_left)) + y1),
            poly_time
        ))
        poly_right.append((
            (((xo - x1) * math.cos(alpha_right)) - ((y1 - y1) * math.sin(alpha_right)) + x1),
            (((xo - x1) * math.sin(alpha_right)) + ((y1 - y1) * math.cos(alpha_right)) + y1),
            poly_time
        ))

    poly_right.reverse()
    MISR_poly = poly_left + poly_right
    return MISR_poly


def draw_orbits(date):
    '''
    take in a date and return all the ascending Terra paths on that day seperately, and all the orbit tracks of Calipso on that day
    please give date as a datetime object with only the date (0 hrs, mins, etc)
    '''
    j_date = int(date.strftime('%y%j'))

    def read_TLE_dates(filename):
        tle_list = []
        with open(filename) as tf:
            for line in tf:
                epoch_date = float(line[18:32])
                # if on date requested +- 1 day, save TLE
                if j_date - 1. <= epoch_date <= j_date + 1.:
                    tle_list.append((line, next(tf)))
                # TLE are in desecending order, if after this point, can break from file
                elif epoch_date < j_date - 1.:
                    break
                # only parse `1` lines
                else:
                    next(tf)

        return tle_list

    res = {
        'Terra': {
            'lonlat_arr': [],
            'ascending': [],
            'last_nearest': False,
            'tle': read_TLE_dates('/fs/home/fs1/eccc/aq/arqi/mas001/computer_programs/python/MINX/data/Terra_TLE.txt'),
            'orbit': None,
            'lonlat_paths': [],
        },
        'Calipso': {
            'lonlat_arr': [],
            'last_nearest': False,
            'tle': read_TLE_dates('/fs/home/fs1/eccc/aq/arqi/mas001/computer_programs/python/MINX/data/Calipso_TLE.txt'),
            'orbit': None
        }
    }

    iterdate = date
    index = 0
    while date + datetime.timedelta(1) > iterdate:

        time_frac = ((iterdate.hour * 60. + iterdate.minute) * 60. + iterdate.second) / (24. * 60. * 60.)

        for sat in res:

            nearest_tle = min(res[sat]['tle'], key=lambda line: abs(float(line[0][18:32]) - float(iterdate.strftime('%y%j.')) + time_frac))

            # either create a new orbit object or use the previous
            if nearest_tle != res[sat]['last_nearest']:
                res[sat]['orbit'] = porb.Orbital(sat, line1=nearest_tle[0], line2=nearest_tle[1])

            lon, lat, alt = res[sat]['orbit'].get_lonlatalt(iterdate)
            res[sat]['lonlat_arr'].append((lon, lat, iterdate))

            # Only keep descending (daytime) orbits of Terra Satellite
            # Day and night for CALIOP is usable (and pretty much the same) (J. R. Campbell et al., 2012)
            try:
                if sat == 'Terra' and res[sat]['lonlat_arr'][-2][1] < lat:
                    res[sat]['ascending'].append(index)
            except IndexError:
                pass

            res[sat]['last_nearest'] = nearest_tle

        index += 1
        iterdate += datetime.timedelta(minutes=0.5)

    for i in reversed(res['Terra']['ascending']):
        del res['Terra']['lonlat_arr'][i]

    # TODO make much better
    # Split Terra data into seperate paths, to avoid polar/meridian problems
    path_index = 0

    for ind, (lon, lat, time)in enumerate(res['Terra']['lonlat_arr']):

        try:
            res['Terra']['lonlat_paths'][path_index].append((lon, lat, time))
        except IndexError:
            res['Terra']['lonlat_paths'].append([])
            res['Terra']['lonlat_paths'][path_index].append((lon, lat, time))

        try:
            # funny thing is we actually exploit the meridian/polar problems here to avoid them later
            if (lon < 0 and res['Terra']['lonlat_arr'][ind + 1][0] > 0) or (lat < 0 and res['Terra']['lonlat_arr'][ind + 1][1] > 0):
                path_index += 1
        except IndexError:
            pass

    return (res['Terra']['lonlat_paths'], res['Calipso']['lonlat_arr'])


def orbit_intercept(t_path, c_track, delta_t=datetime.timedelta(hours=9)):
    '''
    take in a single (why?) Terra path and a day of Calipso orbits and see if they intercept
    will account for Calipso crossing the antimeridian
    will not accept interecepts in polar circles

    will return a list of dictionaries, for each interecept, containing relevant info

    checks if lines between each calipso point cross over the Terra line +- a swath width outwards, within a certain time limit of MISR
    '''

    # make a polygon around the terra path, and polygon check each c line segment (can gather return info all from calipso segment)

    # t_path.reverse()  # oops did all the math for a ascending orbit
    MISR_poly = draw_swath(t_path, 3.)

    if not MISR_poly:
        raise IndexError('Not enough points, run simulation for shorter time interval')

    MISR_poly_coords = [(x[0], x[1]) for x in MISR_poly]
    intercept_results = []
    checked_inds = []
    for Cind, Cpoint in enumerate(c_track):

        if not (Cpoint[1] > 70. or Cpoint[1] < -70.) and Cind not in checked_inds:
            # skip iterations which cross the antimeridian, and skip the last point, to avoid later indexError
            try:
                if abs(Cpoint[0] - c_track[Cind + 1][0]) > 180.:
                    continue
            except IndexError:
                continue

            if Cpoly.polygon_check(MISR_poly_coords, [(Cpoint[0], Cpoint[1]), (c_track[Cind + 1][0], c_track[Cind + 1][1])], ccwsorted=True):

                # This Calipso segment overlaps with the MISR path, check if time is okay
                dist_arr = [distance(MISR_point[0], MISR_point[1], Cpoint[0], Cpoint[1]) for MISR_point in MISR_poly]
                int_index, _ = min(enumerate(dist_arr), key=lambda dist: dist[1])  # get index of closest point to MISR

                C_time = Cpoint[2]
                *M_coord, M_time = MISR_poly[int_index]

                if abs(C_time - M_time).total_seconds() <= 18000.:  # diff less than 9 hours
                    coordinates = [tuple(M_coord)]

                    # start iterating over the next calipso points until the lines leave the swath again, and then combine all
                    # these intercepts into one and skip ahead
                    # TODO actually probably only need the first and last one (when enter swath and exit swath)
                    for cont_ind, cont_point in enumerate(c_track[Cind + 1:]):
                        cont_Cind = Cind + 1 + cont_ind
                        checked_inds.append(cont_Cind)
                        cont_pnt = (cont_point[0], cont_point[1])
                        cont_pnt_next = (c_track[cont_Cind][0], c_track[cont_Cind][1])

                        if Cpoly.polygon_check(MISR_poly_coords, [cont_pnt, cont_pnt_next], ccwsorted=True):
                            coordinates.append(cont_pnt)

                        else:
                            # TODO this should be the last intecepting, not the first one *not* interrceptign, also find better way to get coordinate
                            dist_arr = [distance(MISR_point[0], MISR_point[1], cont_pnt[0], cont_pnt[1]) for MISR_point in MISR_poly]
                            cont_int_index, _ = min(enumerate(dist_arr), key=lambda dist: dist[1])  # get index of closest point to MISR
                            *M_final_coord, M_final_time = MISR_poly[cont_int_index]

                            coordinates.append(tuple(M_final_coord))
                            break

                    intercept_results.append({
                        'CALIOP_time': C_time,
                        'CALIOP_path': c_track[Cind - 10: cont_Cind + 10],
                        'MISR_time': M_time,
                        'MISR_end_time': M_final_time,
                        'MISR_path': t_path,
                        'coordinates': coordinates
                    })

    return intercept_results


def fire_pixel_check(intercepts, FRP_threshhold=1000.):  # TODO test with a lower frp?
    '''
    Given a list of intercept times and locations, find out if there may have been a fire on that day, in that location

    will return a list of intercepts which should be investigated, along with the name of the MOD14 granule to download
    '''
    # TODO optimize, still the slowest function of the program by minutes
    # TODO reduce false positives
    import time
    res = []
    grans = []

    http_url = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD14.006/'
    folder_date = intercepts[0]['MISR_time'].strftime('%Y.%m.%d')  # TODO is it fine to assume all on same day folder?

    folder_content = requests.get(http_url + folder_date).text
    file_list = [file for file in re.findall('(?<=<a href=")[^"]*', folder_content) if 'xml' in file]  # list of all xml files, for finding the production id

    for ind, inter in enumerate(intercepts):
        t0 = time.time()

        rounded_time = inter['MISR_time'] - datetime.timedelta(minutes=inter['MISR_time'].minute % 5)
        gran_time = rounded_time.strftime('%Y%j.%H%M')

        # find identifier from folder file list
        for file in file_list:
            if 'A{0}'.format(gran_time) in file:
                identifier = file[24:37]
                break

        try:
            gran_filename = 'MOD14.A{date}.006.{id}.hdf'.format(date=gran_time, id=identifier)
        except UnboundLocalError:
            mssg = 'no MOD14 granule found with such a time: {0}'.format(gran_time)
            sys.stdout.write('Warning: ' + mssg + '; skipping intercept\n')
            continue
            # raise UnboundLocalError(mssg)

        if gran_filename in grans:
            continue
        else:
            grans.append(gran_filename)

        # Read L3 fire product file
        L3_path = '/fs/site2/dev/eccc/aq/r1/nod001/MINX/Comp/workdir/MCD14ML/MCD14ML.{Yr}.006.01.nc'.format(Yr=inter['MISR_time'].strftime('%Y'))
        L3_file = netCDF4.Dataset(L3_path)

        MISR_start_srch_time = float(inter['MISR_time'].strftime('%H%M'))
        MISR_end_srch_time = float(inter['MISR_end_time'].strftime('%H%M'))
        day_ind_rng = np.where(L3_file['YYYYMMDD'][0] == float(inter['MISR_time'].strftime('%Y%m%d')))[0][[0, -1]]

        # TODO 5 or ten extra minutes?
        time_inds = np.where(
            (MISR_start_srch_time - 5 <= L3_file['HHMM'][0][day_ind_rng[0]:day_ind_rng[-1]]) &
            (MISR_end_srch_time + 5 >= L3_file['HHMM'][0][day_ind_rng[0]:day_ind_rng[-1]])
        )[0]

        tot_inds = time_inds + day_ind_rng[0]

        MISR_poly_plus = draw_swath(inter['MISR_path'], 3.3)  # TODO I think the size of this may be leading to false positives

        FRP_vals = L3_file['FRP'][0]

        min_int_lat = min(inter['coordinates'], key=lambda coord: coord[1])[1]
        max_int_lat = max(inter['coordinates'], key=lambda coord: coord[1])[1]
        polypoints = [(x[0], x[1]) for x in MISR_poly_plus if min_int_lat < x[1] < max_int_lat]

        total_FRP = 0
        FRP_coords = []
        for lind, (lon, lat) in enumerate(zip(L3_file['lon'][0][tot_inds], L3_file['lat'][0][tot_inds])):

            if Cpoly.polygon_check_single(polypoints, (lon, lat)):
                total_FRP += FRP_vals[tot_inds[lind]]
                FRP_coords.append((lon, lat))

                # TODO this wont grab all the pixels, only those until we reach the threshold
                # if this isnt the bottleneck of this func., should probably do all pixels so we can see them all on plot

                # <--

        if total_FRP > FRP_threshhold:
            inter['gran_filename'] = gran_filename
            inter['pixels'] = FRP_coords
            res.append(inter)
            # break

    return res


def reduce_intercepts(path_ints, basemap, cut_africa=False):
    # attempt to reduce false positives TODO make more reduction methods and improve the quality of these, and speed
    intercepts = []
    for inter in path_ints:
        if args.reduce:
            if not any([basemap.is_land(*coord) for coord in inter['coordinates']]):  # might lose a few fires this way, near coasts, but whatever
                continue
            # TODO a polygon check may be a bit extra for just a constant box
            if cut_africa and Cpoly.polygon_check_single([(-10., 10.), (50., 10.), (50., -30.), (-10., -30.)], inter['coordinates'][0]):
                continue

        intercepts.append(inter)
    return intercepts


def user_check(path_num, orb_num, inter, terr_file):
    '''Take an intercept and MISR file (AF) and display it to the user, who can decide if there is a valid plume'''

    # first block is northernmost block, because blocks descend from the north in number
    first_block, *_ = mu.LonLattoBLS(path_num, *max(inter['coordinates'], key=lambda coord: coord[1]))
    last_block, *_ = mu.LonLattoBLS(path_num, *min(inter['coordinates'], key=lambda coord: coord[1]))

    try:
        savepath = '{dir}/O{orb:06d}_B{sb:03d}_B{eb:03d}.png'.format(dir=directory, orb=orb_num, sb=first_block, eb=last_block)
        BRF, Bulc, Blrc, offset = mu.drawmisr(directory + '/' + terr_file, first_block - 1, last_block + 1, filepath=savepath)
    except NameError:
        mssg = 'No blocks on this path ({path}) match the intercept near {coord}'.format(path=path_num, coord=inter['coordinates'][0])
        raise RuntimeError(mssg)

    # TODO these things refuse to line up right, sometimes the sample is just plain off, I think maybe to do with offset? latitude?

    # convert caliop path and fire pixel coordinates to BLS
    int_B, int_L, int_S = zip(*map(lambda coord: mu.LonLattoBLS(path_num, coord[0], coord[1]), inter['coordinates']))
    cal_B, cal_L, cal_S = zip(*map(lambda coord: mu.LonLattoBLS(path_num, coord[0], coord[1]), inter['CALIOP_path']))
    frp_B, frp_L, frp_S = zip(*map(lambda coord: mu.LonLattoBLS(path_num, coord[0], coord[1]), inter['pixels']))

    mis_B, mis_L, mis_S = zip(*map(lambda coord: mu.LonLattoBLS(path_num, coord[0], coord[1]), inter['MISR_path']))
    pol_B, pol_L, pol_S = zip(*map(lambda coord: mu.LonLattoBLS(path_num, coord[0], coord[1]), draw_swath(inter['MISR_path'], 3.3)))

    # plotting over N blocks, therefore y-axis = Line + (blocknum[i] - blocknum[0]) * 512, x-axis = Sample
    int_BL = [L + (B - (first_block - 1)) * 512. for B, L in zip(int_B, int_L)]
    cal_BL = [L + (B - (first_block - 1)) * 512. for B, L in zip(cal_B, cal_L)]
    frp_BL = [L + (B - (first_block - 1)) * 512. for B, L in zip(frp_B, frp_L)]

    mis_BL = [L + (B - (first_block - 1)) * 512. for B, L in zip(mis_B, mis_L)]
    pol_BL = [L + (B - (first_block - 1)) * 512. for B, L in zip(pol_B, pol_L)]

    # TODO redo this entirely, using same method as is done in MISRutils, stop returning offsets from there
    # adjust for offset from image, then recorrect for relative offset of blocks
    int_S = [(S - offset) + mu.rel_offset_table[last_block - int_B[Sind] - 1] for Sind, S in enumerate(int_S)]
    cal_S = [(S - offset) + mu.rel_offset_table[last_block - cal_B[Sind] - 1] for Sind, S in enumerate(cal_S)]
    frp_S = [(S - offset) + mu.rel_offset_table[last_block - frp_B[Sind] - 1] for Sind, S in enumerate(frp_S)]

    mis_S = [(S - offset) + mu.rel_offset_table[last_block - mis_B[Sind] - 1] for Sind, S in enumerate(mis_S)]
    pol_S = [(S - offset) + mu.rel_offset_table[last_block - pol_B[Sind] - 1] for Sind, S in enumerate(pol_S)]

    # get relative offset from table and roll data to eliminate
    # rel_offset = 0
    # bb = None
    # frp_S = list(frp_S)
    # for sind, sample in enumerate(frp_S):
    #     if frp_B[sind] != bb:
    #         bb = frp_B[sind]
    #         rel_offset += mu.rel_offset_table[frp_B[sind] - 1]

    #     frp_S[sind] = sample - (rel_offset * 4)

    plt.plot(cal_S, cal_BL, 'bo-')
    plt.plot(int_S, int_BL, 'cs')
    plt.plot(frp_S, frp_BL, 'ro')

    plt.plot(mis_S, mis_BL, 'ko-')
    plt.plot(pol_S, pol_BL, 'go-')

    plt.imshow(np.ma.dstack([BRF['Red'], BRF['Green'], BRF['Blue'], ~BRF['NIR'].mask]))

    plt.title(Bulc)

    plt.show()

    resp = input('Intercept in plume?  Y / [N] ')

    return True if 'y' in resp.lower() else False


def download_intercepts(terr_url, terr_dir, intercepts, loader=False):
    larc = ftplib.FTP(terr_url)
    larc.login()
    larc.cwd(terr_dir)

    for ind, inter in enumerate(intercepts):
        # download the misr AF terrain file and display the real colour image with calipso path overlaid, then ask the user if its good

        # TODO sometimes the day is off, ie orbit doesnt exist in day
        orb_num = mu.TimetoOrbit(inter['MISR_time'])
        path_num = mu.OrbittoPath(orb_num)

        terr_file = 'MISR_AM1_GRP_TERRAIN_GM_P{path:03d}_O{orbit:06d}_AF_F03_0024.hdf'.format(path=path_num, orbit=orb_num)

        with open('{dir}/{file}'.format(dir=directory, file=terr_file), 'wb') as f:
            larc.retrbinary('RETR ' + terr_file, f.write)

        try:
            loader.pause()
        except AttributeError:
            pass

        # draw image of intercept and wait for user visual verification
        if not user_check(path_num, orb_num, inter, terr_file):
            sys.stdout.write('    Skipping intercept\n')
            os.remove('{dir}/{file}'.format(dir=directory, file=terr_file))

            try:
                loader.restart()
            except AttributeError:
                pass

            continue

        try:
            loader.restart()
        except AttributeError:
            pass

        inter_url = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD14.006/{date}/{gran}'.format(date=day.strftime('%Y.%m.%d'), gran=inter['gran_filename'])
        redirected_r = requests.request('get', inter_url)
        MODIS_r = requests.request('get', redirected_r.url, stream=True, auth=('DicksonN', 'SSHSBench2015!'))

        with open("{dir}/granules/{gran}".format(dir=directory, gran=inter['gran_filename']), 'wb') as f:
            for chunk in MODIS_r.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)

    larc.close()

    with open('{dir}/fire_intercepts_{date}.p'.format(dir=directory, date=day.strftime('%Y%m%d')), 'wb') as inter_file:
        pickle.dump(intercepts, inter_file)

    preMINX_bng = "#!/bin/sh\n"
    preMINX_cmd = "preMINX -L -180 180 -90 90 -d {date} ".format(date=day.strftime('%Y.%m.%d'))
    preMINX_opt = "-u DicksonN -p SSHSBench2015! --projdir {dir}/MINX/ --relaxed --localgran --grandir {dir}/granules/".format(dir=directory)

    with open('{dir}/preMINX_cmd_{date}.sh'.format(dir=directory, date=day.strftime('%Y%m%d')), 'w') as preMINX_file:
        preMINX_file.write(preMINX_bng + preMINX_cmd + preMINX_opt)
    os.chmod(preMINX_file.name, 0o775)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Find CALIOP/MISR/Fire Granule intercepts")

    parser.add_argument('-d', '--date', required=True, nargs='*', metavar='YYYY.MM.DD',
                        help="Date, must have format YYYY.MM.DD. Can also search over range of dates, \
                        simply by providing a second, later date.")
    parser.add_argument('--reduce', action='store_true',
                        help="Attempt to reduce false positives. Will eliminate over-ocean intersects and cut africa if specified.\
                        May cause the loss of real fires occasionally.")
    parser.add_argument('--cut-africa', action='store_true',
                        help="Don't accept intercepts near central Africa, in order to avoid the large amounts of plume-less biomass burnings.\
                        This is simply an attempt to reduce false positives.")
    args = parser.parse_args()

    argdates = [datetime.datetime.strptime(date, '%Y.%m.%d') for date in args.date]
    s = argdates[0]
    dates = []
    while s <= argdates[-1]:
        dates.append(s)
        s += datetime.timedelta(days=1)

    for day in dates:
        directory = "/fs/site2/dev/eccc/aq/r1/nod001/MINX/Comp/workdir/{date}".format(date=day.strftime('%Y/%m/%d'))

        if not os.path.exists(directory):
            os.makedirs(directory)
        if not os.path.exists(directory + '/granules/'):
            os.makedirs(directory + '/granules/')

        sys.stdout.write('Date: {0}\n'.format(day.strftime('%Y, %m, %d')))

        # timing and loading #
        load = pyloader.Load_message('-> drawing orbit')

        T_paths, C_arr = draw_orbits(day)

        load.update('-> finding intercepts')

        intercepts = []
        if args.reduce:
            bm = Basemap(llcrnrlat=-90, urcrnrlat=90, llcrnrlon=-180, urcrnrlon=180)  # this is just an attempt to reduce false positives
        # refract this to do orbit_intercept on a list of paths so we can move the bm definition into reduce_int
        for path in T_paths:
            try:
                path_ints = orbit_intercept(path, C_arr)
                intercepts += reduce_intercepts(path_ints, bm, args.cut_africa) if args.reduce else path_ints
            except IndexError:
                continue

        if not intercepts:
            load.stop()
            sys.stdout.write("No potential intercepts found, skipping day.\n")

        # TODO fancier printing, maybe some kind of table with more info
        sys.stdout.write("{0} potential CALIOP/MISR intercepts found.\n".format(len(intercepts)))
        load.update('-> testing for FRP')

        fire_intercepts = fire_pixel_check(intercepts)

        sys.stdout.write("{0} intercepts flagged for investigation.\n".format(len(fire_intercepts)))
        load.update('-> downloading files')

        terr_url, terr_dir = 'l5ftl01.larc.nasa.gov', 'MISR/MI1B2T.003/{date}/'.format(date=day.strftime('%Y.%m.%d'))
        if fire_intercepts:
            # TODO move AF to MISRFiles after accepted
            download_intercepts(terr_url, terr_dir, fire_intercepts, load)

        load.stop()

        sys.stdout.write("Total runtime: {t} seconds\n".format(t=load.delta_t))