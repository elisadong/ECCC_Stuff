#! /usr/bin/env python

# function to use FIRECAL intercept results
# NOTE:
# Possible to have up to two overpasses for a fire on one day, more if delta_t is a few days
# Only retain overpasses that are during the day time

# Requirements:
## Must have intercept file from FIRECAL_intersect.py results
## NOTE: FIRECAL_intersect must be run in python 3 environment
### Looks like: intercepts_s20170901_e20170930.p
# Suggestions:
## FIRECAL_intersect should be run with:
### padding = 10000 #10km padding around fire cluster
### frp = 300 #300MWatts per fire cluster at the minimum to be considered a 'fire'
### delat_t = 12 #default value, search +/- 12 hours of fire cluster to find overpass
### radius = 3000 #default value, fire pixels within radius of each other will be considred part of the same 'fire'

# example for use (in interactive python):
# import buildCALFires as calFires
# fireDictionary = calFires.cal_intersects()
# fireDictionary2 = calFires.cal_intersects(savePlts = True, showPlot = False, firePath = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/CALIOP/intercepts_s20170401_e20170430.p')
## Example for september fires
# import os
# fireDictionary3 = calFires.cal_intersects(savePlts=True, showPlot = True, fst_dir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/CFFEPSLowResMant4165", ctrl_dir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/GEMMACHOPLowRes', firePath = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/CALIOP/intercepts_s20170901_e20170930.p', cal_path = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/CALIOP', savePath = os.path.join(cal_path, 'Images'), threshold = 5, filestart = 00)

# example for running FIRECAL_intersect:
# python3.4 FIRECAL_intersect.py -d 2017.07.01 2017.07.31 -L -142.0 -50.7 42.0 83.0 -f 300 -p 10000

import os
import datetime
import numpy as np
import niceUtilities as nu
import rpnpy.librmn.all as rmn
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry import Polygon

try:
    import cPickle as pickle
except ImportError:
    import pickle

# VARIABLES
savePlts = True #set to True to save images and pkl files
showPlot = True #set to True to show the image and save as you like
#fst_dir='/fs/site2/dev/eccc/aq/r1/eld001/MINX/fireworkOPLowRes' #FIREWORK
fst_dir= "/fs/site2/dev/eccc/aq/r1/eld001/MINX/CFFEPSLowResMant4165" #CFFEPS directory
ctrl_dir='/fs/site2/dev/eccc/aq/r1/eld001/MINX/GEMMACHOPLowRes'
firePath = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/CALIOP/intercepts_s20170901_e20170930.p'
cal_path = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/CALIOP'
savePath = os.path.join(cal_path, 'Images')
threshold= 5
filestart = 00 #this is the FST file start time
dayTimeRange= [8, 20] # want to keep this to a 12 hour range to represent day light hours
# NOTE: CALIOP overpasses can potentially occur at the same location twice a day.
# dayTimeRange ensures that we pick the 'daylight' version and only have one overpass per day

# NOTE: Include function to get the fire pickle files
# include examples

# CALIOP stuff
def get_cal(fire, ftype = 3, cal_path=os.getcwd(),auth=('elisadong', 'elisaPass1')):
    '''
    Retrieve CALIOP data for a fire in the form of 'overpass' dictionaries.

    It is possible to have more than one overpass per fire (depends on the delta_t setting when running FIRECAL_intersect.py).
    '''

    import ftplib
    import pyhdf.SD as SD

    # deal with potentially more than one overpass per fire
    overpasses = []
    for overTime in fire.overpasses:
        # Pick out the overpasses that happen during the day time.
        if overTime.hour >= dayTimeRange[0] and overTime.hour < dayTimeRange[1]:
            overpasses += [overTime]

    overpassDict = {}
    for overpass in overpasses:
        ftp_url = 'ftp.icare.univ-lille1.fr'
        ftp_dir = 'SPACEBORNE/CALIOP/VFM.v4.10/{yr}/{date}/'.format(yr=overpass.year, date=overpass.strftime('%Y_%m_%d'))
        icare = ftplib.FTP(ftp_url, user=auth[0], passwd=auth[1])
        icare.cwd(ftp_dir)
        file_list = icare.nlst()

        # find the best fitting hdf file, should already have been downloaded
        vind, vfile = min(enumerate(file_list),
                          key=lambda filename: abs(datetime.datetime.strptime(filename[1][30:49], '%Y-%m-%dT%H-%M-%S') - overpass
                                                   if datetime.datetime.strptime(filename[1][30:49], '%Y-%m-%dT%H-%M-%S') < overpass
                                                   else datetime.timedelta(999)))
        delta_t = datetime.datetime.strptime(vfile[30:49], '%Y-%m-%dT%H-%M-%S') - overpass
        if 0 < delta_t.seconds < 300:
            inter_files = [vfile, file_list[vind + 1]]
        elif -300 < delta_t.seconds < 0:
            inter_files = [file_list[vind - 1], vfile]
        else:
            inter_files = [vfile]

        # Look for the file in the directory, otherwise download the file
        for find, ftp_VFM_file in enumerate(inter_files):
            local_VFM_path = '{dir}/{file}'.format(dir=cal_path, file=ftp_VFM_file)
            if not os.path.exists(local_VFM_path):  # if dont have locally, download
                with open(local_VFM_path, 'wb') as f:
                    icare.retrbinary('RETR ' + ftp_VFM_file, f.write)
            vfile = SD.SD(local_VFM_path, SD.SDC.READ)

        icare.close() #close the FTP connection

        # Load in the lon/lat data for the overpass
        lon_arr = vfile.select('ssLongitude')[:, 0]
        lat_arr = vfile.select('ssLatitude')[:, 0]

        polypoints = fire.boundaries
        polyBounds = Polygon(polypoints)
        # polyBounds is a pretty small area, can increase size with a buffer
        # NOTE: can also increase the size of polyBounds using the 'padding' argument from the FIRECAL_intersect function
        ## This will result in larger boundaries for the fire objects inside the pickle file
        newPolyBounds = Polygon(polyBounds.buffer(1.0).exterior)

        # Keep indices of the lon/lat arrays that are within the boundaries of the newPolyBounds area
        keepInds = []
        for pind, (lon, lat) in enumerate(zip(lon_arr, lat_arr)):
            llPoint = Point(lon,lat)
            if newPolyBounds.contains(llPoint):
                keepInds += [pind]

        keepInds = np.array(keepInds)

        # Load in the feature classification flag data
        fcf_arr = vfile.select('Feature_Classification_Flags').get()

        # Break up the fcf_arr into 3 layers corresponding to the CALIOP retrieval breakdown
        layer2 = fcf_arr[:, 0:165]
        layer1 = fcf_arr[:, 165:1165]
        layer0 = fcf_arr[:, 1165:5515]

        sslayer2 = np.flip(layer2.reshape(layer2.shape[0] * 3, 55), axis=1)
        sslayer1 = np.flip(layer1.reshape(layer1.shape[0] * 5, 200), axis=1)
        sslayer0 = np.flip(layer0.reshape(layer0.shape[0] * 15, 290), axis=1)

        # store the fcf data into array that has reasonable dimensions
        ## ie. each layer now has the same number of columns
        ss_fcf_arr = np.empty((lon_arr.shape[0], 545), dtype=np.uint16)
        ss_fcf_arr[:, 490:545] = sslayer2.repeat(5, axis=0)
        ss_fcf_arr[:, 290:490] = sslayer1.repeat(3, axis=0)
        ss_fcf_arr[:, :290] = sslayer0

        # Only keep columns of data that are within the newPolyBounds area
        ss_fcf_arr = ss_fcf_arr[keepInds, :]

        # Get the feature types and indices
        type_arr = ss_fcf_arr & 7
        feat_inds = np.where((type_arr == ftype))
        inds_1D = np.unique(feat_inds[0])

        # Mask indices without desired feature
        mask = np.ones_like(ss_fcf_arr)
        mask[feat_inds] = 0
        mask = mask[inds_1D, :]
        fcf_ones = np.ones_like(ss_fcf_arr[inds_1D, :])
        feat_mask = np.ma.masked_array(fcf_ones, mask=mask, dtype=np.float32)

        # Build altitude array using the inds_1D
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

        # apply mask to altitudes.
        ## Mask = True hides the corresponding altitudes without the desired features
        alt_mask = feat_mask * ss_alt * 1000.
        lon_arr = lon_arr[keepInds][inds_1D]
        lat_arr = lat_arr[keepInds][inds_1D]

        # Get the maximum altitudes where we see our feature
        max_alt = alt_mask.max(axis=1)

        # save results to our dictionary
        overpassDict[overpass] = {'lons': lon_arr, 'lats': lat_arr, 'calHeights': max_alt.data}
    return overpassDict

# Model stuff
def data_from_key(key, file_id, lonarr, latarr):
    meta = rmn.fstprm(key)
    meta['iunit'] = file_id
    data = rmn.fstluk(key)['d']
    grid = rmn.ezqkdef(meta)
    xypos = rmn.gdxyfll(grid, latarr, lonarr)
    val = rmn.gdxysval(grid, xypos['x'], xypos['y'], data)
    return meta, grid, xypos, val

def get_model(overpass, lonArray,latArray,fst_dir, ctrl_dir=None, var='AF', threshold=4, filestart=00):
    if filestart == 12:
        file_time = (overpass if overpass.minute < 30 else overpass+ datetime.timedelta(hours=1)) - datetime.timedelta(hours=12)
        file_path = "{0}/{1}".format(os.path.abspath(fst_dir), file_time.strftime('%Y%m%d12_0%H'))
    elif filestart == 0:
        file_time = (overpass if overpass.minute < 30 else overpass + datetime.timedelta(hours=1))
        file_path = "{0}/{1}".format(os.path.abspath(fst_dir), file_time.strftime('%Y%m%d00_0%H'))
    rmn.fstopt(rmn.FSTOP_MSGLVL, rmn.FSTOPI_MSG_CATAST)
    rmn.ezsetopt(rmn.EZ_OPT_INTERP_DEGREE, rmn.EZ_INTERP_NEAREST)
    fid = rmn.fstopenall(file_path, rmn.FST_RO)
    lonarr = np.array(lonArray,dtype='float32')
    latarr = np.array(latArray,dtype='float32')
    keylist = rmn.fstinl(fid, nomvar=var)
    dir = os.path.abspath(fst_dir)
    ctrl_dir = os.path.abspath(ctrl_dir) if ctrl_dir else None
    # assuming this for now
    ref_lvl = 'sea'
    height = [float('-inf')] * len(lonarr)
    value = [float('-inf')] * len(lonarr)
    points = [{} for _ in lonarr]
    iplist = []
    for key in keylist:
        meta = rmn.fstprm(key)
        iplist.append(rmn.DecodeIp(meta['ip1'], meta['ip2'], meta['ip3'])[0].v1)
    sorted_keylist = (x for _, x in sorted(zip(iplist, keylist), reverse=True))
    next(sorted_keylist, None)
    before_val = [float('inf')] * len(lonarr)
    cur_meta, cur_grid, cur_xypos, cur_val = data_from_key(next(sorted_keylist), fid, lonarr,latarr)
    if ctrl_dir:
        if filestart == 12:
            ctrl_path = "{0}/{1}".format(os.path.abspath(ctrl_dir), file_time.strftime('%Y%m%d12_0%H'))
        elif filestart == 0:
            ctrl_path = "{0}/{1}".format(os.path.abspath(ctrl_dir), file_time.strftime('%Y%m%d00_0%H'))
        ctrl_fid = rmn.fstopenall(ctrl_path, rmn.FST_RO)
        ctrl_keylist = rmn.fstinl(ctrl_fid, nomvar=var)
        sorted_ctrl_keylist = (x for _, x in sorted(zip(iplist, ctrl_keylist), reverse=True))
        next(sorted_ctrl_keylist, None)
        _, _, _, ctrl_val = data_from_key(next(sorted_ctrl_keylist), ctrl_fid,lonarr,latarr)
        cur_val -= ctrl_val
    for progress_ind, after_key in enumerate(sorted_keylist):
        after_meta, after_grid, after_xypos, after_val = data_from_key(after_key, fid,lonarr,latarr)
        if ctrl_dir:
            after_ctrl_key = next(sorted_ctrl_keylist)
            _, _, _, after_ctrl_val = data_from_key(after_ctrl_key, ctrl_fid,lonarr,latarr)
            after_val -= after_ctrl_val
        for ind, val in enumerate(cur_val):
            if ((val > before_val[ind]) and (val > after_val[ind]) and (val >= threshold) and (val > value[ind])):
                try:
                    if int(ind) <= 20:
                        print('Updating GZ, val: {}, existing val: {}'.format(val, value[ind]), ind, cur_meta['ip1'])
                    gzkey = rmn.fstinf(fid, nomvar='GZ', ip1=cur_meta['ip1'])['key']
                    gzdata = rmn.fstluk(gzkey)['d']
                    meta = rmn.fstprm(gzkey)
                    meta['iunit'] = fid
                    gz_grid = rmn.ezqkdef(meta)
                    heightList = rmn.gdxysval(gz_grid, cur_xypos['x'], cur_xypos['y'], gzdata) * 10
                    height[ind] = float(heightList[ind])
                    value[ind] = float(val)
                    #print (height[ind], ind)
                    print (height, value)
                except TypeError:
                    continue
        before_val = cur_val
        cur_meta, cur_grid, cur_xypos, cur_val = after_meta, after_grid, after_xypos, after_val
        print height
    rmn.fstcloseall(fid)
    if ctrl_dir:
        rmn.fstcloseall(ctrl_fid)
    #print(height)
    return height, value

# Plots for reference
def scatter_3D(lons,lats, calHeights, modHeights, fireBoundaries, savePath = os.getcwd(), title = '', showPlot = showPlot):
    from mpl_toolkits.mplot3d import Axes3D
    import copy

    #fireBoundaries is a polygon
    if os.path.exists(savePath) == False:
        nu.makeDir(savePath)
        print('Adding new save path {}'.format(savePath))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.autoscale()
    ax.plot(lons,lats,calHeights, 'o',  label='CALIOP')
    ax.plot(lons,lats,modHeights, 'o', label = 'Model')
    x,y = fireBoundaries.exterior.xy
    ax.plot(x,y,color='r',label='FireBound')
    ax.legend()
    ax.set_title(title)
    try:
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    except:
        pass

    if showPlot == True:
        plt.show()
    #name and save the image to png and pkl file
    filename =  "".join(title.split())
    filename = '_'.join(filename.split(':'))
    pklfig = copy.copy(fig)
    fig.savefig(os.path.join(savePath, filename + '.png'))
    pklfile = '{}.pkl'.format(filename)
    with open(os.path.join(savePath, pklfile), 'wb') as f:
        pickle.dump(pklfig,  f, protocol = 2)
    return

# Actual function
# NOTE: potentially load fires from a FIRECAL_intersect result
def cal_intersects(cal_path = cal_path, firePath = firePath, fst_dir=fst_dir,ctrl_dir=ctrl_dir, threshold=threshold,filestart=filestart, savePlts = savePlts, showPlot = showPlot):
    fires = pickle.load(open(firePath,'rb'))
    keptFires = []

    for f in fires:
        if f.overpasses != []:
            keptFires += [f]

    fireDict = {}
    for k in keptFires:
        try:
            fireDict[k] = get_cal(k, cal_path = cal_path)
            passNum = 1
            for overpass in fireDict[k]:
                lons = fireDict[k][overpass]['lons']
                lats = fireDict[k][overpass]['lats']
                modHeights, modVals = get_model(overpass, lons,lats,fst_dir, ctrl_dir=ctrl_dir, threshold=threshold, filestart=filestart, var='AF')
                print(modHeights)
                fireDict[k][overpass]['modHeights'] = modHeights
                fireDict[k][overpass]['modVals'] = modVals
                if savePlts == True:
                    fireBounds = k.boundaries
                    calHeights = fireDict[k][overpass]['calHeights']
                    title = 'FirePass: {}\nOverpass Date: {}'.format(passNum, overpass)
                    scatter_3D(lons,lats, calHeights, modHeights, Polygon(fireBounds), savePath, title)
                passNum += 1
        except:
            continue
    return fireDict
