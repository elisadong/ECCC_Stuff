# Adjustment to the monthly means program
# Will attempt to borrow the from the original program as much as possible
# NOTE: I've done this so it deals with only one year at a time
# This has been tested in the python 2.7 environment using anaconda

# Example for use:
#import seasonal_means as seas
#season = ['Dec', 'Jan','Feb']
#outputName = 'DJF'
#spcsFiles = ["/fs/site2/dev/eccc/aq/r1/eld001/testFSTs/speciesNetCDFs/GO3/2010.nc", "/fs/site2/dev/eccc/aq/r1/eld001/testFSTs/speciesNetCDFs/CO/2010.nc"]
#spcs = ['go3', 'co']
#seas.plotFSTs(season=season, outputName = outputName, spcsFiles=spcsFiles, spcs=spcs)

#imports

import pygeode as pyg # need to get this, conda install -c aph42 pygeode
import numpy as np
import rpnpy.librmn.all as rmn
import time, glob
import matplotlib.pyplot as plt
import rpnpy.vgd.all as vgd
from mpl_toolkits.basemap import Basemap
from levels import * #this is from a python program that has a series of stored variables
import niceUtilities as nu
import os

#VARIABLES/Defaults:
# List of months in season
season = ['Jun','Jul','Aug']
outputName = 'JJA'
saveDir = '/fs/site2/dev/eccc/aq/r1/eld001/testFSTs/outputFSTs'

#netCDF files, one per year
# only one lnsp file, should point to the appropriate year
lnsp_file = "/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/LNSP/2010.nc"
#number of species files must match up to number of spcs
spcsFiles = ["/fs/site2/dev/eccc/aq/r1/eld001/testFSTs/speciesNetCDFs/GO3/2010.nc",
            "/fs/site2/dev/eccc/aq/r1/eld001/testFSTs/speciesNetCDFs/CO/2010.nc",
            "/fs/site2/dev/eccc/aq/r1/eld001/testFSTs/speciesNetCDFs/NOX/2010.nc",
            "/fs/site2/dev/eccc/aq/r1/eld001/testFSTs/speciesNetCDFs/HCHO/2010.nc",
            "/fs/site2/dev/eccc/aq/r1/eld001/testFSTs/speciesNetCDFs/SO2/2010.nc"]
spcs = ['go3', 'co', 'nox', 'hcho', 'so2']

# Conversions for various species.
# Feel free to add to the list (add name to allSpcs, and weight to molecWeights).

allSpcs = ['go3', 'co', 'nox', 'hcho', 'so2']
molecWeights = [48, 28.01, 30, 30.031, 64.066] # for each species, g/mol

#Don't touch this. Left the function in this location so you can track how the conversion works
scaleSpcs = list((1e9 * 28.97/weight) for weight in molecWeights) # now in ppmv

# Should probably not touch these hardcoded values
# Some useful ways of breaking down the nc file (and indices for dates)
monthList = ['Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Oct','Sep','Nov','Dec']
indRangeLeap = [[0,124], [124,240], [240,364], [364,484],[484,608], [608,728], [728,852], [852,976],[976,1096], [1096,1220], [1220,1340], [1340,1464]]
indRangeNoLeap = [[0,124], [124,236], [236,360], [360,480],[480,604], [604,724], [724,848], [848,972],[972,1092], [1092,1216], [1216,1336], [1336,1460]]

# get species(variable) data from netcdf
def get_var(open_file, nomvar, mInds):
    repr_var = "<Var '{0}'>".format(nomvar)
    print("Getting " + repr_var)
    var_data = splice_month([data for data in open_file.vars if repr(data) == repr_var][0], mInds)
    return var_data

# pass in two 4D arrays, modify the second dimension (which stores the vertical levels)
def vert_interp(pressure, org):
    global const_pressure
    print "Conducting vertical interpolation..."
    lon         = len(org[0, 0, 0])
    lat         = len(org[0, 0])
    lev         = len(org[0])
    timesteps   = len(org)
    y_interp    = np.zeros((timesteps, len(const_pressure), lat, lon))
    x_interp    = const_pressure
    for timestep in range(timesteps):
        for j in xrange(lon):
            for i in xrange(lat):
                try:
                    x_initial = pressure[timestep, :, i, j]
                    y_initial = org[timestep, :, i, j]
                    y_interp[timestep, :, i, j] = np.interp(x_interp, x_initial, y_initial)
                except:
                    print ('Issue with vertical interpretation. Pressure shape : {}'.format(pressure.shape))
                    raise
    return y_interp

# Cut up the netCDF data(now in arrays) by month indices
# Refers to the indRanges that are provided for leap and non leap years
def splice_month(open_var, mInds):
    print "Splicing..."
    if len(open_var.time) == 1464:
        leap_year = True
    elif len(open_var.time) == 1460:
        leap_year = False
    else:
        print("File must contain all 4 time values for each day of the year")
    #NOTE: MODIFY THE STORED ARRAY. Check to see which ints are wanted, then
    ## create the entire array
    if leap_year:
        nc_data = []
        for mInd in mInds:
            print('Month Index is: {}, spliceRange is{},{}: {}:{}'.format(mInd, indRangeLeap[mInd][0], indRangeLeap[mInd][1]))
            start,end = indRangeLeap[mInd][0], indRangeLeap[mInd][1]
            mthData = open_var[start:end]
            nc_data += [mthData]
        seasonData = np.concatenate(nc_data)
    else:
        nc_data = []
        for mInd in mInds:
            print('Month Index is: {}, spliceRange is:{}.{}'.format(mInd, indRangeNoLeap[mInd][0], indRangeNoLeap[mInd][1]))
            start,end = indRangeNoLeap[mInd][0], indRangeNoLeap[mInd][1]
            mthData = open_var[start:end]
            nc_data += [mthData]
        seasonData = np.concatenate(nc_data)
    print('Data for Season is now in shape: {}'.format(np.shape(seasonData)))
    return seasonData

#NOTE: check to see if this is just one nc file, or the merged one
# Gets pressure, 4D array?
def get_pressures(open_nc, mInds):
    # TODO: check the units
    global Ak, Bk, const_pressure
    print "Getting pressures..."
    lnsp = open_nc.lnsp
    lnsp_array = splice_month(lnsp, mInds)
    print lnsp_array.shape
    shape = list(lnsp_array.shape); shape.insert(1, 60)  # creates a dimension for levels
    print "PRESSURE SHAPE IS {0}".format(str(shape))
    pressure = np.zeros(shape)
    try:
        for lev in xrange(len(pressure[0])):
            pressure[:, lev, :, :] = ((Ak[lev] + Bk[lev]*np.exp(lnsp_array)) / 100.)
    except:
        debug_tuple = (lev)
        print "LEV: {}".format(lev)
        raise
    return pressure

#NOTE: This function appears to be used for basemap plotting
def shift_lon(array):
    global timelen, latlen, lonlen
    new_array = np.zeros(array.shape)
    for k in xrange(60):
        tmp = np.zeros((timelen, latlen, lonlen))
        tmp[:, :, :lonlen/2], tmp[:, :, lonlen/2:] = array[:, k, :, lonlen/2:], array[:, k, :, :lonlen/2]
        new_array[:,k] = tmp
    return new_array

def plotFSTs(season=season, spcs=spcs, spcsFiles=spcsFiles, outputName = outputName, saveDir=saveDir):
    # print minimum outputs
    rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)

    mInds = []
    for m in season:
        mInds += [monthList.index(m)]

    if os.path.exists(saveDir) == False:
        nu.makeDir(saveDir)

    for spcInd, nomvar in enumerate(spcs):
        try:
            filename = os.path.join(saveDir, 'output_file_{0}_{1}.fst'.format(outputName, nomvar))
            print('Creating and saving to {}'.format(filename))
            tmp = open(filename, 'w+'); tmp.close()

            output_file = filename
            file_id = rmn.fnom(output_file)
            open_fst = rmn.fstouv(file_id, rmn.FST_RW)
            open_file = spcsFiles[spcInd]
            print "Parameter: " + nomvar
            seaSpcData = get_var(pyg.open(open_file), nomvar, mInds)
            nc_lnsp = pyg.open(lnsp_file)
            pressures = get_pressures(nc_lnsp, mInds)

            timelen, levlen, latlen, lonlen = seaSpcData.shape
            #NOTE: uncomment the following three lines to prep data for basemap use
            #lonShiftSSData = shift_lon(seaSpcData)
            #vertInterpSSData = vert_interp(pressures, lonShiftSSData)
            #meanSSData = np.mean(vertInterpSSData, axis=0)
            #NOTE: uncommment the following four liness to use for fst plotting
            vertInterpSSData = vert_interp(pressures, seaSpcData)
            meanSSData = np.mean(vertInterpSSData, axis=0)  # temp
            for lvl, ray in enumerate(meanSSData):
                meanSSData[lvl] = np.flipud(ray)
            scaleFac = scaleSpcs[allSpcs.index(nomvar)]
            scaledSSData = meanSSData*scaleFac

            #define grid for this file - note that the MACC grid in the file is
            #defined for lons -180 to 180, but the python defGrid_L can't deal
            #with that and defines the grid from 0 to 360 so will have to reorder
            #the MACC fields a bit, or they end up 180 deg out of phase
            # Also necessary to add one more longitude to wrap around
            dlatlon = 360./lonlen   # this is equal to the resolution of the grid

            params0 = {
                    'grtyp' : 'Z',
                    'grref' : 'L',
                    'nj'    : latlen,
                    'ni'    : lonlen,
                    'lat0'  : -90.,
                    'lon0'  : 0.,
                    'dlat'  : dlatlon,
                    'dlon'  : dlatlon
                    }

            MACC_grid= rmn.encodeGrid(params0)
            print("Grids created.")
            print 'Grid Shape:' + str(MACC_grid['shape'])

            # copies the default record
            new_record = rmn.FST_RDE_META_DEFAULT.copy()
            tic_record = rmn.FST_RDE_META_DEFAULT.copy()
            tac_record = rmn.FST_RDE_META_DEFAULT.copy()

            try:
                rmn.writeGrid(file_id, MACC_grid)

                tac = rmn.fstinl(file_id, nomvar='>>')[0]
                tic = rmn.fstinl(file_id, nomvar='^^')[0]

                tic_record.update(rmn.fstprm(tic))
                tac_record.update(rmn.fstprm(tac))

                tic_record.update({'datyp' : rmn.FST_DATYP_LIST['float']})
                tac_record.update({'datyp' : rmn.FST_DATYP_LIST['float']})

                rmn.fsteff(tic)
                rmn.fsteff(tac)

                tic_record.update({'d': MACC_grid['ay']})
                tac_record.update({'d': MACC_grid['ax']})
                toc_record = vgd.vgd_new_pres(const_pressure, ip1=MACC_grid['ig1'], ip2=MACC_grid['ig2'])

                rmn.fstecr(file_id, tic_record)  # write the dictionary record to the file as a new record
                rmn.fstecr(file_id, tac_record)  # write the dictionary record to the file as a new record
                vgd.vgd_write(toc_record, file_id)

            except:
                raise

            for rp1 in xrange(len(const_pressure)):  # writes a record for every level (as a different ip1)
                try:
                    # converts rp1 into a ip1 with pressure kind
                    ip1 = rmn.convertIp(rmn.CONVIP_ENCODE, const_pressure[rp1], rmn.KIND_PRESSURE)
                    new_record.update(MACC_grid)
                    new_record.update({  # Update with specific meta
                        'nomvar': nomvar,
                        'typvar': 'C',
                        'etiket': 'MACCRean',
                        'ni'    : MACC_grid['ni'],
                        'nj'    : MACC_grid['nj'],
                        'ig1'   : tic_record['ip1'],
                        'ig2'   : tic_record['ip2'],
                        'ig3'   : tic_record['ip3'],
                        'ig4'   : tic_record['ig4'],
                        'dateo' : rmn.newdate(rmn.NEWDATE_PRINT2STAMP, 20120101, 0000000),
                        'deet'  : 0,  # Timestep in sec
                        'ip1'   : ip1
                        })

                    #tmp_nparray = np.asfortranarray(monthly_mean[rp1])
                    tmp = scaledSSData[rp1]
                    tmp = np.transpose(tmp)
                    # data array is structured as tmp = monthly_mean[level] where monthly_mean is [level, lat, lon]
                    new_record.update({'d': tmp.astype(np.float32)}) # Updates with data array in the form (lon x lat)

                    print "Defined a new record with dimensions ({0}, {1})".format(new_record['ni'], new_record['nj'])
                    rmn.fstecr(file_id, new_record)  # write the dictionary record to the file as a new record

                except:
                    #rmn.closeall(file_id)
                    rmn.fstfrm(file_id)
                    rmn.fclos(file_id)
                    raise
            rmn.fstfrm(file_id)
            rmn.fclos(file_id)
            print('{} complete~'.format(filename))
        except:
            rmn.fstfrm(file_id)
            rmn.fclos(file_id)
            raise
    print('Finished plotting all FSTs. ')
