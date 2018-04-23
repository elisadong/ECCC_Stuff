# Monthly Mean of Chemical Species
import pygeode as pyg
import numpy as np
import rpnpy.librmn.all as rmn
import time, glob
import matplotlib.pyplot as plt
import rpnpy.vgd.all as vgd
from mpl_toolkits.basemap import Basemap
from levels import *

##### FXNS #####

def get_var(open_file, nomvar, m_int):
    repr_var = "<Var '{0}'>".format(nomvar)
    print("Getting " + repr_var)
    
    var_data = splice_month([data for data in open_file.vars if repr(data) == repr_var][0], m_int)

    return var_data


def get_num_params(nomvar):
    '''
    returns a tuple of time, levels, lat, lon values.
    '''
    
    print "Now obtaining dimensions: {0}".format(str(nomvar.shape))
    return (len(nomvar), len(nomvar[0]), len(nomvar[0,0]), len(nomvar[0,0,0])) 


def vert_interp(pressure, org):
    '''
    org should be a 4d array, pressure should be a 4d array
    of pressures.
    '''
    global const_pressure    
    start_time = time.time()
    print "Conducting vertical interpolation..."

    lon         = len(org[0, 0, 0])
    lat         = len(org[0, 0])
    lev         = len(org[0])
    timesteps   = len(org)
    y_interp    = np.zeros((timesteps, len(const_pressure), lat, lon))
    x_interp    = const_pressure
    
    for timestep in xrange(timesteps):  # the array for np.interp must be 1d - this is slow
        for j in xrange(lon):
            for i in xrange(lat):
                try:
                    x_initial = pressure[timestep, :, i, j]
                    y_initial = org[timestep, :, i, j]
                    y_interp[timestep, :, i, j] = np.interp(x_interp, x_initial, y_initial)
                except:
                    print pressure.shape
                    #debug_tuple = (len(y_initial), len(x_initial), len(x_interp), len(y_interp[timestep, :, i, j]), lev)
                    #print "LEV: {4}, Y: {0}, X: {1}, XI: {2}, YI: {3}".format(debug_tuple)
                    raise
                
    print "That took " + str(time.time() - start_time) + " seconds"
    return y_interp


def splice_month(open_var, m_int):
    '''
    takes an open var, splices out the timesteps that dont regard
    to that month (with m_int)
    '''
    print "Splicing..."
    if len(open_var.time) == 1464:
        leap_year = True
    elif len(open_var.time) == 1460:
        leap_year = False
    else:
        print("File must contain all 4 time values for each day of the year")
        exit()
    
    if leap_year:
        stored_array = [
                        open_var[:124], open_var[124:240], open_var[240:364], open_var[364:484],
                        open_var[484:608], open_var[608:728], open_var[728:852], open_var[852:976],
                        open_var[976:1096], open_var[1096:1220], open_var[1220:1340], open_var[1340:]
                        ]
    else:
        stored_array = [
                        open_var[:124], open_var[124:236], open_var[236:360], open_var[360:480],
                        open_var[480:604], open_var[604:724], open_var[724:848], open_var[848:972],
                        open_var[972:1092], open_var[1092:1216], open_var[1216:1336], open_var[1336:]
                        ]
    return stored_array[m_int]


def get_pressures(open_nc, m_int):
    '''
    (open .nc, int) -> np.ndarray
    gets pressures from open_file, based on m_int as an integer
    representing the month you wish to obtain.
    '''
    # TODO: check the units
    global Ak, Bk, const_pressure
                                 
    start_time = time.time()
    print "Getting pressures..."
    
    lnsp = open_nc.lnsp  # open file containing every lnsp for that year
    lnsp_array = splice_month(lnsp, m_int)
    print lnsp_array.shape

    shape = list(lnsp_array.shape); shape.insert(1, 60)  # creates a dimension for levels
    print "PRESSURE SHAPE IS {0}".format(str(shape))
    pressure = np.zeros(shape)
    try:
        for lev in xrange(len(pressure[0])):
            pressure[:, lev, :, :] = ((Ak[lev] + Bk[lev]*np.exp(lnsp_array)) / 100.)
    
    except:
        debug_tuple = (lev)
        print "LEV: {1}, TS: {0}".format(timestep, lev)
        raise
    
#    print pressure[5, :, 12, 33]
#    print pressure[15, :, 12, 23]
#    print pressure[45, :, 12, 53]
#    print pressure[75, :, 2, 39]
    
    print "That took " + str(time.time() - start_time) + " seconds"
    return pressure 


def shift_lon(array):
    '''
    shifts the longitude from -180 - 180 to 0 - 360
    '''
    global timelen, latlen, lonlen
    new_array = np.zeros(array.shape)
    for k in xrange(60):
        # all this does is move the indexes corresponding from -180 to 180 to indexes corresponding to 0 - 360
        tmp = np.zeros((timelen, latlen, lonlen))
        tmp[:, :, :lonlen/2], tmp[:, :, lonlen/2:] = array[:, k, :, lonlen/2:], array[:, k, :, :lonlen/2]
        new_array[:,k] = tmp

    return new_array


def GEMMACCify(array):  # not used.. 
    '''
    np.array -> np.array
    takes an array, horizontally interpolates it to GM specs
    '''
    global GM_grid, MACC_grid, latlen, lonlen

    for k in xrange(60):
        # all this does is move the indexes corresponding from -180 to 180 to indexes corresponding to 0 - 360
        tmp = np.zeros((latlen, lonlen))
        tmp[:, 0:lonlen/2 - 1], tmp[:, lonlen/2 - 1:lonlen - 1] = array[k, :, lonlen/2 - 1:lonlen - 1], array[k, :, 0:lonlen/2 - 1]
        array[k] = tmp 

    cubed_array = np.zeros((60, GM_grid['ni'], GM_grid['nj']))
    array = array.transpose(0,2,1)  # flips lon and lat axes
    for lev in xrange(len(cubed_array)):
        cubed_array[lev] = rmn.ezsint(GM_grid, MACC_grid, array[lev])
    print cubed_array.shape
    return cubed_array


def merge_nc(nc1, nc2, nc3, nomvar, m_int):
    '''
    (open .nc, open .nc) -> np.array
    merges northern hemisphere (part 1+2) file data with tropics 
    and southern hemisphere (part 3) data. 
    '''
    global monthly_mean

    start_time = time.time()
    print "Now performing .nc merge..."

    hcho1 = get_var(nc1, nomvar).get()
    hcho2 = splice_month(get_var(nc2, nomvar), m_int)
    hcho3 = get_var(nc3, nomvar).get()

    param_lens = list(get_num_params(hcho3))  # makes a list of params (time, lev, lat, lon)
    
    # shrinks lat by 2+1 (redundancy), adds in hcho1 and hcho2
    param_lens[2] = param_lens[2] + get_num_params(hcho1)[2] - 2 + get_num_params(hcho2)[2] - 1 
    hcho = np.zeros(param_lens)
    
    hcho3.setflags(write=1)  # sets the np array to be writeable. it's read only, must be writeable to splice out the repeat data

    print "the shape is " + str(param_lens)
    print len(hcho[3])
    
    temp_shape = list(hcho2.shape)
    temp_shape[2] = temp_shape[2] - 1
    fixed_hcho2 = np.zeros(temp_shape)
    
    print "DEBUG: fixed_hcho2 shape = " + str(fixed_hcho2.shape)
    
    temp_shape = list(hcho3.shape)
    temp_shape[2] = temp_shape[2] - 2
    
    # smooths out the two latitudes
    param_lens[2] = 2
    smoothing = np.zeros((param_lens))
    smoothing[:,:,0,:] = (hcho1[:,:,0,:] + hcho3[:,:,-2,:])/ 2.
    smoothing[:,:,1,:] = (hcho1[:,:,1,:] + hcho3[:,:,-1,:])/ 2.

#    for t in xrange(smoothing.shape[0]):
#        for k in xrange(smoothing.shape[1]):
#            for j in xrange(smoothing.shape[2]):
#                for i in xrange(smoothing.shape[3]):
#                    smoothing[t,k,j,i] = np.mean([hcho1[t,k,-j-1,i], hcho3[t,k,j,i]])
    
    hcho[:,:,0:108,:]=hcho3#[:,:,2:108,:]
    #hcho[:,:,106:108,:] = smoothing
    hcho[:,:,108:159,:]=hcho1[:,:,:-2,:]
    hcho[:,:,159:161,:]=hcho2[:,:,0:2,:]

    check1, check2 = True, True #hcho2[2,2,0,4] == hcho[2,2,-1,4], hcho1[2,2,0,4] == hcho[2,2,-3,4]
    if check1 and check2:
        print "MERGE CONFIRMED SUCCESSFUL"
        print "That took " + str(time.time() - start_time) + " seconds"
        return hcho
    else:
        print check1, check2
        print "VALUES DO NOT MATCH"
        exit()


#def merge_nc(nc1, nc2, nc3, nomvar, m_int):
#    '''
#    (open .nc, open .nc) -> np.array
#    merges northern hemisphere (part 1+2) file data with tropics 
#    and southern hemisphere (part 3) data. 
#    '''
#
#    start_time = time.time()
#    print "Now performing .nc merge..."
#
#    #hcho1 = nc1.go3.get()
#    #hcho3 = nc3.go3.get()
#
#    hcho1 = get_var(nc1, nomvar).get()
#    hcho2 = splice_month(nc2.go3, m_int)
#    hcho3 = get_var(nc3, nomvar).get()
#
#    param_lens = list(get_num_params(hcho3))  # makes a list of params (time, lev, lat, lon)
#    # shrinks lat by 2+1 (redundancy), adds in hcho1 and hcho2
#    param_lens[2] = param_lens[2] + get_num_params(hcho1)[2] - 2 + get_num_params(hcho2)[2] - 1 
#    hcho = np.zeros((param_lens[0], param_lens[1], param_lens[2], param_lens[3]))
#    
#    hcho3.setflags(write=1)  # sets the np array to be writeable. it's read only, must be writeable to splice out the repeat data
#
#    print "the shape is " + str(param_lens)
#    print len(hcho[3])
#    
#    temp_shape = list(hcho2.shape)
#    temp_shape[2] = temp_shape[2] - 1
#    fixed_hcho2 = np.zeros(temp_shape)
#    print "DEBUG: fixed_hcho2 shape = " + str(fixed_hcho2.shape)
#    temp_shape = list(hcho3.shape)
#    temp_shape[2] = temp_shape[2] - 2
#    fixed_hcho3 = np.zeros(temp_shape)
#    print "DEBUG: fixed_hcho3 shape = " + str(fixed_hcho3.shape)
#    print "DEBUG: Broadcasting {0} into {1}".format(str(hcho2[:,:,:2].shape), fixed_hcho2[:,:].shape)
#    fixed_hcho2[:,:] = hcho2[:,:,1::-1]  # splices the first two lat indexes INTO fixed_hcho2
#    print "DEBUG: Broadcasting {0} into {1}".format(str(hcho3[:,:,2:].shape), fixed_hcho3[:,:].shape)
#    fixed_hcho3[:,:] = hcho3[:,:,-3::-1]  # splices first two indexes of hcho3 out (shared lats)
#
#    hcho = np.concatenate((fixed_hcho3, hcho1[:,:,::-1], fixed_hcho2), axis=2)
#    #hcho[1] = np.mean(hcho[0], axis=0)
#    print "the shapes are ORG: {0}".format(str(hcho[0].shape))
#   
#    check1, check2 = hcho2[2,2,0,4] == hcho[2,2,-1,4], hcho1[2,2,0,4] == hcho[2,2,-3,4]
#    if check1 and check2:
#        print "MERGE CONFIRMED SUCCESSFUL"
#        print "That took " + str(time.time() - start_time) + " seconds"
#        return hcho
#    else:
#        print check1, check2
#        print "VALUES DO NOT MATCH"
#        exit()


##### MAIN CODE #####

start_time = time.time()

lnsp_files = ["/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/LNSP/2008.nc", 
              "/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/LNSP/2009.nc", 
              "/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/LNSP/2010.nc", 
              "/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/LNSP/2011.nc", 
              "/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/LNSP/2012.nc"]  # in this list, place the directory of the year
#lnsp_files = [lnsp_files[1]]  # this is temporary, just to work with one file
file_dir = ["/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/2008", 
            "/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/2009", 
            "/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/2010", 
            "/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/2011", 
            "/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/2012"]  # in this list, place the directory of the year
#file_dir = file_dir[:4]  # this is temporary, just to work with files i have
#file_dir = [file_dir[1]]  # this is temporary, just to work with one file
arc_dir = "/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/PART4/*"  # directory for the arctic files
file_dir = ["/space/hall1/sitestore/eccc/aq/r1/alh002/NCDF/SPECIES/CO"]
#lnsp_files = ["/home/ords/aq/alh002/NCDF/LNSP/2008.nc", 
#              "/home/ords/aq/alh002/NCDF/LNSP/2009.nc", 
#              "/home/ords/aq/alh002/NCDF/LNSP/2010.nc", 
#              "/home/ords/aq/alh002/NCDF/LNSP/2011.nc", 
#              "/home/ords/aq/alh002/NCDF/LNSP/2012.nc"]  # in this list, place the directory of the year
##lnsp_files = [lnsp_files[1]]  # this is temporary, just to work with one file
#file_dir = ["/home/ords/aq/alh002/NCDF/2008", 
#            "/home/ords/aq/alh002/NCDF/2009", 
#            "/home/ords/aq/alh002/NCDF/2010", 
#            "/home/ords/aq/alh002/NCDF/2011", 
#            "/home/ords/aq/alh002/NCDF/2012"]  # in this list, place the directory of the year
##file_dir = file_dir[:4]  # this is temporary, just to work with files i have
##file_dir = [file_dir[1]]  # this is temporary, just to work with one file
#arc_dir = "/home/ords/aq/alh002/NCDF/PART4/*"  # directory for the arctic files

params = ['go3', 'hcho', 'ch4', 'co', 'nox']
params = [params[3]]

month_list = ['/01JAN/','/02FEB/', '/03MAR/',  
              '/04APR/', '/05MAY/', '/06JUN/', 
              '/07JLY/', '/08AUG/', '/09SEP/', 
              '/10OCT/', '/11NOV/', '/12DEC/']

try:  # opens the file
#    fileid = rmn.fstopenall(output_files[0], rmn.FST_RW)        
    # makes an empty .fst file
    for m_int, month in enumerate(month_list):
        temp = month.strip('/')
    #    day = time.gmtime()
    #    temp = '' 
    #    for x in day: temp += str(x)
        filename = '/home/jac001/sitestore_ppp1+hare/r2/model_output/GM-g_HTAP2010_analysis/MACC_analysis/FST_STORAGE/output_file_{0}.fst'.format(temp) 
        #filename = '/home/ords/aq/alh002/pyscripts/workdir/output_file_{0}.fst'.format(temp) 
        #filename = '/home/ords/aq/alh002/pyscripts/workdir/output_file_{0}.fst'.format(temp) 
        tmp = open(filename, 'w+'); tmp.close()

        output_file = filename
        output_files = [output_file]

        file_id = rmn.fnom(output_file)
        open_fst = rmn.fstouv(file_id, rmn.FST_RW)
        print(file_id, open_fst)

        for nomvar in params:
            #monthly_mean = np.zeros((len(file_dir), 60, 161, 320))  # use this when not interpolating
            monthly_mean = np.zeros((len(file_dir), len(const_pressure), 161, 320))  # stores the np arrays for the mean of each year
            # start at first year.
            # its just to crawl through the monthly mean array and input the mean for that year
            year = 0
            
            # goes thru each param
            for i, open_file in enumerate(file_dir):
                print "Now processing:", open_file[-4:] + month
                print "Parameter: " + nomvar

                # this whole mess is just me opening the files
                # without causing errors
                param = get_var(pyg.openall(file_dir[0]+"/*"), nomvar, m_int)
#                files = glob.glob(open_file + "{0}*.nc".format(month))
#                p2file = glob.glob(arc_dir)
#                p2file.sort()
#                p2file = p2file[i]  # change to i instead of 1 when you quit debugging
#                p3file = glob.glob(open_file + "{0}*_PART3.nc".format(month))
#                files.remove(p3file[0])
#                nc = pyg.openall(files)
#                nc2 = pyg.open(p2file)
                nc_lnsp = pyg.open(lnsp_files[i])
                pressures = get_pressures(nc_lnsp, 0)
#                nc3 = pyg.openall(p3file)
            
                ''' for each nc variable, depending on the month, they are
                4D np arrays of (time, level, lat, lon). time is 4*days,
                level is from 1-60, lat and lon are usually 53+2+106,320 '''
                
                #param = merge_nc(nc, nc2, nc3, nomvar, m_int)  # change 0 to a proper m_int
                print param.shape
                timelen, levlen, latlen, lonlen = param.shape  # get the number of datapoints of each
                param = shift_lon(param)
                param = vert_interp(pressures, param)  # temp
                print param.shape
                print "Obtaining monthly mean..."
                mean_param = np.mean(param, axis=0)  # temp
                #mean_param = param[50]
                print mean_param.shape
                monthly_mean[year] = mean_param
                year += 1

            ''' the first set of .nc files (part 1 and 2) are limited to northern
            hemisphere, whereas the rest are tropics and southern hemisphere.
            they have colliding vars (species), thus they must be processed individually and combined
            '''

            print "{0} months of data were analyzed".format(monthly_mean.shape[0])
            monthly_mean = np.mean(monthly_mean, axis=0)
    #        tmp = list(monthly_mean.shape)
    #        tmp[2] = tmp[2]+1
    #        tmp = np.zeros(tmp)
    #        tmp[:,:,:len(monthly_mean[2,2])] = monthly_mean
    #        monthly_mean = tmp
    #        monthly_mean[:,:,320] = monthly_mean[:,:,319]
            print "The shape of the data is " + str(monthly_mean.shape)

            # everything above is the "slow step" of the "reaction"
            print "That took {0} seconds!".format(str((time.time() - start_time)))

            # because all the grids within a month have the same values, i just pulled the most recent
            #timelen, levlen, latlen, lonlen = get_num_params(hcho)  # get the number of datapoints of each

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
                #tic_record.update(MACC_grid)
                #tac_record.update(MACC_grid)
                
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
                #rmn.closeall(file_id)
    #            rmn.fstfrm(file_id)
    #            rmn.fclos(file_id) 
                raise
                exit()

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
                    tmp = monthly_mean[rp1]
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
                    exit()

        rmn.fstfrm(file_id)
        rmn.fclos(file_id) 
        print "File Closed"
        print "That took {0} seconds!".format(str((time.time() - start_time)))
except:
    rmn.fstfrm(file_id)
    rmn.fclos(file_id)
    raise
    exit()

