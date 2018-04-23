## ConcentrationPlotter returns a map of the difference between two concentration plots for the same area, latitude and longitude
## Can be run standalone (use main to pass arguments if so desired) or use a python script to pass through values to generate multiple plots

import numpy as np
import os
import rpnpy.librmn.all as rmn
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from math import cos, radians
import errno
import logging

## constants
currentDir = os.path.dirname(os.path.abspath(__file__))
defaultPathStnd = os.path.join(currentDir,'testFiles/Global_401x200/')
defaultPathAlt = os.path.join(currentDir,'testFiles/adomKPPKPPB_adomYB/20100710000000/model/')
defaultSave = currentDir
# defaultFile is the name of the FST files to be compared
defaultFile = '2010011400_024'
defaultIp1=76696048
# defaultSpc not in list format since this function only takes one species
defaultSpc='TO3'
defaultDiffCmap = 'RdBu'
defaultConcCmap = 'spectral'
defaultRev = True
defaultExt = 'neither'
#defaultProjection = 'stere'
defaultProjection = 'robin'
defaultBins = None
# defaultClear removes a set of x/y coordinates from the data
defaultClear = 0
# buff indicates the degress of extra buffer space to be added to the map dimensions. May be needed if the map does not display all of the data
buff = 2
R=6400E3 #radius of earth
# cmaps is a list of all possible default colormaps. Currently the list is incomplete, though maps not in this list may still be plotted
cmaps = ['RdBu', 'Greys','cubehelix','jet','spectral']
# setting constraints on total data plotted. Note that to use lon lat, need to pass through getFuzzyNiNj
defaultLL = [0,0]
defaultUR = [400,100]

def getFuzzyNiNj(lonVal, latVal, gridData):
    xypos = rmn.gdxyfll(gridData, latVal, lonVal)
    ni = int(round(xypos['x'])) - 1 #-1 to use 0 as list index
    nj = int(round(xypos['y'])) - 1
    print ('Closest lonlat values: ({}, {})'.format(lonData[ni][nj], latData[ni][nj]))
    return ni, nj

def getConc(filePath, level, spc, niRange, njRange):
    """ Get the concentration data for an FST file.

    Returns the data in an array, the datakey as a key and the fileID of the FST file.
    """

    try:
       fileID = rmn.fstopenall(filePath,rmn.FST_RO)
       dataKey = rmn.fstinf(fileID,nomvar=spc,ip1=level)['key']
       dataRec = rmn.fstluk(dataKey)
       tempConc = dataRec['d']
       print(len(tempConc),len(tempConc[0]))
       print(tempConc)
       concData = dataRec['d'][niRange[0]:niRange[1]+1, njRange[0]:njRange[1]+1]
       print ('File {} recorded'.format(filePath))
       return {'concData':concData, 'dataKey':dataKey, 'fileID':fileID}
    except TypeError:
       print('Unable to record file {}. Please see log for details'.format(filePath))
       # log an error into the log file
       logging.warning('nomvar {} and ip1 {} could not be found for file {}.'.format(fileID, spc, level, filePath))
       pass

def concDiff(concDataStnd, concDataAlt):
    """ Get the differences between two sets of concentration data.

    Returns the difference between two concentration data sets in an array.
    """

    # Get the dimensions of the difference data, assume both concentration data sets are the same size
    dim1 = concDataStnd.shape[0]
    dim2 = concDataStnd.shape[1]
    # Initialize an array of zeros.
    diffData = np.zeros((dim1,dim2))
    # Replace the zeros with differences between the concentrations.
    # Note that the difference should be Model-Base
    for i in range(dim1):
        for j in range(dim2):
            diffVal = concDataAlt[i][j]-concDataStnd[i][j]
            diffData[i][j] = diffVal
    return diffData

def convertLons(lonArray):
    rowInd = 0
    for row in lonArray:
        lonInd = 0
        for lon in row:
            if lon > 180:
                lonArray[rowInd][lonInd] = lon-360
            lonInd += 1
        rowInd += 1
    return lonArray

def getGrid(dataKey, fileID, niRange, njRange):
    """ Get the grid details in form of lon/lat of the FST file. """

    # Get the file metadata, then add on the keypair for 'iunit'
    fileMeta = rmn.fstprm(dataKey)
    fileMeta['iunit'] = fileID
    # Get the grid data and decode it
    gridData = rmn.ezqkdef(fileMeta)
    gridDecode = rmn.decodeGrid(gridData)
    llGridData = rmn.gdll(gridData)
    latData = llGridData['lat'][niRange[0]:niRange[1]+1, njRange[0]:njRange[1]+1]
    lonData = convertLons(llGridData['lon'])[niRange[0]:niRange[1]+1, njRange[0]:njRange[1]+1]
    return {'gridll':llGridData,'gridID':gridData, 'latData':latData,'lonData':lonData}

def closeFST(fileID):
    """ Closes the FST file once relevant data has been saved. """

    rmn.fstcloseall(fileID)
    print ('File has been closed.')

def gridCheck(gridStnd, gridAlt):
    """ Does a sanity check in case the grids are different for the FST files. """

    if (np.array_equiv(gridStnd['lat'],gridAlt['lat'])) and (np.array_equiv(gridStnd['lon'],gridAlt['lon'])):
        True
    else:
        print('The area of the files do not match each other. Program may not function as expected. See log for more details.')
        logging.warning('Grids for current run do not match up. The differences are : \n' + str(np.setdiff1d(gridStnd,gridAlt)))

def maxDimensions(lowLon,highLon,lowLat,highLat,buff):
    """ Get the maximum dimensions of the map.

    Note that this function is still underdevelopment and may not represent the ideal dimensions for the resulting map.
    """

    lonDiff = radians(highLon - lowLon)
    latDiff = radians(highLat - lowLat)
    radBuffer = 2*(radians(buff))
    yDist = R*(latDiff + radBuffer)
    Rad1 = R*cos(lowLat)
    Rad2 = R*cos(highLat)
    xDist1 = Rad1*(lonDiff + radBuffer)
    xDist2 = Rad2*(lonDiff + radBuffer)
    xDist = (xDist1+xDist2)/2
    return {'xDist':xDist, 'yDist':yDist}

def midLatLon(lonData, latData):
    """ Provides an alternate manner in which to calculate the central lon/lat values. """

    ni = len(lonData)
    nj = len(lonData[0])
    #find the middle value between these
    if ni%2==0:
        if nj%2==0:
            midLon=(lonData[ni//2-1,nj//2]+lonData[ni//2,nj//2])/2
            midLat=(latData[ni//2,nj//2-1]+latData[ni//2,nj//2])/2
        else:
            midLon=(lonData[ni//2-1,nj//2]+lonData[ni//2,nj//2])/2
            midlat=(latData[ni//2,nj//2-1]+latData[ni//2,nj//2])/2
    else:
        if nj%2==0:
            midLon=(lonData[ni//2,nj//2]+lonData[ni//2,nj//2+1])/2
            midLat=(latData[ni//2,nj//2]+latData[ni//2,nj//2+1])/2
        else:
            midLon=lonData[ni//2,nj//2]
            midLat=latData[ni//2,nj//2]
    return {'midLat':midLat, 'midLon':midLon}

def reverseName(mapType, reverse):
    """ Returns the reverse mapType if reverse=True, otherwise returns mapType. """

    if reverse == True:
        mapType = mapType + '_r'
        return mapType
    else:
        return mapType

def isCmap(mapType):
    """ Checks to see if the mapType is in the list of default colormaps. """

    if mapType in cmaps:
        return True
    else:
        print('Warning: Type mapType is not in the list of cmaps, may resort to default cmap.')
         # TODO: If cmap is invalid, set cmap = default cmap
         # note that this behaviour might be default! May set to cmap 'jet' though

def cmapType(cmap,reverse):
    """ Returns map details. """

    isCmap(cmap)
    mapType = reverseName(cmap,reverse)
    return mapType

def makeDir(path):
    """ Creates a directory if it doesn't exist.

    Existence errors will be ignored. All other errors will be raised.
    """

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def plotConc(concData, lonData, latData, saveLoc, modelRun, prtype, ip1, spc, cmapType, bins, minVal, maxVal, extension, name, removed,buff):
    """ Plots and saves concentration data.

    Will plot difference data if difference data was passed.
    """

    # useful things for determining projection details
    lowLon = np.amin(lonData)
    lowLat = np.amin(latData)
    highLon = np.amax(lonData)
    highLat = np.amax(latData)
    maxDim = maxDimensions(lowLon,highLon,lowLat,highLat,buff)
    midLon = (highLon+lowLon)/2
    midLat = (highLat + lowLat)/2
    # Uncomment the following lines for an alternative midLatLon calculation
    #mids = midLatLon(lonData,latData)
    #midLon = mids['midLon']
    #midLat = mids['midLat']
    max_width = abs(maxDim['xDist'])
    max_height = abs(maxDim['yDist'])
    # Initialize the figure
    fig = plt.figure(figsize=(8,8))

    # Create the map based on projection type
    print('Maybe I only worked up until I tried to make a map...')
    if (prtype=='ortho') or (prtype == 'nsper') or (prtype == 'laea') or (prtype == 'aeqd') or (prtype == 'gnom') or (prtype == 'lcc'):
        concMap = Basemap(projection=prtype, resolution = 'c', lon_0=midLon, lat_0=midLat, width=max_width, height=max_height)
    elif prtype == 'stere' or prtype == 'robin':
        concMap = Basemap(projection=prtype, lon_0=midLon,lat_0=midLat,width=max_width,height=max_height)
    elif (prtype == 'cyl') or (prtype == 'merc'):
        concMap = Basemap(projection=prtype, resolution='c', llcrnrlat=lowLat,  urcrnrlat=highLat, llcrnrlon=lowLon, urcrnrlon=highLon)
    elif (prtype == 'aea') or (prtype == 'eqdc'):
        concMap = Basemap(projection = prtype, lon_0=midLon,  lat_0=midLat, llcrnrlat=lowLat, urcrnrlat=highLat, llcrnrlon=lowLon, urcrnrlon=highLon)
    else:
        print('Error: Could not generate map. Try a different projection.')
        # Can check the available basemap types and add to existing if statements

    # Add in other map details
    mapColor = cmapType
    mapBins = bins
    x, y = concMap(lonData, latData)
    print('I worked up until here~.')
    concMap.pcolormesh(x, y, concData, cmap=plt.cm.get_cmap(mapColor, mapBins))
    concMap.drawcoastlines(color='lightgray')
    concMap.drawcountries(color='gray')
    concMap.drawstates(color='gray')
    # Comment out the following to remove lonlat lines
    concMap.drawmeridians(np.arange(0,360,10), labels=[0,0,0,1], fontsize=6)
    concMap.drawparallels(np.arange(-180,180,10), labels=[1,0,0,0],fontsize=6)

    # Add colorbar and details
    #TODO: Fix the set label option for the colorbar, right now the colorbar doesn't have a title
    cbar = plt.colorbar(extend = extension, shrink=0.5)
    cbar.set_label=('Concentration: '+spc)
    plt.clim(minVal, maxVal)

    # Name and save the figure
    hy = ((os.popen('r.ip1 {}'.format(ip1))).read()).lstrip()
    plt.title('{}, {} \n  hy: {}, Spc: {}'.format(name,modelRun, hy,spc))
    fig.savefig(os.path.join(saveLoc, modelRun) + '_' + spc + '_' + name + '.png', dpi=300, bbox_inches='tight', pad_inches=0.3)

    plt.close('all')

def diffPlot(fileRun=defaultFile, baseFile=defaultPathStnd+defaultFile, modelFile=defaultPathStnd+defaultFile, savePath=defaultSave, level=defaultIp1, species=defaultSpc, buffering = buff, mapType=defaultDiffCmap, reverse=False, extension='neither', projType='stere', vmin=None, vmax=None, totalBins=defaultBins,partName='test',removed=defaultClear, niRange=[defaultLL[0], defaultUR[0]], njRange=[defaultLL[1], defaultUR[1]]):
    """ Execute difference plot creation. """

    logging.basicConfig(filename='ConcentrationPlotter.log', level=logging.DEBUG,format='%(asctime)s %(message)s')

    # print the min number of messages produced by librmn
    rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)
    # get the data for the standard first, close file to continue
    concDataStnd = getConc(baseFile,level,species, niRange, njRange)
    if concDataStnd == None:
       return None
    else:
       gridStnd = getGrid(concDataStnd['dataKey'],concDataStnd['fileID'], niRange, njRange)
       closeFST(concDataStnd['fileID'])
       # get the data for the alternative, close file to continue
       concDataAlt = getConc(modelFile,level,species, niRange, njRange)
       if concDataAlt == None:
          return None
       else:
          gridAlt = getGrid(concDataAlt['dataKey'],concDataAlt['fileID'], niRange, njRange)
          closeFST(concDataAlt['fileID'])
          # get the difference between the datasets
          diffData = concDiff(concDataStnd['concData'],concDataAlt['concData'], niRange, njRange)
          # check that the grid sizes make sense
          gridCheck(gridStnd['gridll'],gridAlt['gridll'])
          # minor cleanup
          cmapDetails = cmapType(mapType,reverse)
          makeDir(savePath)
          # plot and save the figure
          plotConc(diffData, gridStnd['lonData'], gridStnd['latData'], savePath, fileRun, projType, level, species,cmapDetails, totalBins, vmin, vmax, extension, partName,removed,buffering)

def concPlot(fileRun=defaultFile, modelFile=defaultPathStnd+defaultFile, savePath=defaultSave, level=defaultIp1, species=defaultSpc, buffering=buff, mapType=defaultConcCmap, reverse=False, extension='neither', projType='stere', vmin=None, vmax=None, totalBins=defaultBins, partName='test', removed=defaultClear, niRange = [defaultLL[0],defaultUR[0]], njRange = [defaultLL[1], defaultUR[1]]):
    """ Execute the concentration plot creation. """

    logging.basicConfig(filename='ConcentrationPlotter.log', level=logging.DEBUG,format='%(asctime)s %(message)s')

    # print the min number of messages produced by librmn
    rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)
    # get concentration and grid data, then close file
    concData = getConc(modelFile,level,species, niRange, njRange)
    if concData == None:
       return None
    else:
       gridData = getGrid(concData['dataKey'],concData['fileID'], niRange, njRange)
       closeFST(concData['fileID'])
       # minor cleanup
       cmapDetails =cmapType(mapType, reverse)
       makeDir(savePath)
       # plot and save the figure
       plotConc(concData['concData'], gridData['lonData'], gridData['latData'], savePath, fileRun, projType, level, species, cmapDetails, totalBins, vmin, vmax, extension, partName, removed, buffering)

# To run this program alone, uncomment one of the following and include desired parameters
# Otherwise the defaults will run (and I can't guarantee it'll work if the files have been moved!)

# diffPlot()
concPlot()
