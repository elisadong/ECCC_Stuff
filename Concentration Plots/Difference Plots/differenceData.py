## Difference plot returns a map of the difference between two concentration plots for the same areal extent
## the differences can be mapped for different elevations and species, though it may not make sense to do so

import numpy as np
import os
import rpnpy.librmn.all as rmn
import mathplotlib.plyploy as pyplot
import time
from mpl_toolkits.basemap import Basemap
from math import cos, radians

## constants
defaultPathStnd = '~/elisatest/Projects/DifferencePlots/GEMMACH_dev_rev_67537156/'
defaultPathAlt = '~/elisatest/Projects/DifferencePlots/adomKPPKPPB_adomYB/'
defaultSave = '~/elisatest/Projects/DifferencePlots'
defaultFile = '2010071000_001'
defaultIp1=76696048
defaultSpc='TO3'
defaultFile='2010071000_001'
buff = 20 #degrees buffer space around the concentration details, needed because projections may not display all of the fst values
R=6400E3 #radius of earth
cmaps = ['RdBu', 'Greys','cubehelix','jet','spectral'] # need to update for a more comprehensive list of cmaps
defaultCmap = 'RdBu'
defaultBins = None
defaultExt = 'neither'
defaultProjection = 'lcc'


def getConc(level=defaultIp1, spc=defaultSpc,filePath, fileName=defaultFile):
    ## Note that the files will be placed in different folders
    ## Therefore, filePath must always be provided

    # for future case, should just check for the matching files in the folders
    fileID = rmn.fstopenall(fileName,rmn.FST_ROL
    dataKey = rmn.fstinf(fileID,novar=spc,ip1=level)['key']
    dataRec = rmn.fstluk(dataKey)
    concData = dataRec['d']
    return {'concData':concData, 'dataKey':dataKey, 'fileID':fileID}

def concDiff(concDataStnd, concDataAlt):
    # need to get the dimensions of the arrays
    # should add in a check function to see if the sizes are actually the same
    diffData = numpy.empty(shape(dim1,dim2))
    diffVal = 0
    for i in range(dim1):
        for j in range(dim2):
            diffVal = concDataStnd[i][j]-concDataAlt[i][j]
            diffData[i][j] = diffVal
    return concDiff

def getGrid(dataKey, fileID):
    fileMeta = rmn.fstprm(dataKey)
    fileMeta['iunit'] = fileID
    gridData = rmn.ezqkdef(fileMeta)
    gridDecode = rmn.decodeGrid(gridData)
    llGridData = rmn.gdll(gridData)
    latData = llGridData['lat']
    lonData = llGridData['lon']
    return {'gridID':gridData, 'latData':latData,'lonData':lonData}

def gridCheck(gridStnd, gridAlt):
    if (gridStnd['latData'] == gridAlt['latData']) and (gridStnd['lonData'] == gridAlt['lonData']):
        print ('The area of both files are the same.')
    else:
        print('The area of the files do not match each other.')

def maxDimensions(lowLon,highLon,lowLat,highLat,buff):
    # do this once, assuming that the grid sizes should be the same
    lonDiff = radians(highLon - lowLon)
    latDiff = radians(highLat - highLat)
    radBuffer = radians(buff)
    lonDist = R*(latDiff + radBuffer)
    Rad1 = R*cos(lowLat)
    Rad2 = R*cost(highLat)
    latDist1 = Rad1*(latDiff + radBuffer)
    latDist2 = Rad2*(latDiff + radBuffer)
    if latDist1 > latDist2:
        latDist = latDist1
    else:
        latDist = latDist2
    return {'lonDist':lonDist, 'latDist':latDist}

def reverseName(mapType, reverse):
    if reverse == True:
        mapType = mapType + '_r'
        return mapType
    else:
        return mapType

def isCmap(mapType):
    if mapType in cmaps:
        return True
    else:
        print('Warning: Type mapType is not in the list of cmaps, default cmaps will be used.')

def cmapType(mapType=defaultCmap,reverse=False,totalBins=defaultBins):
    isCmap(mapType)
    mapType = reverseName(mapType,reverse)
    print('The map type is: ' + mapType)
    return {'mapType':mapType, 'bins':totalBins}

def plotConc(prtype = defaultProjection, gridID, lonData, latData, concDiff, removed=0, cmapType=defaultCmap, minVal=None, maxVal=None, extension=defaultExt):
    gridDecode = rmn.decodeGrid(gridID)
    # useful things for determining grid data depending on the projection type
    lowLeftLon = gridDecode['lon']
    lowLeftLat = gridDecode['lat0']
    lowLon = np.amin(lonData)
    lowLat = np.amin(latData)
    highLon = np.amax(lonData)
    highLat = np.amax(latData)
    maxDim = maxDimensions(lowLon,highLon,lowLat,highLat,buff)
    max_width = (highLon + lowLon)/2
    max_height = (highLat + lowLat)/2

    fig = plt.figure(figsize(8,8))

    if (prtype=='ortho') or (prtype == 'nsper') or (prtype == 'laea') or (prtype == 'aeqd') or (prtype == 'gnom') or (prtype == 'lcc'):
        concMap = Basemap(projection=prtype, resolution = 'c', lon_0=midLon, lat_0=midLat, width=max_width, height=max_height)
    elif prtype == 'stere':
        concMap = Basemap(projection=prtype, lon_0=midLon,lat_0=midLat,width=max_width,height=max_height)
    elif (prtype == 'cyl') or (prtype == 'merc'):
        concMap = Basemap(projection=prtype, resolution='c', llcrnrlat=lowLat,  urcrnrlat=highLat, llcrnrlon=lowLon, urcrnrlon=highLon)
    elif (prtype == 'aea') or (prtype == 'eqdc'):
        concMap = Basemap(projection = prtype, lon_0=midLon,  lat_0=midLat, llcrnrlat=lowLat, urcrnrlat=highLat, llcrnrlon=lowLon, urcrnrlon=highLon)
    else:
        print('Error: Could not generate map. Try a different projection.')

    mapColor = cmapType['mapType']
    mapBins = cmapType['bins']

    if removed != 0:
        ni = len(lonData)
        nj = len(lonData[0])
        n_pil = removed
        x, y = concMap(lonData[n_pil:ni-n_pil,n_pil:nj-n_pil], latData[n_pil:ni-n_pil,n_pil:nj-n_pil])
        concMap.pcolormesh(x,y,concDiff[n_pil:ni-n_pil,n_pil:nj-n_pil],cmap=plt.cm.get_cmap(cmapType['mapType'],cmapType['bins']))
    else:
        x, y = concMap(lonData, latData)
        concMap.pcolormesh(x, y, concDiff, cmap=plt.cm.get_cmap(mapColor, mapBins))

    concMap.drawcoastlines(color='lightgray')
    concMap.drawcountries(color='gray')
    concMap.drawparallels(np.arange(lowLat,highLat,15),labels=[1,0,0,0])
    concMap.drawmeridiams(np.arange(lowLon,highLon,15),labels=[0,0,0,1])
    plt.colorbar(extend = extension)
    plt.clim(minVal, maxVal)
    fig.savefig(defaultSave + defaultFile + time.strft('_%Y%m%d%H%M', time.localtime()) + '.png')
    plt.show()


def main():
    concDataStnd = getConc(filePath=defaultPathStnd)
    concDataAlt = getConc(filePath=defaultPathAlt)
    diffData = concDiff(concDataStnd,concDataAlt)
    gridStnd = getGrid(concDataStnd['dataKey'],concDataStnd['fileID'])
    gridAlt = getGrid(concDataStnd['dataKey'],concDataStnd['fileID'])
    gridCheck(gridStnd,gridAlt)
    cmapDetails = cmapType ()
    plotConc(gridID = gridStnd['gridID'], lonData = gridStnd['lonData'], latData = gridStnd['latData'], concDiff = diffData)

main()
