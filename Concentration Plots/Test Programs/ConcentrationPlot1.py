
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
import rpnpy.librmn.all as rmn
defaultIp1=76696048
defaultSpc='TO3'
defaultFile='2010071000_001'
level=defaultIp1
spc=defaultSpc
fileName = defaultFile
fileID = rmn.fstopenall(fileName, rmn.FST_RO)
dataKey = rmn.fstinf(fileID, nomvar=spc, ip1=level)['key']
dataRec = rmn.fstluk(dataKey)
concData = dataRec['d']
fileMeta = rmn.fstprm(dataKey)
fileMeta['iunit'] = fileID
gridData = rmn.ezqkdef(fileMeta)
llGridData = rmn.gdll(gridData)
latData = llGridData['lat']
lonData = llGridData['lon']

ni = dataRec['ni']
nj = dataRec['nj']
n_pil=0 #crop the number of points from top/bot

## ToDo:
# find way to calculate width, height based on data available
# Adjust so there are projection options

# plt.figure(figsize=(8,8))

# get: max lat,max lon, min lat, min lon
# ID lon/lat of corners, and adjust by how many degrees
gridDecode = rmn.decodeGrid(gridData)

## Test values

prtype='lcc'

##note that the lon value is always positive, from E
## lat is not always positive
lowLon=gridDecode['lon0']
lowLat=gridDecode['lat0']
highLon=np.amax(lonData)
highLat=np.amax(latData)
midLon=(highLon-lowLon)/2
midLat=(highLat-lowLat)/2
circumEquator=40070E3
circumMerids=39931E3
max_width=(highLon-lowLon)/360*circumEquator #used to calculate the maximum width(m) assuming is at the equator
max_height=(highLat-lowLat)/180*circumMerids #used to calculate the maximum hiehgt(m) along the meridian

if prtype =='lcc':
    concMap = Basemap(projection=prtype,resolution = 'c', lon_0=midLon,lat_0=midLat,width=max_width,height=max_height)

if prtype == 'ortho':
    cornerAdjust = 0
    lowLeftLon= gridDecode['lon0']-cornerAdjust
    lowLeftLat= gridDecode['lat0']-cornerAdjust
    upRightLon= gridDecode['lon0']+gridDecode['dlon']*ni+cornerAdjust
    upRightLat= gridDecode['lat0']+gridDecode['dlat']*ni+cornerAdjust


    #concMap = Basemap(projection='lcc', resolution = 'c', width=8E6, height=8E6, lat_0=45, lon_0=-100)
    concMap = Basemap(projection=prtype, resolution = 'c', llcrnrlon=lowLeftLon, llcrnrlat=lowLeftLat, urcrnrlon=upRightLon,urcrnrlat=upRightLat)

## adjustments to x,y are based on n_pil for reduction (ie.rewrite the arrays)
# else x,y = storing the x,y data in arrays as lonData and latData
if n_pil != 0:
    x, y = concMap(lonData[n_pil:ni-n_pil,n_pil:nj-n_pil], latData[n_pil:ni-n_pil,n_pil:nj-n_pil])
    concMap.pcolormesh(x,y,concData[n_pil:ni-n_pil,n_pil:nj-n_pil])
else:
    x,y = concMap(lonData,latData)
    concMap.pcolormesh(x,y,concData)

concMap.drawcoastlines(color='lightgray')
plt.show()
