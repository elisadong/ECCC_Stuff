#! /usr/bin/env python

# used to filter out a folder(fPixDir) containing firepixel files (MODIS hotspots) by removing all all spots with fuel values that are not provided by the canfuels map (ie. remove all hot spots that are not inside of canada or do not have an associated fuel value).
# The newly modified firepixel files get saved into saveDir. Recommended quick check is to see the total number of tiles is the same as in fPixDir as firepixel files without hotspots will still have their headers saved in the file.
# also plots a map of all the filtered out hotspots so you can visually check that you haven't removed any hotspots in canada by accident.
# fuel file is used to get an array of fuel values from the canfuels map
# fuel types are all possible fuel types that are not '255' in the fuel array
# run 'gdalinfo vat.adf' in command line for more information on the adf file 

import os
import numpy as np
import sys
import pandas as pd
from mpl_toolkits.basemap import Basemap
import gdal

fuelFile = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/canfuels/canfuels/bns2013/vat.adf'
lonInd = 0
latInd = 1
saveDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/fpFuelFilt/"
fPixDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/allFirePix/"
fuelTypes = [101, 102, 103, 104, 105, 106, 107,108,109, 113, 114, 116, 118, 119, 120, 121, 122]

def toFloats(row):
    floats = map(float, row.split())
    return floats

def readFile(path):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    return lines

def writeFile(path, lines, writeInds):
    f = open(path, 'w')
    for rowInd, string in enumerate(lines):
        if rowInd in writeInds:
            f.write(lines[rowInd])
    f.close()

a = 6378137
f = 1/298.257222101
b = a*(1-f)
w = 5351000
m = Basemap(projection='lcc', rsphere = (a,b), lat_1=49, lat_2=77, lat_0=49, lon_0=-95, llcrnrlon = -121.112, llcrnrlat = 38.35766, urcrnrlon =-12.02121, urcrnrlat = 63.05306, resolution = 'i')
fuelData = gdal.Open(fuelFile, gdal.GA_ReadOnly)
fuelArray = fuelData.ReadAsArray()

def getFuel(lon,lat):
    h = 4558250
    xPix, yPix = m(lon, lat)
    yInd = int((h-yPix)/250)
    xInd = int(xPix/250)
    try:
        fuel = fuelArray[yInd][xInd]
        return int(fuel), xPix, yPix
    except:
        return np.nan, np.nan, np.nan

txtList = os.listdir(fPixDir)
for txt in txtList:
    if os.path.isfile(os.path.join(saveDir, txt)) == False:
        tLines =readFile(os.path.join(fPixDir,txt))
        writeInds = [0,1,2,3]
        fpList =[]
        for rowInd, row in enumerate(tLines, 4):
            if rowInd < len(tLines):
                rowList = toFloats(tLines[rowInd])
                #add fuel value
                rowList += [getFuel(rowList[0], rowList[1])[0]]
                fpList += [rowList]
        try:
            tdf = pd.DataFrame(fpList)
            tdf['index'] = tdf.index + 4
            tdf['isFuel'] = tdf[12].map(lambda x: x in fuelTypes)
            print (tdf)
            for r in tdf.index:
                if tdf['isFuel'][r] == True:
                    writeInds += [tdf['index'][r]]
        except:
            continue
        writeFile(os.path.join(saveDir, txt), tLines, writeInds)
        print('File {} written to {}'.format(txt, saveDir))

# plot a map of all the non 255 fuels
xList = []
yList = []
for txt in txtList:
    tLines = readFile(os.path.join(fPixDir,txt))
    for rowInd, row in enumerate(tLines, 4):
        if rowInd < len(tLines):
            rowList = toFloats(tLines[rowInd])
            fuel, x, y = getFuel(rowList[0],rowList[1])
            if fuel not in fuelTypes:
                xList += [x]
                yList += [y]

fig = plt.figure()
m.scatter(xList, yList, 1, marker = 'o', color = 'k')
m.drawcountries(linewidth =0.5, color = 'gray')
m.drawcoastlines(linewidth = 0.5, color = 'gray')
m.drawstates(linewidth = 0.5, color = 'gray')
plt.title('Hotspots with no fuel')
fig.savefig('255s.png', dpi = 300)
