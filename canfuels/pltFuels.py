# Nice function used to plot the canfuels map
# Added functionality for plotting points

import gdal
import os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

fuelFile = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/canfuels/canfuels/bns2013/vat.adf'

def plotFuelMap(fuelDir, plmOrigins=None, pltFuels = False, pltPoints = True, modifier = ''):
    a = 6378137
    f = 1/298.257222101
    b = a*(1-f)
    w = 5351000
    cellWidth = 250

    fuelData = gdal.Open(fuelFile, gdal.GA_ReadOnly)
    fuelArray = fuelData.ReadAsArray()
    yList = []
    count = 0
    for y in range(len(fuelArray)):
        yRow = [count]*len(fuelArray[1])
        yList += [yRow]
        count += cellWidth
    yArray = np.array(yList)
    yArray = np.flipud(yArray)
    xList = []
    for x in range(len(fuelArray)):
        xRow = range(0,w,250)
        xList += [xRow]
    xArray = np.array(xList)

    fig = plt.figure(figsize = (6,6))

    m = Basemap(projection='lcc', rsphere = (a,b), lat_1=49, lat_2=77, lat_0=49, lon_0=-95, llcrnrlon = -121.112, llcrnrlat = 38.35766, urcrnrlon =-12.02121, urcrnrlat = 63.05306, resolution = 'i')
    if pltPoints:
        xoList = []
        yoList = []
        for o in plmOrigins:
            xoList += [o[0]]
            yoList += [o[1]]
        m.scatter(xoList, yoList, 1, marker = '^',color = 'k', zorder=4) #,zorder = 10)

    if pltFuels:
        m.pcolormesh(xArray,yArray,fuelArray, vmin = 100, vmax = 123, cmap=plt.cm.get_cmap('RdBu',22)) #, zorder=10)
    #m.drawrivers()
    m.drawcountries(linewidth =0.5, color = 'gray')
    m.drawcoastlines(linewidth = 0.5, color = 'gray')
    m.drawstates(linewidth = 0.5, color = 'gray')
    #m.drawmeridians(np.arange(-170, -20,10), linewidth = 0.2, labels=[0,0,0,1], fontsize=3)
    #m.drawparallels(np.arange(40,60,10), linewidth = 0.2, labels=[1,0,0,0],fontsize=3)
    #plt.colorbar(shrink=0.5)
    #assume passed in the xy origins

    fig.savefig(os.path.join(fuelDir, 'fuel{}_distribution.png'.format(modifier)), dpi = 300)
    plt.close('all')
