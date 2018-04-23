## plot concentration values for each latlon coordinate


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

need: gridValues['lat'], gridValues['lon'], concValues['concData']

latDic = llGridData['lat']
lonDic = llGridData['lon']
concValDic = concData

def plotConcMap(latDic,lonDic,concValDic):
    plt.figure(figsize=(8,8))
    concMap = Basemap(projection='lcc', resolution = 'c', width=8E6, height=8E6, lat_0=45, lon_0=-100) ## beed to accomodate for different lon/lat sources
    latList = latDic.tolist() ##note that these individual values are still list format, need to process more
    lats = [] # initialize storage of all lat values
    i = 0
    for i in range(len(latList)):
        if i == 0:
            lats = latList[i]
        else:
            lats += latList[i]
    ##check
    # lats = np.array(lats)
    lonList = lonDic.tolist() #Note that there are as many arrays of lons as there are lat values
    lons = []
    j = 0
    for j in range(len(lonList)):
        if j == 0:
            lons =lonList[j]
        else:
            lons += lonList[j]
    # lons = np.array(lons)
    lons,lats = np.meshgrid(lons,lats)
    concMap.pcolormesh(lons,lats,concValDic,latlon=True)

    ##need to plot some smaller areas

    concMap.drawcoastlines()
    concMap.drawcountries()

    plt.colorbar()
    plt.clim()



## plot concentration values for each latlon coordinate


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

need: gridValues['lat'], gridValues['lon'], concValues['concData']

latDic = llGridData['lat']
lonDic = llGridData['lon']
concValDic = concData

def plotConcMap(latDic,lonDic,concValDic):
    plt.figure(figsize=(8,8))
    concMap = Basemap(projection='lcc', resolution = 'c', width=8E6, height=8E6, lat_0=45, lon_0=-100)
    latList = latDic.tolist()
    lonList = lonDic.tolist()

    ##should loop over for each lon value and construct opposing lat set,
    ## demo will only have one instance

    lon0 = lonList[0]
    ##check length of lon0, should be same as len(latList)
    i = 0
    lat0 = []
    for i in range(len(latList)):
        lat0 += [latList[i][0]]
    ## check len lat0, should be the same as lon0

    lon0,lat0 = np.meshgrid(lon0,lat0)
