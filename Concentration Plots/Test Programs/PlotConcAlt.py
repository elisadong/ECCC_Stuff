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
    latList = latDic.tolist(
