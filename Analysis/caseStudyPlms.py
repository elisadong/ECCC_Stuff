# general plume analysis

import os
import postMINX as pm
import matplotlib.pyplot as plt
import numpy as np
import plumeAnalysis as pa

try:
    import cPickle as pickle
except ImportError:
    import pickle

import gdal
from mpl_toolkits.basemap import Basemap

fuelFile = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/canfuels/canfuels/bns2013/vat.adf'
fuelData = gdal.Open(fuelFile, gdal.GA_ReadOnly)
fuelArray = fuelData.ReadAsArray()

a = 6378137
f = 1/298.257222101
b = a*(1-f)
w = 5351000
h = 4558250
cellWidth = 250
m = Basemap(projection='lcc', rsphere = (a,b), lat_1=49, lat_2=77, lat_0=49, lon_0=-95, llcrnrlon = -121.112, llcrnrlat = 38.35766, urcrnrlon =-12.02121, urcrnrlat = 63.05306)

CFpklPath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/CF-GM_start00_pklPlmsOps/"
FWpklPath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/FW-GM_start00_pklPlmsOps/"

julPltPlms = ['Plume_O093359-B003-SPWR01', 'Plume_O093359-B003-SPWR02']
augPlms = ['Plume_O093906-B037-SPWR01', 'Plume_O093906-B037-SPWR02','Plume_O093906-B037-SPWR03','Plume_O093906-B037-SPWR04']
septPlms = ['Plume_O094211-B051-SPWR02','Plume_O094211-B051-SPWR03','Plume_O094211-B051-SPWR04','Plume_O094211-B051-SPWR05','Plume_O094211-B051-SPWR06','Plume_O094211-B051-SPWR07','Plume_O094211-B051-SPWR08','Plume_O094211-B052-SPWR01','Plume_O094211-B052-SPWR02','Plume_O094211-B052-SPWR03','Plume_O094211-B052-SPWR04','Plume_O094211-B052-SPWR05']

def getFuel(lon,lat, fuelArray):
    h = 4558250
    xPix, yPix = m(lon,lat)
    #NOTE: fuel array assumes topleft starts at 0,0, but map uses lowerleft as 0,0
    yInd = int((h-yPix)/250) #manual checking shows that you don't need to index down
    xInd = int(xPix/250)
    try:
        fuel = fuelArray[yInd][xInd]
        print('Fuel {} for origin {}, {} get'.format(fuel, lon, lat))
    except IndexError:
        print('Origin of fire out of bounds')
        return np.nan
    return fuel, [xPix, yPix]

plmCFs = []
plmFWs = []
savePlms = []
for p in septPlms:
    pklFile = p + '.pkl'
    try:
        CFplm = pickle.load(open(os.path.join(CFpklPath, pklFile), 'rb'))
        FWplm = pickle.load(open(os.path.join(FWpklPath, pklFile), 'rb'))
        plmCFs += [CFplm]
        plmFWs += [FWplm]
        savePlms += [p]
    except:
        print('did not add {}'.format(pklFile))

for pInd, (cPlm, fPlm, plm) in enumerate(zip(plmCFs, plmFWs, savePlms)):
    cName = 'CF-' + plm
    fName = 'FW-' + plm
    pm.Plume.scatter_3D(cPlm, option=['filtered_height', 'model_height'], save=True, filename='scatter_' + cName)
    plt.close('all')
    pm.Plume.scatter_3D(fPlm, option=['filtered_height', 'model_height'], save=True, filename='scatter_' + fName)
    plt.close('all')
    pm.Plume.hist(cPlm, option=['filtered_height', 'model_height'], save=True, filename='hist_' + cName)
    plt.close('all')
    pm.Plume.hist(fPlm, option=['filtered_height', 'model_height'], save=True, filename='hist_' + fName)
    plt.close('all')

plmHeightsCF = []
plmDistsCF = []
plmModHeightsCF = []
plmBiomesCF = {}
plmFuelsCF = {}
plumesCF = list(septPlms) #to prevent changing the original list
keptPlumesCF = {}
droppedPlumesCF = []

for pInd, (plume, plm) in enumerate(zip(septPlms, plmCFs)):
    keptPlumesCF[plume] = {'plume':plm}
    plmOrigin = plm['origin']
    xPix, yPix = m(plmOrigin[0], plmOrigin[1])
    yInd = int((h-yPix)/250)
    xInd = int(xPix/250)
    fuel = fuelArray[yInd][xInd]
    keptPlumesCF[plume]['xyOrigin'] = [xPix, yPix]
    keptPlumesCF[plume]['fuel'] = fuel

for kp in keptPlumesCF:
    kPlm = keptPlumesCF[kp]['plume']
    plmHeight, plmModHeight, plmDist, plmBiome, plmFuel, plmXYOrigin = kPlm['filtered_height'], kPlm.Model.height, kPlm['distance'], kPlm['biome'][0], keptPlumesCF[kp]['fuel'], keptPlumesCF[kp]['xyOrigin']
    plmHeightsCF += plmHeight
    plmModHeightsCF += plmModHeight
    plmDistsCF += plmDist
    if plmBiome not in plmBiomesCF.keys():
        plmBiomesCF[plmBiome] = {'plmHeights':[], 'plmModHeights':[], 'plmDists':[]}
    plmBiomesCF[plmBiome]['plmHeights'] += plmHeight
    plmBiomesCF[plmBiome]['plmModHeights'] += plmModHeight
    plmBiomesCF[plmBiome]['plmDists'] += plmDist
    if plmFuel not in plmFuelsCF.keys():
        plmFuelsCF[plmFuel] = {'plmHeights':[], 'plmModHeights':[], 'plmDists':[], 'plmXYOrigins':[], 'plmCount': 0}
    plmFuelsCF[plmFuel]['plmHeights'] += plmHeight
    plmFuelsCF[plmFuel]['plmModHeights'] += plmModHeight
    plmFuelsCF[plmFuel]['plmDists'] += plmDist
    plmFuelsCF[plmFuel]['plmXYOrigins'] += [plmXYOrigin]
    plmFuelsCF[plmFuel]['plmCount'] += 1

saveDir = os.getcwd()
pa.saveStats(saveDir, keptPlumesCF, plmHeightsCF, plmModHeightsCF, minFRP = 0, modifier = 'CFFEPS')
pa.overallImages(plmDistsCF, plmHeightsCF, plmModHeightsCF, saveDir)
plmClean, modClean = pa.clean_inv(plmHeightsCF, plmModHeightsCF)
pa.pltPlmHeightComp('CFFEPS','points', plmClean, modClean, saveDir)
pa.distanceStats(plmHeightsCF, plmModHeightsCF, plmDistsCF, saveDir,distBins=5, modifier='CFFEPS')
pa.plotMonthFuels(keptPlumesCF, saveDir)

plmHeightsFW = []
plmDistsFW = []
plmModHeightsFW = []
plmBiomesFW = {}
plmFuelsFW = {}
plumesFW = list(septPlms) #to prevent changing the original list
keptPlumesFW = {}
droppedPlumesFW = []

for pInd, (plume, plm) in enumerate(zip(septPlms, plmFWs)):
    keptPlumesFW[plume] = {'plume':plm}
    plmOrigin = plm['origin']
    xPix, yPix = m(plmOrigin[0], plmOrigin[1])
    yInd = int((h-yPix)/250)
    xInd = int(xPix/250)
    fuel = fuelArray[yInd][xInd]
    keptPlumesFW[plume]['xyOrigin'] = [xPix, yPix]
    keptPlumesFW[plume]['fuel'] = fuel

for kp in keptPlumesFW:
    kPlm = keptPlumesFW[kp]['plume']
    plmHeight, plmModHeight, plmDist, plmBiome, plmFuel, plmXYOrigin = kPlm['filtered_height'], kPlm.Model.height, kPlm['distance'], kPlm['biome'][0], keptPlumesFW[kp]['fuel'], keptPlumesFW[kp]['xyOrigin']
    plmHeightsFW += plmHeight
    plmModHeightsFW += plmModHeight
    plmDistsFW += plmDist
    if plmBiome not in plmBiomesFW.keys():
        plmBiomesFW[plmBiome] = {'plmHeights':[], 'plmModHeights':[], 'plmDists':[]}
    plmBiomesFW[plmBiome]['plmHeights'] += plmHeight
    plmBiomesFW[plmBiome]['plmModHeights'] += plmModHeight
    plmBiomesFW[plmBiome]['plmDists'] += plmDist
    if plmFuel not in plmFuelsFW.keys():
        plmFuelsFW[plmFuel] = {'plmHeights':[], 'plmModHeights':[], 'plmDists':[], 'plmXYOrigins':[], 'plmCount': 0}
    plmFuelsFW[plmFuel]['plmHeights'] += plmHeight
    plmFuelsFW[plmFuel]['plmModHeights'] += plmModHeight
    plmFuelsFW[plmFuel]['plmDists'] += plmDist
    plmFuelsFW[plmFuel]['plmXYOrigins'] += [plmXYOrigin]
    plmFuelsFW[plmFuel]['plmCount'] += 1

saveDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/CF-GM_start00/ops_postCanFilt/plumeApr_Oct/filt20k/20kplumes/plmTests/Sept_094211/FireWork/"
pa.saveStats(saveDir, keptPlumesFW, plmHeightsFW, plmModHeightsFW, minFRP = 0, modifier = 'FireWork')
pa.overallImages(plmDistsFW, plmHeightsFW, plmModHeightsFW, saveDir)
plmClean, modClean = pa.clean_inv(plmHeightsFW, plmModHeightsFW)
pa.pltPlmHeightComp('FireWork','points', plmClean, modClean, saveDir)
pa.distanceStats(plmHeightsFW, plmModHeightsFW, plmDistsFW, saveDir,distBins=5, modifier='FireWork')
pa.plotMonthFuels(keptPlumesFW, saveDir)
