# check for plumes with model height points > 20000
# checking for CFFEPS, remember to change sam and sm to point at CFFEPS pkls

import os
import numpy as np
import loopSpeedMINX as sm

saveDirName = 'CF-GM_start00/ops_postCanFilt/plumeApr_Oct' #UPDATE THIS, FW-GM is firework, CF-GM is CFFEPS
curDir = os.getcwd()
saveDir = os.path.join(curDir, saveDirName)

plumeDir = os.path.join(curDir, 'plumesP50C60/plumesCan')
pixDir = os.path.join(curDir, 'fpFuelFilt')
timeStart = 00

def toFloats(row):
    floats = map(float, row.split())
    return floats

def readFile(path):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    return lines

plumePath = os.path.join(curDir, plumeDir)
tempList = os.listdir(plumePath)
plumeList= []
plmBandInd = 23
plmBandList = []
# make sure is only plumes, if maybePlume is plume, add to list
for ind, mPlume in enumerate(tempList):
    if mPlume.endswith('txt'):
        plumeList += [mPlume]
        if mPlume[plmBandInd] == 'B':
            plmBandList += 'B'
        else:
            plmBandList += 'R'


pwrInd = 5
pixList = os.listdir(pixDir)
pwrs = []
for t in pixList:
    tLines = readFile(os.path.join(pixDir,t))
    for rowInd, row in enumerate(tLines, 4):
        if rowInd < len(tLines):
            pwrs += [toFloats(tLines[rowInd])[pwrInd]]

def getplms20k(plumes):
    plms20k = {}
    plms20kLess = {}
    for p in plumes:
        modHeights = plumes[p]['plume'].Model.height
        if any(x > 20000 for x in modHeights):
            plms20k[p] = plumes[p]
        else:
            plms20kLess[p] = plumes[p]
    plmHeights = []
    plmModHeights = []
    plmDists = []
    plmBiomes = {}
    plmFuels = {}
    for p20L in plms20kLess:
        pl = plms20kLess[p20L]['plume']
        try:
            plmHeight, plmModHeight, plmDist, plmBiome, plmFuel, plmXYOrigin = pl['filtered_height'], pl.Model.height, pl['distance'], pl['biome'][0], plms20kLess[p20L]['fuel'], plms20kLess[p20L]['xyOrigin']
            plmHeights += plmHeight
            plmModHeights += plmModHeight
            plmDists += plmDist
            if plmBiome not in plmBiomes.keys():
                plmBiomes[plmBiome] = {'plmHeights':[], 'plmModHeights':[], 'plmDists':[]}
            plmBiomes[plmBiome]['plmHeights'] += plmHeight
            plmBiomes[plmBiome]['plmModHeights'] += plmModHeight
            plmBiomes[plmBiome]['plmDists'] += plmDist
            if plmFuel not in plmFuels.keys():
                plmFuels[plmFuel] = {'plmHeights':[], 'plmModHeights':[], 'plmDists':[], 'plmXYOrigins':[], 'plmCount': 0}
            plmFuels[plmFuel]['plmHeights'] += plmHeight
            plmFuels[plmFuel]['plmModHeights'] += plmModHeight
            plmFuels[plmFuel]['plmDists'] += plmDist
            plmFuels[plmFuel]['plmXYOrigins'] += [plmXYOrigin]
            plmFuels[plmFuel]['plmCount'] += 1
        except:
            continue
    return plms20k, plms20kLess, plmHeights, plmModHeights, plmDists, plmBiomes, plmFuels

save20Dir = os.path.join(saveDir, 'filt20k')
sm.makeDir(save20Dir)
percentiles = [50,95,98]
distRanges = [[0,10],[10,30],[30,50],[50,500]]
distBins = [0,5,10,20,30,50,100]
plms20kList = [] #look at these after
# pass percentiles as the modifier
for pt in percentiles:
    savePtDir = os.path.join(save20Dir, '{}_percentile'.format(pt))
    sm.makeDir(savePtDir)
    biomeDir = os.path.join(savePtDir, 'biomeBreakdown')
    fuelDir = os.path.join(savePtDir,'fuelBreakdown')
    pwrMax = np.percentile(pwrs, pt)
    print('pwrMax: {}'.format(pwrMax))
    plumes, _, _, _, _, _, threshVal = sm.processPlumes(plumeDir, list(plumeList), plmBandList, pt=pt, minFRP=pwrMax, startTime = timeStart)
    plms20k, plms20kLess, plmHeights, plmModHeights, plmDists, plmBiomes, plmFuels = getplms20k(plumes)
    plms20kList += [plms20k]
    print('Saved plumes for percentile: {}'.format(pt))
    try:
        sm.saveStats(savePtDir, plms20kLess, plmHeights, plmModHeights, minFRP = pwrMax, threshVal=threshVal, modifier = '{}perc_'.format(pt))
        print('Saved stats for percentile: {}'.format(pt))
    except:
        print('Insufficient number of plumes to get stats.')
    try:
        sm.overallImages(plmDists, plmHeights, plmModHeights, pt, threshVal, savePtDir)
        print('Plotted basic images for percentile: {}'.format(pt))
        plmClean, modClean = sm.clean_inv(plmHeights, plmModHeights)
        sm.pltPlmHeightComp('pt',pt, plmClean, modClean, savePtDir, threshVal)
    except:
        print('Insufficient number of points to plot.')
    try:
        sm.distanceStats(plmHeights, plmModHeights, plmDists, distBins, savePtDir, modifier='_prct{}'.format(pt))
        print('Saved distance stats for percentile: {}'.format(pt))
    except:
        print('Insufficient number of points to run distance stats.')
    try:
        #TODO: Can potentially move this out, only really need the 98percentile
        sm.makeDir(biomeDir)
        sm.plotBiomes(biomeDir, plmBiomes, modifier='_prct{}'.format(pt))
        print('Plotted plumes by biomes.')
    except:
        print('Distances by biomes were not plotted.')
    try:
        sm.makeDir(fuelDir)
        sm.plotFuels(fuelDir, plmFuels, threshVal, pt, modifier='_prct{}'.format(pt))
    except:
        print('Distances by fuel type were not plotted.')
    try:
        sm.makeDir(fuelDir)
        sm.plotMonthFuels(plms20kLess, fuelDir)
        print('Monthly fuel maps have been plotted')
    except:
        print('Failed to save monthly fuel maps')
    distanceDir = os.path.join(savePtDir, 'distancePlots')
    sm.makeDir(distanceDir)
    for dRange in distRanges:
        try:
            dPlmHeights, dPlmModHeights, dPlmDists = sm.filterPlmRange(dRange,plmHeights, plmModHeights, plmDists)
            sm.overallImages(dPlmDists, dPlmHeights, dPlmModHeights, threshVal, pt, distanceDir, modifier='_dist{}-{}'.format(dRange[0], dRange[1]))
        except:
            print('Did not save distance breakdown plot for dRange: {}'.format(dRange))
    print('All done! Please check to see plots and stats for percentile: {} have been generated correctly!'.format(pt))

toCheck = plms20kList[0]
plms20ktxt = toCheck.keys()


import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

a = 6378137
f = 1/298.257222101
b = a*(1-f)
h = 4558250
m = Basemap(projection='lcc', rsphere = (a,b), lat_1=49, lat_2=77, lat_0=49, lon_0=-95, llcrnrlon = -121.112, llcrnrlat = 38.35766, urcrnrlon =-12.02121, urcrnrlat = 63.05306)

xList = []
yList = []
lonlatList = []
fuelList = []
filtHeightList = []
modHeightList = []
dateList = []
for p20 in toCheck:
    x = toCheck[p20]['xyOrigin'][0]
    y = toCheck[p20]['xyOrigin'][1]
    fuel = toCheck[p20]['fuel']
    plmHeights = toCheck[p20]['plume']['filtered_height']
    modHeights = toCheck[p20]['plume'].Model.height
    date = str(toCheck[p20]['plume']['datetime'])
    xList += [x]
    yList += [y]
    fuelList += [fuel]
    filtHeightList += [plmHeights]
    modHeightList += [modHeights]
    lonlatList += [m(x,y,inverse=True)]
    dateList += [date]

t = open(os.path.join(save20Dir, 'plumes_over_20k.txt'),'w')
t.write('plume, x, y, lon/lat, date, fuel\n')
for pInd, p20 in enumerate(toCheck):
    tStr = '{}, {}, {}, {}, {}, {}\n'.format(p20, xList[pInd], yList[pInd],lonlatList[pInd],dateList[pInd],fuelList[pInd])
    t.write(tStr)
t.close()



fig = plt.figure()
fuelTypes = [255, 101, 102, 103, 104, 105, 106, 107,108,109, 113, 114, 116, 118, 119, 120, 121, 122]
fuelColours = ['black', 'darkolivegreen', 'blanchedalmond', 'tan', 'gold', 'olive', 'c', 'navy', 'darkviolet', 'lime', 'peru', 'gainsboro','silver', 'aquamarine', 'orange', 'saddlebrown', 'y', 'g']

fuelBreak = {}
for flInd, f in enumerate(fuelList):
    if f not in fuelBreak.keys():
        fuelBreak[f] = {'xoList': [], 'yoList': []}
    fuelBreak[f]['xoList']+= [xList[flInd]]
    fuelBreak[f]['yoList']+=[yList[flInd]]

for fb in fuelBreak:
    fInd = fuelTypes.index(fb)
    fCol = fuelColours[fInd]
    m.scatter(xList, yList, 1, marker = '^', color = fCol, label =fb)

m.drawcountries(linewidth =0.5, color = 'gray')
m.drawcoastlines(linewidth = 0.5, color = 'gray')
m.drawstates(linewidth = 0.5, color = 'gray')
plt.legend()
plt.title('Plumes with heights over 20k')
fig.savefig(os.path.join(save20Dir, 'plumes over 20k.png'), dpi=300)
