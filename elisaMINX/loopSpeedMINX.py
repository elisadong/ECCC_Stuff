#! /usr/bin/env python

# writes a loop for each plume text file in the directory provided
# utilizes speedyAltMINXdemo.py

import os
import speedyAltMINXdemo as sam
import numpy as np
import math
import matplotlib.pyplot as plt
import errno
import pandas as pd
from mpl_toolkits.basemap import Basemap
import gdal

# constants, change as needed, some of these are not used if passing through from another function

plumeDir = 'plumesP50C60/plumesCan'
firePixDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/fpFuelFilt'
threshVal = 0
saveDirName = 'plumeAnalysis'
curDir = os.getcwd()
saveDir = os.path.join(curDir, saveDirName)
threshVal = 0
a = 6378137
f = 1/298.257222101
b = a*(1-f)
h = 4558250
#w = 5351000

# uncomment the following for the appropriate directories for models with and without fire

# directories for CFFEPS High Res (sept2017 only)
#fireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/septFSTs/ModelInput/HighRes2p5km'
#noFireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/CFFEPSHighRes_nofire'

# directories for CFFEPS low res (sept2017 only)
#fireDir ='/fs/site2/dev/eccc/aq/r1/eld001/MINX/septFSTs/ModelInput/LowRes10km/processing'

# directories for FireWork Low Res (note that the nofire is just GEMMACH), note that both point to the 'science' folder
#fireDir= '/fs/site2/dev/eccc/aq/r1/eld001/MINX/fireworkLowRes/mach'
#noFireDir ='/fs/site2/dev/eccc/aq/r1/eld001/MINX/fireworkLowRes_nofire'
#operational directories
fireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/fireworkOPLowRes'
#fireDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/CFFEPSLowResMant4165"
noFireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/GEMMACHOPLowRes'

# files for fuel type
fuelFile = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/canfuels/canfuels/bns2013/vat.adf'

def makeDir(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def clean_inv(*listN, **kwargs):
    if not kwargs:
        kwargs = {'nans': True, 'infs': True, 'negatives': False, 'fills': (False, [-9999])}
    listN = [list_x for list_x in listN if isinstance(list_x, list)]
    if not listN:
        raise TypeError("arguments supplied must be lists, with an optional **kwargs dict.")
    try:
        if kwargs['nans']:
            listN = [[x for i, x in enumerate(sublistN) if not any([math.isnan(chk_sublistN[i]) for chk_sublistN in listN])] for sublistN in listN]
        if kwargs['infs']:
            listN = [[x for i, x in enumerate(sublistN) if not any([math.isinf(chk_sublistN[i]) for chk_sublistN in listN])] for sublistN in listN]
        if kwargs['negatives']:
            listN = [[x for i, x in enumerate(sublistN) if not any([chk_sublistN[i] < 0. for chk_sublistN in listN])] for sublistN in listN]
        if kwargs['fills'][0]:
            for fill_val in kwargs['fills'][1]:
                listN = [[x for i, x in enumerate(sublistN) if not any([chk_sublistN[i] == fill_val for chk_sublistN in listN])] for sublistN in listN]
    except KeyError:
        raise KeyError("kwargs must be in same format as defaults, see `help(clean_inv)`.")
    except IndexError:
        raise IndexError("all supplied lists must be of same size.")
    return listN if len(listN) > 1 else listN[0]

def getFuel(lon,lat, fuelArray):
    ''' Takes a give lon/lat and retrieves the corresponding fuel at the location.

    Note that this function may still be slightly broken (ie. not retrieving the correct fuel values).
    '''
    #offset values are only if you're checking against a GIS display with coordinates. In that case, subtract the offets from the GIS returned coordinates to find pixel location on the map
    #yOff = -724826.023
    #xOff = -2341250.000
    # NOTE: Indexing tends to be off by one yInd for pixels further north
    h = 4558250
    #Using x/y for visual reference. x is the horizontal axis and longitude, array first axis is latitude
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

    #NOTE: This method has now been updatd. For any given point in the qgis map, for gis (xVal, yVal) coordinates
    #xPix = xVal - xOff
    #yPix = yVal - yOff

def RMSD(optVal, modVal):
    RMS = []
    optval, modelval = clean_inv(optVal, modVal)
    numer = sum([(oval - modelval[i])**2 for i, oval in enumerate(optval)])
    RMS.append((numer / len(optval))**0.5)
    return RMS[0] if len(RMS) == 1 else RMS

def linRegress(x,y):
    from scipy import stats
    arg1, arg2 = clean_inv(x,y)
    slope, intercept, rVal, pVal, stndErr = stats.linregress(arg1, arg2)
    return (rVal**2, pVal, 'y = {}x + {}'.format(slope, intercept))

def corr_coeff(optVal, modVal):
    from scipy import stats
    arg1, arg2 = clean_inv(optVal, modVal)
    return stats.pearsonr(arg1, arg2)

def mape(x,y):
    # mean absolute percent error,
    # forecast value is model, true value is MINX plume
    arg1, arg2 = clean_inv(x,y)
    errors = []
    for ind, (xVal,yVal) in enumerate(zip(arg1,arg2)):
        try:
            err = np.abs((xVal-yVal)/xVal)
            errors.append(err)
        except:
            continue
    return np.mean(errors)* 100

def mpe(x,y):
    # mean percent error
    err = 0
    arg1, arg2 =clean_inv(x,y)
    for ind, (xVal, yVal) in enumerate(zip(arg1, arg2)):
        try:
            err += (xVal-yVal)/xVal
        except:
            continue
    return err* 100/len(arg1)

def mean_bias(optVal, modVal, error=False, normalized=True):
    bias = []
    optval, modelval = clean_inv(optVal, modVal)
    if not error:
        numer = sum([(oval - modelval[i])**2 for i, oval in enumerate(optval)])
    else:
        numer = sum([abs(oval - modelval[i])**2 for i, oval in enumerate(optval)])
    denom = len(optVal) if not normalized else sum(modelval)
    bias.append(numer / denom)
    return bias[0] if len(bias) == 1 else bias

def pltPlmHeightCompBR(plmCleanB, modCleanB,plmCleanR, modCleanR, plmMax, modMax, saveDir, threshVal, modifier=''):
    # plot comparison of heights
    # Probably don't need this one anymore, since not a big difference between B/R plumes
    fig = plt.figure()
    plt.plot(plmCleanB, modCleanB, 'bo')
    plt.plot(plmCleanR, modCleanR, 'ro')
    plt.plot([1, max(plmMax, modMax)], [1, max(plmMax,modMax)], 'k-')
    plt.axis('equal')
    plt.xlabel('MINX Plume Height')
    plt.ylabel('Model Plume Height')
    plt.title('MINX, Model Plume Ratio')
    fig.savefig(os.path.join(saveDir, '{}PlumeHeightComp_thresh{}.png'.format(modifier, threshVal)))
    plt.close('all')

def pltHistFRP(plmFRPs, saveDir, threshVal, modifier=''):
    # Also not likely to be used, since FRP is not representative of plume size/height given that it takes a cluster of fire pixels within the digitized plume
    fig2 = plt.figure()
    plt.hist(clean_inv(plmFRPs))
    plt.title('Plume FRPs')
    fig2.savefig(os.path.join(saveDir, '{}PlumeFRPs_thresh{}.png'.format(modifier, threshVal)))
    plt.close('all')

def pltHistModValues(plmModValues, saveDir, threshVal, modifier=''):
    fig3 = plt.figure()
    plt.hist(clean_inv(plmModValues))
    plt.title('Model Plume Values (Local Max of PM2.5)')
    fig3.savefig(os.path.join(saveDir, '{}PlumeModValuesHist_thresh{}.png'.format(modifier, threshVal)))
    plt.close('all')

def pltDistHeight(plmDist, plmHeights,plmModHeights,saveDir,threshVal,pt,modifier=''):
    figDist = plt.figure()
    plt.plot(plmDist, plmHeights, 'ko', label = 'MINX Plume')
    plt.plot(plmDist, plmModHeights, 'go', label = 'Model Plume')
    plt.xlabel('Distance (km)')
    plt.ylabel('Plume Height (m)')
    plt.title('Plume Heights by Distance from Plume Origin')
    plt.legend()
    figDist.savefig(os.path.join(saveDir, 'DistancePlot_thresh{}_prct{}{}.png'.format(threshVal, pt, modifier)))
    plt.close('all')

def pltHistHeight(plmHeights, plmModHeights,saveDir, threshVal, pt, modifier=''):
    figHist = plt.figure()
    clnPlmHeights, clnPlmModHeights = clean_inv(plmHeights, plmModHeights)
    plt.hist(clnPlmHeights, alpha = 0.5)
    plt.hist(clnPlmModHeights, alpha = 0.5)
    plt.legend(['MINX Plume Heights', 'Model Plume Heights'])
    plt.xlabel('Plume Heights')
    plt.title('Plume Height Histogram')
    figHist.savefig(os.path.join(saveDir, 'DistanceHist_thresh{}_prct{}{}.png'.format(threshVal, pt, modifier)))
    plt.close('all')

def pltResids(plmHeights, plmModHeights, saveDir, threshVal,pt):
    # likely not to be used, mostly to see if there is linear trend or not
    cleanPlm, cleanMod = clean_inv(plmHeights, plmModHeights)
    residuals =[]
    for point, (pl, md) in enumerate(zip(cleanPlm, cleanMod)):
        residuals.append(pl-md)
    figResid = plt.figure()
    plt.scatter(cleanPlm, residuals)
    plt.xlabel('MINX Plume Height')
    plt.ylabel('Plume Height (MINX - Model)')
    plt.title('Plume Comparison Residual Plot')
    figResid.savefig(os.path.join(saveDir, 'ResidualPlot_thresh{}_prct{}.png'.format(threshVal, pt)))
    plt.close('all')

def pltHistDistances(cleanDist, distBins, saveDir, modifier=''):
    #NOTE: May be more useful to only plot the furthest points
    figDistH = plt.figure()
    #plt.hist(dist, bins = distBins, alpha=0.5)
    #plt.hist(cleanDist, bins=distBins,alpha=0.5)
    #plt.hist([distances,cleanDist], bins=distBins, label=['Distances', 'Plotted Distances'])
    plt.hist(cleanDist, bins = distBins)
    plt.xlim([np.min(distBins), np.max(distBins)])
    plt.xlabel('Distance from Origin (km)')
    plt.title('Histogram of Points, at Distance from Plume Origin')
    figDistH.savefig(os.path.join(saveDir,'distanceHistogram{}.png'.format(modifier)))
    plt.close('all')

def pltDistRatio(cleanDist, cleanModPlmRatio, saveDir, mod=''):
    figDistRatio = plt.figure()
    plt.scatter(cleanDist, cleanModPlmRatio)
    plt.axhline(y=1, color='k', linestyle='-')
    plt.xlim(xmin=0)
    plt.xlabel('Distance from Origin (km)')
    plt.ylabel('Model/MINX Plume Height')
    plt.title('Plume Height Ratios by Distance from Plume Origin')
    figDistRatio.savefig(os.path.join(saveDir, 'distanceRatio{}.png'.format(mod)))
    plt.close('all')

def pltPlmHeightComp(type, subtype, plmClean, modClean, saveDir, threshVal, modifier=''):
    compfig = plt.figure()
    plmMax, modMax = max(plmClean), max(modClean)
    plt.plot(plmClean, modClean, 'bo')
    plt.plot([0.1, max(plmMax, modMax)], [0.1, max(plmMax, modMax)], 'k-')
    plt.axis('equal')
    plt.xlabel('MINX Plume Height')
    plt.ylabel('Model Plume Height')
    plt.title('Plume Height Ratio, {}: {}'.format(type, subtype)) #this needs to be updated, can be fuel or biome
    compfig.savefig(os.path.join(saveDir, 'PlumeHeightComp_{}_{}{}thresh{}.png'.format(type, subtype, modifier, threshVal)))
    plt.close('all')

def saveStats(saveDir, plumeDict, plmHeights, plmModHeights,minFRP, threshVal=threshVal, modifier = ''):
    #TODO: modify stats to deal with non BR split
    statsList = ['Minimum FRP (MWatts)', 'Plumes Processed', 'Number of Valid Pairs', 'Number of Dropped Pairs', 'MINX Plume Max', 'Model Plume Max', 'MINX Plume Min', 'Model Plume Min', 'MINX Plume Mean', 'Model Plume Mean', 'RMSD', 'Pearson Correlation Coefficient, pValue', 'R Squared, pValue', 'Mean Absolute Percent Error', 'Mean Percentage Error' ]

    # may need to adjust the min/max functions if there are inf vals
    plmClean, modClean = clean_inv(plmHeights, plmModHeights)
    plmMax, plmMin, plmMean = np.nanmax(plmClean), np.nanmin(plmClean), np.nanmean(plmClean)
    modMax, modMin, modMean = np.nanmax(modClean), np.nanmin(modClean), np.nanmean(modClean)
    rmsdVal = RMSD(plmClean, modClean)
    corrVal = corr_coeff(plmClean, modClean)
    rSquareVal = linRegress(plmClean, modClean)
    mapeVal = mape(plmClean, modClean)
    mpeVal = mpe(plmClean, modClean)
    mbVal = mean_bias(plmClean, modClean)
    numPlms = len(plumeDict.keys())
    pairVal = len(plmClean)
    unPairVal = len(plmHeights)-len(plmClean)
    statsVals = [minFRP, numPlms, pairVal, unPairVal, plmMax, modMax, plmMin, modMin, plmMean, modMean, rmsdVal, corrVal, rSquareVal, mapeVal, mpeVal]

    f = open(os.path.join(saveDir, '{}PlumeStats_thresh{}.txt'.format(modifier, threshVal)), 'w')
    for statInd, stat in enumerate(statsList):
        f.write('{}, {} \n'.format(stat,str(statsVals[statInd])))
    f.close()

    # plot FRP from plumes, used to see what kind of plumes were retained
    #pltHistFRP(plmFRPs, saveDir, threshVal, modifier)

    # plot the AF values from models, used to see what kind of AF max values were picked out
    # may help in determining future threshold values
    #pltHistModValues(plmModValues, saveDir, threshVal, modifier)

def saveStatsLite(saveDir, plmCount, plmHeights, plmModHeights, threshVal = threshVal, modifier = ''):
    statsList = ['Plumes Processed', 'Number of Valid Pairs', 'Number of Dropped Pairs', 'MINX Plume Max', 'Model Plume Max', 'MINX Plume Min', 'Model Plume Min', 'MINX Plume Mean', 'Model Plume Mean', 'RMSD', 'Pearson Correlation Coefficient, pValue', 'R Squared, pValue', 'Mean Absolute Percent Error', 'Mean Percentage Error' ]
    plmClean, modClean = clean_inv(plmHeights, plmModHeights)
    plmMax, plmMin, plmMean = np.nanmax(plmClean), np.nanmin(plmClean), np.nanmean(plmClean)
    modMax, modMin, modMean = np.nanmax(modClean), np.nanmin(modClean), np.nanmean(modClean)
    rmsdVal = RMSD(plmClean, modClean)
    corrVal = corr_coeff(plmClean, modClean)
    rSquareVal = linRegress(plmClean, modClean)
    mapeVal = mape(plmClean, modClean)
    mpeVal = mpe(plmClean, modClean)
    mbVal = mean_bias(plmClean, modClean)
    numPlms = plmCount
    pairVal = len(plmClean)
    unPairVal = len(plmHeights)-len(plmClean)
    statsVals = [numPlms, pairVal, unPairVal, plmMax, modMax, plmMin, modMin, plmMean, modMean, rmsdVal, corrVal, rSquareVal, mapeVal, mpeVal]
    f = open(os.path.join(saveDir, '{}PlumeStats_thresh{}.txt'.format(modifier, threshVal)), 'w')
    for statInd, stat in enumerate(statsList):
        f.write('{}, {} \n'.format(stat,str(statsVals[statInd])))
    f.close()

def processPlumes(plumePath,plumeList, plmBandList, plmBandInd = 23, pt=0, minFRP=0, startTime = 00):
    # can save plm and plmModHeights as two separate lists if distinguishing between red and blue band recorded data
    plmHeights = []
    plmDists = []
    plmModHeights = []
    plmFRPs = []
    plmBiomes = {}
    plmFuels = {}
    plumes = list(plumeList) #to prevent changing the original list
    keptPlumes = {}
    droppedPlumes = []
    m = Basemap(projection='lcc', rsphere = (a,b), lat_1=49, lat_2=77, lat_0=49, lon_0=-95, llcrnrlon = -121.112, llcrnrlat = 38.35766, urcrnrlon =-12.02121, urcrnrlat = 63.05306)
    fuelData = gdal.Open(fuelFile, gdal.GA_ReadOnly)
    fuelArray = fuelData.ReadAsArray()

    for pInd, plume in enumerate(plumes):
        try:
            plm = sam.pltPlume(plumePath, plume, minFRP, mod_fire=fireDir, mod_nofire=noFireDir, startTime = startTime)
            # shifted plmFRP check to go inside sam, though can check here again anyways, should get skipped over
            if plm['frp'] >= minFRP:
                keptPlumes[plume] = {'plume':plm}
                keptPlumes[plume]['band'] = plmBandList[pInd]
                plmOrigin = plm['origin']
                #print('origin get! {}'.format(plmOrigin))
                xPix, yPix = m(plmOrigin[0], plmOrigin[1])
                yInd = int((h-yPix)/250) #manual checking shows that you don't need to index down
                xInd = int(xPix/250)
                #print xInd, yInd
                fuel = fuelArray[yInd][xInd]
                #print('Fuel {} for origin {}, {} get'.format(fuel, plmOrigin[0], plmOrigin[1]))
                keptPlumes[plume]['xyOrigin'] = [xPix, yPix]
                keptPlumes[plume]['fuel'] = fuel
        except:
            droppedPlumes += [plume]
            #print('fuel not found, xyorigin not added')
            pass

    for kp in keptPlumes:
        kPlm = keptPlumes[kp]['plume']
        try:
            plmHeight, plmModHeight, plmDist, plmBiome, plmFuel, plmXYOrigin = kPlm['filtered_height'], kPlm.Model.height, kPlm['distance'], kPlm['biome'][0], keptPlumes[kp]['fuel'], keptPlumes[kp]['xyOrigin']
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

    return keptPlumes, plmHeights, plmModHeights, plmDists, plmBiomes, plmFuels, threshVal

def overallImages(plmDist, plmHeights, plmModHeights, pt = '',threshVal=threshVal, saveDir = os.getcwd(), modifier=''):
    pltDistHeight(plmDist, plmHeights,plmModHeights,saveDir,threshVal,pt,modifier)
    pltHistHeight(plmHeights, plmModHeights, saveDir, threshVal, pt, modifier)
    #plmClean, modClean = clean_inv(plmHeights, plmModHeights)
    #pltPlmHeightComp('all','', plmClean, modClean, saveDir, threshVal)

    # adding residual plot
    pltResid = False
    if pltResid == True:
        pltResids(plmHeights, plmModHeights, saveDir, threshVal,pt)

def distanceStats(plumeHeights, modHeights, distances, distBins=5, saveDir = os.getcwd(), modifier=''):
    #TODO: add functions processing/binning stuff into distances
    #Histogram for all points dealt with, have overlay of used points ontop of total points
    #plumeHeights = plmHeightsR + plmHeightsB
    #modHeights = plmModHeightsR + plmModHeightsB
    cleanPlm, cleanMod, cleanDist = clean_inv(plumeHeights, modHeights, distances)
    # histogram needs to be updated so that the final edge is undefined
    pltHistDistances(cleanDist, distBins, saveDir, modifier)

    cleanModPlmRatio = []
    for pInd, (plm, modPlm) in enumerate(zip(cleanPlm, cleanMod)):
        modPlmRatio = modPlm/plm # is this proper division?
        cleanModPlmRatio += [modPlmRatio]

    pltDistRatio(cleanDist, cleanModPlmRatio, saveDir, modifier)
    #use saveStats to get the same stats, but save them all to the same file into table format?
    # ^ need to change saveStats to return stats, and have option for save in that case
    #plot for model plume over MINX plume ie. model height/minx height
    # have a 1 line showing where the MINX plume is
    # have dashed lines for the various intervals

def filterPlmRange(distRange, plmHeights, plmModHeights, plmDists):
    distDF = pd.DataFrame(columns=['plmHeight', 'plmModHeight', 'plmDists'])
    distDF['plmHeight']=plmHeights
    distDF['plmModHeight']=plmModHeights
    distDF['plmDists']=plmDists
    distDF = distDF.drop(distDF[distDF['plmDists'] < distRange[0]].index)
    distDF = distDF.drop(distDF[distDF['plmDists'] > distRange[1]].index)
    dPlmHeights =list(distDF['plmHeight'])
    dPlmModHeights =list(distDF['plmModHeight'])
    dPlmDists = list(distDF['plmDists'])
    return dPlmHeights, dPlmModHeights, dPlmDists

def plotBiomes(biomeDir, plmBiomes, modifier=''):
    for biome in plmBiomes:
        try:
            cleanPlm, cleanMod, cleanDist = clean_inv(plmBiomes[biome]['plmHeights'], plmBiomes[biome]['plmModHeights'], plmBiomes[biome]['plmDists'])
            cleanModPlmRatio = []
            for pInd, (plm, modPlm) in enumerate(zip(cleanPlm, cleanMod)):
                modPlmRatio = modPlm/plm # is this proper division?
                cleanModPlmRatio += [modPlmRatio]
            pltDistRatio(cleanDist, cleanModPlmRatio, biomeDir, mod = modifier + '_biome{}'.format(biome))
            pltPlmHeightComp('biome', biome, cleanPlm, cleanMod, biomeDir, threshVal, modifier)
        except:
            pass


def plotFuels(fuelDir, plmFuels,threshVal, pt, modifier=''):
    makeDir(fuelDir)
    distRanges = [[0,10],[10,30],[30,60],[60,500]]
    distDir = os.path.join(fuelDir, 'distances')
    for fuel in plmFuels:
        try:
            plmHeights = plmFuels[fuel]['plmHeights']
            plmModHeights = plmFuels[fuel]['plmModHeights']
            plmDists = plmFuels[fuel]['plmDists']
            plmCount = plmFuels[fuel]['plmCount']
            cleanPlm, cleanMod, cleanDist = clean_inv(plmHeights, plmModHeights, plmDists)
            overallImages(plmDists,plmHeights,plmModHeights, threshVal = threshVal, pt = pt, saveDir = fuelDir, modifier = '_fuel{}'.format(fuel))
            saveStatsLite(fuelDir, plmCount, plmHeights, plmModHeights, threshVal = threshVal, modifier = '_fuel{}'.format(fuel))
            for dRange in distRanges:
                try:
                    dPlmHeights, dPlmModHeights, dPlmDists = sm.filterPlmRange(dRange,plmHeights, plmModHeights, plmDists)
                    overallImages(dPlmDists, dPlmHeights, dPlmModHeights, threshVal, pt, distDir, modifier='_{}_dist{}-{}'.format(fuel, dRange[0], dRange[1]))
                    saveStatsLite(distDir, 'N/A', dPlmHeights, dPlmModHeights, threshVal=threshVal, modifier='{}_dist{}-{}'.format(fuel,dRange[0], dRange[1]))
                except:
                    pass
            cleanModPlmRatio = []
            for pInd, (plm, modPlm) in enumerate(zip(cleanPlm, cleanMod)):
                modPlmRatio = modPlm/plm
                cleanModPlmRatio += [modPlmRatio]
            pltDistRatio(cleanDist, cleanModPlmRatio, fuelDir, mod = modifier + '_fuel{}'.format(fuel))
            pltPlmHeightComp('fuel', fuel, cleanPlm, cleanMod, fuelDir, threshVal, modifier)
            sam.plotFuelMap(fuelDir, plmOrigins = plmFuels[fuel]['plmXYOrigins'], modifier = '{}_pt{}_'.format(fuel, pt))
            print('Fuel {} map plotted'.format(fuel))
        except:
            pass
        # plot a fuel map for each

def plotMonthFuels(plumeDict, saveDir = os.getcwd()):
    import datetime
    import gdal
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    a = 6378137
    f = 1/298.257222101
    b = a*(1-f)
    w = 5351000
    cellWidth = 250
    fuelTypes = [255, 101, 102, 103, 104, 105, 106, 107,108,109, 113, 114, 116, 118, 119, 120, 121, 122]
    fuelColours = ['black', 'darkolivegreen', 'blanchedalmond', 'tan', 'gold', 'olive', 'c', 'navy', 'darkviolet', 'lime', 'peru', 'gainsboro','silver', 'aquamarine', 'orange', 'saddlebrown', 'y', 'g']
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
    m = Basemap(projection='lcc', rsphere = (a,b), lat_1=49, lat_2=77, lat_0=49, lon_0=-95, llcrnrlon = -121.112, llcrnrlat = 38.35766, urcrnrlon =-12.02121, urcrnrlat = 63.05306, resolution = 'i')

    plmMonths = {}
    for p in plumeDict:
        dateObj = plumeDict[p]['plume']['datetime']
        mth = int(dateObj.month)
        if mth not in plmMonths.keys():
            plmMonths[mth] ={'fuels': {}}
        try:
            fuel = plumeDict[p]['fuel']
            xyOrigin = plumeDict[p]['xyOrigin']
            if fuel not in plmMonths[mth]['fuels'].keys():
                plmMonths[mth]['fuels'][fuel] = {'xOrigin':[], 'yOrigin':[]}
            plmMonths[mth]['fuels'][fuel]['xOrigin'] += [xyOrigin[0]]
            plmMonths[mth]['fuels'][fuel]['yOrigin'] += [xyOrigin[1]]
        except:
            pass

    for mth in plmMonths:
        fig = plt.figure()
        for f in plmMonths[mth]['fuels']:
            fInd = fuelTypes.index(f)
            fCol = fuelColours[fInd]
            xoList = plmMonths[mth]['fuels'][f]['xOrigin']
            yoList = plmMonths[mth]['fuels'][f]['yOrigin']
            m.scatter(xoList, yoList, 2, marker = '^',color = fCol, label=f)
        m.drawcountries(linewidth =0.5, color = 'gray')
        m.drawcoastlines(linewidth = 0.5, color = 'gray')
        m.drawstates(linewidth = 0.5, color = 'gray')
        plt.legend()
        plt.title('Plumes in Month: {} '.format(mth))
        fig.savefig(os.path.join(saveDir, '{}_plmMapbyFuel.png'.format(mth)), dpi=300)
        plt.close('all')

    # for each biome
    # plot with modifier biome name, just need the height over height?
    # update the modifier for each biome


    #generate comparisons for distances > 2km and <10km
