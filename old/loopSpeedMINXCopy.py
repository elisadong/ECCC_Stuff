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

# constants, change as needed, some of these are not used if passing through from another function

plumeDir = 'BRPlumes'
firePixDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/plumesP50C60'
threshVal = 0
saveDirName = 'plumeAnalysis'
curDir = os.getcwd()
saveDir = os.path.join(curDir, saveDirName)

# uncomment the following for the appropriate directories for models with and without fire

# directories for CFFEPS High Res (sept2017 only)
#fireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/septFSTs/ModelInput/HighRes2p5km'
#noFireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/CFFEPSHighRes_nofire'

# directories for CFFEPS low res (sept2017 only)
#fireDir ='/fs/site2/dev/eccc/aq/r1/eld001/MINX/septFSTs/ModelInput/LowRes10km/processing'

# directories for FireWork Low Res (note that the nofire is just GEMMACH), note that both point to the 'science' folder
fireDir= '/fs/site2/dev/eccc/aq/r1/eld001/MINX/fireworkLowRes/mach'
noFireDir ='/fs/site2/dev/eccc/aq/r1/eld001/MINX/fireworkLowRes_nofire'

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

def pltFRPHist(plmFRPs, saveDir, threshVal, modifier=''):
    fig2 = plt.figure()
    plt.hist(clean_inv(plmFRPs))
    plt.title('Plume FRPs')
    fig2.savefig(os.path.join(saveDir, '{}PlumeFRPs_thresh{}.png'.format(modifier, threshVal)))
    plt.close('all')

def pltModValues(plmModValues, saveDir, threshVal, modifier=''):
    fig3 = plt.figure()
    plt.hist(clean_inv(plmModValues))
    plt.title('Model Plume Values (Local Max of PM2.5)')
    fig3.savefig(os.path.join(saveDir, '{}PlumeModValues_thresh{}.png'.format(modifier, threshVal)))
    plt.close('all')

def pltDistHist(plmDist, plmHeights,plmModHeights,saveDir,threshVal,pt,modifier=''):
    figDist = plt.figure()
    plt.plot(plmDist, plmHeights, 'ko', label = 'MINX Plume')
    plt.plot(plmDist, plmModHeights, 'go', label = 'Model Plume')
    plt.xlabel('Distance (km)')
    plt.ylabel('Plume Height (m)')
    plt.legend()
    figDist.savefig(os.path.join(saveDir, 'DistancePlot_thresh{}_prct{}{}.png'.format(threshVal, pt, modifier)))
    plt.close('all')

def pltHeightHist(plmHeights, plmModHeights,saveDir, threshVal, pt, modifier=''):
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
    cleanPlm, cleanMod = clean_inv(plmHeights, plmModHeights)
    residuals =[]
    for point, (pl, md) in enumerate(zip(cleanPlm, cleanMod)):
        residuals.append(pl-md)
    figResid = plt.figure()
    plt.scatter(cleanPlm, residuals)
    plt.xlabel('MINX Plume Height')
    plt.ylabel('Plume Height (MINX - Model)')
    figResid.savefig(os.path.join(saveDir, 'ResidualPlot_thresh{}_prct{}.png'.format(threshVal, pt)))
    plt.close('all')

def pltDistBins(cleanDist, distBins, saveDir, modifier=''):
    figDistH = plt.figure()
    #plt.hist(dist, bins = distBins, alpha=0.5)
    #plt.hist(cleanDist, bins=distBins,alpha=0.5)
    #plt.hist([distances,cleanDist], bins=distBins, label=['Distances', 'Plotted Distances'])
    plt.hist(cleanDist, bins = distBins)
    plt.xlim([np.min(distBins), np.max(distBins)])
    plt.xlabel('Distance from Origin (km)')
    figDistH.savefig(os.path.join(saveDir,'distanceHistogram{}.png'.format(modifier)))
    plt.close('all')

def pltDistRatio(cleanDist, cleanModPlmRatio, saveDir, mod=''):
    figDistRatio = plt.figure()
    plt.scatter(cleanDist, cleanModPlmRatio)
    plt.axhline(y=1, color='k', linestyle='-')
    plt.xlim(xmin=0)
    plt.xlabel('Distance from Origin (km)')
    plt.ylabel('Model/MINX Plume Height')
    figDistRatio.savefig(os.path.join(saveDir, 'distanceRatio{}.png'.format(mod)))
    plt.close('all')

def pltPlmHeightComp(biome, plmClean, modClean, saveDir, threshVal, modifier=''):
    compfig = plt.figure()
    plmMax, modMax = max(plmClean), max(modClean)
    plt.plot(plmClean, modClean, 'bo')
    plt.plot([0.1, max(plmMax, modMax)], [0.1, max(plmMax, modMax)], 'k-')
    plt.axis('equal')
    plt.xlabel('MINX Plume Height')
    plt.ylabel('Model Plume Height')
    plt.title('Plume Height Ratio, Biome: {}'.format(biome)) #this needs to be updated, can be fuel or biome
    compfig.savefig(os.path.join(saveDir, '{}{}_PlumeHeightComp_thresh{}.png'.format(biome, modifier, threshVal)))
    plt.close('all')

def saveStats(saveDir, plumeList, plmHeights, plmModHeights, minFRP, threshVal=threshVal, modifier = ''):
    statsList = ['Minimum FRP (MWatts)', 'MINX Plume Max', 'Model Plume Max', 'MINX Plume Min', 'Model Plume Min', 'MINX Plume Mean', 'Model Plume Mean', 'RMSD', 'Pearson Correlation Coefficient, pValue', 'R Squared, pValue', 'Mean Absolute Percent Error', 'Mean Percentage Error', 'Plumes Processed', 'Number of Valid Pairs', 'Number of Dropped Pairs']

    # may need to adjust the min/max functions if there are inf vals
    plmClean, modClean = clean_inv(plmHeights, plmModHeights)
    plmCleanB, modCleanB = clean_inv(plmHeightsB, plmModHeightsB)
    plmCleanR, modCleanR = clean_inv(plmHeightsR, plmModHeightsR)
    plmMax, plmMin, plmMean = np.nanmax(plmClean), np.nanmin(plmClean), np.nanmean(plmClean)
    modMax, modMin, modMean = np.nanmax(modClean), np.nanmin(modClean), np.nanmean(modClean)
    rmsdVal = RMSD(plmClean, modClean)
    corrVal = corr_coeff(plmClean, modClean)
    rSquareVal = linRegress(plmClean, modClean)
    mapeVal = mape(plmClean, modClean)
    mpeVal = mpe(plmClean, modClean)
    mbVal = mean_bias(plmClean, modClean)
    numPlms = len(plumeList)
    pairVal = len(plmClean)
    unPairVal = len(plmHeights)-len(plmClean)
    statsVals = [minFRP, plmMax, modMax, plmMin, modMin, plmMean, modMean, rmsdVal, corrVal, rSquareVal, mapeVal, mpeVal, numPlms, pairVal, unPairVal]

    f = open(os.path.join(saveDir, '{}PlumeStats_thresh{}.txt'.format(modifier, threshVal)), 'w')
    for statInd, stat in enumerate(statsList):
        f.write('{}, {} \n'.format(stat,str(statsVals[statInd])))
    f.close()

    # plot methods
    #pltPlmHeightCompBR(plmCleanB, modCleanB,plmCleanR, modCleanR, plmMax, modMax, saveDir, threshVal, modifier)

    # plot FRP from plumes, used to see what kind of plumes were retained
    #pltFRPHist(plmFRPs, saveDir, threshVal, modifier)

    # plot the AF values from models, used to see what kind of AF max values were picked out
    # may help in determining future threshold values
    #pltModValues(plmModValues, saveDir, threshVal, modifier)

def processPlumes(plumePath,plumeList, plmBandList,plmBandInd = 23, pt=0, minFRP=0, startTime = 00):
    # can save plm and plmModHeights as two separate lists if distinguishing between red and blue band recorded data
    plmHeightsB = []
    plmHeightsR = []
    plmDistB = []
    plmDistR = []
    plmModHeightsB = []
    plmModHeightsR = []
    plmModValuesB = []
    plmModValuesR = []
    plmFRPs = []
    plmBiomes = {}
    plmFuels = {}
    plumes = list(plumeList) #to prevent changing the original list
    keptPlumes = []

    for pInd, plume in enumerate(plumes):
        try:
            plmHeight, modHeight, modVal, plmFRP, plmDist, plmBiome,plmOrigin = sam.pltPlume(plumePath, plume, minFRP=minFRP, mod_fire=fireDir, mod_nofire=noFireDir, startTime = startTime)
            # shifted plmFRP check to go inside sam, though can check here again anyways, should get skipped over
            if plmFRP >= minFRP:
                plmFRPs +=[plmFRP]
                keptPlumes += [plume]
                if plmBiome not in plmBiomes.keys():
                    plmBiomes[plmBiome] = {'plmHeights':[], 'plmModHeights':[], 'plmDists':[]}
                plmBiomes[plmBiome]['plmHeights'] += plmHeight
                plmBiomes[plmBiome]['plmModHeights'] += modHeight
                plmBiomes[plmBiome]['plmDists'] += plmDist
                plmFuel = getFuel(plmOrigin[0], plmOrigin[1])
                if plmFuel not in plmFuels.keys():
                    plmFuels[plmFuel] = {'plmHeights':[], 'plmModHeights':[], 'plmDists':[]}
                plmFuels[plmFuel]['plmHeights'] += plmHeight
                plmFuels[plmFuel]['plmModHeights'] += modHeight
                plmFuels[plmFuel]['plmDists'] += plmDist
                # print(plmFRPs)
                if plmBandList[pInd] =='B':
                    plmHeightsB += plmHeight
                    plmModHeightsB += modHeight
                    plmDistB += plmDist
                    plmModValuesB += modVal
                else:
                    plmHeightsR += plmHeight
                    plmModHeightsR += modHeight
                    plmDistR += plmDist
                    plmModValuesR += modVal
        except:
            pass

    #plmModValues = plmModValuesB + plmModValuesR
    plmHeights = plmHeightsB + plmHeightsR
    plmModHeights = plmModHeightB + plModHeightsR
    plmDists =plmDistB + plmDistR
    return keptPlumes, plmHeights, plmModHeights, plmDist, plmBiomes, plmFuels, threshVal

def overallImages(plmDist, plmHeights, plmModHeights, threshVal, pt, saveDir = os.getcwd(), modifier=''):
    pltDistHist(plmDist, plmHeights,plmModHeights,saveDir,threshVal,pt,modifier)
    pltHeightHist(plmHeights, plmModHeights, saveDir, threshVal, pt, modifier)

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
    pltDistBins(cleanDist, distBins, saveDir, modifier)

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
        cleanPlm, cleanMod, cleanDist = clean_inv(plmBiomes[biome]['plmHeights'], plmBiomes[biome]['plmModHeights'], plmBiomes[biome]['plmDists'])
        cleanModPlmRatio = []
        for pInd, (plm, modPlm) in enumerate(zip(cleanPlm, cleanMod)):
            modPlmRatio = modPlm/plm # is this proper division?
            cleanModPlmRatio += [modPlmRatio]
        pltDistRatio(cleanDist, cleanModPlmRatio, biomeDir, mod = modifier + '_{}'.format(biome))
        pltPlmHeightComp(biome, cleanPlm, cleanMod, biomeDir, threshVal, modifier)

    # for each biome
    # plot with modifier biome name, just need the height over height?
    # update the modifier for each biome


    #generate comparisons for distances > 2km and <10km

#    plmHeightsDB = []
#    plmHeightsDR = []
#    plmModHeightsDB = []
#    plmModHeightsDR = []
#    plmModValuesDB = []
#    plmModValuesDR = []
#    distances = []
#    for dInd, distB in enumerate(plmDistB):
#        if distB >= 2 and distB <= 10:
#            plmHeightsDB += [plmHeightsB[dInd]]
#            plmModHeightsDB += [plmModHeightsB[dInd]]
#            plmModValuesDB += [plmModValuesB[dInd]]
#            distances += [plmDistB[dInd]]
#    for dInd2, distR in enumerate(plmDistR):
#        if distR >= 2 and distR <= 10:
#            plmHeightsDR += [plmHeightsR[dInd2]]
#            plmModHeightsDR += [plmModHeightsR[dInd2]]
#            plmModValuesDR += [plmModValuesR[dInd2]]
#            distances += [plmDistR[dInd2]]

#    plmModValuesD = plmModValuesDB + plmModValuesDR
#    saveStats(plumes, plmFRPs, plmModValuesD, plmHeightsDB, plmHeightsDR, plmModHeightsDB, plmModHeightsDR, threshVal, modifier = 'Dist2_10_')
