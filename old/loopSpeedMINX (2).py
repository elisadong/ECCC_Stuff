# writes a loop for each plume text file in the directory provided
# utilizes speedyAltMINXdemo.py

import os
import speedyAltMINXdemo as sam
import numpy as np
import math
import matplotlib.pyplot as plt

# constants, change as needed

plumeDir = 'BRPlumes'
fireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/septFSTs/ModelInput/HighRes2p5km'
noFireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/CFFEPSHighRes_nofire'
threshVal = 15

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
    
def mean_bias(optVal, modVal, error=False, normalized=False):
    bias = []
    optval, modelval = clean_inv(optVal, modVal)
    if not error:
        numer = sum([(oval - modelval[i])**2 for i, oval in enumerate(optval)])
    else:
        numer = sum([abs(oval - modelval[i])**2 for i, oval in enumerate(optval)])
    denom = len(optVal) if not normalized else sum(modelval)
    bias.append(numer / denom)
    return bias[0] if len(bias) == 1 else bias

def saveStats(plumeList, plmFRPs, plmModValues, plmHeightsB, plmHeightsR, plmModHeightsB, plmModHeightsR, threshVal, modifier = ''):
    plmHeights = plmHeightsB + plmHeightsR
    plmModHeights = plmModHeightsB + plmModHeightsR

    statsList = ['MINX Plume Max', 'Model Plume Max', 'MINX Plume Min', 'Model Plume Min', 'MINX Plume Mean', 'Model Plume Mean', 'RMSD', 'Pearson Correlation Coefficient, pValue', 'R Squared, pValue', 'Mean Absolute Percent Error', 'Mean Percentage Error', 'Mean Bias', 'Plumes Processed', 'Number of Valid Pairs']

    # may need to adjust the min/max functions if there are inf vals
    plmCleanB, modCleanB = clean_inv(plmHeightsB, plmModHeightsB)
    plmCleanR, modCleanR = clean_inv(plmHeightsR, plmModHeightsR)
    plmClean, modClean = clean_inv(plmHeights, plmModHeights)
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
    statsVals = [plmMax, modMax, plmMin, modMin, plmMean, modMean, rmsdVal, corrVal, rSquareVal, mapeVal, mpeVal, mbVal, numPlms, pairVal]

    f = open(os.path.join(curDir, '{}PlumeStats_thresh{}.txt'.format(modifier, threshVal)), 'w')
    for statInd, stat in enumerate(statsList):
        f.write('{}, {} \n'.format(stat,str(statsVals[statInd])))

    f.close()

    # plot methods
    # TODO: actually write these out as functions, not just lying around like this

    # plot comparison of heights
    # consider also outlining clusters from the same plumes
    fig = plt.figure()
    plt.plot(plmCleanB, modCleanB, 'bo')
    plt.plot(plmCleanR, modCleanR, 'ro')
    plt.plot([0.1, max(plmMax, modMax)], [0.1, max(plmMax,modMax)], 'k-')
    plt.axis('equal')
    plt.xlabel('MINX Plume Height')
    plt.ylabel('Model Plume Height')
    fig.savefig(os.path.join(curDir, '{}PlumeHeightComp_thresh{}.png'.format(modifier, threshVal)))
    plt.close('all')

    # plot FRP from plumes, used to see what kind of plumes were retained
    # may need the bins to be log scale if there is a big difference
    fig2 = plt.figure()
    plt.hist(clean_inv(plmFRPs))
    plt.title('Plume FRPs')
    fig2.savefig(os.path.join(curDir, '{}PlumeFRPs_thresh{}.png'.format(modifier, threshVal)))
    plt.close('all')

    # plot the AF values from models, used to see what kind of AF max values were picked out
    # may help in determining future threshold values
    fig3 = plt.figure()
    plt.hist(clean_inv(plmModValues))
    plt.title('Model Plume Values (Local Max of PM2.5)')
    fig3.savefig(os.path.join(curDir, '{}PlumeModValues_thresh{}.png'.format(modifier, threshVal)))
    plt.close('all')

curDir = os.getcwd()
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

print('The list of plumes is: ')
print(plumeList)

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


for pInd, plume in enumerate(plumeList):
    plmHeight, modHeight, modVal, plmFRP, plmDist = sam.pltPlume(plumePath, plume, mod_fire=fireDir, mod_nofire=noFireDir)
    plmFRPs +=[plmFRP]
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

plmModValues = plmModValuesB + plmModValuesR
#generate comparisons for all distances
saveStats(plumeList, plmFRPs, plmModValues, plmHeightsB, plmHeightsR, plmModHeightsB, plmModHeightsR, threshVal, modifier = '')
#generate comparisons for distances > 2km and <10km

plmHeightsDB = []
plmHeightsDR = []
plmModHeightsDB = []
plmModHeightsDR = []
plmModValuesDB = []
plmModValuesDR = []
distances = []
for dInd, distB in enumerate(plmDistB):
    if distB >= 2 and distB <= 10:
        plmHeightsDB += [plmHeightsB[dInd]]
        plmModHeightsDB += [plmModHeightsB[dInd]]
        plmModValuesDB += [plmModValuesB[dInd]]
        distances += [plmDistB[dInd]]
for dInd2, distR in enumerate(plmDistR):
    if distR >= 2 and distR <= 10:
        plmHeightsDR += [plmHeightsR[dInd2]]
        plmModHeightsDR += [plmModHeightsR[dInd2]]
        plmModValuesDR += [plmModValuesR[dInd2]]
        distances += [plmDistR[dInd2]]

plmModValuesD = plmModValuesDB + plmModValuesDR
saveStats(plumeList, plmFRPs, plmModValuesD, plmHeightsDB, plmHeightsDR, plmModHeightsDB, plmModHeightsDR, threshVal, modifier = 'Dist2_10_')

# plots distance vs height of model and MINX plumes
# used to check how many missing model values there are for each MINX plume point

figDist = plt.figure()
dist = plmDistB + plmDistR
plmHeights = plmHeightsB + plmHeightsR
plmModHeights = plmModHeightsB + plmModHeightsR
plt.plot(dist, plmHeights, 'bo', label = 'MINX Plume')
plt.plot(dist, plmModHeights, 'ro', label = 'Model Plume')
plt.xlabel('Distance (km)')
plt.ylabel('Plume Height (m)')
plt.legend()
figDist.savefig(os.path.join(curDir, 'DistancePlot_thresh{}.png'.format(threshVal)))
plt.close('all')

# adding residual plot
cleanPlm, cleanMod = clean_inv(plmHeights, plmModHeights)
residuals =[]
for point, (pl, md) in enumerate(zip(cleanPlm, cleanMod)):
    residuals.append(pl-md)
figResid = plt.figure()
plt.scatter(cleanPlm, residuals)
plt.xlabel('MINX Plume Height')
plt.ylabel('Plume Height (MINX - Model)')
figResid.savefig(os.path.join(curDir, 'ResidualPlot_thresh{}'.format(threshVal)))
print('All done! Please check to see plots and stats have been generated correctly!')
