# get the matching hotspots for all points in plumes
# note that datetime from the MINX output and the time in the hotspot are both in UTC

import postMINX as pm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime
import os
import math

plmDirect='BRPlumes/'
csvDir = '/fs/site2/dev/eccc/aq/r2/jac001/FireWork/CFFEPS/From_Kerry/20171114/outputs.weighted/FST'
curDir = os.getcwd()
plmDir = os.path.join(curDir, plmDirect)

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

def getmin(lat,lon,latList, lonList, zList):
    llDiffs = map(lambda x, y: abs(x-lat)+abs(y-lon), latList,lonList)
    minInd = llDiffs.index(min(llDiffs))
    hsLats = latList[minInd]
    hsLons = lonList[minInd]
    hsZs = zList[minInd]
    return hsLats, hsLons, hsZs

def RMSD(optVal, modVal):
    RMS = []
    optval, modelval = clean_inv(optVal, modVal)
    numer = sum([(oval - modelval[i])**2 for i, oval in enumerate(optval)])
    RMS.append((numer / len(optval))**0.5)
    return RMS[0] if len(RMS) == 1 else RMS

def corr_coeff(optVal, modVal):
    from scipy import stats
    arg1, arg2 = clean_inv(optVal, modVal)
    return stats.pearsonr(arg1, arg2)

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


def saveStats(plmList, plmHeights, zPlumes, modifier='HotSpot_'):
    statsList = ['MINX Plume Max', 'ZPlume Max', 'MINX Plume Min', 'ZPlume Min', 'RMSD', 'Correlation Coefficient', 'Mean Bias', 'Plumes Processed', 'Number of Valid Pairs']

    # may need to adjust the min/max functions if there are inf vals
    plmClean, zClean = clean_inv(plmHeights, zPlumes)
    plmMax, plmMin = np.nanmax(plmClean), np.nanmin(plmClean)
    zMax, zMin = np.nanmax(zClean), np.nanmin(zClean)
    rmsdVal = RMSD(plmClean, zClean)
    corrVal = corr_coeff(plmClean, zClean)
    mbVal = mean_bias(plmClean, zClean)
    numPlms = len(plmList)
    pairVal = len(plmClean)
    statsVals = [plmMax, zMax, plmMin, zMin, rmsdVal, corrVal, mbVal, numPlms, pairVal]

    f = open(os.path.join(curDir, '{}PlumeStats.txt'.format(modifier)), 'w')
    for statInd, stat in enumerate(statsList):
        f.write('{}, {} \n'.format(stat,str(statsVals[statInd])))

    f.close()

    # plot comparison of heights
    # consider also outlining clusters from the same plumes

    fig = plt.figure()
    plt.plot(plmClean, zClean, 'bo')
    plt.plot([0.1, max(plmMax, zMax)], [0.1, max(plmMax,zMax)], 'k-')
    plt.axis('equal')
    plt.xlabel('MINX Plume Height')
    plt.ylabel('Model Plume Height')
    fig.savefig(os.path.join(curDir, '{}PlumeHeightComp.png'.format(modifier)))
    plt.close('all')



plmList = os.listdir(plmDir)
zLats = []
zLons = []
zPlms = []
plmHeights = []
plmLons = []
plmLats = []
dateList = []
f2 = open(os.path.join(curDir, '{}PlumeDetails.txt'.format('HotSpot_')), 'w')
header='Date, MINX Lat, MINX Lon, Filtered Height, HS Lat, HS Lon, ZPlume\n'
f2.write(header)
for p in plmList:
    print(p)
    plm = pm.Plume.from_filename(os.path.join(plmDir, p))
    plmLon = plm['longitude']
    plmLat = plm['latitude']
    plmHeight = plm['filtered_height']
    plmLons += plmLon
    plmLats += plmLat
    plmHeights += plmHeight
    plmDate = plm.datetime
    dateStr = plmDate.strftime('%Y%m%d')
    if plmDate.minute < 30:
        dateTimeStr = plmDate.strftime('%Y%m%d_%H:') +'00'
    else:
        plm.datetime + datetime.timedelta(hours=1)
        dateTimeStr = plmDate.strftime('%Y%m%d_%H:') +'00'
    dateList += [dateTimeStr]*len(plmLat)
    csvFile = os.path.join(csvDir, dateStr, 'CFFEPS_{}.csv'.format(dateStr))
    # import the csv as a dataframe

    cdf = pd.read_csv(csvFile)
    cHeaders = list(cdf.columns.values)
    latVal = 0
    lonVal = 1
    UTCVal = 17
    zPlumeVal = 36

    cdf[cHeaders[UTCVal]] = cdf[cHeaders[UTCVal]].map(lambda x: x.strip())
    cdf = cdf.drop(cdf[cdf[cHeaders[UTCVal]] != dateTimeStr].index)
    latList = cdf[cHeaders[latVal]].tolist()
    lonList = cdf[cHeaders[lonVal]].tolist()
    zList = cdf[cHeaders[zPlumeVal]].tolist()

    for ind, (lon,lat) in enumerate(zip(plmLon, plmLat)):
        zLat, zLon, zPlm = getmin(lat, lon, latList, lonList, zList)
        llDiff = (abs(zLat - lat)+abs(zLon - lon))**2
        if llDiff <=0.0001:
            zLats += [zLat]
            zLons += [zLon]
            zPlms += [zPlm]
        else:
            zLats += [float('nan')]
            zLons += [float('nan')]
            zPlms += [float('nan')]

for pInd, plm in enumerate(plmLats):
    newLine='{}, {}, {}, {}, {}, {}, {}\n'.format(dateList[pInd],plmLats[pInd], plmLons[pInd], plmHeights[pInd], zLats[pInd], zLons[pInd], zPlms[pInd])
    f2.write(newLine)


f2.close()
saveStats(plmList, plmHeights, zPlms, modifier='HotSpot_')
