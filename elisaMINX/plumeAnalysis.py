# collection of all the basic functions used to break down postMINX results
# includes:
## postMINX product retrieval (plume and associated model components)
## basic statistical methods
## plotting methods
## breakdown into fuel, biome and distance components

# some of these imports may be repeated elsewhere. You may find that parts of the code are easier to deal with/copy over if you leave them in
import os
import niceUtilities as nu
import numpy as np
import math
import matplotlib.pyplot as plt
import errno
import pandas as pd
from mpl_toolkits.basemap import Basemap
import gdal

# general constants for retrieving plume data, these can be changed as needed
# these are default values for this function. The plumePath and picklePath can be passed in as arguments through a different function (ex. from exBinPercentiles.py)
plumePath = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/plumesP50C60/plumesCan/'
picklePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/FW-GM_start00_pklPlmsOps/" #FireWork less GEMMACH for time starting at 00
#picklePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/CF-GM_start00_pklPlmsOps/" #CFFEPS less GEMMACH for time starting at 00
currentDir = os.getcwd()
saveDir = os.path.join(currentDir, 'saveTestFolders')
biomeDir = os.path.join(saveDir, 'biomeBreakdown')
fuelDir = os.path.join(saveDir, 'fuelBreakdown')


# some specifics
fireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/fireworkOPLowRes' #FireWork directory
#fireDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/CFFEPSLowResMant4165" #CFFEPS directory
noFireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/GEMMACHOPLowRes' #GEMMACH directory
startTime = 0
threshVal = 0

# More specifics
fuelFile = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/canfuels/canfuels/bns2013/vat.adf'

# individual plume results. Choose to save plots from postMINX
savePlumePlts = False
modelStats = False
imageMethods = ['contour','colour', 'scatter_3D', 'surface_plot', 'direction_plot', 'scatter_2D', 'poly_2D', 'shape', 'hist']
numMethods = ['max', 'min','median','mean','std']
compMethods = ['table', 'abs_diff', 'linear_comp']
compNumMethods = ['RMSD', 'corr_coeff', 'mean_bias']

# import packages that are frequently used
import postMINX as pm
import niceUtilities as nu
from niceUtilities import clean_inv
import matplotlib.pyplot as plt

try:
    import cPickle as pickle
except ImportError:
    import pickle

# dealing with individual plumes
def pltMethods(plumePath, plume, plm, modelStats, saveMethod =True):
    ''' Methods grabbed from postMINX for an individual plume.

    Cycles through imageMethods and numMethods for the MINX retrieved plume,
    has the option to use imageMethods for the model retrieved plume as well.
    Defaults to saveMethod = True to save the images and number methods.

    To view the results without saving, set saveMethod = False and run:
    plt.show()
    Note that this is unadvised as all plots will show up at once.

    To view individual methods, please use postMINX directly.
    '''

    import pandas as pd
    saveDir = os.path.join(currentDir, 'postAltOutputs',plume[7:-4])
    nu.makeDir(saveDir)
    os.chdir(saveDir)
    for im in imageMethods:
        saveName = plume[7:-4] + im + '.png'
        try:
            eval('plm.{}(save={},filename="{}")'.format(im, saveMethod, saveName))
            plt.close('all')
        except:
            plt.close('all')
            continue
        print('Save {} success!'.format(saveName))

    # number methods for minx plume only
    numSave = plume[7:-4] + '.txt'
    option = 'filtered_height'
    numList = []
    for nu in numMethods:
        numVal = eval('plm.{}(option="{}")'.format(nu, option))
        numRow = [nu, numVal]
        numList += [numRow]

    numDF = pd.DataFrame(numList, columns=('method', option))
    numDF.to_csv(numSave, sep = ',', index=False)

    if modelStats:
        options = ['filtered_height', 'model_height']
        for im in imageMethods:
            try:
                saveName = 'Mod_{}{}.png'.format(plume[7:-4], im)
                eval('plm.{}(save={}, option={}, filename="{}")'.format(im,saveMethod, options, saveName))
                plt.close('all')
                print('Save {} success!'.format(saveName))
            except:
                plt.close('all')
                continue
        numSave='Mod_{}.txt'.format(plume[7:-4])
        numList = []
        option1 = 'filtered_height'
        option2 = 'model_height' # note that this won't actually get the model
        opt2 = plm.Model.height

    os.chdir(currentDir)

def pltPlume(plumePath, plume, minFRP, mod_fire, mod_nofire=None, threshVal=0, startTime = 12, toSave=savePlumePlts, modelStats=True, picklePath = picklePath):
    '''
    pltPlume() returns a plume object with FRP, filtered_height, model height, model values(AF at height), and distance from plume origin.

    Additional plotting methods and stats can be turned on by setting toSave=True and modelStats=True.
    See postMINX.py documentation for other parameters.
    '''
    nu.makeDir(picklePath)

    try:
        # check to see if the pickle file exists first
        plm=pickle.load(open(os.path.join(picklePath, plume[0:5] + plume[6:-3]+'pkl'), 'rb'))
        if plm['frp'] >=minFRP:
            if toSave  == True:
                try:
                    pltMethods(path, plume, plm, modelStats=True)
                except:
                    pass
            return plm
        else:
            #print('{} FRP too low, plume will not be used.'.format(plume))
            return None
    except:
        plm = pm.Plume.from_filename(os.path.join(plumePath,plume))
        if plm['frp'] >= minFRP:
            try:
                plm.get_model(fst_dir=fireDir, ctrl_dir=noFireDir, threshold=threshVal, filestart = startTime)
                plm.save_plume(os.path.join(picklePath, 'Plume_{}.pkl'.format(plume[7:26])))
                print('{} has been pickled.'.format(plume))
                if toSave == True:
                    try:
                        pltMethods(plumePath, plume, plm, modelStats=True)
                    except:
                        pass
                try:
                    return plm
                except:
                    #print('Plume and Model height were not returned.')
                    return None
            except:
                return None
        else:


            print('{} FRP too low, plume will not be used.'.format(plume))
            return None

# useful things/methods that help get more information
def getFuel(lon,lat, fuelArray):
    ''' Takes a give lon/lat and retrieves the corresponding fuel at the location. '''

    #offset values are only if you're checking against a GIS display with coordinates. In that case, subtract the offets from the GIS returned coordinates to find pixel location on the map
    #yOff = -724826.023
    #xOff = -2341250.000
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

def filterPlmRange(distRange, plmHeights, plmModHeights, plmDists):
    ''' Filter a given list of distances, plume MINX heights and plume model heights.

    Return the lists within the range provided by distRange.
    '''

    import pandas as pd
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

# statistics
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

# dealing with a list of plumes
def processPlumes(plumePath,plumeList, plmBandList,plmBandInd = 23, pt=0, minFRP=0, startTime = 00):
    ''' Retrieve plume objects from a list of plumes in a directory.

    Also retrieve some extra stuff (plumeFuels and plumeBiomes dictionaries, lists of heights/distances)

    ie. if there is a list of plume.txt files from MINX output, you can put them through here to get plume objects.
    '''

    from mpl_toolkits.basemap import Basemap
    import gdal
    a = 6378137
    f = 1/298.257222101
    b = a*(1-f)
    w = 5351000
    h = 4558250
    cellWidth = 250
    # can save plm and plmModHeights as two separate lists if distinguishing between red and blue band recorded data
    plmHeights = []
    plmDists = []
    plmModHeights = []
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
            plm = pltPlume(plumePath, plume, minFRP, mod_fire=fireDir, mod_nofire=noFireDir, startTime = startTime)
            if plm['frp'] >= minFRP:
                keptPlumes[plume] = {'plume':plm}
                keptPlumes[plume]['band'] = plmBandList[pInd]
                plmOrigin = plm['origin']
                print('origin get! {}'.format(plmOrigin))
                xPix, yPix = m(plmOrigin[0], plmOrigin[1])
                yInd = int((h-yPix)/250) #manual checking shows that you don't need to index down
                xInd = int(xPix/250)
                #print xInd, yInd
                fuel = fuelArray[yInd][xInd] #Note that the fuel data is stored in yInd, xInd sequence
                #print('Fuel {} for origin {}, {} get'.format(fuel, plmOrigin[0], plmOrigin[1]))
                keptPlumes[plume]['xyOrigin'] = [xPix, yPix]
                keptPlumes[plume]['fuel'] = fuel
        except:
            droppedPlumes += [plume]
            #print('fuel not found, xyorigin not added')
            pass
    # could potentially write everything after as a separate function, but it's nice to have it in one big chunk
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

# methods
def pltDistHeight(plmDist, plmHeights,plmModHeights,saveDir, modifier=''):
    ''' Plot plume heights by distance from origin. '''

    figDist = plt.figure()
    plt.plot(plmDist, plmHeights, 'ko', label = 'MINX Plume')
    plt.plot(plmDist, plmModHeights, 'go', label = 'Model Plume')
    plt.xlabel('Distance (km)')
    plt.ylabel('Plume Height (m)')
    plt.title('Plume Heights by Distance from Plume Origin')
    plt.legend()
    figDist.savefig(os.path.join(saveDir, 'DistancePlot{}.png'.format(modifier)))
    plt.close('all')

def pltHistHeight(plmHeights, plmModHeights,saveDir, modifier=''):
    ''' Plot a histogram of the plume MINX and Model heights. '''

    figHist = plt.figure()
    clnPlmHeights, clnPlmModHeights = clean_inv(plmHeights, plmModHeights)
    plt.hist(clnPlmHeights, alpha = 0.5)
    plt.hist(clnPlmModHeights, alpha = 0.5)
    plt.legend(['MINX Plume Heights', 'Model Plume Heights'])
    plt.xlabel('Plume Heights')
    plt.title('Plume Height Histogram')
    figHist.savefig(os.path.join(saveDir, 'DistanceHist{}.png'.format(modifier)))
    plt.close('all')

def pltResids(plmHeights, plmModHeights, saveDir, modifier = ''):
    ''' Plot a residual plot. This is used to see if there is a linear trend in the data. '''

    cleanPlm, cleanMod = clean_inv(plmHeights, plmModHeights)
    residuals =[]
    for point, (pl, md) in enumerate(zip(cleanPlm, cleanMod)):
        residuals.append(pl-md)
    figResid = plt.figure()
    plt.scatter(cleanPlm, residuals)
    plt.xlabel('MINX Plume Height')
    plt.ylabel('Plume Height (MINX - Model)')
    plt.title('Plume Comparison Residual Plot')
    figResid.savefig(os.path.join(saveDir, 'ResidualPlot{}.png'.format(modifier)))
    plt.close('all')

def pltHistDistances(cleanDist, distBins, saveDir, modifier=''):
    ''' Plot historgram of plume distances. '''

    #NOTE: May be more useful to only plot the furthest points
    figDistH = plt.figure()
    plt.hist(cleanDist, bins = distBins)
    plt.xlim([np.min(distBins), np.max(distBins)])
    plt.xlabel('Distance from Origin (km)')
    plt.title('Histogram of Points, at Distance from Plume Origin')
    figDistH.savefig(os.path.join(saveDir,'distanceHistogram{}.png'.format(modifier)))
    plt.close('all')

def pltDistRatio(cleanDist, cleanModPlmRatio, saveDir, modifier=''):
    ''' Plot plume height ratios by distance where cleanModPlmRatio = model/minx plume height. '''

    figDistRatio = plt.figure()
    plt.scatter(cleanDist, cleanModPlmRatio)
    plt.axhline(y=1, color='k', linestyle='-')
    plt.xlim(xmin=0)
    plt.xlabel('Distance from Origin (km)')
    plt.ylabel('Model/MINX Plume Height')
    plt.title('Plume Height Ratios by Distance from Plume Origin')
    figDistRatio.savefig(os.path.join(saveDir, 'distanceRatio{}.png'.format(modifier)))
    plt.close('all')

def pltPlmHeightComp(type, subtype, plmClean, modClean, saveDir, modifier=''):
    ''' Plot a comparison, MINX plume heights on x axis, and model plume heights on y axis. '''

    compfig = plt.figure()
    plmMax, modMax = max(plmClean), max(modClean)
    plt.plot(plmClean, modClean, 'bo')
    plt.plot([0.1, max(plmMax, modMax)], [0.1, max(plmMax, modMax)], 'k-')
    plt.axis('equal')
    plt.xlabel('MINX Plume Height')
    plt.ylabel('Model Plume Height')
    plt.title('Plume Heights , {}: {}'.format(type, subtype)) #this needs to be updated, can be fuel or biome
    compfig.savefig(os.path.join(saveDir, 'PlumeHeightComp_{}_{}{}.png'.format(type, subtype, modifier)))
    plt.close('all')

def saveStats(saveDir, plumeDict, plmHeights, plmModHeights,minFRP, modifier = ''):
    ''' Generate some statistics and save to a txt file. '''

    statsList = ['Minimum FRP (MWatts)', 'Plumes Processed', 'Number of Valid Pairs', 'Number of Dropped Pairs', 'MINX Plume Max', 'Model Plume Max', 'MINX Plume Min', 'Model Plume Min', 'MINX Plume Mean', 'Model Plume Mean', 'RMSD', 'Pearson Correlation Coefficient, pValue', 'R Squared, pValue', 'Mean Absolute Percent Error', 'Mean Percentage Error' ]

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

    f = open(os.path.join(saveDir, 'PlumeStats{}.txt'.format(modifier)), 'w')
    for statInd, stat in enumerate(statsList):
        f.write('{}, {} \n'.format(stat,str(statsVals[statInd])))
    f.close()

def saveStatsLite(saveDir, plmCount, plmHeights, plmModHeights, modifier = ''):
    ''' Generate statistics, Lite version. Used for stats after the overall plume dictionary has been thrown away. '''

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
    f = open(os.path.join(saveDir, 'PlumesStats{}.txt'.format(modifier)), 'w')
    for statInd, stat in enumerate(statsList):
        f.write('{}, {} \n'.format(stat,str(statsVals[statInd])))
    f.close()

# combination methods
def overallImages(plmDist, plmHeights, plmModHeights,saveDir, modifier=''):
    ''' Plot basics for heights (by distance from origin and histogram of heights). '''

    pltDistHeight(plmDist, plmHeights,plmModHeights,saveDir,modifier)
    pltHistHeight(plmHeights, plmModHeights, saveDir, modifier)
    #plmClean, modClean = clean_inv(plmHeights, plmModHeights)
    #pltPlmHeightComp('all','', plmClean, modClean, saveDir, threshVal)

def distanceStats(plumeHeights, modHeights, distances, saveDir,distBins=5,modifier=''):
    ''' This is a misnomer actually.

    Plots a histogram of distances then plots a ratio of model/minc plume points by distance.
    '''

    cleanPlm, cleanMod, cleanDist = clean_inv(plumeHeights, modHeights, distances)
    pltHistDistances(cleanDist, distBins, saveDir, modifier)
    cleanModPlmRatio = []
    for pInd, (plm, modPlm) in enumerate(zip(cleanPlm, cleanMod)):
        modPlmRatio = modPlm/plm # is this proper division?
        cleanModPlmRatio += [modPlmRatio]
    pltDistRatio(cleanDist, cleanModPlmRatio, saveDir, modifier)
    #NOTE: should include the savestats method here as well

    pltResid = False
    if pltResid == True:
        pltResids(plmHeights, plmModHeights, saveDir, threshVal,pt)

def plotMonthFuels(plumeDict, saveDir):
    ''' Plot maps of fuels by month. '''

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

# detailed methods
def plotBiomes(biomeDir, plmBiomes, modifier=''):
    ''' Get results by biome. '''

    nu.makeDir(biomeDir)
    for biome in plmBiomes:
        try:
            cleanPlm, cleanMod, cleanDist = clean_inv(plmBiomes[biome]['plmHeights'], plmBiomes[biome]['plmModHeights'], plmBiomes[biome]['plmDists'])
            cleanModPlmRatio = []
            for pInd, (plm, modPlm) in enumerate(zip(cleanPlm, cleanMod)):
                modPlmRatio = modPlm/plm # is this proper division?
                cleanModPlmRatio += [modPlmRatio]
            pltDistRatio(cleanDist, cleanModPlmRatio, biomeDir, modifier = modifier + '_biome{}'.format(biome))
            pltPlmHeightComp('biome', biome, cleanPlm, cleanMod, biomeDir, modifier)
        except:
            pass

def plotFuels(fuelDir, plmFuels,threshVal, pt, modifier=''):
    ''' Get results by fuel type. '''

    nu.makeDir(fuelDir)
    distRanges = [[0,10],[10,30],[30,60],[60,500]]
    distDir = os.path.join(fuelDir, 'distances')
    for fuel in plmFuels:
        try:
            plmHeights = plmFuels[fuel]['plmHeights']
            plmModHeights = plmFuels[fuel]['plmModHeights']
            plmDists = plmFuels[fuel]['plmDists']
            plmCount = plmFuels[fuel]['plmCount']
            cleanPlm, cleanMod, cleanDist = clean_inv(plmHeights, plmModHeights, plmDists)
            overallImages(plmDists,plmHeights,plmModHeights, saveDir = fuelDir, modifier = modifier +'_fuel{}'.format(fuel))
            saveStatsLite(fuelDir, plmCount, plmHeights, plmModHeights, modifier = modifier +'_fuel{}'.format(fuel))
            for dRange in distRanges:
                try:
                    dPlmHeights, dPlmModHeights, dPlmDists = filterPlmRange(dRange,plmHeights, plmModHeights, plmDists)
                    overallImages(dPlmDists, dPlmHeights, dPlmModHeights, distDir, modifier=modifier +'_{}_dist{}-{}'.format(fuel, dRange[0], dRange[1]))
                    saveStatsLite(distDir, 'N/A', dPlmHeights, dPlmModHeights, modifier=modifier +'{}_dist{}-{}'.format(fuel,dRange[0], dRange[1]))
                except:
                    pass
            cleanModPlmRatio = []
            for pInd, (plm, modPlm) in enumerate(zip(cleanPlm, cleanMod)):
                modPlmRatio = modPlm/plm
                cleanModPlmRatio += [modPlmRatio]
            pltDistRatio(cleanDist, cleanModPlmRatio, fuelDir, modifier = modifier + '_fuel{}'.format(fuel))
            pltPlmHeightComp('fuel', fuel, cleanPlm, cleanMod, fuelDir, modifier)
            nu.plotFuelMap(fuelDir, plmOrigins = plmFuels[fuel]['plmXYOrigins'], modifier = modifier +'_{}'.format(fuel))
            print('Fuel {} map plotted'.format(fuel))
        except:
            pass

# out of use
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
