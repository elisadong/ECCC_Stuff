#!/usr/bin/python

import numpy as np
import os
import sys
import ast
import rpnpy.librmn.all as rmn
import errno

# get the statistics for the concentration data that is passed to this funciton

sample = np.array([[5,3,63246], [4,25,0.4],[-52,5,25]])
sampleData1 = np.array([[5,3,63246], [4,25,0.4],[-52,5,25]])
sampleData2 = np.array([[6,3,63246], [4,35,0.4],[-52,5,2563]])
samples=(sampleData1,sampleData2)
sampleName= 'data'
sampleNames=('Data1', 'Data2')
saveTo = os.getcwd()
fileName = 'test.csv'

def toLit(string):
    """ Converts string to literal representation.

    Useful when taking in arguments from command line which are all passed as strings. In this program, converts to lists and dictionaries.
    """

    nowLit = ast.literal_eval(string)
    return nowLit

def getConc(filePath, level, spc):
    """ Get the concentration data for an FST file.

    Returns the data in an array, the datakey as a key and the fileID of the FST file.
    """

    try:
       fileID = rmn.fstopenall(filePath,rmn.FST_RO)
       dataKey = rmn.fstinf(fileID,nomvar=spc,ip1=level)['key']
       dataRec = rmn.fstluk(dataKey)
       concData = dataRec['d']
       fileMeta = rmn.fstprm(dataKey)
       fileMeta['iunit']=fileID
       gridData = rmn.ezqkdef(fileMeta)
       llGridData = rmn.gdll(gridData)
       print ('File {} recorded'.format(filePath))
       rmn.fstcloseall(fileID)
       return concData,llGridData
    except TypeError:
       # log an error into the log file
       #logging.warning('nomvar {} and ip1 {} could not be found for file {}.'.format(fileID, spc, level, filePath))
       return None


def makeDir(path):
    """ Creates a directory if it doesn't exist.

    Existence errors will be ignored. All other errors will be raised.
    """

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


baseModel = sys.argv[1]
altModels = toLit(sys.argv[2])
fileNames = toLit(sys.argv[3])
params = toLit(sys.argv[4])
ip1 = params['ip1'] # should give an int
spc = params['spc'] # should give a list of species
headers = ['Analysis','Alternate Model', 'Base Model']

if len(altModels) != len(fileNames):
    print('Warning: Number of reference models and save fileNames does not match! \nPlease check your input arguments.')

defaultDir = os.getcwd()
basePath = os.path.join(defaultDir, baseModel)
savePath = os.path.join(defaultDir, 'BasicStatistics/')
print('The basePath and savePath are: {} and {}'.format(basePath,savePath))

statList = ['MIN','MAX','MEAN','MEDIAN','RANGE','VAR','STD', 'P2','P10','P25','P50','P75','P90','P98','maxLatLon']

rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)
# Create folders for each species
spcCounter = 0
for s in spc:
    saveSpcPath = os.path.join(savePath, str(ip1), spc[spcCounter])

    altCounter = 0
    for model in altModels:
        altPath =os.path.join(defaultDir,model)
        baseRunList = os.listdir(basePath)
        altRunList = os.listdir(altPath)
        altIndex = altCounter
        folderName = fileNames[altCounter]
        altCounter += 1

        runCounter = 0
        for runCounter in range(len(baseRunList)):
            run = baseRunList[runCounter]
            if run in altRunList:
                baseFilePath = os.path.join(basePath,run)
                altFilePath = os.path.join(altPath,run)
                baseData,baseGrid = getConc(filePath = baseFilePath, level=ip1, spc=spc[spcCounter])
                altData,altGrid = getConc(filePath = altFilePath, level = ip1, spc=spc[spcCounter])
                if (baseData != None) and (altData != None):
                    makeDir(saveSpcPath)
                    csvName = '{}_{}.csv'.format(run,folderName)
                    f = open(os.path.join(saveSpcPath, csvName), 'w')
                    f.write(','.join(map(str,headers)))
                    runCounter += 1

                    statVals = []
                    dataList=[altList,baseData]
                    for i in range(len(dataList)): #sequece of these has to be the same as the matching headers
                        minStat = np.nanmin(dataList[i])
                        maxStat = np.nanmax(dataList[i])
                        meanStat = np.mean(dataList[i])
                        medianStat = np.median(dataList[i])
                        rangeStat = np.ptp(dataList[i])
                        varStat = np.var(dataList[i])
                        stdStat = np.std(dataList[i])
                        p2Stat = np.percentile(dataList[i],2)
                        p10Stat = np.percentile(dataList[i],10)
                        p25Stat = np.percentile(dataList[i],25)
                        p50Stat = np.percentile(dataList[i],50)
                        p75Stat = np.percentile(dataList[i],75)
                        p90Stat = np.percentile(dataList[i],90)
                        p98Stat = np.percentile(dataList[i],98)
                        maxLoc = np.argmax(dataList[i])
                        if i == 0:
                            maxLon = np.ndarray.flatten(altGrid['lon'])[maxLoc]
                            maxLat = np.ndarray.flatten(altGrid['lat'])[maxLoc]
                            maxLL = '/'.join(map(str,[maxLon,maxLat]))
                            statRow = [minStat,maxStat,meanStat,medianStat,rangeStat,varStat,stdStat,p2Stat,p10Stat,p25Stat,p50Stat,p75Stat,p90Stat,p98Stat,maxLL]
                            statVals += [statRow]
                        else:
                            maxLon = np.ndarray.flatten(baseGrid['lon'])[maxLoc]
                            maxLat = np.ndarray.flatten(baseGrid['lat'])[maxLoc]
                            maxLL = [maxLon,maxLat]
                            maxLL = '/'.join(map(str,[maxLon,maxLat]))
                            statRow = [minStat,maxStat,meanStat,medianStat,rangeStat,varStat,stdStat,p2Stat,p10Stat,p25Stat,p50Stat,p75Stat,p90Stat,p98Stat, maxLL]
                            statVals += [statRow]

                    # statVals should now be array of length statList*samples
                    arrayStatVals = np.array(statVals)
                    swappedStatVals = np.swapaxes(arrayStatVals,0,1)
                    allStats = []
                    statCounter = 0
                    for statCounter in range(len(statList)):
                        statRow = [statList[statCounter]]
                        statRow += swappedStatVals[statCounter]
                        allStats += [statRow] #need to save each stat row as a list
                        statCounter += 1

                    for row in allStats:
                        f.write('\n' + ','.join(map(str,row)))

                    f.close()
                else:
                    print('Stats file was not created for {}, run: {}'.format(folderName,run))
    spcCounter += 1
