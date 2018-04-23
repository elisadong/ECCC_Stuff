# for station data extraction and comparison
# requires the following files/folders to work
## List of stations, IDs, Lon/Lat coordinates in a text file
## For each species, a CSV with station data, times and spc values
## For each model, a folder path to all model files for different run times

import pandas as pd
import rpnpy.librmn.all as rmn
import os
import ast
import time
import errno
import sys


def toLit(string):
    """ Converts string to literal representation.

    Useful when taking in arguments from command line which are all passed as strings. In this program, converts to lists and dictionaries.
    """

    nowLit = ast.literal_eval(string)
    return nowLit

stationFile = sys.argv[1]
spcs = toLit(sys.argv[2])
spcsPaths = toLit(sys.argv[3])
spcHeaderInd = 9
modelPath = sys.argv[4]
modelName = sys.argv[5]
ip1 = 76696048
currentDir = os.getcwd()
saveDir = os.path.join(currentDir, 'SpeciesByStation')
# index the columns from the station data
statIDColInd, statLatColInd, statLonColInd, statNameColInd = 0,1,2,3

rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)

def getID(filePath):
    fileID = rmn.fstopenall(filePath,rmn.FST_RO)
    return fileID

def getConc(fileID, spc):
    dataKey = rmn.fstinf(fileID, nomvar=spc, ip1=ip1)['key']
    dataRec = rmn.fstluk(dataKey)
    concData = dataRec['d']
    return concData, dataKey

def getNiNj(lonVal, latVal, gridID):
    xypos = rmn.gdxyfll(gridID, latVal, lonVal)
    ni = int(round(xypos['x'])) - 1 #-1 to use 0 as list index
    nj = int(round(xypos['y'])) - 1
    gridData = rmn.gdll(gridID)
    if ni >= len(gridData['lon']) or nj >= len(gridData['lon'][0]):
        print('lon lat values of station outside of grid domain. Please check your values.')
        return None
    else:
        return ni, nj

def getGrid(dataKey, fileID):
    fileMeta = rmn.fstprm(dataKey)
    fileMeta['iunit'] = fileID
    gridID = rmn.ezqkdef(fileMeta)
    return gridID

def makeDir(path):
    """ Creates a directory for the provided path if it does not exist. """

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

statInfo = pd.read_table(os.path.join(currentDir,stationFile), sep=',', header=True)
statHeaders = list(statInfo.columns.values)
statIDCol, statLatCol, statLonCol, statNameCol = statHeaders[statIDColInd], statHeaders[statLatColInd], statHeaders[statLonColInd], statHeaders[statNameColInd]

# set up empty species dictionary for each station
stationData = {}
statIDs = statInfo[statIDCol].tolist()
statLats = statInfo[statLatCol].tolist()
statLons = statInfo[statLonCol].tolist()
statNames = statInfo[statNameCol].tolist()

for ID in statIDs:
    stationData[ID] = {}
    statInd = statIDs.index(ID)
    stationData[ID]['name'] = statNames[statInd]
    stationData[ID]['lon'] = statLons[statInd]
    stationData[ID]['lat'] = statLats[statInd]
    for spc in spcs:
        stationData[ID][spc]={}
        stationData[ID][spc]['modelVals'] = []
        stationData[ID][spc]['modelNames'] = []
        stationData[ID][spc]['obsVals'] = []
        stationData[ID][spc]['obsTime'] = []
    # add in name, time, lon, lat, timezone

models = os.listdir(os.path.join(currentDir, modelPath))
models.sort()
# process all model related things

### This section needs testing
## mID get isn't working, can't confirm entire forloop
for model in models:
    try:
        mID = getID(os.path.join(currentDir, modelPath, model))
        if mID == None:
            break
    except:
        pass
    # get all concentration values for a species on a layer
    for spc in spcs:
        spcConcData, mKey = getConc(mID, spc)
        gridID = getGrid(mKey, mID)
        # get the spc val for each station ni nj location
        # add to dictionary that stores values
        for ID in statIDs:
            # need ni, nj value
            mNi, mNj = getNiNj(stationData[ID]['lon'], stationData[ID]['lat'], gridID)
            spcConcVal = spcConcData[mNi][mNj]
            stationData[ID][spc]['modelVals'] += [spcConcVal]
            stationData[ID][spc]['modelNames'] += [model]

# add observed data for all of the species
### Tested on gpsc
for spc in spcs:
    spcInd = spcs.index(spc)
    spcData = pd.read_csv(os.path.join(currentDir, spcsPaths[spcInd]))
    for ID in statIDs:
        cutCSV = spcData.drop(spcData[spcData['stationid'] != ID].index)
        cutCSV['timeStrct']=cutCSV['datetime'].map(lambda x: time.strptime(x, '%Y-%m-%dT%H:%M:%SZ'))
        sortCSV = cutCSV.copy()
        sortCSV = sortCSV.sort('timeStrct')
        spcValList = sortCSV[spcData.columns.values[spcHeaderInd]].tolist()
        timeList = cutCSV['datetime'].tolist()
        for t in timeList:
            stationData[ID][spc]['obsVals'] += [spcValList[timeList.index(t)]]
            stationData[ID][spc]['obsTime'] += [t]


for ID in statIDs:
    xList = []
    for s in spcs:
        xID = ID
        xSpc = s
        xLon = stationData[ID]['lon']
        xLat = stationData[ID]['lat']
        # skilped the modelVals in the test
        modelVals = stationData[ID][s]['modelVals']
        obsVals = stationData[ID][s]['obsVals']
        for mVal in modelVals:
            valInd = modelVals.index(mVal)
            xVal = mVal
            xModelName = stationData[ID][s]['modelNames'][valInd]
            xNote = 'Model Value'
            xRow = [xID, xLon, xLat, xVal, xModelName, xNote]
            xList += [xRow]
        for oVal in obsVals:
            valInd = obsVals.index(oVal)
            xVal = oVal
            xTime = stationData[ID][s]['obsTime'][valInd]
            xNote = 'Observed Value'
            xRow = [xID, xLon, xLat, xVal, xTime, xNote]
            xList += [xRow]
        xdf = pd.DataFrame(xList, columns=('StationID', 'lon','lat', s, 'UTC time','note'))
        xCSVName = '{}_{}_{}.csv'.format(ID, s, modelName)
        xPath = os.path.join(saveDir, xCSVName)
        makeDir(saveDir)
        xdf.to_csv(xPath, sep = ',', index = False)
