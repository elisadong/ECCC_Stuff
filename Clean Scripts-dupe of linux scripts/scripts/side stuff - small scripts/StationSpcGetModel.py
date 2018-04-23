#!/usr/bin/python

# for station data extraction and values using model data
# requires the following files/folders to work
## List of stations, IDs, Lon/Lat coordinates in a text file
## For model, a folder path to all model files for different run times

# Example:
# python StationSpcGetModel.py PM2.5.station_id_name_CFFEPS_lat_lon.txt "['AF']" FST_files/20170901120000/model CFFEPS20170901120000

import pandas as pd
import rpnpy.librmn.all as rmn
import rpnpy.rpndate as rdate
import os
import ast
import datetime
import errno
import sys
import math

def toLit(string):
    """ Converts string to literal representation.

    Useful when taking in arguments from command line which are all passed as strings. In this program, converts to lists and dictionaries.
    """

    nowLit = ast.literal_eval(string)
    return nowLit

# system provided values
stationFile = sys.argv[1]
spcs = toLit(sys.argv[2])
modelPath = sys.argv[3]
modelName = sys.argv[4]

# default values
ip1 = 76696048
currentDir = os.getcwd()
saveDir = os.path.join(currentDir, 'SpeciesByStation_ModelValues')
# index the columns from the station data
statIDColInd, statLatColInd, statLonColInd, statNameColInd = 0,1,2,3

# other helper functions
def getID(filePath):
    """ Get the ID of an FST file. Useful when looking at FST file contents. """

    fileID = rmn.fstopenall(filePath,rmn.FST_RO)
    return fileID

def getConc(fileID, spc):
    """ Get the concentration data for an FST file and a species. IP1 is a default value.

    The dataKey and dataRec is also returned for further use.
    """

    dataKey = rmn.fstinf(fileID, nomvar=spc, ip1=ip1)['key']
    dataRec = rmn.fstluk(dataKey)
    concData = dataRec['d']
    return concData, dataKey, dataRec

def getNiNj(lonVal, latVal, gridID):
    """ Get the ni,nj position in a grid for a given lon/lat value.

    Function returns None if lon/lat value is out of bounds.
    """

    xypos = rmn.gdxyfll(gridID, latVal, lonVal)
    ni = int(math.floor(xypos['x'])) - 1 #-1 to use 0 as list index
    nj = int(math.floor(xypos['y'])) - 1
    gridData = rmn.gdll(gridID)
    if ni >= len(gridData['lon']) or nj >= len(gridData['lon'][0]):
        print('lon lat values of station outside of grid domain. Please check your values.')
        return None
    else:
        return ni, nj

def getGrid(dataKey, fileID):
    """ Get the gridID for an FST file. Can be used to get more grid details. """

    fileMeta = rmn.fstprm(dataKey)
    fileMeta['iunit'] = fileID
    gridID = rmn.ezqkdef(fileMeta)
    return gridID

def getDate(dataRec):
    """ Returns the date time of the data retrieval in YYYY-MM-DDTHH:mm:ssZ format. """

    rpnTime = rdate.RPNDate(dataRec['datev'])
    timeStrct = rpnTime.toDateTime()
    timeStr = datetime.date.strftime(timeStrct,'%Y-%m-%dT%H:%M:%SZ')
    return timeStr

def makeDir(path):
    """ Creates a directory for the provided path if it does not exist. """

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

# get station information, assuming there is a header in the file
statInfo = pd.read_table(os.path.join(currentDir,stationFile), sep=',', header=True)
statHeaders = list(statInfo.columns.values)
# retrieve header names to later reference
statIDCol, statLatCol, statLonCol, statNameCol = statHeaders[statIDColInd], statHeaders[statLatColInd], statHeaders[statLonColInd], statHeaders[statNameColInd]

# set up empty species dictionary for each station
stationData = {}
statIDs = statInfo[statIDCol].tolist()
statLats = statInfo[statLatCol].tolist()
statLons = statInfo[statLonCol].tolist()
statNames = statInfo[statNameCol].tolist()

models = os.listdir(os.path.join(currentDir, modelPath))
models.sort()

for ID in statIDs:
    stationData[ID] = {}
    statInd = statIDs.index(ID)
    stationData[ID]['name'] = statNames[statInd]
    stationData[ID]['lon'] = statLons[statInd]
    stationData[ID]['lat'] = statLats[statInd]
    for spc in spcs:
        stationData[ID][spc]={}
        stationData[ID][spc]['modelVals'] = []
        stationData[ID][spc]['modelTime'] = []

rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)

for model in models:
    try:
        mID = getID(os.path.join(currentDir, modelPath, model))
        if mID == None:
            break
    except:
        pass
    # get all concentration values for a species on a layer
    for spc in spcs:
        spcConcData, mKey, mData = getConc(mID, spc)
        gridID = getGrid(mKey, mID)
        # get the spc val for each station ni nj location
        # add to dictionary that stores values
        for ID in statIDs:
            # need ni, nj value
            mNi, mNj = getNiNj(stationData[ID]['lon'], stationData[ID]['lat'], gridID)
            spcConcVal = spcConcData[mNi][mNj]
            stationData[ID][spc]['modelVals'] += [spcConcVal]
            stationData[ID][spc]['modelTime'] += [getDate(mData)]

for ID in statIDs:
    xList = []
    for s in spcs:
        xID = ID
        xSpc = s
        xLon = stationData[ID]['lon']
        xLat = stationData[ID]['lat']
        # skilped the modelVals in the test
        modelVals = stationData[ID][s]['modelVals']
        for mVal in modelVals:
            valInd = modelVals.index(mVal)
            xVal = mVal
            xModelName = stationData[ID][s]['modelTime'][valInd]
            xNote = 'Model Value'
            xRow = [xID, xLon, xLat, xVal, xModelName, xNote]
            xList += [xRow]
        xdf = pd.DataFrame(xList, columns=('StationID', 'lon','lat', s, 'UTC time','note'))
        xCSVName = '{}_{}_{}.csv'.format(ID, s, modelName)
        xPath = os.path.join(saveDir, xCSVName)
        makeDir(saveDir)
        xdf.to_csv(xPath, sep = ',', index = False)
