#!/usr/bin/python

# for station data extraction and comparison
# requires the following files/folders to work
## List of stations, IDs, Lon/Lat coordinates in a text file
## For each species, a CSV with station data, times and spc values

# Example:
# python StationSpcGetObs.py PM2.5.station_id_name_CFFEPS_lat_lon.txt "['AF']" "['pm25_obs_20170901_20170930_with_station_info.csv']"


import pandas as pd
import rpnpy.librmn.all as rmn
import rpnpy.rpndate as rdate
import os
import ast
import errno
import sys
import time

def toLit(string):
    """ Converts string to literal representation.

    Useful when taking in arguments from command line which are all passed as strings. In this program, converts to lists and dictionaries.
    """

    nowLit = ast.literal_eval(string)
    return nowLit

stationFile = sys.argv[1]
spcs = toLit(sys.argv[2])
spcsPaths = toLit(sys.argv[3])

# defaults
ip1 = 76696048
currentDir = os.getcwd()
saveDir = os.path.join(currentDir, 'SpeciesByStation_ObsVals')
# index the columns from the station data
statIDColInd, statLatColInd, statLonColInd, statNameColInd = 0,1,2,3
# index for relevant headers in the species concentration data
spcStatInd, spcValInd, spcDateInd = 0, 9, 13

def makeDir(path):
    """ Creates a directory for the provided path if it does not exist. """

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

# processing station info for all stations
statInfo = pd.read_table(os.path.join(currentDir,stationFile), sep=',')
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
        stationData[ID][spc]['obsVals'] = []
        stationData[ID][spc]['obsTime'] = []


for spc in spcs:
    spcInd = spcs.index(spc)
    spcData = pd.read_csv(os.path.join(currentDir, spcsPaths[spcInd]))
    # headers for individual species recorded data
    spcDataHeaders = list(spcData.columns.values)
    spcStat, spcVal, spcDate = spcDataHeaders[spcStatInd], spcDataHeaders[spcValInd], spcDataHeaders[spcDateInd]
    for ID in statIDs:
        cutCSV = spcData.drop(spcData[spcData[spcStat] != ID].index)
        cutCSV['timeStrct']=cutCSV[spcDate].map(lambda x: time.strptime(x, '%Y-%m-%dT%H:%M:%SZ'))
        sortCSV = cutCSV.copy()
        sortCSV = sortCSV.sort('timeStrct')
        spcValList = sortCSV[spcData.columns.values[spcValInd]].tolist()
        timeList = cutCSV[spcDate].tolist()
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
        obsVals = stationData[ID][s]['obsVals']
        for oVal in obsVals:
            valInd = obsVals.index(oVal)
            xVal = oVal
            xTime = stationData[ID][s]['obsTime'][valInd]
            xNote = 'Observed Value'
            xRow = [xID, xLon, xLat, xVal, xTime, xNote]
            xList += [xRow]
        xdf = pd.DataFrame(xList, columns=('StationID', 'lon','lat', s, 'UTC time','note'))
        xCSVName = '{}_{}_{}.csv'.format(ID, s, 'obsVals')
        xPath = os.path.join(saveDir, xCSVName)
        makeDir(saveDir)
        xdf.to_csv(xPath, sep = ',', index = False)
