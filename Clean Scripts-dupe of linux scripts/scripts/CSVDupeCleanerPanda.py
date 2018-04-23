#!/usr/bin/python

# CSVDupeCleanerPanda uses the pandas module to clean hotspot csvs.
# Produces two new CSVs,
# one with the earliest hotspot location for the user provided date (clean version)
# second CSV contains all other hotspot data for repeated measurement locations (dupe version)
# Function takes in two user inputs,
# 1: the fileName (or filePath) relative to the current directory
# 2: the date user wants to look at in YYYYMMDD format
# Sample Command: python CSVDupeCleanerPanda.py in20170808.csv 20170808

import os
import pandas as pd
import time
import sys

# name of the file in the current directory (adjust filePath if needed)
fileName = sys.argv[1]
# date must be in the form YYYYMMDD
userDate = sys.argv[2]

currentDir = os.getcwd()
saveCleanPath = os.path.join(currentDir)
saveDupePath = os.path.join(currentDir)
filePath = os.path.join(currentDir, fileName)

# make directories if they don't exist
def makeDir(path):
    """ Creates a directory for the provided path if it does not exist. """

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def nameCSV(baseFile,appendText):
    """ Get the name for the new CSV. """

    base,ext = os.path.splitext(baseFile)
    newCSV = base + appendText + ext
    return newCSV

# main program!

# read the source data
refCSV = pd.read_csv(filePath)
# get the names of the csv headers with hardcoded locations
## (don't have to deal with extra spaces and should cope with removal of spaces)
csvHeaders = list(refCSV.columns.values)
latVal = csvHeaders[0]
lonVal = csvHeaders[1]
dateVal = csvHeaders[2]

# add formatted timestrings
cutCSV = refCSV.copy()
cutCSV['timeStr']= cutCSV[dateVal].map(lambda x: time.strftime('%Y%m%d', time.strptime(x, '%Y-%m-%d %H:%M:%S')))
# drop rows if not user date
cutCSV=cutCSV.drop(cutCSV[cutCSV['timeStr'] != userDate].index)
# drop the string column
cutCSV = cutCSV.drop('timeStr', axis=1)

# create column for time structures
sortCSV = cutCSV.copy()
sortCSV['timeStrct'] = sortCSV[dateVal].map(lambda x: time.strptime(x, '%Y-%m-%d %H:%M:%S'))
# sort sequentially
sortCSV = sortCSV.sort('timeStrct')
# drop the extra column
sortCSV = sortCSV.drop('timeStrct',axis=1)

# Create clean,dupe dataframe
cleanCSV = sortCSV.copy()
dupeCSV = sortCSV.copy()
# drop duplicates from clean, keep first
cleanCSV = cleanCSV.drop_duplicates([latVal,lonVal])
# get list of all duplicates other than first
dupeCSV['isDupe']= dupeCSV.duplicated([latVal,lonVal])
dupeCSV = dupeCSV.drop(dupeCSV[dupeCSV['isDupe'] != True].index)
dupeCSV = dupeCSV.drop(['isDupe'], axis =1)

# write to CSV
csvFile = os.path.basename(filePath)
cleanName = nameCSV(csvFile, '_clean')
dupeName = nameCSV(csvFile, '_dupe')
cleanCSV.to_csv(os.path.join(saveCleanPath,cleanName), index=False)
dupeCSV.to_csv(os.path.join(saveDupePath,dupeName),index=False)
