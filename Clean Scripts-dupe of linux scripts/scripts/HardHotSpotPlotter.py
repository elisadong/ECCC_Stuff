#!/usr/bin/python

## maps csvs to SPI.obs files for CFFEPS csvs
# takes in three arguements from the command lines
# 1: the program name
# 2: The file name from the default directory
# 3: The date time in YYYYMMDD format
# 4: Optional, The range of est Area in hectares in the form of a list
# ex. python HardHotSpotPlotter.py nameof.csv "[3, 30]"

# example command line
# python HardHotSpotPlotter.py in20170912.csv 20170912 "[10,50]"

import os
import csv
import ast
import sys
import errno
import time

# Constants that can be changed
sourcePath = '/home/eld001/elisatest/Projects/HotSpotPlotter/test/'
savePath = '/home/eld001/elisatest/Projects/HotSpotPlotter/obsFiles/10500'
estAreaRange = [50,500]
# Note that LAT and LON are mandatory for the obs file format, so should not be removed
outputCol = ['lat','lon','estarea']

def toLit(string):
    """ Converts string to literal representation.

    Useful when taking in arguments from command line which are all passed as strings. In this program, converts to lists and dictionaries.
    """

    nowLit = ast.literal_eval(string)
    return nowLit

def OBS(fileName):
    """ Converts file type to .obs. """

    base,ext = os.path.splitext(fileName)
    newName = base + '.obs'
    return newName

def toDate(dateString):

    dateTimeStruct = time.strptime(dateString, '%Y-%m-%d %H:%M:%S')
    newDateString = time.strftime('%Y%m%d', dateTimeStruct)
    return newDateString

def OBSHeader(curHead):
    """ Gets a header that falls within the OBS format. """

    leading = ['ID', 'LAT', 'LON']
    newHeader = leading
    headCount = 0
    for head in range(len(curHead)):
        # Turns requested header into uppercase and removes leftside white space
        curHead[headCount] = ((curHead[headCount]).upper()).lstrip()
        # If the header is not one of the mandatory fields in leading
        # Strip the white space and add "DATA." to fit the obs format
        if curHead[headCount] not in leading:
            newHeader += ['DATA.'+ (curHead[headCount]).replace(' ','_') + '.' + userDate]
        headCount += 1
    return newHeader

def makeDir(path):
    """ Creates a directory for the provided path if it does not exist. """

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

fileName = sys.argv[1]
userDate = sys.argv[2]
# If user provides difference range,
if len(sys.argv) == 4:
    estAreaRange = toLit(sys.argv[3])

# Create the save directory
makeDir(savePath)
# Create the .obs file
newFile = OBS(fileName)
newOBS = open(os.path.join(savePath,newFile), 'w')
newOBS.write("Obs 3.1")
newOBS.write("\nThis is the OBS version of the {} for date {}.".format(fileName,userDate))

# Open and read the csv
refPath = os.path.join(sourcePath, fileName)
f =open(refPath)
reader = csv.reader(f)
headers = reader.next()
colKeys={}
for header in headers:
	colKeys[(header).lstrip()]=headers.index(header)

estAreaIndex = colKeys['estarea']
dateIndex = colKeys['rep_date']
colIndex = []
for colKey in outputCol:
    colIndex+=[colKeys[colKey]]

newHeader = OBSHeader(['ID'] + outputCol)
newOBS.write('\n'+" ".join(map(str,newHeader)))

refDate = time.strptime(userDate,'%Y%m%d')
rowCounter = 1 #start row counting in the forloop
for row in reader:
    rowCounter += 1
    try:
        row[estAreaIndex]
        if toDate(row[dateIndex]) == userDate:
            if float(row[estAreaIndex]) >= estAreaRange[0] and float(row[estAreaIndex]) <= estAreaRange[1]:
                newID = str(row[colKeys['lat']]) + str(row[colKeys['lon']])
                newRow = [newID.replace(' ','')]

                for c in colIndex:
                     # print(c)
                     newRow += [(row[c]).replace(' ','')] #row should be in list format
                newString = " ".join(map(str, newRow))
                newOBS.write('\n'+newString)
            else:
                print('Row {}: The estarea of this row is not within the range. estarea(hec): {}'.format(rowCounter,float(row[estAreaIndex])))
    except:
        pass

f.close()
newOBS.close()
print('The file has been converted to obs.')
