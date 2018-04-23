#!/usr/bin/python

## maps csvs to SPI.obs files for CFFEPS csvs
# takes in three arguements from the command lines
# 1: the program name
# 2: The file name from the default directory
# 3: Optional, The range of est Area in hectares in the form of a list
# ex. python HardHotSpotPlotter.py nameof.csv "[3, 30]"

import os
import csv
import ast
import sys
import errno

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
            newHeader += ['DATA.'+ (curHead[headCount]).replace(' ','_')]
        headCount += 1
    return newHeader

def toHec(squareM):
    """ Converts square meters to hectares. """

    hec=squareM/10000
    return hec

def makeDir(path):
    """ Creates a directory for the provided path if it does not exist. """

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
    else:
        raise

# get arguments from command line, else set defaults to 50 to 500 hectares
sourcePath = '/home/eld001/elisatest/Projects/HotSpotPlotter/'
savePath = '/home/eld001/elisatest/Projects/HotSpotPlotter/'
estAreaRange = [50,500]
outputCol = ['lat','lon','estarea']
fileName = sys.argv[1]
# If user provides difference range,
if len(sys.argv) == 3:
    estAreaRange = toLit(sys.argv[2])

# Create the save directory
makeDir(savePath)
# Create the .obs file
newFile = OBS(fileName)
newOBS = open(savePath+newFile, 'w')
newOBS.write("Obs 3.1")
newOBS.write("\nThis is the OBS version of the {}.".format(fileName))

# Open and read the csv
refPath = os.path.join(sourcePath, fileName)
f =open(refPath)
reader = csv.reader(f)
headers = reader.next()
colKeys={}
for header in headers:
	colKeys[(header).lstrip()]=headers.index(header)
# colKeys = {'lat':0, 'lon':1,'rep_date':2,'source':3,'sensor':4,'fwi':5,'fuel':6,'ros':7,'sfc':8,'tfc':9,'bfc':10,'hfi':11,'estarea':12,'area':13,'heat':14,'area fraction':15,'flaming fraction':16,'smoldering fraction':17,'residual fraction':18}
estAreaIndex = colKeys['estarea']
# print(colKeys)
colIndex = []
for colKey in outputCol:
    colIndex+=[colKeys[colKey]]

newHeader = OBSHeader(['ID'] + outputCol)
newOBS.write('\n'+" ".join(map(str,newHeader)))

for row in reader:
    try:
        row[estAreaIndex]
        if toHec(float(row[estAreaIndex])) >= estAreaRange[0] and toHec(float(row[estAreaIndex])) <= estAreaRange[1]:
            newID = str(row[colKeys['lat']]) + str(row[colKeys['lon']])
            newRow = [newID]

            for c in colIndex:
                 # print(c)
                 newRow += [row[c]] #row should be in list format
            newString = " ".join(map(str, newRow))
            newOBS.write('\n'+newString)
        else:
            print('The estarea of this row is not within the range. estarea(hec): {}'.format(toHec(float(row[estAreaIndex]))))
    except:
        pass
#print(colIndex)
f.close()
newOBS.close()
print('The file has been converted to obs.')
