#!/usr/bin/python

## maps csvs to SPI.obs files

import os
import csv
import ast
import sys

def toLit(string):
    """ Converts string to literal representation.

    Useful when taking in arguments from command line which are all passed as strings. In this program, converts to lists and dictionaries.
    """

    nowLit = ast.literal_eval(string)
    return nowLit

def OBS(fileName):
    base,ext = os.path.splitext(fileName)
    newName = base + '.obs'
    return newName

def OBSHeader(curHead):
    leading = ['ID', 'LAT', 'LON']
    newHeader = leading
    headCount = 0
    for head in range(len(curHead)):
        curHead[headCount] = ((curHead[headCount]).upper()).lstrip()
        if curHead[headCount] not in leading:
            newHeader += ['DATA.'+ (curHead[headCount])]
        headCount += 1
    return newHeader

def toHec(squareM):
   hec=squareM/10000
   return hec


sourcePath = sys.argv[1]
estAreaRange = [50,500]
if len(sys.argv) == 3:
    estAreaRange = toList(sys.argv[2])

# constants
currentDir = os.getcwd()
savePath = os.path.join(currentDir, 'obsFiles/')
print 'The savePath is: {}'.format(savePath)

sourceList = os.listdir(sourcePath)
sourceIndex = 0
for csvfile in sourceList:
    sourceFile = sourcePath[0]
    newFile = OBS(sourceFile)
    newOBS = open(savePath+newFile, 'w')
    newOBS.write("Obs 3.1")
    newOBS.write("\nThis is the OBS version of the {}.".format(sourceFile))

    f =open(sourceFile)
    headReader = csv.reader(f)
    headers = headReader.next()
    estAreaHead = headers.index(' estarea')
    f.close()

    tempHeader = ['ID'] + headers
    newHeader = OBSHeader(tempHeader)
    newOBS.write('\n'+" ".join(map(str,newHeader)))

    f = open(defaultFile)
    reader = csv.reader(f)
    reader.next()
    for row in reader:
       if hec(float(row[estAreaHead])) >= estAreaRange[0] and hec(float(row[estAreaHead])) <= estAreaRange[1]:
          newID = str(row[0]) + str(row[1])
          newRow = [newID] + row #row should be in list format
          newString = " ".join(map(str, newRow))
          newOBS.write('\n'+newString)
       else:
          print('The estarea of this row is not within the range. estarea(hec): {}'.format(hec(row[estAreaHead])))
    f.close()
    newOBS.close()
print('The file has been converted to obs.')
