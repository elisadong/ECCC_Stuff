# clean csv
# save for the earliest time frame
# save the duplicates into dupe file
# save all of the components in memory, only write at the end

import os
import csv
import ast
import sys
import errno
import time

# hardcode in the source file path(note that it should be a folder containing the CSVs to process)
# also hardcode in the save file path(should be folders)
filePath = sys.argv[1]
userDate = sys.argv[2]
currentDir = os.getcwd()
originalPath = os.path.join(currentDir, filePath)
saveCleanPath = os.path.join(currentDir)
saveDupePath = os.path.join(currentDir)

# hardcoded locations of the matching indexes
latIndex = 0
lonIndex = 1
repDateIndex = 2

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

def inStorage(lat,lon,storage):
    """ Check to see if the lat lon is in the existing values in the 'storage'. """

    result = []
    for x in storage:
        if lat != x[0] and lon!=x[1]:
            result += [False]
        elif lat == x[0] and lon == x[1]:
            result += [storage.index(x)]
        else:
            result += [False]
    for ans in result:
        if ans != False:
            return ans
    return False

# get all rows for clean and dupe storage
def getLines(readerRows):
    """ Store all of the contents for the clean and dupe CSVs. """

    # initialize clean and dupe lists
    # start reading
    cleanStorage = []
    dupeStorage = []
    print('Getting lines for CSV.')
    for row in readerRows:
        rowDate = time.strptime(row[repDateIndex], '%Y-%m-%d %H:%M:%S')
        rowDateStr = time.strftime('%Y%m%d', rowDate)
        if rowDateStr != userDate:
            continue
        else:
            if cleanStorage == []:
                cleanStorage += [row]
            else:
                latVal = row[latIndex]
                lonVal = row[lonIndex]
                if inStorage(latVal, lonVal, cleanStorage) == False:
                    cleanStorage += [row]
                else:
                    dupeRow =inStorage(latVal,lonVal,cleanStorage)
                    cleanDateVal = time.strptime(cleanStorage[dupeRow][repDateIndex], '%Y-%m-%d %H:%M:%S')
                    if rowDate > cleanDateVal:
                        dupeStorage += [row]
                    else:
                        dupeStorage += [cleanStorage[dupeRow]]
                        cleanStorage[dupeRow]=row


    return cleanStorage, dupeStorage


# This is the main function
makeDir(saveCleanPath)
makeDir(saveDupePath)

# get the names of the save to csvs
csvFile = os.path.basename(originalPath)
print('Getting info on: {}'.format(csvFile))
cleanCSVName = nameCSV(csvFile,'_clean')
dupeCSVName = nameCSV(csvFile,'_dupes')
cleanCSV = open(os.path.join(saveCleanPath,cleanCSVName), 'w')
dupeCSV = open(os.path.join(saveDupePath,dupeCSVName), 'w')
#get header
refCSV = open(os.path.join(originalPath))
reader = csv.reader(refCSV)
headers = reader.next()
cleanCSV.write(','.join(map(str,headers)))
dupeCSV.write(','.join(map(str,headers)))
readerRows = []
for row in reader:
    readerRows += [row]
refCSV.close()
# get the clean and dupe info
cleanRows, dupeRows = getLines(readerRows)
# write rows to csvs
for cRow in cleanRows:
    cleanCSV.write('\n' + ','.join(map(str,cRow)))
for dRow in dupeRows:
    dupeCSV.write('\n' + ','.join(map(str,dRow)))

# end the files
cleanCSV.write('\n')
dupeCSV.write('\n')

cleanCSV.close()
dupeCSV.close()
