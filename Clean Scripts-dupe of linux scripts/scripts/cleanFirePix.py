# take top 10% of all fire pixels, and also filter by confidence level
# run this function from inside a folder with all the text files you want to clean

import os
import numpy as np

pwrInd = 5

def toFloats(row):
    floats = map(float, row.split())
    return floats

def readFile(path):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    return lines

def writeFile(path, lines, writeInds):
    f = open(path, 'w')
    for rowInd, string in enumerate(lines):
        if rowInd in writeInds:
            f.write(lines[rowInd])
    f.close()


# get top 10% in all files
cDir = os.getcwd()
txtList = os.listdir(cDir)
pwrs = []
for txt in txtList:
    tLines = readFile(txt)
    for rowInd, row in enumerate(tLines, 4):
        if rowInd < len(tLines):
            pwrs += [toFloats(tLines[rowInd])[pwrInd]]

pwrMax = np.percentile(pwrs, 90)

# write only the top 10%
# can also add in extra conditions like, confidence >= 60

for txt in txtList:
    writeInds = [0,1,2,3]
    tLines = readFile(txt)
    for rowInd, row in enumerate(tLines, 4):
        if rowInd < len(tLines):
            rowPwr = toFloats(tLines[rowInd])[pwrInd]
            if rowPwr >= pwrMax:
                writeInds += [rowInd]
    writeFile(txt, tLines, writeInds)
