# plume stats by distance
# rerun and update the appropriate pkl and saveDir directories

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
import datetime
import numpy as np
import speedyAltMINXdemo as sam
import loopSpeedMINX as sm

startTime = 0
plumePath ="/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/plumesP50C60/plumesCan/"
fireDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/CFFEPSLowResMant4165"
noFireDir = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/GEMMACHOPLowRes'
pixDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/fpFuelFilt/"

saveDirName = 'CFFEPSbyDistance'
curDir = os.getcwd()
saveDir = os.path.join(curDir, saveDirName)
plmHeights = []
plmDists = []
plmModHeights = []

pt = 95

pwrInd = 5
pixList = os.listdir(pixDir)
pwrs = []

def toFloats(row):
    floats = map(float, row.split())
    return floats

def readFile(path):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    return lines

tempList = os.listdir(plumePath)
plumeList= []
# make sure is only plumes, if maybePlume is plume, add to list
for ind, mPlume in enumerate(tempList):
    if mPlume.endswith('txt'):
        plumeList += [mPlume]

plumes = list(plumeList)

for t in pixList:
    tLines = readFile(os.path.join(pixDir,t))
    for rowInd, row in enumerate(tLines, 4):
        if rowInd < len(tLines):
            pwrs += [toFloats(tLines[rowInd])[pwrInd]]

pwrMax = np.percentile(pwrs, pt)

keptPlumes = {}
sm.makeDir(saveDir)
for pInd, plume in enumerate(plumes):
    plm = sam.pltPlume(plumePath, plume, minFRP=pwrMax, mod_fire=fireDir, mod_nofire=noFireDir, threshVal = 0,  startTime = startTime)
    try:
        if plm['frp'] >= minFRP:
            keptPlumes[plume] = {'plume':plm}
    except:
        continue

for kp in keptPlumes:
    kPlm = keptPlumes[kp]['plume']
    try:
        plmHeight, plmModHeight, plmDist = kPlm['filtered_height'], kPlm.Model.height, kPlm['distance']
        plmHeights += plmHeight
        plmModHeights += plmModHeight
        plmDists += plmDist

distRanges = [[0,10],[10,30],[30,60],[50,500]]
for dRange in distRanges:
    try:
        dPlmHeights, dPlmModHeights, dPlmDists = sm.filterPlmRange(dRange,plmHeights, plmModHeights, plmDists)
        sm.overallImages(dPlmDists, dPlmHeights, dPlmModHeights, threshVal, pt, saveDir, modifier='_dist{}-{}'.format(dRange[0], dRange[1]))
        sm.saveStatsLite(saveDir, 'N/A', dPlmHeights, dPlmModHeights, minFRP = pwrMax, threshVal=threshVal, modifier='_dist{}-{}'.format(dRange[0], dRange[1]))
    except:
        print('Did not save distance breakdown plot for dRange: {}'.format(dRange))
