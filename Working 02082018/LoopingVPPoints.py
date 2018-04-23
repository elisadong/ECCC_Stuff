#!/usr/bin/python

# LoopingVPPoints.py

import sys
import os
import ast
import cleanVerticalProfilePoint as cVPP

# provide a list of base model and name, models, corresponding titles, area and areatype, spc(list)

def toLit(string):
    nowLit = ast.literal_eval(string)
    return nowLit

base = toLit(sys.argv[1])
baseModel = base[0]
baseTitle = base[1]
altModels = toLit(sys.argv[2])
plotTitles = toLit(sys.argv[3])
extraParams=['spc','area', 'type']

cali = [240.6,36.8]
defaultXPVals = {'spc':['TO3'], 'area':cali, 'type':'l'}
updatedXPVals = defaultXPVals
if len(sys.argv) > 4:
    otherParams = toLit(sys.argv[4])
    for p in range(len(extraParams)):
        pName = extraParams[p]
        if pName in otherParams:
            updatedXPVals[pName] = otherParams[pName]

x = updatedXPVals
spc, area, areaType = x['spc'], x['area'], x['type']

currentDir = os.getcwd()
savePath = os.path.join(currentDir, 'VerticalProfiles')
spcCounter = 0
for s in spc:
    saveSpcPath = os.path.join(savePath, spc[spcCounter])
    spcIndex = spcCounter
    spcCounter += 1

    modelCounter = 0
    for model in altModels:
        modelPath = os.path.join(defaultDir, model)
        paramsBase={'filePath':baseModel, 'spc':spc[spcIndex], 'area':area, 'type': areaType, 'title':baseTitle}
        print(paramsBase)
        paramsAlt={'filePath':model, 'spc':spc[spcIndex], 'area':area, 'type': areaType, 'title':plotTitles[modelCounter]}
        cVPP.plotTwo(paramsBase, paramsAlt, saveSpcPath)
        modelCounter += 1

print('Vertical Profiles have been generated. \nPlease check {} for results.'.format(savePath))
