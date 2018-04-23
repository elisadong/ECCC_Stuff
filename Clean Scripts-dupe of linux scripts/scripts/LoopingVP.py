#!/usr/bin/python

# LoopingVP.py

import sys
import os
import ast
import cleanVerticalProfile as cVP

# provide a list of models, corresponding titles, area and areatype, spc(list)

def toLit(string):
    nowLit = ast.literal_eval(string)
    return nowLit

base = toLit(sys.argv[1])
baseModel = base[0]
baseTitle = base[1]
altModels = toLit(sys.argv[2])
plotTitles = toLit(sys.argv[3])
extraParams=['spc','area', 'type']

norCal = {'LL':[238.4,36.7], 'UR':[238.8,41.8]}
BC = {'LL':[233.8, 50.2], 'UR':[239.2,58.8]}
defaultXPVals = {'spc':['TO3'], 'area':BC, 'type':'l'}
updatedXPVals = defaultXPVals
if len(sys.argv) > 4:
    otherParams = toLit(sys.argv[4])
    for p in range(len(extraParams)):
        pName = extraParams[p]
        if pName in otherParams:
            updatedXPVals[pName] = otherParams[pName]

x = updatedXPVals
spc, area, areaType = x['spc'], x['area'], x['type']

defaultDir = os.getcwd()
savePath = os.path.join(defaultDir, 'VerticalProfiles')
spcCounter = 0
for s in spc:
    saveSpcPath = os.path.join(savePath, spc[spcCounter])
    spcIndex = spcCounter
    spcCounter += 1

    modelCounter = 0
    for model in altModels:
        try:
            modelPath = os.path.join(defaultDir, model)
            paramsBase={'filePath':baseModel, 'spc':spc[spcIndex], 'area':area, 'type': areaType, 'title':baseTitle}
            paramsAlt={'filePath':model, 'spc':spc[spcIndex], 'area':area, 'type': areaType, 'title':plotTitles[modelCounter]}
            cVP.plotTwo(paramsBase, paramsAlt, saveSpcPath)
            modelCounter += 1
        except:
            modelCounter += 1
            print('WARNING: Vertical Profile was not generated for {}, {}.'.format(paramsBase, paramsAlt))

print('Vertical Profiles have been generated. \nPlease check {} for results.'.format(savePath))
