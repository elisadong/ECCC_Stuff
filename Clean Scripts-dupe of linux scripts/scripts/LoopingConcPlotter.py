#!/usr/bin/python

## Uses the ConcentrationPlotter program to create concentration plots saved in user determined folders. ConcentrationPlotter.py should be in the same folder this program is run from
## For each FST file found, a concentration plot is created
## Max number of arguments provided should be 4, program name, source models, save destination and modifying parameters
# 1: program name,
## ie. LoopingConcPlotter.py
# 2: list of paths to folders with files different models, relative from current directory
# Please note that if the list has more than one item, you must add additional quotation marks on either end of the list
## ex. ['/FST_files/adomKPPKPPB_adomYB/folderName']
## ex. "['/FST_files/adomKPPKPPB_adomYB/folderName', /FST_files/funModel/folderName]"
# 3: name of the folders where the concentration plots will be stored, one for each comparison model
# Note that these will be stored by default in CurrentDirectory/ConcentrationPlots
## ex. ['adomYBConcPlots']
## ex. "['adomYBConcPlots', 'funPlots']"
# 4: Optional. Include a dictionary of extra parameters with selected keys from:
# 'ip1', 'spc', 'buff', 'cmap', 'rev', 'extScale', 'projection', 'minVal', 'maxVal', 'bins', 'clear'
# Note that the dictionary must be passed with quotation marks on either side
# Passing some parameters may have no effect (ex. setting rev=True on cmap=spectral)
## ex. "{'bin': 50, 'rev': True, 'cmap': 'RdBu', 'extScale':'min'}"

# example command line
# python LoopingConcPlotter.py "['/FST_files/adomKPPKPPB_adomYB/folderName', /FST_files/funModel/folderName]" "['adomYBConcPlots', 'funPlots']" "{'bin': 50, 'rev': True, 'cmap': 'RdBu', 'extScale':'min'}"

import sys
import ast
import ConcentrationPlotter as cplt
import os

def toLit(string):
    """ Converts string to literal representation.

    Useful when taking in arguments from command line which are all passed as strings. In this program, converts to lists and dictionaries.
    """

    nowLit = ast.literal_eval(string)
    return nowLit


modelList = toLit(sys.argv[1])
saveTo = toLit(sys.argv[2])
if len(modelList) != len(saveTo):
    print('Warning: Number of reference model folders and save to folders does not match! \nPlease check your input arguments.')

defaultDir = os.getcwd()
savePath = os.path.join(defaultDir, 'ConcentrationPlots/')
print('The savePath is: {}'.format(savePath))

extraParams = ['ip1', 'spc', 'buff', 'cmap', 'rev', 'extScale', 'projection', 'minVal', 'maxVal', 'bins', 'clear']

# dictionary of default values
defaultXPVals = {'ip1':76696048, 'spc':['TO3'], 'buff':2, 'cmap':'spectral', 'rev':False,'extScale':'neither', 'projection':'stere', 'minVal':None, 'maxVal':None, 'bins':None, 'clear':0}

# Set the new parameter values to be the same as the current default values.
# This will be updated if the modifiable parameters were passed.
updatedXPVals = defaultXPVals

i = 3 # the modifiers start at the fourth argument
if len(sys.argv) > 3:
    otherParams = toLit(sys.argv[3])
    for p in range(len(extraParams)):
        # for each param in extraParams
        pName = extraParams[p]
        # check to see if user provided an alternate value for the param
        if pName in otherParams:
            updatedXPVals[pName] = otherParams[pName]
    # Uncomment the following lines to get a list of extra parameters passed that were not accepted
    #notParam = []
    #for o in range(len(otherParams)):
    #    if otherParams[p] not in extraParams:
    #        notParam += otherParams[p]
    #print('The following parameters were not added: ' + str(notParam))


print('These are the updated parameters: ',updatedXPVals)

# assign the new values
x = updatedXPVals
ip1, spc, buff, cmap, rev, extScale,projection, minVal,maxVal,bins,clear = (x['ip1'], x['spc'], x['buff'], x['cmap'], x['rev'], x['extScale'], x['projection'], x['minVal'], x['maxVal'], x['bins'], x['clear'])

# Create folders for each species
spcCounter = 0
for s in spc:
    saveSpcPath = os.path.join(savePath, str(ip1), spc[spcCounter])
    spcIndex = spcCounter
    spcCounter += 1

    modelCounter = 0
    for model in modelList:
       modelPath = os.path.join(defaultDir,model)
       modelRunList = os.listdir(modelPath)
       saveCplts = saveSpcPath
       partName = saveTo[modelCounter]
       modelCounter += 1

       runCounter = 0
       for runCounter in range(len(modelRunList)):
          run = modelRunList[runCounter]
          print('The run is: ' + run)
          # Uncomment the following two lines to only create two images per folder
          #if runCounter == 2:
            # print('Warning: Created a max of two plots for this model.')
             #break
          if run in modelRunList:
             modelFilePath = os.path.join(modelPath, run)
             cplt.concPlot(fileRun=run, modelFile=modelFilePath, savePath=saveCplts, level=ip1, species=spc[spcIndex], buffering=buff,
                   mapType=cmap, reverse=rev, extension=extScale, projType=projection, vmin=minVal, vmax=maxVal, totalBins = bins, partName=partName, removed=clear)
             runCounter += 1
          else:
             print('Warning: Unknown error for {} in {}.'.format(modelRunList[runCounter], modelPath))
             runCounter += 1
print('Concentration plots have been generated. \nPlease check {} for results.'.format(savePath))
