#!/usr/bin/python

## loops through all the models that are to be processed
## for each matching file with the base model, a difference plot is created
## plots are saved in a folder provided by the user
## Note that the ConcentrationPlotter.py should be in the same folder as this program
# Parameters are as follows
# 1: program name,
## ie. LoopingDiffPlotter.py
# 2: path to the folder with the files produced by the reference model, relative from current directory,
## ex. '/FST_files/GENMACH/folderName'
# 3: list of paths to folders with files produced by the comparison model, relative from current directory
# Please note that if the list has more than one item, you must add additional quotation marks on either end of the list
## ex. ['/FST_files/adomKPPKPPB_adomYB/folderName']
## ex. "['/FST_files/adomKPPKPPB_adomYB/folderName', /FST_files/funModel/folderName]"
# 4: name of the folders where the comparisons will be stored, one for each comparison model
# Note that these will be stored by default in CurrentDirectory/DifferencePlots

# example command lines
# python LoopingDiffPlotter.py '/here/fancyModel' "['/home/somewhere/GENMACHmodels','/here/NorthAmericamodels']" "['GENMACH', 'NorthAm']"
# python LoopingConcPlotter.py 'FST_files/fancyModel' "['/FST_files/adomKPPKPPB_adomYB/folderName', /FST_files/funModel/folderName]" "['adomYBConcPlots', 'funPlots']" "{'bin': 50, 'rev': True, 'cmap': 'RdBu', 'extScale':'min'}"

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


baseModel = sys.argv[1]
modelList = toLit(sys.argv[2])
saveTo = toLit(sys.argv[3])
if len(modelList) != len(saveTo):
    print('Warning: Number of reference model folders and save to folders does not match! \nPlease check your input arguments.')

defaultDir = os.getcwd()
basePath = os.path.join(defaultDir, baseModel)
savePath = os.path.join(defaultDir, 'DifferencePlots/')
print('The basePath and savePath are: {} and {}'.format(basePath,savePath))

extraParams = ['ip1', 'spc', 'buff', 'cmap', 'rev', 'extScale', 'projection', 'minVal', 'maxVal', 'bins', 'clear']

# dictionary of default values
defaultXPVals = {'ip1':76696048, 'spc':['TO3'], 'buff':2, 'cmap':'RdBu', 'rev':False,'extScale':'neither', 'projection':'stere', 'minVal':None, 'maxVal':None, 'bins':None, 'clear':0}

# Set the new parameter values to be the same as the current default values.
# This will be updated if the modifiable parameters were passed.
updatedXPVals = defaultXPVals

i = 4 # the modifiers start at the fifth argument
if len(sys.argv) > 4:
    otherParams = toLit(sys.argv[4])
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

# assign the values
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
       baseRunList = os.listdir(basePath)
       modelRunList = os.listdir(modelPath)
       saveCplts =saveSpcPath
       modelIndex = modelCounter
       partName = saveTo[modelCounter]
       modelCounter += 1

       runCounter = 0
       for runCounter in range(len(baseRunList)):
          run = baseRunList[runCounter]
          print('The run is: ' + run)
          # Uncomment the following two lines to only create two images per folder
          #if runCounter == 2:
            # print('Warning: Created a max of two plots for this model.')
             #break
          if run in modelRunList:
             baseFilePath = os.path.join(basePath,run)
             modelFilePath = os.path.join(modelPath, run)
             cplt.diffPlot(fileRun=run, baseFile=baseFilePath, modelFile=modelFilePath, savePath=saveCplts, level=ip1, species=spc[spcIndex], buffering=buff,
                   mapType=cmap, reverse=rev, extension=extScale, projType=projection, vmin=minVal,vmax=maxVal, totalBins = bins, partName=partName, removed=clear)
             runCounter += 1
          else:
             print('Warning: Difference plot for {} was not found in {}. Please check to see if the file is missing.'.format(baseRunList[runCounter], modelPath))
             runCounter += 1
print('Difference plots have been generated. \nPlease check {} for results.'.format(savePath))
