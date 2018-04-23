#!/usr/bin/python

import os
import ast
import sys
import rpnpy.librmn.all as rmn
import numpy

# import file location, lon lat, and if want exact values

def toLit(string):
    """ Converts string to literal representation.

    Useful when taking in arguments from command line which are all passed as strings. In this program, converts to lists and dictionaries.
    """

    nowLit = ast.literal_eval(string)
    return nowLit

def getData(filePath, spc, level):
    try:
        fileID = rmn.fstopenall(filePath,rmn.FST_RO)
        print('Opened file!')
        dataKey = rmn.fstinf(fileID,nomvar=spc,ip1=level)['key']
        dataRec = rmn.fstluk(dataKey)
        concData = dataRec['d']
        return concData, dataKey, fileID
    except:
        print('ERROR: Could not find data for {} and {}. Please try again.'.format(spc, level))

def getGrid(dataKey, fileID):
    try:
        fileMeta = rmn.fstprm(dataKey)
        fileMeta['iunit'] = fileID
        gridData = rmn.ezqkdef(fileMeta)
        llGridData = rmn.gdll(gridData)
        latData = llGridData['lat'].tolist()
        lonData = llGridData['lon'].tolist()
        return lonData, latData, gridData
    except:
        print('ERROR: Could not get grid. Please try again.')

def getFuzzyNiNj(lonVal, latVal, gridData):
    xypos = rmn.gdxyfll(gridData, latVal, lonVal)
    ni = int(round(xypos['x'])) - 1 #-1 to use 0 as list index
    nj = int(round(xypos['y'])) - 1
    print ('Closest lonlat values: ({}, {})'.format(lonData[ni][nj], latData[ni][nj]))
    return ni, nj

def getConcentration(ni,nj,concData):
    concVal = concData[ni][nj]
    print('The concentration is: {}'.format(concVal))

#values working with
currentDir = os.getcwd()
filePath = os.path.join(currentDir,sys.argv[1])
LLCoords = toLit(sys.argv[2])
lonVal, latVal = LLCoords[0], LLCoords[1]
spc = sys.argv[3]
level = int(sys.argv[4])

rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)
concData, dataKey, fileID = getData(filePath, spc, level)
lonData, latData, gridData = getGrid(dataKey, fileID)
ni, nj = getFuzzyNiNj(lonVal, latVal, gridData)
print('The ni,nj values are: [{},{}]'.format(ni,nj))
getConcentration(ni, nj, concData)

#TODO:
# include a shebang!
# setup refault values
# convert optional params to dictionary? Or leave as mandatory statements
