#!/usr/bin/python

# WARNING:
# the author recommends that the user make a duplicate of their desired FST file to be processed as the changes are not stored
# the original file is not recoverable if the execution of this program is successful  (I think)

# This program writes AF or AC data to an FST file as run from command line
# The function  will try to create such levels from the sum of the following species
# AF = TSU1 + TOC1 + TNI1 + TAM1 + TCM1 + TEC1 + TPC1
# AC = TSU2 + TOC2 + TNI2 + TAM2 + TCM2 + TEC2 + TPC2
# Value check:
# If checkON = True, the program will check to see if an existing PM level data matches the generated data
# if the data does not match (likely due to different grid sizes), the ip1Check level will not have it's data updated

# Notes:
# **fstecr** function writes new data and metadata to a new or existing FST file
# to ensure tha the FST file is functional and not missing other data, this function writes to an existing FST file

# To read more on fstecr: https://wiki.cmc.ec.gc.ca/wiki/Python-RPN/2.0/rpnpy/librmn/fstd98#fstecr
# **datyp** in metadata is explicitly removed as per suggestion to not force the data type when using the fstecr function
# To read more on datyp: https://wiki.cmc.ec.gc.ca/wiki/Python-RPN/2.1/rpnpy/librmn/fstd98#dtype_fst2numpy

# Examples:
# python addPM.py AF fs/cetus/fs1/aqd/pxarqdr6/arqj/ahl/FST_files/GEMMACH_dev/20100710000000/model/2017091212_007

import os
import rpnpy.librmn.all as rmn
import sys
import numpy as np

# defaults, NOTE: try not to fiddle with the spacing in the defaults section, comments become uncommented in my text editor when that happens.

ip1List = [95528103, 97341399, 95360636, 94517525, 94506747, 95203234, 95177672, 97164386, 97013298, 95348903, 94790638, 94765118, 94740623, 94717051, 96320038, 96279484, 95879133, 95834648, 96886973, 96783497, 95335517, 95151143, 95122572, 94497415, 94489330, 96410389, 96363016, 95092030, 95060058, 95791222, 95750107, 94544377, 94530094, 94650523, 94630473, 95577949, 95551544, 95369037, 96154478, 96109937, 96583883, 98293613, 95304400, 94902189, 94873167, 96064236, 96017812, 94962914, 94932068, 95248115, 95226733, 95639424, 95607291, 94576580, 94560020, 96699938, 96633977, 95320496, 98024801, 97840480, 95287241, 95971125, 95924598, 94694044, 94671769, 97720038, 97644070, 95268471, 94611577, 94593675, 95711238, 95674198, 94844858, 94817232, 95027249, 94994764, 94482342, 94476231, 96239545, 96197821, 76696048]

# ipCheck assumes that the original FST file has AF and AC data at the 1.5m level, can disable the check if otherwise
ip1Check = 76696048
AFspcList = ['TSU1', 'TOC1', 'TNI1', 'TAM1', 'TCM1', 'TEC1', 'TPC1']
ACspcList = ['TSU2', 'TOC2', 'TNI2', 'TAM2', 'TCM2', 'TEC2', 'TPC2']
checkON = False

currentDir = os.getcwd()

# system provided values
pmType = sys.argv[1]
filePath = os.path.join(currentDir, sys.argv[2])
print(filePath)

def getID(filePath):
    """ Gets the fileID and opens the file in read/write mode. """

    rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)
    # open file to process
    fileID = rmn.fstopenall(filePath, rmn.FST_RW)
    return fileID

def getConc(fileID, spc, ip1):
    """ Gets the concentratio and meta data for a fileID, spc and ip1 if it exists. Else prints a warning.

    meta data keys 'd' and 'datyp' are removed for this program. To retain original contents, remove the 'del' commands.
    """

    try:
        dataKey = rmn.fstinf(fileID, nomvar=spc, ip1=ip1)['key']
        dataRec = rmn.fstluk(dataKey)
        concData = dataRec['d']
        meta = dataRec
        del meta['d']
        del meta['datyp']
        return concData, meta
    except:
        print('WARNING: Data does not exist for ip1 {} and spc {}.'.format(ip1, spc))

def sumArrays(dictionary, key1List, key2):
    """ Updates the dictionary passed to this function to read into FST file.

    Format is specific to this program, but the assumption is that key1List is a list of ip1 keys and key2 is the values stored as a list to be summed up. The returned dictionary has the key2 values as a sum of the arrays instead of a list of arrays.
    """

    for key in key1List:
        sumArray = sum(dictionary[key][key2])
        dictionary[key][key2] = sumArray
    return dictionary

def createAFCDict(fileID):
    """ Creates a dictionary with level 1 key: ip1 from ip1List, level 2 keys: 'vals' and 'meta'.

    key 'meta' contains meta data (will save the last species in the spcList meta data) with 'd' and 'datyp' removed
    key 'vals' contains the sum of all concentration data retrieved from spc in spcList using the sumArrays function.
    """

    AFCDict = {}
    if pmType == 'AF':
        spcList = AFspcList
    elif pmType == 'AC':
        spcList = ACspcList
    for ip1 in ip1List:
        AFCDict[ip1] = {}
        AFCDict[ip1]['vals'] = []
        for spc in spcList:
            #print('Getting AFC Data and adding to dictionary...')
            AFCData, AFCMeta = getConc(fileID, spc, ip1)
            try:
                AFCMeta['nomvar'] = pmType
                AFCDict[ip1]['meta'] = AFCMeta
                AFCDict[ip1]['vals'] += [AFCData]
                #print('Adding AFCData for ip1: {}, spc: {}.'.format(ip1, spc))
            except:
                print('Did not add AFCData for ip1: {}, spc: {}.'.format(ip1, spc))
    sumAFCDict = sumArrays(AFCDict, list(AFCDict.keys()), 'vals')
    return sumAFCDict, spcList

def checkPMlvl(AFCDict, ip1Check, pmType):
    """ Checks to see if the original PM data at ip1 level = ip1Check matches the generated PM data.

    This function will not be used if checkON = False.
    Note that this function is still undergoing development as the source data is of different grid size than the generated data.
    """
    AFCVals = AFCDict[ip1Check]['vals']
    pmTypeVals, pmMeta = getConc(fileID, pmType, ip1Check)
    if AFCVals != pmTypeVals:
        print('Generated PM values do not match.')
        return False
    else:
        print('PM values match for ip1: {}.'.format(ip1Check))
        return True

def updateFST(fileID, AFCDict):
    """ Updates the original FST file with the new PM data and corresponding metadata. """

    for ip1 in AFCDict:
        rmn.fstecr(fileID, AFCDict[ip1]['vals'], AFCDict[ip1]['meta'])

def closeFST(fileID):
    """ Closes the FST to prevent file corruption. """

    rmn.fstcloseall(fileID)

def inDict(dictionary, string):
    """ Checks to see if there are any expected keys that are not found in the dictionary.

    Prints a statement if there are missing keys.
    Used in this function to see if everything on the ip1List made it into the AFCDict.
    """

    dictKeys = list(dictionary.keys())
    notPassed = []
    for s in string:
        if s not in dictKeys:
            notPassed += [s]
    if notPassed != []:
        print('The following ip1s were not passed: {}'.format(str(notPassed)))

# Main function
fileID = getID(filePath)
AFCDict, spcList = createAFCDict(fileID)
if checkON == False:
    updateFST(fileID, AFCDict)
    closeFST(fileID)
    print('The FST file has been updated.')
    inDict(AFCDict, ip1List)
else:
    passAFCDict = checkPMlvl(AFCDict, ip1Check, pmType)
    if passAFCDict == True:
        updateFST(fileID, AFCDict)
        closeFST(fileID)
        print('The FST file has been updated.')
        inDict(AFCDict, spcList)
    else:
        del AFCDict[ip1Check]
        updateFST(fileID, AFCDict)
        closeFST(fileID)
        print('The FST file has been updated.')
        print('WARNING: The data for ip1 level: {} has not been updated and the original has been retained.'.format(ip1Check))
