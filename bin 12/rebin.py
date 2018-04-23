import rpnpy.librmn.all as rmn
import os
import numpy as np

filePath = '2017090112_fire_emissions_12bin.fst'

BIN1=0.000005
BIN2=0.001054
BIN3=0.041466
BIN4=0.328735
BIN5=0.490694
BIN6=0.127949
BIN7=0.009841
BIN8=0.000253
BIN9=0.000002
BINA=0.000000001

bins = [BIN1, BIN2, BIN3,BIN4,BIN5,BIN6, BIN7,BIN8, BIN9,BINA]
binApp = [1,2,3,4,5,6,7,8,9,'A']
species = ['ESU', 'ECM', 'EAM','EEC', 'ENT','EPC']
ip1List = [12001]


def getID(filePath):
    """ Gets the fileID and opens the file in read/write mode. """

    rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)
    # open file to process
    fileID = rmn.fstopenall(filePath, rmn.FST_RW)
    return fileID

def delRec(recKey):
    try:
        rmn.fsteff(recKey)
        print('record deleted')
    except:
        pass

def getConc(fileID, spc, ip1):
    try:
        dataKey1 = rmn.fstinf(fileID, nomvar='{}1'.format(spc), ip1=ip1)['key']
        dataRec1 = rmn.fstluk(dataKey1)
        concData1 = dataRec1['d']
        meta = dataRec1
        del meta['d']
        del meta['datyp']
        dataKey2 = rmn.fstinf(fileID, nomvar='{}2'.format(spc), ip1=ip1)['key']
        dataRec2 = rmn.fstluk(dataKey2)
        concData2 = dataRec2['d']
        concData = sum(concData1, concData2)
        return concData, meta, dataKey1, dataKey2
    except:
        print('WARNING: Data does not exist for ip1 {} and spc {}.'.format(ip1, spc))

def scaleArray(FCDict,ip1List, key2, scaleFac, spc, binNum):
    newDict = dict(FCDict)
    for ip1 in ip1List:
        #this should result in an array, not a list
        scaleArray = (np.asarray(newDict[ip1][key2]))*scaleFac#[x*scaleFac for x in newDict[key][key2]]
        newDict[ip1][key2] = scaleArray
        if binNum ==10:
            binNum = 'A'
        newDict[ip1]['meta']['nomvar'] = '{}{}'.format(spc, binNum)
    return newDict

def updateFST(fileID, AFCDict):
    """ Updates the original FST file with the new PM data and corresponding metadata. """

    for ip1 in AFCDict:
        rmn.fstecr(fileID, AFCDict[ip1]['vals'], AFCDict[ip1]['meta'])

def createAFCDict(fileID, spc):
    """ Creates a dictionary with level 1 key: ip1 from ip1List, level 2 keys: 'vals' and 'meta'.

    key 'meta' contains meta data (will save the last species in the spcList meta data) with 'd' and 'datyp' removed
    key 'vals' contains the sum of all concentration data retrieved from spc in spcList using the sumArrays function.
    """

    AFCDict = {}
    for ip1 in ip1List:
        AFCDict[ip1] = {}
        AFCDict[ip1]['vals'] = []
        #print('Getting AFC Data and adding to dictionary...')
        AFCData, AFCMeta, keyF, keyC = getConc(fileID, spc, ip1)
        try:
            AFCMeta['nomvar'] = '{}FC'.format(spc)
            AFCDict[ip1]['meta'] = AFCMeta
            AFCDict[ip1]['vals'] += [AFCData]
            delRec(keyF)
            delRec(keyC)
            #print('Adding AFCData for ip1: {}, spc: {}.'.format(ip1, spc))
        except:
            print('Did not add AFCData for ip1: {}, spc: {}.'.format(ip1, spc))
    #sumAFCDict = sumArrays(AFCDict, list(AFCDict.keys()), 'vals')
    return AFCDict

def closeFST(fileID):
    """ Closes the FST to prevent file corruption. """

    rmn.fstcloseall(fileID)

def getDateO(fileID, s, ip1):
    keyList = rmn.fstinl(fileID, nomvar = '{}1')
fileID = getID(filePath)
for s in species:
    FCDict = createAFCDict(fileID, s)
    dateoList = []

    for bInd, bin in enumerate(bins):
        scaleDict = scaleArray(FCDict, ip1List, 'vals', bin, s, binApp[bInd])
        updateFST(fileID, scaleDict)


closeFST(fileID)

for date in dateoList:
    sumFCDict = sumArrays(FCDict, date, 'vals')
    for bInd, bin in enumerate(bins):
        scaleDict = scaleArray(sumFCDict, date, 'vals', bin, s, binApp[bInd])
        rmn.fstecr(fileID, scaleDict[date]['vals'], scaleDict[date]['meta'])
    print('New records written for spc {} and date {}.'.format(s, date))
