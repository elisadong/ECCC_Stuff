#! /usr/bin/env python

# changes the emission binning in an fst file. See the default bin details for bin breakdown
# to run, either change the filePath below, then run from command line
## ex. python bin12.py
# or, provide the filepath directly from command line
## ex. python bin12.py /home/pathtofile/2017090112_fire_emissions_12bin.fst

# An additional check is built into the createAFCDict function in case there are two records that match the spc/ip1/dateo

import rpnpy.librmn.all as rmn
import os
import numpy as np
import sys

#default file
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
ip1 = 12001

if len(sys.argv) > 1:
    filePath = sys.argv[1]

def getID(filePath):
    ''' Returns the ID of the fst file provided by the filePath. '''

    rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)
    # open file to process
    fileID = rmn.fstopenall(filePath, rmn.FST_RW)
    return fileID

def delRec(recKey):
    '''
    Deletes a record in an fst by a key.
    Used to allow user to use the same nomvar as an existing record,
    where we no longer want to keep the data from the existing record.'''

    try:
        rmn.fsteff(recKey)
        #print('record {} deleted'.format(recKey))
    except:
        pass

def scaleArray(FCDict,key1, key2, scaleFac, spc, binNum):
    '''
    Scales the data/val array in a dictionary by a factor of scaleFac.

    Additionally modifies the meta data in the dictionary containing the data.
    This is where the nomvar is updated for each spc/bin. '''

    newDict = dict(FCDict)
    scaleArray = (np.asarray(newDict[date][key2]))*scaleFac#[x*scaleFac for x in newDict[key][key2]]
    newDict[key1][key2] = scaleArray
    newDict[key1]['meta']['nomvar'] = '{}{}'.format(s, binNum)
    #newDict[key1]['meta']['dateo'] = key1
    #print(key1)
    #print(newDict[key1]['meta']['dateo'])
    return newDict

def getMeta(fileID, key):
    '''
    Gets the meta data from a key.
    Additionally returns the concentration data and the dateo. '''

    meta = rmn.fstluk(key)
    concData = meta['d']
    dateo = meta['dateo']
    del meta['d']
    del meta['datyp']
    return meta, concData, dateo

def getKeys(fileID, ip1, spc, spcApp):
    ''' Gets a list of keys from an FST based on ip1, spc (and any modifications to spc name). '''

    spcs = '{}{}'.format(spc, spcApp)
    keyList = rmn.fstinl(fileID, nomvar = spcs, ip1 = ip1)
    return keyList

def sumArrays(FCDict, key1, key2):
    '''
    Sums the list of arrays indicated by FCDict[key1][key2].
    Normally used for FCDict[date]['vals']. '''

    newDict = dict(FCDict)
    array = newDict[key1][key2] #should be dict[dateo][val]
    sumArray = sum(array)
    newDict[key1][key2] = sumArray
    return newDict

def createAFCDict(fileID, keyList1, keyList2, spc):
    '''
    Creates a dictionary FCDict in the format of:
    FCDict{dateo:{'meta':metadata, 'vals':list([data],[data])}, dateo2:{..}...} '''

    FCDict = {}
    for key1 in keyList1:
        meta, cData, dateo = getMeta(fileID, key1)
        # save this existing information
        meta['nomvar'] = '{}F'.format(spc)
        meta['dateo'] = dateo
        rmn.fstecr(fileID, cData, meta)
        # update the nomvar again for the combo
        meta['nomvar'] = '{}FC'.format(spc)
        if dateo not in FCDict.keys():
            try:
                FCDict[dateo] = {'meta':meta, 'vals':[cData]}
            except:
                print('Did not add data for dateo {} in key1s'.format(dateo))
        else:
            print('there are two keys for dateo {} in key1s'.format(dateo))
        delRec(key1)
    for key2 in keyList2:
        meta2, cData2, dateo2 = getMeta(fileID, key2)
        # save the existing info
        meta2['nomvar'] = '{}C'.format(spc)
        meta2['dateo'] = dateo2
        rmn.fstecr(fileID, cData2, meta2)
        meta2['nomvar'] = '{}FC'.format(spc)
        if dateo2 not in FCDict.keys():
            try:
                FCDict[dateo2] = {'meta':meta2, 'vals': [cData2]}
            except:
                print('Did not add data for dateo {} in keys2'.format(dateo2))
        elif len(FCDict[dateo2]['vals']) < 2:
            try:
                FCDict[dateo2]['vals'] += [cData2]
            except:
                print('Did not add arrays for dateo {} in keys2'.format(dateo2))
        else:
            print('there are two keys for dateo {} in keys2'.format(dateo2))
        delRec(key2)
    return FCDict

def closeFST(fileID):
    ''' Close the FST to prevent file corruption and the like. '''

    rmn.fstcloseall(fileID)


fileID = getID(filePath)
for s in species:
    keyList1 = getKeys(fileID, ip1, s, 1)
    keyList2 = getKeys(fileID, ip1, s, 2)
    #print('keylist1 len is {}'.format(len(keyList1)))
    #print('keyList2 len is {}'.format(len(keyList2)))
    FCDict = createAFCDict(fileID, keyList1, keyList2, s)
    dateoList = FCDict.keys()
    dateoList.sort()
    #print('date list is len {}'.format(len(dateoList)))
    for date in dateoList:
        sumFCDict = sumArrays(FCDict, date, 'vals')
        for bInd, bin in enumerate(bins):
            scaleDict = scaleArray(sumFCDict, date, 'vals', bin, s, binApp[bInd])
            rmn.fstecr(fileID, scaleDict[date]['vals'], scaleDict[date]['meta'])
    #        print('New records written for spc {} and date {}.'.format(scaleDict[date]['meta']['nomvar'], date))
    #print(dateoList)


closeFST(fileID)
