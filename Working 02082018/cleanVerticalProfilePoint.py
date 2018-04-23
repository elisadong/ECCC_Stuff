#!/usr/bin/python

# cleanVerticalProfilePoint
# plots graph for one point, for as many models as required

import sys
import os
import numpy as np
import rpnpy.librmn.all as rmn
import ast
import pandas as pd
import matplotlib.pyplot as plt
import errno

# default values
currentDir = os.getcwd()
saveLoc = os.path.join(currentDir, 'VerticalProfiles')
ip1List = [95528103, 97341399, 95360636, 94517525, 94506747, 95203234, 95177672, 97164386, 97013298, 95348903, 94790638, 94765118, 94740623, 94717051, 96320038, 96279484, 95879133, 95834648, 96886973, 96783497, 95335517, 95151143, 95122572, 94497415, 94489330, 96410389, 96363016, 95092030, 95060058, 95791222, 95750107, 94544377, 94530094, 94650523, 94630473, 95577949, 95551544, 95369037, 96154478, 96109937, 96583883, 98293613, 95304400, 94902189, 94873167, 96064236, 96017812, 94962914, 94932068, 95248115, 95226733, 95639424, 95607291, 94576580, 94560020, 96699938, 96633977, 95320496, 98024801, 97840480, 95287241, 95971125, 95924598, 94694044, 94671769, 97720038, 97644070, 95268471, 94611577, 94593675, 95711238, 95674198, 94844858, 94817232, 95027249, 94994764, 94482342, 94476231, 96239545, 96197821]
defaultSpc = 'TO3'
defaultArea = [241, 36.7783]
defaultType = 'l'
baseModel = os.path.join(currentDir, 'testFiles/adomKPPKPPB_KPP/20100710000000/model/2010071000_048')
altModel = baseModel
paramsBase = {'filePath':baseModel, 'spc':defaultSpc, 'area':defaultArea, 'type': defaultType, 'title':'test'}
paramsAlt = {'filePath':altModel, 'spc':defaultSpc, 'area':defaultArea, 'type': defaultType, 'title':'testagain'}
# constants regarding plot visuals
defaultFontSize = 8
defaultRotation = 10
# set limits on the x axis
xleft = None
xright = None

#convert to literal
def toLit(string):
    nowLit = ast.literal_eval(string)
    return nowLit

# open file
def openFile(filePath):
    fileID = rmn.fstopenall(filePath,rmn.FST_RO)
    return fileID

# get the key for matching spc/ip1
def getKey(fileID, spc, ip1):
    key = rmn.fstinf(fileID,nomvar=spc,ip1=ip1)['key']
    return key

# get the gridData(array),lon(list),lat(list) for one key(one ip1)
def getGridData(fileID, key):
    fileMeta = rmn.fstprm(key)
    fileMeta['iunit'] = fileID
    gridMeta = rmn.ezqkdef(fileMeta)
    gridData = rmn.gdll(gridMeta)
    lonData = gridData['lon'].tolist()
    latData = gridData['lat'].tolist()
    return gridData, lonData, latData

# get the concentration data for one key(one ip1)
def getConcData(key):
    dataRec = rmn.fstluk(key)
    concData = dataRec['d']
    return concData

# get the ni nj for a lon/lat value
def getFuzzyNiNj(lonVal, latVal, gridData):
    xypos = rmn.gdxyfll(gridData, latVal, lonVal)
    ni = int(round(xypos['x'])) - 1 #-1 to use 0 as list index
    nj = int(round(xypos['y'])) - 1
    return ni, nj

# get the concentration for ninj
def getConc(ni,nj,concData):
    concVal = concData[ni][nj]
    return concVal

# get hy and conc lists, (note that for more than one value, should use ni and nj ranges)
def getPlotValues(fileID, ip1List, spc, ni, nj):
    CHdf = pd.DataFrame(columns = ['hy','conc'])
    for ip1 in ip1List:
        key = getKey(fileID, spc, ip1)
        concData = getConcData(key)
        hyVal = ((os.popen('r.ip1 {}'.format(ip1))).read()).rstrip().lstrip()
        concVal = getConc(ni,nj,concData)
        CHdf = CHdf.append({'hy':hyVal, 'conc':concVal},  ignore_index=True)
    CHdf['hy'] = CHdf['hy'].map(lambda x: x[:-3])
    CHdf = CHdf.sort('hy', ascending = False).reset_index(drop = True)
    hyList = CHdf['hy'].tolist()
    concList = CHdf['conc'].tolist()
    return hyList, concList

# plot one plot
def plotVP(concList, hyList, title, spc, savePath=saveLoc):
    fig = plt.figure(figsize=(6,8))
    plt.plot(concList, hyList, 'b')
    plt.xticks(fontsize=defaultFontSize,rotation=defaultRotation)
    plt.yticks(fontsize=defaultFontSize)
    plt.title('Vertical Profile of {}\n{}'.format(spc, title))
    plt.yscale('log')
    ax = plt.gca()
    ax.invert_yaxis()
    ax.set_xlim(left=xleft,right=xright)
    plt.xlabel('Concentration')
    plt.ylabel('Hybrid Level')
    makeDir(savePath)
    fig.savefig(os.path.join(savePath, title+spc+'.png'), dpi = 300)
    plt.close()

# plot stacked plot
def plotStackVP(concList1, concList2, hyList1, hyList2, title1, title2, spc, savePath=saveLoc):
    fig = plt.figure(figsize=(6,8))
    plt.plot(concList1, hyList1, 'b')
    plt.plot(concList2, hyList2, 'ko')
    plt.xticks(fontsize=defaultFontSize,rotation=defaultRotation)
    plt.yticks(fontsize=defaultFontSize)
    plt.title('Vertical Profile of {}\nBase Model: {}, Overlay Model: {}'.format(spc, title1, title2))
    plt.yscale('log')
    ax = plt.gca()
    ax.invert_yaxis()
    ax.set_xlim(left=xleft,right=xright)
    plt.xlabel('Concentration')
    plt.ylabel('Hybrid Level')
    makeDir(savePath)
    fig.savefig(os.path.join(savePath, title1+title2+spc+'.png'), dpi = 300)
    plt.close()
        # plot more than one VP, provide a list of conclist/hy values
    # have sequence of colours and shapes

# make directories if they don't exist
def makeDir(path):
    """ Creates a directory for the provided path if it does not exist. """

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

#for one plot
# note that this can be converted if provide savenames and reference files as list, produces one for each case

def plotOne(params=paramsBase, savePath = saveLoc):
    rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)
    filePath, spc, area, title = params['filePath'], params['spc'], params['area'], params['title']
    fileID = openFile(filePath)
    gridKey = getKey(fileID, spc, 76696048)
    gridData, lonData, latData = getGridData(fileID, gridKey)
    try:
        areaType = params['type']
        if areaType == 'n':
            ni, nj = area[0], area[1]
        else:
            ni, nj = getFuzzyNiNj(area[0], area[1], gridData)
    except: #default to ninj type, only will get if lon/lat failed to get ninj values
        ni, nj = area[0], area[1]
    hyList, concList = getPlotValues(ip1List, spc, ni, nj)
    plotVP(concList, hyList, title, spc, savePath)

def plotTwo(paramsBase=paramsBase, paramsAlt=paramsAlt, savePath = saveLoc):
    rmn.fstopt(rmn.FSTOP_MSGLVL,rmn.FSTOPI_MSG_CATAST)
    filePath1, spc1, area1, title1 = paramsBase['filePath'], paramsBase['spc'], paramsBase['area'], paramsBase['title']
    filePath2, spc2, area2, title2 = paramsAlt['filePath'], paramsAlt['spc'], paramsAlt['area'], paramsAlt['title']
    fileID1= openFile(filePath1)
    print(filePath)
    print(spc1)
    gridKey1 = getKey(fileID1, spc1, 76696048)
    gridData1, lonData1, latData1 = getGridData(fileID1, gridKey1)
    try:
        areaType1 = paramsBase['type']
        if areaType1 == 'n':
            ni1, nj1 = area1[0], area1[1]
        else:
            ni1, nj1 = getFuzzyNiNj(area1[0], area1[1], gridData1)
    except: #default to ninj type, only will get if lon/lat failed to get ninj values
        ni1, nj1 = area1[0], area1[1]

    fileID2 = openFile(filePath2)
    gridKey2 = getKey(fileID2, spc2, 76696048)
    gridData2, lonData2, latData2 = getGridData(fileID2, gridKey2)
    try:
        areaType2 = paramsAlt['type']
        if areaType2 == 'n':
            ni2, nj2 = area2[0], area2[1]
        else:
            ni2, nj2 = getFuzzyNiNj(area2[0], area2[1], gridData2)
    except: #default to ninj type, only will get if lon/lat failed to get ninj values
        ni2, nj2 = area1[0], area1[1]

    hyList1, concList1= getPlotValues(fileID1, ip1List, spc1, ni1, nj1)
    hyList2, concList2= getPlotValues(fileID2, ip1List, spc2, ni2, nj2)
    plotStackVP(concList1, concList2, hyList1, hyList2, title1, title2, spc1, savePath)

# update via sys.argv
try:
    paramsBase = toLit(sys.argv[1])
    paramsAlt = toLit(sys.argv[2])
except:
    pass
plotTwo(paramsBase,paramsAlt,saveLoc)
