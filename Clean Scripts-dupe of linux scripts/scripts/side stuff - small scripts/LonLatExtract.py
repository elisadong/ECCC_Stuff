import os
import rpnpy.librmn.all as rmn
import numpy as np

currentDir = os.path.dirname(os.path.abspath(__file__))
defaultFile = '2010071000_001'
# or should it check through every single file?
defaultSpc='TO3'
defaultIp1=76696048

def getMeta(filePath, level, spc):
    try:
        fileID = rmn.fstopenall(filePath,rmn.FST_RO)
        dataKey = rmn.fstinf(fileID,nomvar=spc,ip1=level)['key']
        dataRec = rmn.fstluk(dataKey)
        # get the concentration data
        concData = dataRec['d']
        fileMeta = rmn.fstprm(dataKey)
        fileMeta['iunit'] = fileID
        gridData = rmn.ezqkdef(fileMeta)
        gridDecode = rmn.decodeGrid(gridData)
        # get the lonlat data
        llGridData = rmn.gdll(gridData)
        latData = llGridData['lat']
        lonData = llGridData['lon']
        return concData, lonData.tolist(), latData.tolist()
    except:
        print('Could not read file.')
        pass

def getConc(ni,nj,concData):
    concVal = concData[ni][nj]
    return concVal

def getNiNj(lonVal,latVal,lonData,latData):
    ninj = []
    for i in range(len(lonData)):
        try:
            lonNjIndex = lonData[i].index(lonVal)
            latNjIndex = latData[i].index(latVal)
            if lonNjIndex == latNjIndex:
                ninj += [i,lonIndex]
        except:
            pass
    return ninj[0], ninj[1]

def printConc (filePath=defaultFilePath, lonVal=defaultLon,latVal=defaultLat):
