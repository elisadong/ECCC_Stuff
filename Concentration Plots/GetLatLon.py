import os
import rpnpy.librmn.all as rmn

## Test Constants:
defaultDir=os.getcwd()
defaultIp1=52454224
defaultSpc='TO3'

def getData(level, spc, fileID):
    keylist = rmn.fstinl(fileID,nomvar=spc,ip1=level)
    rec = rmn.fstluk(keylist)
    data = rec['d']
    return rec, data

def getGrid(directory):
    fileName = os.listdir(directory)
    fileID = rmn.fstopenall(directory+fileName, rmn.FST_RO)
    getKeys(level=defaultIp1,spc=defaultSpc,fileID=fileID)
    #now have data for all applicable fields

    grid=rmn.readGrid(fileID,rec)
    return grid
