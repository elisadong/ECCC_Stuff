#! /usr/bin/env python

import os
import numpy as np
import sys
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry import Polygon
import pp

def toFloats(row):
    floats = map(float, row.split())
    return floats

def readFile(path):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    return lines

def writeFile(path, lines, writeInds):
    f = open(path, 'w')
    for rowInd, string in enumerate(lines):
        if rowInd in writeInds:
            f.write(lines[rowInd])
    f.close()


fPixDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/allFirePix/"
txtList =os.listdir(fPixDir)
shpfilePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/canfuels/gpr_000b11a_e.shp"
saveDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/canFPcleanParallel/"
df = gpd.read_file(shpfilePath)
shapeSeries = df['geometry']
shapeList = list(shapeSeries)
ppservers = ()

def cleanFP(txt,saveDir, fPixDir):
    tLines =readFile(os.path.join(fPixDir,txt))
    fpList =[]
    writeInds = [0,1,2,3]
    for rowInd, row in enumerate(tLines, 4):
        if rowInd < len(tLines):
            rowList = toFloats(tLines[rowInd])
            fpList += [rowList]
    try:
        tdf = pd.DataFrame(fpList)
        lon = list(tdf[0])
        lat = list(tdf[1])
        pts = []
        for x,lonVal in enumerate(lon):
            pts += [Point(lon[x],lat[x])]
        tdf['ptOrg'] = pts
        tdf['inShape'] = tdf['ptOrg'].map(lambda x: any(shp.contains(x) for shp in shapeList))
        tdf['index'] = tdf.index + 4
        for r in tdf.index:
            if tdf['inShape'][r] == True:
                writeInds += [tdf['index'][r]]
        writeFile(os.path.join(saveDir, txt),tLines,writeInds)
    except:
        pass

jobs = []
for txt in txtList:
    if os.path.isfile(os.path.join(saveDir,txt)) == False:
        job_server = pp.Server(ppservers=ppservers)
        jobs.append(job_server.submit(cleanFP(txt,saveDir,fPixDir)))
