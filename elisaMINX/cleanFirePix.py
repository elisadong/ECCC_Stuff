# take top 10% of all fire pixels, and also filter by confidence level
# run this function from inside a folder with all the text files you want to clean

import os
import numpy as np
import sys

pwrInd = 5
pwrPrct = 80
confInd = 11
confMin = 60

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

#TODO:
# rewrite these options using data frames rather than looping through by line
# reminder to ignore the first 4 rows

def cleanPrct():
    # get top 10% in all files
    cDir = sys.argv[1]
    txtList = os.listdir(cDir) #should probably add in some sorting
    pwrs = []
    for txt in txtList:
        tLines = readFile(os.path.join(cDir,txt))
        for rowInd, row in enumerate(tLines, 4):
            if rowInd < len(tLines):
                pwrs += [toFloats(tLines[rowInd])[pwrInd]]
    pwrMax = np.percentile(pwrs, pwrPrct)
    # write only the top 10%
    # can also add in extra conditions like, confidence >= 60
    for txt in txtList:
        writeInds = [0,1,2,3]
        tLines = readFile(txt)
        for rowInd, row in enumerate(tLines, 4):
            if rowInd < len(tLines):
                rowList = toFloats(tLines[rowInd])
                rowPwr = rowList[pwrInd]
                rowConf = rowList[confInd]
                if rowPwr >= pwrMax and rowConf >= confMin:
                    writeInds += [rowInd]
        writeFile(txt, tLines, writeInds)

# added in, clean by country

def cleanCtryLoop():
    lonInd = 0
    latInd = 1
    import geopandas as gpd
    from shapely.geometry import Point
    fPixDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/canFirePix/"
    txtList =os.listdir(fPixDir)
    shpfilePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/canfuels/gpr_000b11a_e.shp"

    df = gpd.read_file(shpfilePath)
    shapeSeries = df['geometry']
    shapeList = list(shapeSeries)

    for txt in txtList:
        writeInds = [0,1,2,3]
        tLines =readFile(os.path.join(fPixDir,txt))
        print(txt)
        for rowInd, row in enumerate(tLines, 4):
            if rowInd < len(tLines):
                rowList = toFloats(tLines[rowInd])
                rowLon = rowList[lonInd]
                rowLat = rowList[latInd]
                rowOrg = Point(rowLon, rowLat)
                ptInShape = False
                for shp in shapeList:
                    ptInShape =shp.contains(rowOrg)
                    if ptInShape == True:
                        print('FirePix in Canada')
                        writeInds += [rowInd]
                        break
                    else:
                        continue
            if ptInShape == False:
                print('FirePix in US')
        writeFile(txt, tLines, writeInds)

def cleanCtryPnd():
    import pandas as pd
    import geopandas as gpd
    from shapely.geometry import Point
    fPixDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/canFirePix/"
    txtList =os.listdir(fPixDir)
    shpfilePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/canfuels/gpr_000b11a_e.shp"
    saveDir = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/newCanFirePix/"
    df = gpd.read_file(shpfilePath)
    shapeSeries = df['geometry']
    shapeList = list(shapeSeries)
    for txt in txtList:
        if os.path.isfile(os.path.join(saveDir,txt)) == False:
            print(txt)
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
                print(tdf)
                for r in tdf.index:
                    if tdf['inShape'][r] == True:
                        writeInds += [tdf['index'][r]]
                    #tdf = tdf.drop('ptOrg', axis=1)
                    #tdf = tdf.drop('inShape', axis=1)
                    #tdf = tdf.drop('index', axis=1)
                writeFile(os.path.join(saveDir, txt),tLines,writeInds)
                print('saved {}'.format(txt))
            except:
                pass
