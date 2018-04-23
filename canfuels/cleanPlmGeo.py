import geopandas as gpd
from shapely.geometry import Point
import os
import postMINX as pm

# note, when reading in the shapefile, need to also have the .shx in the same folder (possibly the .dbf and the .prj as well)

shpfilePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/canfuels/gpr_000b11a_e.shp"
plumePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/plumesP50C60/"
canPlumes = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/plumesP50C60/plumesCan/"
usPlumes = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/plumesP50C60/plumesUS/"

def moveFile(file, sourceLoc, newLoc):
    oldFileLoc = os.path.join(sourceLoc, file)
    newFileLoc = os.path.join(newLoc, file)
    os.rename(oldFileLoc, newFileLoc)
    print('File moved to {}'.format(newFileLoc))

df = gpd.read_file(shpfilePath)
shapeSeries = df['geometry']
shapeList = list(shapeSeries)

tempList = os.listdir(plumePath)
plumeList =[]
for i in tempList:
    if '.txt' in i:
        plumeList +=[i]

for p in plumeList:
    plm = pm.Plume.from_filename(os.path.join(plumePath, p))
    origin = plm['origin']
    ptOrigin = Point(origin[0], origin[1])
    print(ptOrigin)
    ptInShape = False
    for shp in shapeList:
        ptInShape = shp.contains(ptOrigin)
        if ptInShape == True:
            print('plume is in Canada')
            moveFile(p, plumePath, canPlumes)
            break
        else:
            continue
    if ptInShape == False:
        print('plume is in US')
        moveFile(p,plumePath, usPlumes)
