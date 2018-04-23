import os
import sys
import rpnpy.librmn.all as rmn
from mpl_tookits.basemap import Basemap
import matplotlib.pyplot as plt


## concentration_plot takes the following inputs
## directory, path to the folder that contains the files to be reader
## spc, species of interest (will be a list later if needed)
## hy, hy of interest (will be a list later if needed)

// something that intakes arguments from command line
programName = sys.argv[0]
arguments = sys.argv[1:] //get all arguments that are provided
// need to ID the inputs, can they be lists?, first one is the argument, second is the value?
// alternatively, save the inputs into a list

## default values if not provided by commandline
defaultPath = os.path.direname(os.path.realpath(__file__)) //writes the full path, can maybe use current working directory: os.getcwd()
defaultSpc = 'TO3'
deafaultHy = [1.5M] //double check that this is right
defaultIncrLatLon = 0.1

def getSpcData (direct=defaultPath, spc=defaultSpc, hy=deafaultHy, incrLatLon=defaultIncrLatLon): //, dlat, dlong, date - don't think I need these
//use sys arguments if provided
# open file found in directory, for the concentration_plot function, loop through each one
# for the demo file, assuming only one file and one hy
    fileName = os. something //get file name
    fileID = rmn.fstopenall(fileName, rmn.FST_RO) //open file in read only mode
    keylist = rmn.fstinl(fileName, nomvar=spc, ip1=hy)
    rec = rmn.fstluk(keylist) //reads metadata from keylist
    data = rec['d'] //saves metadata in a list
    // add something to record the data value, x and y position
    // plot end results

latRange = (-90,90,defaultIncrLatLon) /all possible ranges of latitude
lonRange = (0,180,defaultIncrLatLon) /all possible ranges of longitude
//should actually take lat and lon as inputs, otherwise generating over very large area
dataList = []

    for i in data:
        #all rows of data
        // maybe need to get lat lon to get correct data?
        subList = [data[i]] //correct syntax for creating a list?
        //note,can't use tuples because need to modify later
        for lat in range latRange:
            subList.append[lat]
            for lon in range lonRange:
                subList.append[lon]
                dataList.append(subList)
    return dataList

# set up orthographic map projection with
# perspective of satellite looking down at 50N, 100W.
# use low resolution coastlines.
defaultLat=45
defaultLon=-100
defaultProjection='ortho'
defaultWidth = 5E6
defaultHeight = 5E6
//otherwise takes sys provided arguments, latlon seems to be perspective, not range of info
//note that example varies the values of lon_0 and lat_0 based on the ni,nj values
//values of width and height calculations to be added
def generateMap(data):
    //need to load lat,lon,dataVal
    map.pcolormesh(lon,lat,dataVal,latlon=True)//plot the data values, can deal with colour set later
    map = Basemap(projection=defaultProjection, lat_0=defaultLat, lon_0=defaultLon, resolution='l',width=defaultWidth,height=defaultHeight)
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.drawmapboundary
    map.drawmeridians
    map.drawparallels
    plt.colorbar(label='something for now')//get to read the data type for this
    plt.clim(minVal,maxVal)//values will have to be grabbed from the results in the data





## concentration_plot //don't actually need to write definition
