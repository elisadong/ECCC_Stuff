import os
import rpnpy.librmn.all as rmn

## Test Constants:
defaultDir=os.getcwd()
defaultIp1=76696048 ##this is the 1.5M level for GEMMACH_dev_rev_67537156/20100710000000/model/ 2010071000_013
## note that this appears to be consistent, can use ip1 instead of hy values to look for
## should also make an IP1 and level matching index
## can create a script to look through all files and find IP1's. if doesn't exist in index, add to
## else return hy/ip1 as needed
defaultSpc='TO3'

def getData(level, spc, fileID):
    keylist = rmn.fstinl(fileID,nomvar=spc,ip1=level)
    rec = rmn.fstluk(keylist)
    data = rec['d']
    return rec, data
    ## returning just data assumes that the data also contains xy information in association with the dataVals

def getGrid(directory):
    fileName = os.listdir(directory)
    fileID = rmn.fstopenall(directory+fileName, rmn.FST_RO)
    getKeys(level=defaultIp1,spc=defaultSpc,fileID=fileID)
    #now have data for all applicable fields

    grid=rmn.readGrid(fileID,rec)
    print ('Grid Type is: ' + grid[grtyp]) ##grid type  (one of 'Z', '#', 'Y', 'U')
    ## note that if there are variable gridtypes, would be problematic when coding in conversion to lat/lon or x/y
    ## if plotted in latlon, then is not a concern, and can convert to xy for new grid layout to match the basemap
    return grid

    ##things wanted out of grid:
    ## shape, (ni,nj) to resize to


##   Returns:
       {
           'id'    : grid id, same as input arg
           'shape'  : (ni, nj) # dimensions of the grid
           'ni'     : first dimension of the grid
           'nj'     : second dimension of the grid
           'grtyp'  : type of geographical projection
                      (one of 'Z', '#', 'Y', 'U')
           'ig1'    : first grid descriptor
           'ig2'    : second grid descriptor
           'ig3'    : third grid descriptor
           'ig4'    : fourth grid descriptor
           'grref'  : grid ref type (one of 'A', 'B', 'E', 'G', 'L', 'N', 'S')
           'ig1ref' : first grid descriptor of grid ref
           'ig2ref' : second grid descriptor of grid ref
           'ig3ref' : third grid descriptor of grid ref
           'ig4ref' : fourth grid descriptor of grid ref
           ...
           list of other parameters is grtyp dependent,
           See defGrid_* specific function for details
       }
