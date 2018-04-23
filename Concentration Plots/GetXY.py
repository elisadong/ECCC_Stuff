# Convert latlon to xy
# Use results to add to lat and lon list of coordinates for given concentration

import os
import rpnpy.librmn.all as rmn
import numpy as np

defaultLat=0
defaultLon=0

getGrid(directory)
xypost = rmn.gdxyfll(grid,lat=defaultLat,lon=defaultLon)
gridx = np.int
