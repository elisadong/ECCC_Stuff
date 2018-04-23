
import numpy as np
import os
import rpnpy.librmn.all as rmn
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
f = open('test.txt')
h = 4558250
w = 4682500
lines = f.readlines()


def toInts(row):
    ints = map(int, row.split(','))
    return ints



m = Basemap(projection='lcc', width = w, height=h, lat_1=49, lat_2=77, lat_0=49, lon_0=-95)

def conv(x,y):
    lon,lat=m(x,y,inverse=True)
    return lon,lat

w = open('fuelRef.txt','w+')
w.write('Lon,Lat,Fuel\n')

for x, line in enumerate(lines):
    row = toInts(line)
    for y, val in enumerate(row):
        lon, lat = m((x+1)*250, (y+1)*250, inverse=True)
        w.write('{},{},{}\n'.format(lon,lat,val))

w.close()
