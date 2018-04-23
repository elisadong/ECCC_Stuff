## function to calculate the distances between two points based off of lon/lat coordinates
## intent is to get the max height and width for concentration plot
## Note that this can be done with geopy, but can't be installed

from math import sin,cos,sqrt,acos,radians,pi

## radius of earth in m
R = 6400E3
lowLon=179.08134
highLon=347.04001
lowLat=13.749236
highLat=82.181656

def maxDimensions(lowLon,highLon,lowLat,highLat,buff):
    ## convert to radians
    lonDiff = radians(highLon - lowLon)
    latDiff = radians(highLat - lowLat)
    radBuffer = radians(buff)
    ## calculate distance for longitude, radius of earth
    lonDist = R*(latDiff+radBuffer)

    ## calculate distance for latitude, based on longitude
    Rad1=R*cos(lowLat)
    Rad2=R*cos(highLat)
    latDist1=Rad1*(latDiff+radBuffer)
    latDist2=Rad2*(latDiff+radBuffer)
    if latDist1>latDist2:
        latDist=latDist1
    else:
        latDist=latDist2
