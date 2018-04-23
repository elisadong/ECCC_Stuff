def standardizeLon(lon):
    if lon > 180:
        lon=lon-180
        return lon
    else:
        return lon

def standardizeLonLat(lon,lat):
    return {'lon':standardizeLon(lon), 'lat':standardizeLat(lat)}
