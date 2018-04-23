#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2018-03-09 16:49:28
# Contact:    Nolan.Dickson@Canada.ca

'''
Utils functions for converting from MISR imagery block, line, sample indices to
longitude, latitude degrees, creating a real-colour imagery given a MISR terrain file
and calculating orbital paths from satellite TLE's

Note
----

All conversions are only accurate to a floating point rounding error (ish) and you should
not rely on BLStoLonLat == LonLattoBLS being true, as with all floating point arithmetic

Also, this straight up doesn't work yet.. so maybe definitely don't rely on it

Adapted from:
MISR-Toolkit [https://github.com/nasa/MISR-Toolkit]
GCTP [https://github.com/OkoSanto/GCTP]
MINX [https://github.com/nasa/MINX]
Space Oblique Mercator Projection Mathematical Development, by John P. Snyder
Numerical Recipes in C, Press et al.

Examples
--------

>>> import MISRutils as mu
>>>
>>> mu.BLStoLonLat(201, 51, 168, 567)
(-1.304313176667397, 49.49334168048807)
>>>
>>> mu.LonLattoBLS(201, -1.3043, 49.4933)
(51, 168, 567)

>>> import MISRutils as mu
>>> BRF_dict, ulc, lrc = mu.drawmisr('MISR_AM1_GRP_TERRAIN_GM_P042_O094342_AN_F03_0024.hdf', 45, 50)

>>> import MISRutils as mu
>>> import datetime
>>> orbit_paths = mu.draworbit(datetime.datetime(2017, 6, 15), direction='descending')
>>>
>>> for path in orbit_paths:
>>>    path_swath = mu.drawswath(path, 360 * 1000)  # MISR swath width = 360 km

'''


# ----------------------
# Data tables
# ----------------------


# relative offsets of all MISR blocks  #TODO whats up with the missing initial 0.?
rel_offset_table = [
    0.0, 16.0, 0.0, 16.0, 0.0, 0.0, 0.0, 16.0, 0.0, 0.0,
    0.0, 0.0, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0,
    -16.0, 0.0, 0.0, -16.0, 0.0, 0.0, -16.0, 0.0, -16.0,
    0.0, -16.0, 0.0, -16.0, -16.0, 0.0, -16.0, 0.0, -16.0,
    -16.0, 0.0, -16.0, -16.0, -16.0, 0.0, -16.0, -16.0, -16.0,
    -16.0, 0.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -32.0, -16.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -32.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, -16.0, -16.0, 0.0, -16.0,
    -16.0, -16.0, -16.0, -16.0, 0.0, -16.0, -16.0, -16.0, 0.0,
    -16.0, -16.0, 0.0, -16.0, 0.0, -16.0, -16.0, 0.0, -16.0,
    0.0, -16.0, 0.0, 0.0, -16.0, 0.0, -16.0, 0.0, 0.0, -16.0,
    0.0, 0.0, 0.0, 0.0, -16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    16.0, 0.0, 0.0, 16.0, 0.0, 0.0, 16.0, 0.0
]

# reference table for MISR orbit times, to account for drift
orbit_table = [   # julian numb    #  orb         calendar date   |  orb       calendar date
    2451530.84792, 2451599.51745,  # 0               extrapolated | 1000   25 Feb 2000  01:14:34
    2451668.18783, 2451736.85696,  # 2000    3 May 2000  17:19:55 | 3000   11 Jul 2000  09:23:28
    2451805.52616, 2451874.19467,  # 4000   18 Sep 2000  01:27:07 | 5000   25 Nov 2000  17:29:46
    2451942.86333, 2452011.53161,  # 6000    2 Feb 2001  09:32:38 | 7000   12 Apr 2001  01:34:58
    2452080.19981, 2452148.86869,  # 8000   19 Jun 2001  17:37:10 | 9000   27 Aug 2001  09:40:21
    2452217.53632, 2452286.20435,  # 10000   4 Nov 2001  01:41:45 | 11000  11 Jan 2002  17:43:42
    2452354.87290, 2452423.54264,  # 12000  21 Mar 2002  09:46:25 | 13000  29 May 2002  01:50:51
    2452492.21218, 2452560.88142,  # 14000   5 Aug 2002  17:54:59 | 15000  13 Oct 2002  09:58:41
    2452629.55123, 2452698.22027,  # 16000  21 Dec 2002  02:03:13 | 17000  27 Feb 2003  18:06:38
    2452766.88943, 2452835.55849,  # 18000   7 May 2003  10:10:13 | 19000  15 Jul 2003  02:13:40
    2452904.22810, 2452972.89852,  # 20000  21 Sep 2003  18:17:54 | 21000  29 Nov 2003  10:23:19
    2453041.56793, 2453110.23683,  # 22000   6 Feb 2004  02:27:16 | 23000  14 Apr 2004  18:30:29
    2453178.90627, 2453247.57547,  # 24000  22 Jun 2004  10:34:28 | 25000  30 Aug 2004  02:38:07
    2453316.24466, 2453384.91458,  # 26000   6 Nov 2004  18:41:45 | 27000  14 Jan 2005  10:46:26
    2453453.58422, 2453522.25365,  # 28000  24 Mar 2005  02:50:43 | 29000  31 May 2005  18:54:42
    2453590.92318, 2453659.59252,  # 30000   8 Aug 2005  10:58:49 | 31000  16 Oct 2005  03:02:40
    2453728.26194, 2453796.93158,  # 32000  23 Dec 2005  19:06:38 | 33000   2 Mar 2006  11:10:55
    2453865.60132, 2453934.27109,  # 34000  10 May 2006  03:15:21 | 35000  17 Jul 2006  19:19:49
    2454002.94042, 2454071.60990,  # 36000  24 Sep 2006  11:23:39 | 37000   2 Dec 2006  03:27:42
    2454140.28003, 2454208.94952,  # 38000   8 Feb 2007  19:32:41 | 39000  18 Apr 2007  11:36:45
    2454277.61895, 2454346.28864,  # 40000  26 Jun 2007  03:40:44 | 41000   2 Sep 2007  19:45:05
    2454414.95767, 2454483.62756,  # 42000  10 Nov 2007  11:48:29 | 43000  18 Jan 2008  03:53:08
    2454552.29679, 2454620.96679,  # 44000  26 Mar 2008  19:56:49 | 45000   3 Jun 2008  12:01:37
    2454689.63629, 2454758.30577,  # 46000  11 Aug 2008  04:05:42 | 47000  18 Oct 2008  20:09:45
    2454826.97564, 2454895.64529,  # 48000  26 Dec 2008  12:14:22 | 49000   5 Mar 2009  04:18:40
    2454964.31481, 2455032.98389,  # 50000  12 May 2009  20:22:46 | 51000  20 Jul 2009  12:26:15
    2455101.65328, 2455170.32292,  # 52000  27 Sep 2009  04:30:10 | 53000   4 Dec 2009  20:34:27
    2455238.99233, 2455307.66179,  # 54000  11 Feb 2010  12:38:24 | 55000  21 Apr 2010  04:42:25
    2455376.33139, 2455445.00076,  # 56000  28 Jun 2010  20:46:39 | 57000   5 Sep 2010  12:50:32
    2455513.67020, 2455582.33953,  # 58000  13 Nov 2010  04:54:32 | 59000  20 Jan 2011  20:58:22
    2455651.00923, 2455719.67873,  # 60000  30 Mar 2011  13:02:44 | 61000   7 Jun 2011  05:06:49
    2455788.34839, 2455857.01773,  # 62000  14 Aug 2011  21:11:07 | 63000  22 Oct 2011  13:14:58
    2455925.68751, 2455994.35708,  # 64000  30 Dec 2011  05:19:27 | 65000   7 Mar 2012  21:23:38
    2456063.02668, 2456131.69605,  # 66000  15 May 2012  13:27:52 | 67000  23 Jul 2012  05:31:45
    2456200.36543, 2456269.03520,  # 68000  29 Sep 2012  21:35:40 | 69000   7 Dec 2012  13:40:08
    2456337.70461, 2456406.37452,  # 70000  14 Feb 2013  05:44:05 | 71000  23 Apr 2013  21:48:45
    2456475.04394, 2456543.71351,  # 72000   1 Jul 2013  13:52:43 | 73000   8 Sep 2013  05:56:54
    2456612.38295, 2456681.05237,  # 74000  15 Nov 2013  22:00:53 | 75000  23 Jan 2014  14:04:51
    2456749.72210, 2456818.39189,  # 76000   2 Apr 2014  06:09:16 | 77000   9 Jun 2014  22:13:46
    2456887.06129, 2456955.73079,  # 78000  17 Aug 2014  14:17:42 | 79000  25 Oct 2014  06:21:47
    2457024.40029, 2457093.06968,  # 80000   1 Jan 2015  22:25:52 | 81000  11 Mar 2015  14:29:47
    2457161.73915, 2457230.40861,  # 82000           extrapolated | 83000           extrapolated
    2457299.07808, 2457367.74754,  # 84000           extrapolated | 85000           extrapolated
    2457436.41701, 2457505.08647,  # 86000           extrapolated | 87000           extrapolated
    2457573.75594, 2457642.42540,  # 88000           extrapolated | 89000           extrapolated
    2457711.09487, 2457779.76434,  # 90000           extrapolated | 91000           extrapolated
    2457848.43380, 2457917.10327,  # 92000           extrapolated | 93000           extrapolated
    2457985.77273, 2458054.44220,  # 94000           extrapolated | 95000           extrapolated
    2458123.11166, 2458191.78113,  # 96000           extrapolated | 97000           extrapolated
    2458260.45059, 2458329.12006,  # 98000           extrapolated | 99000           extrapolated
    2458396.41623                  # 99999           extrapolated
]


# ----------------------
# Internal Functions
# ----------------------


def _som_init(path_num):
    '''
    function for defining mathematical constants necessary in the forward and reverse SOM to geographic coordinate
    '''
    import math

    def _som_series(dlam):
        '''
        function representing mathematical series in inverse SOM equations (see John Snyder, 1981)

        takes in dlam and uses parameters setup initially in _som_init to return fa2, fa4, fb, fc1, fc3
        '''
        import math

        dlam = dlam * 0.0174532925      # Convert dlam to radians

        sd = math.sin(dlam)
        sdsq = sd * sd

        s = p21 * sa * math.cos(dlam) * math.sqrt((1.0 + t * sdsq) / ((1.0 + w * sdsq) * (1.0 + q * sdsq)))
        h = math.sqrt((1.0 + q * sdsq) / (1.0 + w * sdsq)) * (((1.0 + w * sdsq) / ((1.0 + q * sdsq) * (1.0 + q * sdsq))) - p21 * ca)

        sq = math.sqrt(xj * xj + s * s)

        fb = (h * xj - s * s) / sq
        fa2 = fb * math.cos(2.0 * dlam)
        fa4 = fb * math.cos(4.0 * dlam)

        fc = s * (h + xj) / sq
        fc1 = fc * math.cos(dlam)
        fc3 = fc * math.cos(3.0 * dlam)

        return {'a2': fa2, 'a4': fa4, 'b': fb, 'c1': fc1, 'c3': fc3}

    # Define Terra and MISR physical constants

    r_maj = 6378137.                    # use WGS 84 spheroid
    r_min = 6356752.314245              # r_major & r_minor, respectively

    S2R = 4.848136811095359e-6
    inc_ang = 98018013.7520             # Incidence angle
    asc_lon = 127045037.928240340       # Ascending longitude
    time = 98.88                        # Orbit time in minutes (i.e. P2)

    # Compute longitude of ascending node for this path from first path

    asc_dd = _DMStoDD(asc_lon)                                    # into dd
    asc_rad = asc_dd * (math.pi / 180.0)                          # into rads
    asc_path = asc_rad - (2.0 * math.pi / 233) * (path_num - 1)   # calc in rads
    asc_lon = asc_path / (math.pi / 180.0)                        # into dd

    lon_center = asc_lon * 3600 * S2R
    alf = _DMStoDD(inc_ang) * 3600 * S2R

    es = 1 - (r_min / r_maj)
    p21 = time / 1440.0
    ca, sa = math.cos(alf) if abs(math.cos(alf)) > 1.e-9 else 1.e-9, math.sin(alf)

    e2c = es * ca * ca
    e2s = es * sa * sa
    q = e2s / (1. - es)
    u = e2c / (1. - es)

    w = ((1.0 - e2c) / (1.0 - es)) * ((1.0 - e2c) / (1.0 - es)) - 1
    t = (e2s * (2.0 - es)) / ((1. - es) * (1. - es))
    xj = (1. - es) * (1. - es) * (1. - es)

    # compute sums to get a, b and c coefficients

    sumdict = _som_series(0.)  # initial a2, a4, b, c1, c3 coefficients

    for i in range(9, 81 + 1, 18):
        fdict = _som_series(i)
        for coeff, val in fdict.items():
            sumdict[coeff] += 4. * val

    for i in range(18, 72 + 1, 18):
        fdict = _som_series(i)
        for coeff, val in fdict.items():
            sumdict[coeff] += 2. * val

    dlam = 90.0
    fdict = _som_series(dlam)
    for coeff, val in fdict.items():
        sumdict[coeff] += val

        if coeff == 'c1':
            sumdict[coeff] /= 15.
        elif coeff == 'c3':
            sumdict[coeff] /= 45.
        elif coeff == 'a4':
            sumdict[coeff] /= 60.
        else:
            sumdict[coeff] /= 30.

    # TODO package these constants better?
    return r_maj, sumdict['b'], es, p21, t, w, q, xj, sumdict['a2'], sumdict['a4'], sumdict['c1'], sumdict['c3'], ca, sa, u, lon_center


def _DMStoDD(angle):
    '''converts a packed DMS angle to dd'''

    fac = (-1 if (angle < 0.0) else 1)

    # The degrees are separated out
    angle = abs(angle)
    deg_factor = 1000000.0

    if angle / deg_factor < 360:
        deg = angle / deg_factor

        # The minutes are separated out
        min_sep = angle - deg * deg_factor
        mn_factor = 1000.0

        if min_sep / mn_factor < 360:
            mn = min_sep / mn_factor

            # The seconds are computed
            sec_sep = min_sep - mn * mn_factor
            if sec_sep < 60:

                # The total angle in seconds is computed
                sec = (fac * (deg * 3600.0 + mn * 60.0 + sec_sep)) / 3600.
                return sec

    raise ValueError('Illegal DMS field')


def _mod_lon(angle):
    '''Function to adjust a longitude angle to range from -pi to pi radians'''
    import math
    tau = 2. * math.pi
    while not -math.pi < angle < math.pi:
        angle -= tau * (-1. if angle < 0. else 1.)
    return angle


def _cal_to_jd(cal):
    '''
    Function to convert a given datetime date to a Julian day number
    Algorithm from (Press et al.)
    '''
    import math

    year, month, day = cal.year, cal.month, cal.day
    hour, minute, second = cal.hour, cal.minute, cal.second

    jy = year if month > 2 else year - 1
    jm = month + 1 if month > 2 else month + 13

    jdint = math.floor(math.floor(365.25 * jy) + math.floor(30.6001 * jm) + day + 1720995)

    ja = math.floor(0.01 * jy)
    jdint += 2 - ja + math.floor(0.25 * ja)

    # correct for half-day offset and set fraction of day
    frac = hour / 24.0 - 0.5
    if frac < 0.:
        frac += 1.
        jdint -= 1.

    jdfrac = frac + (minute + second / 60.0) / 60.0 / 24.0

    # round to nearest second
    jd0 = (jdint + jdfrac) * 100000
    jd = math.floor(jd0)
    if jd0 - jd > 0.5:
        jd += 1.

    return jd / 100000


def _jd_to_cal(jd):
    '''
    Function to convert a given Julian day number to a calendar date
    Algorithm from (Press et al.)
    '''
    import math

    jdint = math.floor(jd)
    jdfrac = jd % 1

    temp = math.floor(((jdint - 1867216) - 0.25) / 36524.25)
    j1 = jdint + 1 + temp - math.floor(0.25 * temp)

    j2 = j1 + 1524
    j3 = math.floor(6680. + ((j2 - 2439870) - 122.1) / 365.25)
    j4 = math.floor(j3 * 365.25)
    j5 = math.floor((j2 - j4) / 30.6001)

    dd = int(math.floor(j2 - j4 - math.floor(j5 * 30.6001)))
    mm = int(math.floor(j5 - 1) if ((j5 - 1) < 12) else math.floor(j5 - 1) - 12)
    yr = int(math.floor(j3 - 4715) if mm <= 2 else math.floor(j3 - 4716))

    maxdd = 31 if mm % 2 else 30 if mm != 2 else 28
    dd = dd if jdfrac < 0.5 else dd + 1
    if dd > maxdd:
        mm += 1
        dd = 1

    return (yr, mm, dd)


def _read_TLE_dates(filename, date):
    '''From TLE json file, return all TLE which surround this date'''
    import os
    import json
    import datetime

    if not os.path.exists(filename):
        directory, file = os.path.split(filename)
        get_TLE(filename[:-9], directory, 'Nolan.Dickson@canada.ca', 'trackthissatplease')

    tle_results = []
    with open(filename) as tf:
        try:
            tle_list = json.load(tf)
        except (ValueError, IOError):
            mssg = 'invalid TLE json file: {}\nPlease ensure the path is correct and the file was retrieved using the Space-Track api'.format(filename)
            raise ValueError(mssg)

        for entry in tle_list:

            epoch_date = datetime.datetime.strptime(entry['EPOCH'], '%Y-%m-%d %H:%M:%S')

            # if on date requested +- 1 day, save TLE
            if date - datetime.timedelta(days=1) <= epoch_date <= date + datetime.timedelta(days=1):
                tle_results.append((entry['TLE_LINE1'], entry['TLE_LINE2']))

            # TLE are in desecending order, if after this point, can break from file
            elif epoch_date < date - datetime.timedelta(days=1):
                break
    return tle_results


# ----------------------
# Orbit Calculations
# ----------------------

# TODO docstrings, Test, use in other scripts
# TODO switch time_frac to ['epoch'] read and datetime compare,
# TODO make all extensible (not use my paths and etc)


def get_TLE(sat, directory, username, password):
    '''get a json file containing all the TLE information for this satellite, and store it in the directory'''
    import datetime
    import requests
    import json
    import os

    with open('{ppp}/platforms.txt'.format(ppp=os.environ['PPP_CONFIG_DIR']), 'r') as f:
        satid = [line.split() for line in f if line.split()[0].lower() == sat.lower()]

    if not satid:
        import os
        mssg = 'satellite {sat} not included in platforms list at {ppp}/platforms.txt'.format(sat=sat, ppp=os.environ['PPP_CONFIG_DIR'])
        raise ValueError(mssg)

    filename = '{dir}/{sat}_TLE.json'.format(dir=directory, sat=satid[0])

    login_data = {'identity': username, 'password': password}
    login_url = 'https://www.space-track.org/ajaxauth/login'

    s = requests.Session()
    s.post(login_url, data=login_data)

    json_url = 'https://www.space-track.org/basicspacedata/query/class/tle/EPOCH/%3C{d:%Y-%m-%d}/NORAD_CAT_ID/{satid}/orderby/EPOCH%20DESC/format/json'.format(
        d=datetime.datetime.now().date(),
        satid=satid[1]
    )

    r = s.get(json_url)
    r.raise_for_status()

    with open(filename, 'w') as savefile:
        json.dump(r.json(), savefile)

    return


def haver_dist(long1, lat1, long2, lat2):
    '''
    Calculate distance in metres between
    two coordinates using the Haversine formula
    '''
    import math

    R = 6371000  # Radius of Earth, in metres
    phi_c = lat2 * math.pi / 180.
    phi_p = lat1 * math.pi / 180.
    del_phi = (lat1 - lat2) * math.pi / 180.
    del_lam = (long1 - long2) * math.pi / 180.

    a = math.sin(del_phi / 2.)**2 + math.cos(phi_p) \
        * math.cos(phi_c) * math.sin(del_lam / 2.)**2
    dist = R * 2 * math.atan2(a**0.5, (1 - a)**0.5)

    return dist


def width_to_degrees(distance, latitude):
    '''Given a width, in metres, and a latitude, convert to degrees of horizontal longitude'''
    import math
    return distance / (math.cos(latitude * math.pi / 180.) * 111320)


def drawswath(path, swath_width, swath_function=width_to_degrees):
    '''Given an orbit path, return a polygon of a determined width in metres around it. Swath function must accept (swath_width, latitude)'''
    import math
    # TODO the misr swath seemed a bit large
    # TODO this doesnt seem to work well for small widths

    # --------------------------------------------------------------------------
    # If swath width is None, create a swath which follows the orbit exactly,
    # but still passes all the way up and back down the path, for better use
    # in the polygon functions
    # --------------------------------------------------------------------------

    if swath_width is None:
        return path + list(reversed(path))

    # --------------------------------------------------------------------------
    # Loop through every point on the given path and grab the proceeding
    # point, to be used in calculations
    # --------------------------------------------------------------------------

    poly_left, poly_right = [], []
    for pind, point in enumerate(path[:-1]):
        x1, y1, time1 = point
        x2, y2, time2 = path[pind + 1]

        # --------------------------------------------------------------------------
        # Calculate distances and times between points
        # --------------------------------------------------------------------------

        dx = x2 - x1
        dy = y2 - y1
        d = math.sqrt((dx * dx) + (dy * dy)) / 2.

        poly_time = time1 + (time2 - time1) / 2

        # --------------------------------------------------------------------------
        # Calculate width of the swath to be drawn at this latitude
        # --------------------------------------------------------------------------

        # if being exact would use WGS84 ellipsoid equation but whatever
        try:
            r = swath_function(swath_width, y1 + dy)
        except TypeError:
            r = swath_width

        # --------------------------------------------------------------------------
        # Calculate the horizontal distance from the point to where the swath
        # point would sit if it were on the x-axis. Then calculate the angle from
        # the x-axis to the path line (phi) and the angle from the path line to
        # the location of the swath point (theta), which sits at a distance r
        # perpendicular from the middle of the path line
        # --------------------------------------------------------------------------

        xo = x1 + math.sqrt((r * r) + (d * d))
        theta = math.atan(r / d)
        phi = math.atan(dy / dx)

        # --------------------------------------------------------------------------
        # Calculate the angle that the swath point needs to be rotated off of the
        # x-axis by, for both the left and right swath points
        # --------------------------------------------------------------------------

        alpha_left = phi + theta if phi > 0. else phi - theta
        alpha_right = phi - theta if phi > 0. else phi + theta

        # --------------------------------------------------------------------------
        # Apply the alpha rotations, using a translated rotation matrix, to both
        # the left and right points
        # --------------------------------------------------------------------------

        poly_left.append((
            (((xo - x1) * math.cos(alpha_left)) - ((y1 - y1) * math.sin(alpha_left)) + x1),
            (((xo - x1) * math.sin(alpha_left)) + ((y1 - y1) * math.cos(alpha_left)) + y1),
            poly_time
        ))
        poly_right.append((
            (((xo - x1) * math.cos(alpha_right)) - ((y1 - y1) * math.sin(alpha_right)) + x1),
            (((xo - x1) * math.sin(alpha_right)) + ((y1 - y1) * math.cos(alpha_right)) + y1),
            poly_time
        ))

    # --------------------------------------------------------------------------
    # Flip the right swath array and append it to the end of the left swath
    # array, resulting in a polygon surrounding the path
    # --------------------------------------------------------------------------

    poly_right.reverse()
    swath_poly = poly_left + poly_right

    return swath_poly


def draworbit(date, sat='Terra', delta_t=30, direction='all', split_paths=True, tlefile=''):
    '''
    take in a date and return all the ascending Terra paths on that day seperately, and all the orbit tracks of Calipso on that day
    please give date as a datetime object with only the date (0 hrs, mins, etc)
    '''
    import os
    import datetime
    import pyorbital.orbital as porb

    # --------------------------------------------------------------------------
    # Argument checking
    # --------------------------------------------------------------------------

    if not isinstance(date, datetime.datetime):
        mssg = 'invalid date format, must be given as datetime'
        raise TypeError(mssg)

    if direction not in ['all', 'ascending', 'descending']:
        mssg = "invalid direction option, must be either 'all', 'ascending' or 'descending'"
        raise ValueError(mssg)

    tlefile = '{ppp}/{sat}_TLE.json'.format(ppp=os.environ['PPP_CONFIG_DIR'], sat=sat) if not tlefile else tlefile

    # --------------------------------------------------------------------------
    # Initialize containers and flags, read TLE's for this date from file
    # --------------------------------------------------------------------------

    lonlat_arr = []
    lonlat_paths = [[]]

    ascending = []
    descending = []

    last_nearest = False
    tle = _read_TLE_dates(tlefile, date)
    orbit = None

    # --------------------------------------------------------------------------
    # Simulate the orbit at the requested interval size for an entire day
    # --------------------------------------------------------------------------

    iterdate = date
    index = 0
    while date + datetime.timedelta(1) > iterdate:

        # --------------------------------------------------------------------------
        # From the list of TLE's on this date, find the one with the nearest
        # epoch to the iteration datetime. The TLE epoch is always at
        # position 18 to 32 of the first TLE line
        # --------------------------------------------------------------------------

        time_frac = ((iterdate.hour * 60. + iterdate.minute) * 60. + iterdate.second) / (24. * 60. * 60.)

        nearest_tle = min(tle, key=lambda line: abs(float(line[0][18:32]) - float(iterdate.strftime('%y%j.')) + time_frac))

        # --------------------------------------------------------------------------
        # If this tle is different from that used in the last iteration, create
        # a new orbital object using pyorbital. If it is the same, continue to
        # reuse the old orbit
        # --------------------------------------------------------------------------

        if nearest_tle != last_nearest:
            orbit = porb.Orbital(sat, line1=nearest_tle[0], line2=nearest_tle[1])

        # --------------------------------------------------------------------------
        # Use the SGP4 formulae embedded in pyorbital to calculate the lat/long
        # position of the satellite at this iteration datetime, and save it to
        # the full orbit array
        # --------------------------------------------------------------------------

        lon, lat, alt = orbit.get_lonlatalt(iterdate)
        lonlat_arr.append((lon, lat, iterdate))

        # --------------------------------------------------------------------------
        # If requested, determine if the point is on the ascending or descending
        # portion of the orbit. Must use seperate indices arrays and a later
        # deletion so that all points are available to be used in the next
        # iteration of this check.
        # `-2` index gives preceeding coordinate
        # --------------------------------------------------------------------------

        try:
            if direction == 'ascending' and lonlat_arr[-2][1] > lat:
                descending.append(index)
            elif direction == 'descending' and lonlat_arr[-2][1] < lat:
                ascending.append(index)
        except IndexError:
            pass

        # --------------------------------------------------------------------------
        # Iterate the orbit
        # --------------------------------------------------------------------------

        last_nearest = nearest_tle

        index += 1
        iterdate += datetime.timedelta(seconds=delta_t)

    # --------------------------------------------------------------------------
    # Delete all points whose direction opposes that specified in the
    # direction argument
    # --------------------------------------------------------------------------

    if direction == 'ascending':
        # delete all descending points
        for i in reversed(descending):
            del lonlat_arr[i]

    elif direction == 'descending':
        # delete all ascending points
        for i in reversed(ascending):
            del lonlat_arr[i]

    # --------------------------------------------------------------------------
    # If requested, split the entire orbit data into seperate paths, based on
    # the crossing of the antimeridian and poles. This is useful for avoiding
    # any problems associated with crossing those areas in polygon
    # checking, plotting, etc
    # --------------------------------------------------------------------------

    # TODO need to figure out a better way, I think this relies too much on terra's orbit orientation to be general

    if split_paths:

        path_index = 0
        lonlat_paths[path_index].append(lonlat_arr[0])

        for ind, (lon, lat, time) in enumerate(lonlat_arr[1:], 1):

            lonlat_paths[path_index].append((lon, lat, time))

            try:
                # if the next crosses the antimeridian/poles or switches directions, split it into a new 'path'

                if (lonlat_arr[ind - 1][1] > lat and  # descending and crossing outer border
                   ((lon < 0 and lonlat_arr[ind + 1][0] > 0) or (lat < 0 and lonlat_arr[ind + 1][1] > 0))):
                        path_index += 1
                        lonlat_paths.append([])

                elif (lonlat_arr[ind - 1][1] < lat and  # ascending and crossing outer border
                      ((lon < 0 and lonlat_arr[ind + 1][0] > 0) or (lat > 0 and lonlat_arr[ind + 1][1] < 0))):
                        path_index += 1
                        lonlat_paths.append([])

                # if this iteration and the one before were *not* in the same direction  # TODO -2 will break this a bit?
                elif not (((lonlat_arr[ind - 1][1] > lat) == (lonlat_arr[ind - 2][1] > lonlat_arr[ind - 1][1]))):
                    path_index += 1
                    lonlat_paths.append([])

            except IndexError:
                pass

        return lonlat_paths

    return lonlat_arr


def find_overpass(coord, daterange, swath_width=None, sat='Terra', direction='all', tlefile='', show=False):
    '''
    tuple of lon, lat; tuple of start datetime, end datetime
    return datetime of all coord/swath intercepts
    to a resolution of 30 seconds
    '''
    # TODO its a bit slow, allow coord to bbox
    import datetime
    from Cpoly.Cpoly import polygon_check_single, polygon_check

    start = daterange[0].date()
    end = daterange[-1].date()

    if start > end:
        raise ValueError('invalid date range')

    # TODO check that tuple coord isnt a tuple of bbox coordinates
    if not (isinstance(coord, tuple) or (isinstance(coord, list) and all([isinstance(pnt, tuple) for pnt in coord]))):
        raise TypeError('invalid coordinate values')

    if show:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if isinstance(coord, tuple):
            ax.plot(coord[0], coord[1], 'x', mew=4, ms=8)
        elif isinstance(coord, list):
            cx, cy = zip(*coord)
            ax.plot(cx, cy, 'x-')

    day = start
    overpasses = []
    while day <= end:
        path_list = draworbit(datetime.datetime(day.year, day.month, day.day), sat=sat, direction=direction, tlefile=tlefile)

        for path in path_list:
            swath = drawswath(path, swath_width)

            if swath:
                swath_coords = [(x[0], x[1]) for x in swath]

                if ((isinstance(coord, tuple) and polygon_check_single(swath_coords, coord, ordered=True)) or
                        (isinstance(coord, list) and polygon_check(swath_coords, coord, ordered=True))):

                    # TODO this way or the block split way?
                    # TODO maybe switch from center to a list of range of overpasses if bbox (this is a temp fix)
                    cntr = coord if isinstance(coord, tuple) else ((sum([p[0] for p in coord]) / len(coord), sum([p[1] for p in coord]) / len(coord)))

                    dist_arr = [haver_dist(sat_point[0], sat_point[1], cntr[0], cntr[1]) for sat_point in path]
                    near_index, _ = min(enumerate(dist_arr), key=lambda dist: dist[1])

                    if show:
                        ox, oy, _ = zip(*path)
                        sx, sy, _ = zip(*swath)
                        ax.plot(ox, oy, '-')
                        ax.fill(sx, sy, alpha=0.5)
                        ax.text(ox[near_index], oy[near_index], path[near_index][2].strftime('%H:%M:%S'))

                    overpasses.append(path[near_index][2])

        day += datetime.timedelta(days=1)

    if show:
        from mpl_toolkits.basemap import Basemap

        m = Basemap(projection='cyl', resolution='l')
        m.drawcoastlines(linewidth=0.5)
        m.drawparallels(list(range(-90, 90, 10)), labels=[True, False, False, True])
        m.drawmeridians(list(range(-180, 180, 10)), labels=[True, False, False, True])
        m.drawcountries(linewidth=0.5)

        plt.show()

    return overpasses


def find_position(srch_time, sat='Terra', tlefile='', show=False):
    '''datetime should be in UTC, i think'''
    # TODO show
    import os
    import datetime
    import pyorbital.orbital as porb

    tlefile = '{ppp}/{sat}_TLE.json'.format(ppp=os.environ['PPP_CONFIG_DIR'], sat=sat) if not tlefile else tlefile

    tle = _read_TLE_dates(tlefile, datetime.datetime(srch_time.year, srch_time.month, srch_time.day))

    time_frac = ((srch_time.hour * 60. + srch_time.minute) * 60. + srch_time.second) / (24. * 60. * 60.)
    nearest_tle = min(tle, key=lambda line: abs(float(line[0][18:32]) - float(srch_time.strftime('%y%j.')) + time_frac))

    orbit = porb.Orbital(sat, line1=nearest_tle[0], line2=nearest_tle[1])

    lon, lat, alt = orbit.get_lonlatalt(srch_time)

    return (lon, lat)


def find_intercept(daterange, sat1, sat2, swath_width1=None, swath_width2=None, direction1='all', direction2='all', latrange=(-70, 70), show=False):
    import datetime
    from Cpoly.Cpoly import polygon_check

    raise NotImplementedError()

    start = daterange[0].date()
    end = daterange[-1].date()

    if start > end:
        raise ValueError('invalid date range')

    if show:
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(111)

    day = start
    intercepts = []
    while day <= end:

        orbit1 = draworbit(datetime.datetime(day.year, day.month, day.day), sat=sat1, direction=direction1)
        swath_list1 = [drawswath(path, swath_width1) for path in orbit1] if swath_width1 else orbit1

        orbit2 = draworbit(datetime.datetime(day.year, day.month, day.day), sat=sat2, direction=direction2)
        swath_list2 = [drawswath(path, swath_width2) for path in orbit2] if swath_width2 else orbit2

        # TODO there must be a faster way, I dont think the checking will work well for two lines (will try to connect ends?)
        # TODO probably have to do the block split way? that will be V slow
        for path1 in swath_list1:
            for path2 in swath_list2:
                if polygon_check(path1, path2, ordered=True):
                    # TODO how to implement latrange (else they will all intercept)
                    intercepts.append()

        day += datetime.timedelta(days=1)

    return intercepts

# ----------------------
# Reverse Conversions
# ----------------------

# TODO make sure these conversions are all correct, they are slightly off, based on CALIOP_intersect images


def BLStoSOM(block_num, line_num, sample_num, res=275.):
    '''*Convert a Block Line Sample index to a Space Oblique Mercator x and y coordinate*

    Given the 3D index of a MISR datapoint as it's block, line and sample indices, convert
    to the x and y coordinates of a SOM projection. This method is valid for all MISR terrain
    files, as long as the correct resolution is given, in metres (see DPS table 6-22).

    This conversion is purely geometric, based on the size of the data and its position in the orbit.
    The number of BLS indices can be seen as `data.shape -> (nblocks, nlines, nsamples)`

    Note
    ----

    SOM x and y axes are inverted from what you might normally expect, that is: N-S -> x & E-W -> y

    Args
    ----
    block_num : int
        0th axis index representing the block number.
    line_num : int
        1st axis index representing the line number.
    sample_num : int
        2nd axis index representing the sample number.
    res : float, optional
        The resolution of the dataset being used. For MISR terrain files, this is 275 for all bands of the AN file.
        See DPS table 6-22 for other camera and band combinations.

    Returns
    -------
    som_x : float
        x-coordinate of the requested index on the SOM grid.
    som_y : float
        y-coordinate of the requested index on the SOM grid.

    '''

    # TODO convert to doing for an entire array, not a point?
    # TODO test to make sure its not off, seems slightly pushed towards meridian

    res_factor = 1100. / res

    # project parameters
    ulc_coord = [7460750.0, 1090650.0]  # from pp general init
    lrc_coord = [7601550.0, 527450.0]

    ulc = [ulc_coord[0], lrc_coord[1]]  # ulc y and lrc y are reversed in the structural metadata
    lrc = [lrc_coord[0], ulc_coord[1]]

    nlines = 128 * res_factor
    nsamples = 512 * res_factor
    rel_offset = rel_offset_table

    # size of pixels
    size_x = (lrc[0] - ulc[0]) / nlines
    size_y = (lrc[1] - ulc[1]) / nsamples

    # adjust to center of pixels
    x_cent = ulc[0] + size_x / 2.0
    y_cent = ulc[1] + size_y / 2.0

    # calculate absolue offsets for each block (keep calcing whole array for use after I change to doing on array)
    abs_offset = [0.]
    for rind, rel in enumerate(rel_offset):
        abs_offset.append(abs_offset[rind - 1] + rel)

    # calculate SOM x and y coordinates
    n = (block_num - 1) * nlines * size_x
    som_x = x_cent + n + (line_num * size_x)
    som_y = y_cent + (sample_num + abs_offset[block_num - 1]) * size_y

    return (som_x, som_y)


def SOMtoLonLat(path_num, x, y):
    '''*Convert a Space Oblique Mercator x and y coordinate to a geographic longitude and latitude*

    Inverse SOM conversion function for converting between a coordinate point on the SOM grid of
    a certain MISR path to a geographic coordinate point, in degrees.

    Math seen here adapted from GCTP and Space Oblique Mercator Projection Mathematical Development (John Snyder, 1981)
    These equations require certain physical constants be known, which are thus hard-coded to match MISR and Terra, and will
    not function for any other satellite or instrument.

    Args
    ----
    path_num : int
        The MISR path number of this coordinate.
    x : float
        x-coordinate of the requested point on the SOM grid.
    y : float
        y-coordinate of the requested point on the SOM grid.

    Returns
    -------
    lon : float
        Longitude of the requested point, in degrees.
    lat : float
        Latitude of the requested point, in degrees.

    Todo
    ----

    Documentation explaining the math.

    '''
    import math

    # get initial parameters for MISR SOM conversions

    a, b, es, p21, t, w, q, xj, a2, a4, c1, c3, ca, sa, u, lon_center = _som_init(path_num)

    # SOMinv

    # Inverse equations. Begin inverse computation with approximation for tlon.
    # Solve for transformed longitude (tlon)

    tlon = x / (a * b)

    conv_check = False

    for inumb in range(50):
        sd = math.sin(tlon)
        sdsq = sd * sd

        s = p21 * sa * math.cos(tlon) * math.sqrt((1.0 + t * sdsq) / ((1.0 + w * sdsq) * (1.0 + q * sdsq)))
        blon = (x / a) + (y / a) * s / xj - a2 * math.sin(2.0 * tlon) - a4 * math.sin(4.0 * tlon) - (s / xj) * (c1 * math.sin(tlon) + c3 * math.sin(3.0 * tlon))

        next_tlon = blon / b

        if abs(next_tlon - tlon) < 1.e-9:
            conv_check = True
            break

        tlon = next_tlon

    if not conv_check:
        raise ValueError("50 iterations without convergence")

    # Compute transformed latitude (tlat)

    st = math.sin(tlon)
    defac = math.exp(math.sqrt(1.0 + s * s / xj / xj) * (y / a - c1 * st - c3 * math.sin(3.0 * tlon)))
    actan = math.atan(defac)
    tlat = 2.0 * (actan - (math.pi / 4.0))

    # Compute geodetic longitude (xlam_t)

    if abs(math.cos(tlon) < 1.e-7):
        tlon = tlon - 1.e-7

    dd = st * st
    bigk = math.sin(tlat)
    bigk2 = bigk * bigk

    xlam_t = math.atan(((1.0 - bigk2 / (1.0 - es)) * math.tan(tlon) * ca - bigk * sa *
                        math.sqrt((1.0 + q * dd) * (1.0 - bigk2) - bigk2 * u) / math.cos(tlon)) / (1.0 - bigk2 * (1.0 + u)))

    # Correct inverse quadrant

    sl = 1. if xlam_t >= 0. else -1.
    scl = 1. if math.cos(tlon) >= 0. else -1.

    xlam_t -= ((math.pi / 2.0) * (1.0 - scl) * sl)
    dlon = xlam_t - p21 * tlon

    # Compute geodetic latitude

    if abs(sa) < 1.e-7:
        dlat = math.asin(bigk / math.sqrt((1.0 - es) * (1.0 - es) + es * bigk2))
    else:
        dlat = math.atan((math.tan(tlon) * math.cos(xlam_t) - ca * math.sin(xlam_t)) / ((1.0 - es) * sa))

    # adjust and convert to degrees
    lon = _mod_lon(dlon + lon_center) * (180. / math.pi)
    lat = dlat * (180. / math.pi)

    return (lon, lat)


def BLStoLonLat(path_num, block_num, line_num, sample_num, res=275.):
    '''*Convert a Block Line Sample index to a geographic longitude and latitude*

    Note
    ----

    Conversion is first made from BLStoSOM and then SOMtoLonLat. See the respective function documentation for more details.

    '''
    # TODO its slightly off towards the meridian I believe
    som_x, som_y = BLStoSOM(block_num, line_num, sample_num, res)
    lon, lat = SOMtoLonLat(path_num, som_x, som_y)
    return lon, lat


# ----------------------
# Forward Conversions
# ----------------------

# TODO documentation, testing, fix conversions

def SOMtoBLS(som_x, som_y, res=275.):
    import math

    res_factor = 1100. / res

    # project parameters
    # ulc_coord = [7460750.0, 1090650.0]  # from pp general init
    # lrc_coord = [7601550.0, 527450.0]

    # ulc = [ulc_coord[0], lrc_coord[1]]  # ulc y and lrc y are reversed in the structural metadata
    # lrc = [lrc_coord[0], ulc_coord[1]]
    ulc = [7460750.0, 527450.0]

    nlines = 128 * res_factor
    # nsamples = 512 * res_factor
    rel_offset = rel_offset_table

    # size of pixels
    # size_x = (lrc[0] - ulc[0]) / nlines
    # size_y = (lrc[1] - ulc[1]) / nsamples
    size_x = res
    size_y = res  # TODO

    # adjust to center of pixels
    x_cent = ulc[0] + size_x / 2.0
    y_cent = ulc[1] + size_y / 2.0

    # calculate absolue offsets for each block (keep calcing whole array for use after I change to doing on array)?
    abs_offset = [0.]
    for rind, rel in enumerate(rel_offset[1:]):
        abs_offset.append(abs_offset[rind - 1] + rel)

    # calculate SOM x and y coordinates
    x_int = (som_x - x_cent) / size_x
    y_int = (som_y - y_cent) / size_y

    block_num = int(math.floor((x_int + 0.5) / nlines) + 1)
    line_num = int(x_int - ((block_num - 1) * (nlines)))
    sample_num = int(y_int - abs_offset[block_num - 1])

    return (block_num, line_num, sample_num)


def LonLattoSOM(path_num, lon, lat):
    import math

    # get initial parameters for MISR SOM conversions

    a, b, es, p21, t, w, q, xj, a2, a4, c1, c3, ca, sa, u, lon_center = _som_init(path_num)

    # convert from degrees to rad?

    radlon = lon * (math.pi / 180.) - lon_center
    radlat = lat * (math.pi / 180.)

    tlamp = 1.5 * math.pi if lat < 0. else math.pi / 2.

    n = 0
    while True:

        # L230

        sav = tlamp
        xlamp = radlon + p21 * tlamp
        ab1 = math.cos(xlamp)
        if abs(ab1) < 1.e-7:
            xlamp -= 1.e-7
        scl = -1. if ab1 < 0. else 1.
        ab2 = tlamp - (scl) * math.sin(tlamp) * (math.pi / 2.)

        # L240

        conv_check = False

        for inumb in range(50):
            xlamt = radlon + p21 * sav
            c = math.cos(xlamt)
            if abs(c) < 1.e-7:
                xlamt -= 1.e-7
            xlam = (((1.0 - es) * math.tan(radlat) * sa) + math.sin(xlamt) * ca) / c
            tlam = math.atan(xlam)
            tlam = tlam + ab2
            tabs = abs(sav) - abs(tlam)

            if abs(tabs) < 1.e-7:
                conv_check = True
                break

            sav = tlam

        if not conv_check:
            raise ValueError("50 iterations without convergence")

        # L250

        rlm = math.pi * 0.5201613
        rlm2 = rlm + 2. * math.pi
        n += 1

        if n >= 3 or rlm < tlam < rlm2:
            break
        elif rlm < tlam:
            tlamp = 2.5 * math.pi
        elif tlam >= rlm2:
            tlamp = math.pi / 2.

    dp = math.sin(radlat)
    tphi = math.asin(((1.0 - es) * ca * dp - sa * math.cos(radlat) * math.sin(xlamt)) / math.sqrt(1.0 - es * dp * dp))

    xtan = (math.pi / 4.) + (tphi / 2.)
    tanlg = math.log(math.tan(xtan))
    sd = math.sin(tlam)
    sdsq = sd * sd
    s = p21 * sa * math.cos(tlam) * math.sqrt((1.0 + t * sdsq) / ((1.0 + w * sdsq) * (1.0 + q * sdsq)))
    d = math.sqrt(xj * xj + s * s)

    x = b * tlam + a2 * math.sin(2. * tlam) + a4 * math.sin(4. * tlam) - tanlg * s / d
    x = a * x
    y = c1 * sd + c3 * math.sin(3. * tlam) + tanlg * xj / d
    y = a * y

    # x, y = y, x  # make sure this is right
    return (x, y)


def LonLattoBLS(path_num, lon, lat, res=275.):
    som_x, som_y = LonLattoSOM(path_num, lon, lat)
    B, L, S = SOMtoBLS(som_x, som_y, res=res)

    return (B, L, S)


# ----------------------
# Imaging
# ----------------------


def drawmisr(MISRfile, startblock, endblock=None, save=True, filepath='', band='all'):
    '''*Given a MISR terrain file, create a true-colour image of the requested blocks*

    Given the path to a MISR terrain file, will, if requested, create and save a true colour image of the requested MISR blocks.
    A dictionary containing the reflectance of the requested bands will be returned for use in other functions.
    Will also return the geographic lat/long coordinates of the upper left hand and lower right hand corners of the image.

    The radiance scale factor cannot be read properly from python, but as it is consistenly ~0.047203224, this value is used.
    This means that this is *only* for creating images, not for analysis. *These radiance numbers could be wrong.*

    If you wish to create your own real colour image from the returned BRF, you can do so as follows:

    >>> BRF, *_ = MISRutils.drawmisr(filename, 50, 55, save=False)
    >>> rgba = np.ma.dstack([BRF['Red'], BRF['Green'], BRF['Blue'], ~BRF['NIR'].mask])
    >>> plt.imshow(rgba)

    Note
    ----

    There is a limit to the size of image matplotlib can create, which works out to be about 80 blocks. Do not request more than this at one time.

    No image will be created for blocks over the ocean, as MISR terrain files do not store any such information.

    Args
    ----
    MISRfile : str
        Path of the MISR terrain file.
    startblock : int
        Block number of first block in requested range.
    endblock : int, optional
        Block number of last block in requested range. If not specified, only the startblock will be used.
    save : bool, optional
        If specified, will save a true-colour image.
    filepath : str, optional
        Name or path to which you wish to save the image. `save` arg must be specified in order to use `filepath`.
    band : str or list of str, optional
        The name of the band(s) you wish to create an image using. The default 'all' will result in an RGBA image.

    Returns
    -------
    BRF : dict
        Dictionary containing the masked Bidirectional Radiance Factor arrays for each colour layer.
    ulc : tuple
        Tuple containing the longitude and latitude, in degrees, of the upper left corner of the first block
    lrc : tuple
        Tuple containing the longitude and latitude, in degrees, of the lower right corner of the last block
    (Will save image if requested)

    '''
    import os
    import numpy as np
    import pyhdf.SD as SD
    import matplotlib.image as mpimg

    # TODO fix bright spots, seems kinda overly red?
    # TODO extend to work with bbox
    # TODO default filepath seems weird
    # TODO hand picking bands creates *alot* of bugs

    # arg checking
    if not os.path.exists(MISRfile):
        mssg = 'no such MISR file: {file}'.format(file=MISRfile)
        raise IOError(mssg)
    if filepath and not os.path.exists(os.path.dirname(filepath)):
        mssg = 'no such output directory: {dir}'.format(dir=os.path.dirname(filepath))
        raise IOError(mssg)
    if endblock and endblock + 1 < startblock:
        mssg = 'startblock must be greater than endblock'
        raise ValueError(mssg)
    if (isinstance(band, str) and band not in ['all', 'red', 'blue', 'green', 'nir']) or \
            (isinstance(band, list) and not all([b.lower() in ['all', 'red', 'blue', 'green', 'nir'] for b in band])):
        mssg = 'invalid band option: {band}'.format(band=band)
        raise ValueError(mssg)

    endblock = endblock + 1 if endblock else startblock + 1
    filepath = filepath if filepath else '{dir}/{hdf}_B{sb:03d}_B{eb:03d}.png'.format(dir=os.getcwd(), hdf=MISRfile[:-13], sb=startblock, eb=endblock - 1)

    # open terrain file
    hdf = SD.SD(MISRfile, SD.SDC.READ)
    path_num = hdf.attributes()['Path_number']

    # Get lat and lon of corners of the blocks requested
    # can use the extent of a 275m res file (512, 2048) no matter what, as long as we keep the default res kwarg at 275
    ulc = BLStoLonLat(path_num, startblock, 0, 0)
    lrc = BLStoLonLat(path_num, endblock, 512, 2048)

    dataset_names = ['Blue Radiance/RDQI', 'Green Radiance/RDQI', 'Red Radiance/RDQI', 'NIR Radiance/RDQI']

    BRF = {}

    for dataset in dataset_names:

        colour = dataset.split()[0]

        # Check if this band was requested, need 'isinstance' to avoid error on band.lower
        # TODO this lower() stuff is getting excessice and breaking
        if (isinstance(band, str) and (band.lower() != 'all' and band.lower() != colour.lower())) or (isinstance(band, list) and colour not in band):
            continue

        # Read dataset. slice out the blocks requested
        RDQI_data = hdf.select(dataset)[startblock - 1:endblock - 1, :, :]

        # Set scale factor and fill value manually
        scale_factor = 0.047203224
        fill_value = 16377

        # We need to shift bits for "RDQI" to get "Blue Band" only, which is scaled up.
        scaled_rad_data = np.right_shift(RDQI_data, 2).astype(np.double)

        # Apply all the fill values. (in MISR terrain files all shifted values above 16377 are fill, see DPS table 6-23 for specifics)
        scaled_rad_data[scaled_rad_data >= fill_value] = np.nan

        scaled_rad_mask = np.ma.masked_array(scaled_rad_data, mask=np.isnan(scaled_rad_data))

        # Apply scale factor to get actual radiances.
        rad_mask = scale_factor * scaled_rad_mask

        # Apply BRF conversion factors to get Bidirectional Radiance Factor
        # (BRF factors are on a 17.6 km resolution, therefore 1 for each 64 or 16 MISR pixels)

        BRF_factors = hdf.select('{clr}ConversionFactor'.format(clr=colour))[startblock - 1:endblock - 1, :, :]

        BRF_factors = BRF_factors.repeat(rad_mask.shape[1] / BRF_factors.shape[1], 1).repeat(rad_mask.shape[2] / BRF_factors.shape[2], 2)

        BRF_data = rad_mask * BRF_factors

        # adjust any 1.1km resolution bands to the 275m resolution used by the red band in all cameras
        if BRF_data.shape[1] < 512:
            # these fractions should work out to 4 each time
            BRF_data = BRF_data.repeat(512 / BRF_data.shape[1], 1).repeat(2048 / BRF_data.shape[2], 2)

        # adjust for the offset between blocks
        # extra space must be created to absorb the adjustment (not ideal)
        tot_offset = abs(int(sum(rel_offset_table[startblock - 2: endblock - 2]))) * 4
        BRF_data = np.ma.concatenate((np.ma.masked_all((BRF_data.shape[0], BRF_data.shape[1], tot_offset)), BRF_data), axis=2)

        # get relative offset from table and roll data to eliminate
        rel_offset = 0
        for bind, block in enumerate(BRF_data):

            rel_offset += rel_offset_table[startblock - 2 + bind]
            BRF_data[bind] = np.roll(block, int(4 * rel_offset), axis=1)  # TODO not always 4, need to check res_factor

        # reshape to align all blocks in the same array (i.e. squish the 3D stacked blocks into one long 2D array)
        BRF_data = BRF_data.reshape(BRF_data.shape[0] * BRF_data.shape[1], BRF_data.shape[2])

        # save the gamma corrected (brightened) brf values
        BRF[colour] = np.ma.sqrt(BRF_data)

    # stack all bands with each other (RGB) and plot
    # TODO alpha with other band choices
    alpha = ~BRF['NIR'].mask

    rgb_arr = [BRF['Red'], BRF['Green'], BRF['Blue'], alpha] if band == 'all' else \
              [BRF[bn.capitalize()] for bn in band] if isinstance(band, list) else [BRF[band.capitalize()]]

    rgb = np.ma.dstack(rgb_arr)

    if save:
        try:
            mpimg.imsave(filepath, rgb.astype(np.float64))
        except ValueError:
            mssg = "matplotlib image creation size reached, please use a smaller number of blocks (<70)"
            raise ValueError(mssg)

    # TODO are ulc, lrc even useful? its too shifted to be a horizontal/vertical rectangle on a lat long grid
    # TODO figure out best way to convey that the data is offset and that needs to be accounted for, and how to do so
    return BRF, ulc, lrc, rel_offset * 4


def drawbands(MISRfile, startblock, endblock=None):
    '''Draw all four MISR bands on seperate plots for visual comparison'''
    import matplotlib.pyplot as plt

    brf, ulc, lrc, off = drawmisr(MISRfile, startblock, endblock, save=False)

    fig, axes = plt.subplots(2, 2, sharex='col', sharey='row')

    for bind, band in enumerate(['Red', 'Blue', 'Green', 'NIR']):
        axes[bind % 2, bind // 2].imshow(brf[band])
        axes[bind % 2, bind // 2].set_title(band)

    plt.show()


# ----------------------
# Orbit attributes
# ----------------------


def TimetoOrbit(dtime):
    '''given a datetime, return a corresponding MISR orbit number'''

    orb_incr = 1000

    caltime = _cal_to_jd(dtime)

    ref_num = int(((233.0 / 16.0 * (caltime - orbit_table[1])) + orb_incr) / orb_incr)
    num_refs = len(orbit_table)
    orb_ind = ref_num if ref_num < num_refs else num_refs - 1

    orbit = ((233.0 / 16.0 * (caltime - orbit_table[orb_ind])) + ((orb_ind) * orb_incr))

    return int(orbit)


def OrbittoTime(orb_num):
    '''given an orbit number, return a corresponding datetime'''
    # TODO return time not just date
    import datetime

    orb_incr = 1000

    # Select nearest orbit in table. This table ends at orbit 100 000 and thus
    # a new table will be necessary after this time to account for drift

    ref_num = int(orb_num / orb_incr)
    num_refs = len(orbit_table)
    orb_ind = ref_num if ref_num < num_refs else num_refs - 1

    # Calculate days per orbit using the delta time and delta orbit, only using non-extrapolated orbit tables

    del_time = orbit_table[80] - orbit_table[1]
    del_orbit = (80 - 1) * orb_incr
    days_per_orbit = (16.0 / 233.0 + del_time / del_orbit) / 2.0

    # Calculate the julian time from the table reference + the days since
    # the last table increment and convert to a calendar date

    juliantime = orbit_table[orb_ind] + days_per_orbit * (orb_num - orb_ind * orb_incr)

    date = datetime.datetime.strptime('{0:02d}.{1:02d}.{2:02d}'.format(*_jd_to_cal(juliantime)), '%Y.%m.%d')

    return date


def OrbittoPath(orb_num):
    '''given an orbit number, return the corresponding path'''
    return int((176 + orb_num * 16) % 233 + 1)
