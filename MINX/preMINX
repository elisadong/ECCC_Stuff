#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2018-02-19 13:54:21
# Contact:    Nolan.Dickson@Canada.ca


import re
import os
import sys
import csv
import math
import glob
import getpass
import requests
import argparse
import tempfile
import datetime
import subprocess

import pyloader
import postMINX as pm
import MISRutils as mu


# ******************************************************************************
def polygon_check_single(polypoints, point, ccwsorted=False):
# ******************************************************************************
# Function for checking if a single point lies within a given N-sided
# convex polygon. (Concave polygons not possible without prior info,
# please sort concave polygons before passing points to function)
# Input: list of N tuples (lon, lat) of polygon corners
#        Tuple (lon, lat) of point to check
#        Bool, wether polygon has been pre-sorted counterclockwise, from pi rad
# Output: Bool True if point lies within polygon else False
# ------------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Ensure that arguments are all floats. try/except to check that `point`
    # is a single point, not a list of points.
    # --------------------------------------------------------------------------

    try:
        polypoints = [tuple(float(elem) for elem in coord) for coord in polypoints]
        point = tuple(float(elem) for elem in point)
    except TypeError:
        mssg = "`point` argument must be a single coordinate tuple. See polygon_check() for use of a bounding box."
        raise TypeError(mssg)

    # --------------------------------------------------------------------------
    # Crop area obviously outside of polygon and check if point is directly on
    # a corner of the polygon
    # --------------------------------------------------------------------------

    leftP, rightP = sorted(polypoints, key=lambda tup: tup[0])[::len(polypoints) - 1]
    bottomP, topP = sorted(polypoints, key=lambda tup: tup[1])[::len(polypoints) - 1]

    if not (leftP[0] < point[0] < rightP[0]) or not (bottomP[1] < point[1] < topP[1]):
        return False
    if point in polypoints:
        return True

    # --------------------------------------------------------------------------
    # Sort polygon points ccw, starting at first point after pi rads, if not
    # already sorted, and create a linear function (y=mx+b) for each
    # polygon edge.
    # --------------------------------------------------------------------------

    if not ccwsorted:
        cent = (sum([p[0] for p in polypoints]) / len(polypoints), sum([p[1] for p in polypoints]) / len(polypoints))
        polypoints.sort(key=lambda p: math.atan2(p[1] - cent[1], p[0] - cent[0]))

    polylines = [{
        'm': (p[1] - polypoints[i - 1][1]) / (p[0] - polypoints[i - 1][0]) if (p[0] - polypoints[i - 1][0]) else None,
        'b': p[1] - ((p[1] - polypoints[i - 1][1]) / (p[0] - polypoints[i - 1][0])) * p[0] if (p[0] - polypoints[i - 1][0]) else None,
        'xrange': sorted([p[0], polypoints[i - 1][0]]),
        'yrange': None if (p[0] - polypoints[i - 1][0]) else sorted([p[1], polypoints[i - 1][1]])  # Only care about for vertical lines
    } for i, p in enumerate(polypoints)]

    # --------------------------------------------------------------------------
    # Find value of polygon lines at x value of point, if that point lies
    # within x range of the line.
    # Vertical lines are not relevant for final x-int check, so must check
    # explicitly here
    # --------------------------------------------------------------------------

    results = []
    for line in polylines:
        if line['yrange']:
            if point[0] == line['xrange'][0] and line['yrange'][0] < point[1] < line['yrange'][1]:
                return True
            else:
                continue

        elif line['xrange'][0] <= point[0] <= line['xrange'][1]:
            results.append(line['m'] * point[0] + line['b'])

    # --------------------------------------------------------------------------
    # Return False if either both or no values are above the point
    # If only one value is above, point must be between lines, and in polygon
    # --------------------------------------------------------------------------

    check = [result > point[1] for result in results]

    return False if (all(check) or not any(check)) else True


# ******************************************************************************
def polygon_check(polypoints, bbox, ccwsorted=False):
# ******************************************************************************
# Function for checking if any point within a box with given boundaries lies
# also within a given N-sided, convex polygon. (Concave polygons not possible
# without prior info, please sort concave polygons before passing
# points to function (also works with lines apparently)
# Input: list of N tuples (lon, lat) of polygon corners
#        list of N tuples (lon, lat) of box corners
#        Bool, wether polygon has been pre-sorted counterclockwise, from pi rad
# Output: Bool True if any point in box lies within polygon else False
# ------------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Ensure that arguments are all floats. try/except to check that `bbox`
    # is a list of points, not a single point
    # --------------------------------------------------------------------------

    try:
        polypoints = [tuple(float(elem) for elem in coord) for coord in polypoints]
        bbox = [tuple(float(elem) for elem in coord) for coord in bbox]
    except TypeError:
        mssg = "`bbox` argument must be a list of coordinate tuples. See polygon_check_single() for use of a single point."
        raise TypeError(mssg)

    # --------------------------------------------------------------------------
    # Crop area obviously outside of polygon, check if any box corner is
    # directly on a corner of the polygon and check if box completely
    # surrounds polygon
    # --------------------------------------------------------------------------

    leftP, rightP = sorted(polypoints, key=lambda tup: tup[0])[::len(polypoints) - 1]
    bottomP, topP = sorted(polypoints, key=lambda tup: tup[1])[::len(polypoints) - 1]

    leftBox, rightBox = sorted(bbox, key=lambda tup: tup[0])[::len(bbox) - 1]
    bottomBox, topBox = sorted(bbox, key=lambda tup: tup[1])[::len(bbox) - 1]

    if (topBox[1] < bottomP[1]) or (bottomBox[1] > topP[1]) or (rightBox[0] < leftP[0]) or (leftBox[0] > rightP[0]):
        return False
    elif any(b in polypoints for b in bbox):
        return True
    elif (topBox[1] > topP[1]) and (bottomBox[1] < bottomP[1]) and (rightBox[0] > rightP[0]) and (leftBox[0] < leftP[0]):
        return True

    else:
        # --------------------------------------------------------------------------
        # Sort polygon points ccw, if not already sorted, and create a linear
        # function (y=mx+b) for each polygon and box edge.
        # If vertical line, must check yrange, cannot rely on intercept like normal
        # --------------------------------------------------------------------------

        # TODO stop repetition?, I dont think it NEEDS to be CCW sorted per se anymore, jsut sorted

        if not ccwsorted:
            cent = (sum([p[0] for p in polypoints]) / len(polypoints), sum([p[1] for p in polypoints]) / len(polypoints))
            polypoints.sort(key=lambda p: math.atan2(p[1] - cent[1], p[0] - cent[0]))

            centBox = (sum([b[0] for b in bbox]) / len(bbox), sum([b[1] for b in bbox]) / len(bbox))
            bbox.sort(key=lambda b: math.atan2(b[1] - centBox[1], b[0] - centBox[0]))

        polylines = [{
            'm': (p[1] - polypoints[i - 1][1]) / (p[0] - polypoints[i - 1][0]) if (p[0] - polypoints[i - 1][0]) else None,
            'b': p[1] - ((p[1] - polypoints[i - 1][1]) / (p[0] - polypoints[i - 1][0])) * p[0] if (p[0] - polypoints[i - 1][0]) else None,
            'xrange': sorted([p[0], polypoints[i - 1][0]]),
            'yrange': None if (p[0] - polypoints[i - 1][0]) else sorted([p[1], polypoints[i - 1][1]])  # Only care about for vertical lines
        } for i, p in enumerate(polypoints)]

        boxlines = [{
            'm': (b[1] - bbox[i - 1][1]) / (b[0] - bbox[i - 1][0]) if (b[0] - bbox[i - 1][0]) else None,
            'b': b[1] - ((b[1] - bbox[i - 1][1]) / (b[0] - bbox[i - 1][0])) * b[0] if (b[0] - bbox[i - 1][0]) else None,
            'xrange': sorted([b[0], bbox[i - 1][0]]),
            'yrange': None if (b[0] - bbox[i - 1][0]) else sorted([b[1], bbox[i - 1][1]])  # Only care about for vertical lines
        } for i, b in enumerate(bbox)]

        # --------------------------------------------------------------------------
        # Check the x intercept between each box and polygon point.
        # Return True if any x-int lies within x range of the line.
        # If line is vertical, return true if y-value of box at x is within y
        # range of the line
        # --------------------------------------------------------------------------

        for poly in polylines:
            for box in boxlines:

                if poly['m'] == box['m']:

                    if poly['m'] is None and poly['xrange'][0] == box['xrange'][0]:
                        return True  # both vertical, with same x
                    else:
                        continue  # parallel, don't intersect

                else:
                    # TODO compress this vertical stuff
                    if poly['yrange']:
                        # vertical polygon line
                        x_int = poly['xrange'][0]

                        if (poly['yrange'][0] <= box['m'] * x_int + box['b'] <= poly['yrange'][1]) and (box['xrange'][0] <= x_int <= box['xrange'][1]):
                            return True

                    elif box['yrange']:
                        # vertical box line
                        x_int = box['xrange'][0]

                        if (box['yrange'][0] <= poly['m'] * x_int + poly['b'] <= box['yrange'][1]) and (poly['xrange'][0] <= x_int <= poly['xrange'][1]):
                            return True

                    else:
                        x_int = (box['b'] - poly['b']) / (poly['m'] - box['m'])

                        if (poly['xrange'][0] <= x_int <= poly['xrange'][1]) and (box['xrange'][0] <= x_int <= box['xrange'][1]):
                            return True

    # --------------------------------------------------------------------------
    # Check if a point on the box (top chosen randomly) lies within the
    # polygon, using single point check.
    # Due to lines not intersecting, if any point is within, all points must
    # lie within
    # --------------------------------------------------------------------------

    return polygon_check_single(polypoints, topBox, ccwsorted)


# --------------------------------------------------------------------------
# If using python3.x, import faster polygon checking from Cpoly instead
# --------------------------------------------------------------------------

try:
    from Cpoly import polygon_check, polygon_check_single
except ImportError:
    # if running on python2, use python implementation of polygon checking from preMINX
    pass


# ******************************************************************************
def jd_to_cal(jd):
# ******************************************************************************
# Function to convert a given Julian day number to a calendar date
# Input: Float Julian day number
# Output: Tuple (yr, mm, dd) of calendar date year, month and day
# ------------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Julian number decimals only relevant for time, work only with integers
    # and check if decimals is > 0.5. If so, add a day, due to julian numbers
    # beginning at noon UTC
    # --------------------------------------------------------------------------

    jdint = math.floor(jd)
    jddec = jd % 1

    # --------------------------------------------------------------------------
    # Correction for switch to Gregorian calendar in 1582
    # --------------------------------------------------------------------------

    if jdint >= 2299161.:
        temp = math.floor(((jdint - 1867216) - 0.25) / 36524.25)
        j1 = jdint + 1 + temp - math.floor(0.25 * temp)
    else:
        j1 = jdint

    # --------------------------------------------------------------------------
    # Julian conversion calculations
    # --------------------------------------------------------------------------

    j2 = j1 + 1524
    j3 = math.floor(6680. + ((j2 - 2439870) - 122.1) / 365.25)
    j4 = math.floor(j3 * 365.25)
    j5 = math.floor((j2 - j4) / 30.6001)

    dd = int(math.floor(j2 - j4 - math.floor(j5 * 30.6001)))
    mm = int(math.floor(j5 - 1) if ((j5 - 1) < 12) else math.floor(j5 - 1) - 12)
    yr = int(math.floor(j3 - 4715) if mm <= 2 else math.floor(j3 - 4716))

    return (yr, mm, dd if jddec < 0.5 else dd + 1)


# ******************************************************************************
def read_TLE_dates(filename, j_date):
# ******************************************************************************
# Function to read all Two Line Element blocks for a certain date from a file
# containing multiple elements
# Input: String path to TLE file
#        Int (YYDDD) of date to search for
# Output: list of TLE elements
# -------------------------------------------------------------------------------

    tle_list = []
    with open(filename) as tf:
        for line in tf:
            epoch_date = float(line[18:32])
            # if on date requested +- 1 day, save TLE
            if j_date - 1. <= epoch_date <= j_date + 1.:
                tle_list.append((line, next(tf)))
            # TLE are in desecending order, if after this point, can break from file
            elif epoch_date < j_date - 1.:
                break
            # only parse `1` lines
            else:
                next(tf)

    return tle_list


# ******************************************************************************
def MISR_fetch(orb_num, path_num, orb_date, projdir):
# ******************************************************************************
# Function for downloading all necessary MISR files. Will download
# MISR 1B2 Terrain Data, Ancillary Geographic Product, MISR Aerosol Parameters
# and MISR Geometric Parameters from NASA ftp site using wget.
# Input: Int requested orbit number
#        Int corresponding path number
#        String (YYYY.MM.DD) of corresponding calendar date
#        String path of project directory
# Output: None (will place 12 files in projdir)
# ------------------------------------------------------------------------------

    # --------------------------------------------------------------------------
    # Ensure the correct amount of digits in path and orbit numbers, and
    # create MISR directory to store downloaded files
    # --------------------------------------------------------------------------

    # orb_num = ('0' * (6 - len(str(orb_num)))) + str(orb_num) if len(str(orb_num)) < 6 else str(orb_num)
    # path_num = ('0' * (3 - len(str(path_num)))) + str(path_num) if len(str(path_num)) < 3 else str(path_num)
    cameras = ['DF', 'CF', 'BF', 'AF', 'AN', 'AA', 'BA', 'CA', 'DA']

    if not os.path.exists(projdir + '/MISRFiles'):
        os.makedirs(projdir + '/MISRFiles')

    # --------------------------------------------------------------------------
    # Use ftplib to download MI1B2T (9), MIANCAGP (1), MIB2GEOP (1)
    # and MIL2ASAE (1) files from NASA ftp site
    # --------------------------------------------------------------------------

    import ftplib
    # TODO loop this

    larc = ftplib.FTP('l5ftl01.larc.nasa.gov')
    larc.login()

    # terrain
    larc.cwd('MISR/MI1B2T.003/{date}/'.format(date=orb_date))
    for cam in cameras:
        GRP = 'MISR_AM1_GRP_TERRAIN_GM_P{p:03d}_O{o:06d}_{c}_F03_0024.hdf'.format(p=path_num, o=orb_num, c=cam)

        with open('{dir}/MISRFiles/{file}'.format(dir=projdir, file=GRP), 'wb') as f:
            larc.retrbinary('RETR ' + GRP, f.write)

    # geo params
    larc.cwd('../../MIB2GEOP.002/{date}/'.format(date=orb_date))
    GP_file = 'MISR_AM1_GP_GMP_P{p:03d}_O{o:06d}_F03_0013.hdf'.format(p=path_num, o=orb_num)

    with open('{dir}/MISRFiles/{file}'.format(dir=projdir, file=GP_file), 'wb') as f:
        larc.retrbinary('RETR ' + GP_file, f.write)

    # aerosols
    try:
        larc.cwd('../../MIL2ASAE.002/{date}/'.format(date=orb_date))  # TODO catch error if doesnt exist

        AS_file = 'MISR_AM1_AS_AEROSOL_P{p:03d}_O{o:06d}_F03_0013.hdf'.format(p=path_num, o=orb_num)

        with open('{dir}/MISRFiles/{file}'.format(dir=projdir, file=AS_file), 'wb') as f:
            larc.retrbinary('RETR ' + AS_file, f.write)
    except ftplib.error_perm:
        # will pass if the .002 aerosol files dont exist for this date, they were phased out on 2017.06.01
        pass

    # ancillary geo product
    larc.cwd('../../MIANCAGP.001/1999.11.07/')
    AGP_file = 'MISR_AM1_AGP_P{p:03d}_F01_24.hdf'.format(p=path_num)

    with open('{dir}/MISRFiles/{file}'.format(dir=projdir, file=AGP_file), 'wb') as f:
        # the product is stored in one of two days, with no way of telling where beforehand
        try:
            larc.retrbinary('RETR ' + AGP_file, f.write)
        except ftplib.error_perm:
            larc.cwd('../1999.11.08/')
            larc.retrbinary('RETR ' + AGP_file, f.write)

    larc.close()

    return


# ******************************************************************************
if __name__ == "__main__":
# ******************************************************************************
# Main procedure for preparing for MINX plume digitization.
# This script (1) Downloads relevant MODIS MOD14 Fire Granules. (2) Processes
# MOD14 Fire Granules through MINX. (3) Downloads relevant MISR data files.
# Upon the completion of this program, plume digitization should be completely
# prepped for.
# This program may take some time to complete, depending on download speeds.
# Please be patient and do not cancel the script while it is running. If you do
# so, you may be left with corrupted data files, and leftover metadata files.
# To fix this, you should clear your project directory and begin again.
# This script may also fail in certain circumstances. Some unlikely examples
# include: for orbits greater than 100 000, for coordinates near the equator/
# date line intersect and for coordinates near the planetary poles. If any
# errors do occur, MINX can still be used manually. See the appendices of the
# wiki page for manual instructions.
# More information and a usage guide can be found at:
# https://wiki.cmc.ec.gc.ca/wiki/Dicksonn/MINX4
# Input: Longitude, Latitude coordinates to search
#        Date to search
#        EarthData Username and Password
#        [[MOD14 and project directories]]
# Output: None (will download and process files)
# ------------------------------------------------------------------------------

# --------------------------------------------------------------------------
# (0) Parse and validate command line arguments
# --------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description="Main procedure for preparing for MINX plume digitization.\
     This script (1) Downloads relevant MODIS MOD14 Fire Granules. (2) Processes MOD14 Fire Granules through MINX. \
     (3) Downloads relevant MISR data files. Upon the completion of this program, plume digitization should be completely prepped for. \
     This program may take some time to complete, depending on download speeds. Please be patient and do not cancel the script while it is running. \
     If you do so, you may be left with corrupted data files, and leftover metadata files. To fix this, you should clear your project directory and begin \
     again. This script may also fail in certain circumstances. Some unlikely examples include: for orbits greater than 100 000, for coordinates near the \
     equator/dateline intersect and for coordinates near the planetary poles. If any  errors do occur, MINX can still be used manually. \
     See the appendices of the wiki page for manual instructions.  \
     More information and a usage guide can be found at: https://wiki.cmc.ec.gc.ca/wiki/Dicksonn/MINX4")

    parser.add_argument('-L', '--coords', type=float, nargs='*', required=True, metavar='FLOAT',
                        help='Longitude, latitude to search for hdf files around, or range to search within. minLON maxLON minLAT maxLAT')
    parser.add_argument('-d', '--date', required=True, nargs='*', metavar='YYYY.MM.DD',
                        help='Date you wish to search for hdf files on. Must have format YYYY.MM.DD. Can also search over range of dates, \
                        simply by providing a second, later date.')
    parser.add_argument('-u', '--user', required=True,
                        help='Your EarthData username. You must have a NASA EarthData account to download these files.')
    parser.add_argument('-p', '--password', default=None,
                        help='Your EarthData password. If not provided, you will be asked to enter it before the script may run.')
    parser.add_argument('--grandir', default='/space/hall2/sitestore/eccc/aq/r1/nod001/MINX/Working/MOD14_Granules',
                        help='Directory where MOD14 Fire Granules will be stored. Must provide full path. Will default to  \
                        /space/hall2/sitestore/eccc/aq/r1/nod001/MINX/Working/MOD14_Granules')
    parser.add_argument('--projdir', default=os.getcwd(),
                        help='Project Directory. Will be assumed as current directory, if not given.')
    parser.add_argument('--gui', action='store_true',
                        help='Fire pixel processing gui. Include if you wish to set specific settings for the MODIS fire pixel processing, \
                        a GUI will appear at the time of processing requesting your settings. Default setting will be applied otherwise.')
    parser.add_argument('--relaxed', action='store_true',
                        help='By default, only fire pixels within your requested boundaries will be proccesed, if you requested a range of coordinates. \
                        Include this tag to search all downloaded MODIS Granules. You will receive processed fires slightly outside your range.')
    parser.add_argument('--localgran', action='store_true',
                        help='Include flag if you already have all necessary MODIS MOD14 Fire Granules stored locally. Please ensure that they are stored \
                        in the correct directory specified in the --grandir option.')
    args = parser.parse_args()

    # --------------------------------------------------------------------------
    # EarthData Login
    # --------------------------------------------------------------------------

    username = args.user
    password = args.password

    while not password:
        password = getpass.getpass("EarthData Password: ")

    # --------------------------------------------------------------------------
    # Dates
    # --------------------------------------------------------------------------

    if len(args.date) > 2:
        parser.error("Invalid date format. Please provide either a single date or a range of the form YYYY.MM.DD YYYY.MM.DD")
    if not all([re.compile('\d{4}\.\d{2}\.\d{2}').match(date) for date in args.date]):
        parser.error("Invalid date format. Dates must be given as YYYY.MM.DD")
    argdates = [datetime.datetime.strptime(date, '%Y.%m.%d') for date in args.date]
    s = argdates[0]
    dates = []
    while s <= argdates[-1]:
        dates.append(s)
        s += datetime.timedelta(days=1)

    # --------------------------------------------------------------------------
    # Coordinates
    # --------------------------------------------------------------------------

    if len(args.coords) == 2:
        coords = (args.coords[0], args.coords[1])  # (Lon, Lat)
        args.relaxed = True
    elif len(args.coords) == 4:
        minLON, maxLON, minLAT, maxLAT = args.coords
        coords = [(minLON, minLAT), (minLON, maxLAT), (maxLON, minLAT), (maxLON, maxLAT)]

        if maxLON < minLON or maxLAT < minLAT or -90 > maxLAT > 90 or -180 > maxLON > 180:
            parser.error("Invalid coordinate range. Please provide a range of the form: minLON, maxLON, minLAT, maxLAT. If you wish to wrap around the antimeridian, \
                you must split your region and run this program twice, once on either side of 180 degrees longitude.")
    else:
        parser.error("Invalid coordinates. Please provide either a single point or a range of the form: minLON maxLON minLAT maxLAT")

    # --------------------------------------------------------------------------
    # Directories and files
    # --------------------------------------------------------------------------

    directory = os.path.abspath(args.grandir)       # TODO careful with this, if no slash at beginning appends current dir
    projdir = os.path.abspath(args.projdir)

    if not os.path.exists(directory):
        os.makedirs(directory)

    if not os.path.exists(projdir):
        os.makedirs(projdir)

    if args.gui:
        import Tkinter
        import tkFileDialog
        rt = Tkinter.Tk()
        rt.withdraw()
        directory = tkFileDialog.askdirectory(parent=rt, initialdir=directory, title="Select MOD14 Granules directory")
        projdir = tkFileDialog.askdirectory(parent=rt, initialdir=projdir, title="Select project directory")


# --------------------------------------------------------------------------
# (1) Download algorithm for MODIS MOD14 Fire Granules
# --------------------------------------------------------------------------

    load = pyloader.Load_message('Beginning')

    if not args.localgran:
        for date in dates:

            load.update('> {date} | Simulating Terra orbit'.format(date=date.strftime('%Y/%m/%d')))

            # --------------------------------------------------------------------------
            # Get orbital path of Terra satellite for this date using the most relevant
            # two line element, with a Fsimulated time interval of 30 seconds.
            # These TLE files must be kept updated, for use with more modern dates, as
            # the accuracy of this calculation will decrease exponentially otherwise
            # --------------------------------------------------------------------------

            import pyorbital.orbital as porb

            j_date = int(date.strftime('%y%j'))

            res = {
                'lonlat_arr': [],
                'last_nearest': False,
                'tle': read_TLE_dates('/fs/home/fs1/eccc/aq/arqi/mas001/computer_programs/python/MINX/data/Terra_TLE.txt', j_date),
                'orbit': None,
                'lonlat_paths': [[]],
            }

            iterdate = date
            index = 0
            while date + datetime.timedelta(1) > iterdate:

                time_frac = ((iterdate.hour * 60. + iterdate.minute) * 60. + iterdate.second) / (24. * 60. * 60.)

                nearest_tle = min(res['tle'], key=lambda line: abs(float(line[0][18:32]) - float(iterdate.strftime('%y%j.')) + time_frac))

                # either create a new orbit object or use the previous
                if nearest_tle != res['last_nearest']:
                    res['orbit'] = porb.Orbital('Terra', line1=nearest_tle[0], line2=nearest_tle[1])

                lon, lat, alt = res['orbit'].get_lonlatalt(iterdate)
                res['lonlat_arr'].append((lon, lat, iterdate))

                res['last_nearest'] = nearest_tle

                index += 1
                iterdate += datetime.timedelta(minutes=0.5)

            # --------------------------------------------------------------------------
            # Split the Terra orbit into seperate path sections for each 5 minute
            # interval, or when crossing the antimeridian
            # --------------------------------------------------------------------------

            path_index = 0
            for ind, (lon, lat, time) in enumerate(res['lonlat_arr']):

                if not float(time.strftime('%H%M.%S')) % 5 and ind:
                    path_index += 1
                    res['lonlat_paths'].append([])

                res['lonlat_paths'][path_index].append((lon, lat, time))

                try:
                    # if the next crosses the antimeridian/poles, split it into a new 'path'
                    if ((lon < 0 and res['lonlat_arr'][ind + 1][0] > 0) or (lat < 0 and res['lonlat_arr'][ind + 1][1] > 0)):
                        path_index += 1
                        res['lonlat_paths'].append([])
                except IndexError:
                    pass

            # --------------------------------------------------------------------------
            # Create a larger polygon around each path section to approximate the
            # swath coverage of the MODIS camera
            # --------------------------------------------------------------------------

            hdf_timelist = []
            for path_section in res['lonlat_paths']:
                if path_section:
                    poly_left, poly_right = [], []
                    for pind, point in enumerate(path_section[:-1]):
                        x1, y1, time1 = point
                        x2, y2, time2 = path_section[pind + 1]

                        dx = x2 - x1
                        dy = y2 - y1
                        d = math.sqrt((dx * dx) + (dy * dy)) / 2.

                        swath_width = 2330000  # MODIS swath width in metres
                        r = pm.width_to_degrees(swath_width, y1 + dy)

                        xo = x1 + math.sqrt((r * r) + (d * d))
                        theta = math.atan(r / d)
                        phi = math.atan(dy / dx)

                        alpha_left = phi + theta if phi > 0. else phi - theta
                        alpha_right = phi - theta if phi > 0. else phi + theta

                        poly_time = time1 + (time2 - time1) / 2.

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

                    poly_right.reverse()
                    MODIS_poly = poly_left + poly_right

                    # --------------------------------------------------------------------------
                    # Check if the requested coordinates lie within any of the path sections.
                    # If so, grab the 5 minute time interval for use in downloading said
                    # granule
                    # --------------------------------------------------------------------------

                    rounded_time = path_section[0][2] - datetime.timedelta(minutes=path_section[0][2].minute % 5)

                    if isinstance(coords, list):
                        if polygon_check(coords, MODIS_poly, ccwsorted=True):  # TODO is coords sorted here?
                            hdf_timelist.append(rounded_time)
                    else:
                        if polygon_check_single(MODIS_poly, point=coords, ccwsorted=True):
                            hdf_timelist.append(rounded_time)

            # --------------------------------------------------------------------------
            # Get a list of all granules for this date and find the relevant
            # ones based on the times found above
            # --------------------------------------------------------------------------

            load.update('> {date} | Downloading MODIS fire granules'.format(date=date.strftime('%Y/%m/%d')))

            http_url = 'https://e4ftl01.cr.usgs.gov/MOLT/MOD14.006/{0}'.format(date.strftime('%Y.%m.%d'))
            folder_content = requests.get(http_url).text

            # MOD14.A2003223.0230.006.2015076112402.hdf
            # TODO this regex may need to be altered slightly to make sure its only finding these times in the time slot and not anywhere else in the filename

            srch_times = r'|'.join([tm.strftime('%H%M') for tm in hdf_timelist])
            http_filelist = re.findall(r'(?<=<a href=")[^"]*.(?:{times}).*.hdf(?=")'.format(times=srch_times), folder_content)

            # --------------------------------------------------------------------------
            # Download all relevant granules. The request is redirected, when
            # downloading a file from this site, upon the first attempt. Must supply
            # authorization details to the redirected url to gain access to files
            # --------------------------------------------------------------------------

            for gran in http_filelist:
                gran_url = '{url}/{gran}'.format(url=http_url, gran=gran)
                redirected_r = requests.request('get', gran_url)
                MODIS_r = requests.request('get', redirected_r.url, stream=True, auth=(username, password))

                with open("{dir}/{gran}".format(dir=directory, gran=gran), 'wb') as f:
                    for chunk in MODIS_r.iter_content(chunk_size=1024):
                        if chunk:
                            f.write(chunk)

    sys.stdout.write("All required MODIS MOD14 files downloaded, proceeding to process fire pixels through MINX.\n")


# --------------------------------------------------------------------------
# (2) Process MODIS MOD14 Fire Granules
# --------------------------------------------------------------------------

    load.update('> Processing relevant fire granules')

    # --------------------------------------------------------------------------
    # Create and run temporary IDL batch file containing commands to restore
    # the preMINX.sav file, initialize parameters, process fire pixels and
    # return the MISR orbit reference table. (No better way in std library)
    # The save file should be stored in the same directory as this script, and
    # IDL must be installed.
    # The save file contains only the procedures necessary for fire pixel
    # processing, which are taken directly from the MINX source code, with
    # slight alterations to ignore gui elements, if desired.
    # More information can be found in said code, retrievable at:
    # https://github.com/nasa/MINX
    # --------------------------------------------------------------------------

    idl_args = {
        "ProjName": os.path.basename(projdir).replace('_', '-').replace(' ', '-'),
        "DefaultProjDir": projdir + '/',
        "gui": 1 if args.gui else 0,
        "DefaultMODDir": directory + '/',
        "DefaultBeg": str(argdates[0].year) + str(argdates[0].timetuple().tm_yday),
        "DefaultEnd": str(argdates[-1].year) + str(argdates[-1].timetuple().tm_yday)
    }

    f = tempfile.NamedTemporaryFile()
    f.write(
        "RESTORE, '{0}/preMINX.sav'\n".format(os.path.dirname(os.path.abspath(__file__))).encode() +
        "INITMINXPARAMS()\n".encode() +
        "PROCESSMODISFIREPIXELS, '{ProjName}', '{DefaultProjDir}', {gui}, '{DefaultMODDir}', '{DefaultBeg}', '{DefaultEnd}'\n".format(**idl_args).encode() +
        "WritePreferencesFile\n".encode() +
        "PRINT, !KON.MISR_ORBIT_REF\n".encode() +
        "exit".encode()
    )
    f.flush()
    fout = subprocess.Popen(['idl', f.name], stdout=subprocess.PIPE).communicate()[0]
    f.close()

    orbit_ref_table = [float(x) for x in fout.split()]  # TODO switch to grabbing table from misrutils

    # --------------------------------------------------------------------------
    # If INITMINXPARAMS() looks for a .tbl file, and cannot find it, the first
    # index of orbit_ref_table will incorrectly become 0 (from error code)  # TODO I moved the tbl file to /python/ need to recompile the .sav file
    # --------------------------------------------------------------------------

    if not orbit_ref_table[0]:
        orbit_ref_table = orbit_ref_table[1:]

    sys.stdout.write("All MODIS MOD14 granules processed, proceeding to download necessary MISR data files.\n")


# --------------------------------------------------------------------------
# (3) Download algorithm for MISR data files
# --------------------------------------------------------------------------

    load.update('> Beginning MISR download')

    # --------------------------------------------------------------------------
    # Download MISR files for each orbit in MisrOrderList text file
    # --------------------------------------------------------------------------

    if args.gui:
        orderfile = tkFileDialog.askopenfilename(parent=rt, initialdir=projdir, title="Select the file containing the Misr Order List")
        with open(orderfile, 'r') as f:
            orbits = list(map(int, f.read()[:-1].split(',')))
    else:
        glb_srch = glob.glob('%s/MisrOrderList_MOD14_*.txt' % projdir)
        if len(glb_srch) > 1:
            mssg = "Too many projects in given project directory. Please split different projects into seperate directories."
            raise RuntimeError(mssg)
        else:
            try:
                with open(glb_srch[0], 'r') as f:
                    orbits = list(map(int, f.read()[:-1].split(',')))
            except IndexError as e:
                mssg = "No MISR Process List found in given project directory. \
                        Please ensure you are in the right directory and have processed the MODIS MOD14 Fire Granules."
                raise IndexError(mssg)
            except ValueError as e:
                mssg = "No valid orbits provided. Please ensure the MisrOrderList text file is correct."
                raise ValueError(mssg)

    del_list = []
    for orb_num in orbits:

        load.update('> Orbit {orb} | Downloading MISR files'.format(orb=orb_num))

        if not args.relaxed:

            # --------------------------------------------------------------------------
            # Check if any fire pixels from this orbit are within the requested box.
            # If not, skip download of this orbit and remove orbit from Processing list
            # --------------------------------------------------------------------------

            file_orb_num = ('0' * (6 - len(str(orb_num)))) + str(orb_num) if len(str(orb_num)) < 6 else str(orb_num)
            if args.gui:
                ini = projdir + ('/FirePixels_MOD14_' + idl_args['ProjName']
                                 if os.path.exists('{0}/FirePixels_MOD14_{1}/'.format(projdir, idl_args['ProjName'])) else '')
                fire_file = tkFileDialog.askopenfilename(
                    parent=rt, initialdir=ini,
                    title="Select the file containing the processed Fire Pixel data",
                    filetypes=[('orb', 'FirePixels_MOD14_{0}_{1}.txt'.format(file_orb_num, idl_args['ProjName']))]
                )
            else:
                fire_file = '{0}/FirePixels_MOD14_{1}/FirePixels_MOD14_{2}_{1}.txt'.format(projdir, idl_args['ProjName'], file_orb_num)

            with open(fire_file, 'r') as f:
                for _ in xrange(4):
                    next(f)
                fire_box = [polygon_check_single(coords, point=(float(row[0]), float(row[1]))) for row in csv.reader(f, delimiter=' ', skipinitialspace=True)]

            if not any(fire_box):
                del_list.append(orb_num)
                continue

        # --------------------------------------------------------------------------
        # Determine the dates and path numbers of each orbit and download the
        # relevant files from the NASA ftp server
        # --------------------------------------------------------------------------

        path_num = mu.OrbittoPath(orb_num)
        date = mu.OrbittoTime(orb_num).strftime('%Y.%m.%d')

        MISR_fetch(orb_num, path_num, date, projdir)

    if args.gui:
        procfile = tkFileDialog.askopenfilename(parent=rt, initialdir=projdir, title="Select the file containing the Misr Process List")
        os.rename(procfile, "%s/PlumeProjOrbitList.txt" % projdir)
    else:
        glb_srch = glob.glob('%s/MisrProcessList_MOD14_*.txt' % projdir)
        if glb_srch:
            if len(glb_srch) > 1:
                mssg = "Too many projects in given project directory. Please split different projects into seperate directories."
                raise RuntimeError(mssg)
            else:
                os.rename(glb_srch[0], "%s/PlumeProjOrbitList.txt" % projdir)

        elif not glob.glob("%s/PlumeProjOrbitList.txt" % projdir):
            mssg = "No MISR Process List found in given project directory. \
                    Please ensure you are in the right directory and have processed the MODIS MOD14 Fire Granules. "
            raise RuntimeError(mssg)

    # --------------------------------------------------------------------------
    # Edit file to include new location of MISR data files and remove unwanted
    # orbits
    # --------------------------------------------------------------------------

    with open('%s/PlumeProjOrbitList.txt' % projdir, 'r+') as f:
        data = f.read().split('\n')
        MISR_path = projdir + '/MISRFiles/'
        read_data, write_data = data[1:], [MISR_path]

        for row in read_data:
            try:
                if int(row.split()[0]) not in del_list:
                    write_data.append(row)
            except ValueError:
                write_data.append(row)
            except IndexError:
                continue

        f.seek(0)
        f.write('\n'.join(write_data))
        f.truncate()

    # --------------------------------------------------------------------------
    # Ensure that not all of the orbits were excluded
    # --------------------------------------------------------------------------

    if len(orbits) == len(del_list):
        mssg = "No fire pixels were found that satisfy criteria.\n\
                If this is unexpected, double-check the coordinates and dates requested, and consider the --relaxed option."
        raise RuntimeError(mssg)

    load.stop()

    sys.stdout.write("All required MISR files downloaded, proceed to processing smoke plumes through MINX.\n")

    sys.stdout.write('Total time taken: {0} seconds'.format(load.delta_t))
