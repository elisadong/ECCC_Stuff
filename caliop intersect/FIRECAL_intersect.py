#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2018-03-08 15:12:55
# Contact:    Nolan.Dickson@Canada.ca

# TODO this isnt the fastest thing in the world

import netCDF4
import argparse
import datetime
import numpy as np
import MISRutils as mu
from Cpoly.Cpoly import polygon_check_single

try:
    import cPickle as pickle  # 2/3 compatibility
except ImportError:
    import pickle


class Fire(object):
    """
    Fire object which is created from a close bundle of MODIS fire pixels,
    which surpass a total FRP threshold
    """

    def __init__(self, pixels, FRP):
        super(Fire, self).__init__()

        self.pixels = pixels  # coordinates of each pixel
        self.FRP = FRP        # total FRP of fire

        sum_t = sum([p[3].timestamp() for p in pixels])
        self.datetime = datetime.datetime.fromtimestamp(sum_t / len(pixels))

    def create_boundaries(self, padding):
        '''
        From the fire pixels in this fire, create a surrounding polygon,
        which will constitute an estimate of the near plume
        '''
        import scipy.spatial as spat
        # TODO include wind direction and do this only in general wind direction, right now just creates a padded polygon around fires

        def _rot90(x, y):
            return (y, -x)

        def _scaled_norm(x, y, s):
            r = x * x + y * y
            # TODO width to degrees is only meant for horizontal conversions I believe
            scale = mu.width_to_degrees(s, y)
            return ((x / r) * scale, (y / r) * scale)

        # -----------------------------------------------------------------
        # Get the coordinates of each pixel in this fire and use them,
        # along with the SciPy function, to compute the Convex Hull of
        # the fire
        # -----------------------------------------------------------------

        coords = [(p[0], p[1]) for p in self.pixels]

        self.hull = spat.ConvexHull(coords)
        hull_coords = [coords[i] for i in self.hull.vertices]

        # TEMP FIX #

        leftP, rightP = sorted(hull_coords,
                               key=lambda t: t[0])[::len(hull_coords) - 1]
        bottomP, topP = sorted(hull_coords,
                               key=lambda t: t[1])[::len(hull_coords) - 1]

        top_pad = mu.width_to_degrees(padding, topP[1])
        bot_pad = mu.width_to_degrees(padding, bottomP[1])

        self.boundaries = [
            (leftP[0] - top_pad, topP[1] + top_pad),
            (rightP[0] + top_pad, topP[1] + top_pad),
            (rightP[0] + bot_pad, bottomP[1] - bot_pad),
            (leftP[0] - bot_pad, bottomP[1] - bot_pad)
        ]

        # END TEMP FIX #

        # TODO fix the polygon expansion code, gives weird extended shape
        #       may be due to sharp vertices? use square as approx till then
        # !! # -----------------------------------------------------------------
        # !! # Calculate the vector lines which run parallel to the convex
        # !! # hull lines by finding the CW normal and scaling each vertex
        # !! # along it
        # !! # -----------------------------------------------------------------

        # !! ext_lines = []
        # !! for pind, pnt1 in enumerate(hull_coords):
        # !!     pnt0 = hull_coords[pind - 1]

        # !!     V_edge = (pnt1[0] - pnt0[0], pnt1[1] - pnt0[1])
        # !!     rv = _rot90(V_edge[0], V_edge[1])
        # !!     d = _scaled_norm(rv[0], rv[1], padding)

        # !!     pnt0_ext = (pnt0[0] + d[0], pnt0[1] + d[1])
        # !!     pnt1_ext = (pnt1[0] + d[0], pnt1[1] + d[1])

        # !!     ext_lines.append((pnt0_ext, pnt1_ext))

        # !! # -----------------------------------------------------------------
        # !! # Find the intersects of these expanded lines by using the eqn's:
        # !! #   P = P1 + t * (P2 - P1)
        # !! #   P = P3 + u * (P4 - P3)
        # !! #   a1*t + b1*u = c1
        # !! #   a2*t + b2*u = c2
        # !! # and solving for the value t, to find the intersection.
        # !! # These intersections represent the vertices of the expanded
        # !! # polygon
        # !! # -----------------------------------------------------------------

        # !! # find intersect of lines
        # !! ext_vert = []
        # !! for lind, line1 in enumerate(ext_lines):
        # !!     line0 = ext_lines[lind - 1]

        # !!     a1 = line0[1][0] - line0[0][0]
        # !!     b1 = line1[0][0] - line1[1][0]
        # !!     c1 = line1[0][0] - line0[0][0]

        # !!     a2 = line0[1][1] - line0[0][1]
        # !!     b2 = line1[0][1] - line1[1][1]
        # !!     c2 = line1[0][1] - line0[0][1]

        # !!     t = (b1 * c2 - b2 * c1) / (a2 * b1 - a1 * b2)

        # !!     ext_vert.append((
        # !!         line0[0][0] + t * (line0[1][0] - line0[0][0]),
        # !!         line0[0][1] + t * (line0[1][1] - line0[0][1])
        # !!     ))

        # !! self.boundaries = ext_vert
        return

    def get_overpasses(self, delta_t):
        '''
        From the created boundary find all calipso overpasses within some time
        '''

        # TODO INTERSECT OPTIONS:
        # (i) find all Caliop overpasses of these ranges and polycheck all the fires (may be fastest, doesnt redraw orbit)
        # (ii) draw the orbit of caliop at the time of each fire +- some hours and poly check the fire (may be sketchy)
        # (iii) pass the bounds (bbox) of each fire to find_overpass, along with date range <<<

        td = datetime.timedelta(hours=delta_t)

        self.overpasses = mu.find_overpass(
            coord=self.boundaries,
            daterange=[self.datetime - td, self.datetime + td],
            swath_width=None,
            sat='Calipso'
        )


def read_fire_pixels(coords, startdate, enddate):
    '''read the fire pixels from the L3 file which match the date/coord range'''

    # ---------------------------------------------------------------------
    # Open MCD14ML fire pixel file and find the indices of all points
    # which lie within the specified date range
    # ---------------------------------------------------------------------

    L3_path = '/fs/site2/dev/eccc/aq/r1/nod001/MINX/Comp/workdir/MCD14ML/' \
              'MCD14ML.{Yr:%Y}.006.01.nc'.format(Yr=startdate)
    L3_file = netCDF4.Dataset(L3_path)
    YYYYMMDD_vals = L3_file['YYYYMMDD'][0]

    day_inds = np.where(
        (YYYYMMDD_vals >= float(startdate.strftime('%Y%m%d'))) &
        (YYYYMMDD_vals <= float(enddate.strftime('%Y%m%d')))
    )[0]

    # ---------------------------------------------------------------------
    # Read in the relevant data for all necessary date indices.
    # Date and time arrays must be converted to relevant strings
    # ---------------------------------------------------------------------

    lon_arr = L3_file['lon'][0][day_inds[0]:day_inds[-1]]
    lat_arr = L3_file['lat'][0][day_inds[0]:day_inds[-1]]

    FRP_arr = L3_file['FRP'][0][day_inds[0]:day_inds[-1]]

    # TODO better string conversions
    d_arr = YYYYMMDD_vals[day_inds[0]:day_inds[-1]].astype('|U8')

    tfloat_arr = L3_file['HHMM'][0][day_inds[0]:day_inds[-1]]
    t_arr = np.fromiter(map('{:04g}'.format, tfloat_arr), dtype='|U4')

    # ---------------------------------------------------------------------
    # Check if each point lies within the coordinate range, and return
    # a list of its coordinates, frp value and datetime for each day
    # if it does
    # ---------------------------------------------------------------------

    fire_pixels = {}
    for ind, (lon, lat) in enumerate(zip(lon_arr, lat_arr)):

        if d_arr[ind] not in fire_pixels:
            fire_pixels[d_arr[ind]] = []

        dt = datetime.datetime.strptime(d_arr[ind] + t_arr[ind], '%Y%m%d%H%M')
        FRP = FRP_arr[ind]

        if polygon_check_single(coords, (lon, lat), ordered=True):
            fire_pixels[d_arr[ind]].append((lon, lat, FRP, dt))

    return list(fire_pixels.values())


def find_fires(pixel_list, dist=3000., threshold=500., delta_t=12, pad=5000.):  # TODO fiddle with default dist, think MODIS is 1k res
    '''
    From a list of fire pixels, find those grouped together with a total FRP
    greater than the threshold, and make a Fire out of them
    '''
    def _dfs(node, index):
        '''recursive depth-first search'''
        taken[index] = True
        for gind, item in enumerate(groups):
            if not taken[gind] and not node.isdisjoint(item):
                node.update(_dfs(item, gind))
        return node

    # ---------------------------------------------------------------------
    # Group together all fire pixels which lie within `dist` of
    # one another, using the haversine formula
    # ---------------------------------------------------------------------

    groups = []
    for pind, pnt in enumerate(pixel_list):
        grp = {pnt}
        for testpnt in pixel_list[pind:]:
            testdist = mu.haver_dist(pnt[0], pnt[1], testpnt[0], testpnt[1])

            if testdist <= dist:
                grp.update([testpnt])

        groups.append(grp)

    # ---------------------------------------------------------------------
    # Use a recursive depth-first search to solve the connected-component
    # problem of combining all sets which contain common elements
    # ---------------------------------------------------------------------

    fire_list = []
    taken = [False for _ in groups]

    for ind, node in enumerate(groups):
        if not taken[ind]:
            fire_list.append(list(_dfs(node, ind)))

    # ---------------------------------------------------------------------
    # Determine which groupings of fire pixels contain a total FRP which
    # exceeds the threshold and create Fire objects from them
    # ---------------------------------------------------------------------

    fire_objects = []
    for fire_grp in fire_list:
        total_FRP = sum(p[2] for p in fire_grp)

        if len(fire_grp) > 2 and total_FRP >= threshold:
            fire = Fire(fire_grp, total_FRP)
            fire.create_boundaries(padding=pad)
            fire.get_overpasses(delta_t)
            fire_objects.append(fire)

    return fire_objects


if __name__ == "__main__":
    from GEMCAL_intersect import Fire  # to correct pickle issues

    # ---------------------------------------------------------------------
    # Parse and check all arguments
    # ---------------------------------------------------------------------

    parser = argparse.ArgumentParser(
        description="Find CALIOP/Fire Granule intercepts"
    )
    parser.add_argument(
        '-d', '--daterange', required=True, nargs='*',
        help="Date, must have format YYYY.MM.DD. Can also "
             "search over range of dates, "
             "simply by providing a second, later date. Please only "
             "do a date range inside the same year."
    )
    parser.add_argument(
        '-L', '--coords', type=float, nargs='*', required=True,
        help="Longitude, latitude range to search for hdf files within: "
             "minLON maxLON minLAT maxLAT"
    )
    parser.add_argument(
        '-t', '--deltat', type=int, required=False, default=12,
        help="Max amount of time allowed between a fire's observation by MODIS "
             "to the overpass of Calipso, given in hours."
    )
    parser.add_argument(
        '-f', '--frp', type=float, required=False, default=500.,
        help="Minimum total FRP used to define a fire."
    )
    parser.add_argument(
        '-r', '--radius', type=float, required=False, default=3000.,
        help="Maximum distance between fire pixels used to group them together "
        "into a larger 'Fire'."
    )
    parser.add_argument(
        '-p', '--padding', type=float, required=False, default=5000.,
        help="Padded distance around fire boundaries to look for a Calipso "
        "overpass above."
    )
    parser.add_argument(
        '--debug', action='store_true',
        help='Create a log file with information for use in debugging.'
    )
    args = parser.parse_args()

    debug_flag = args.debug  # TODO add debug info

    dates = [datetime.datetime.strptime(d, '%Y.%m.%d') for d in args.daterange]

    if len(dates) != 2:
        parser.error("Invalid daterange. Please provide a range of two dates" +
                     " in the form: YYYY.MM.DD YYYY.MM.DD")

    if len(args.coords) == 4:
        minLON, maxLON, minLAT, maxLAT = args.coords
        coords = [(minLON, minLAT), (minLON, maxLAT),
                  (maxLON, minLAT), (maxLON, maxLAT)]
    else:
        parser.error("Invalid coordinates. Please provide a range of the " +
                     "form: minLON maxLON minLAT maxLAT")

    # ---------------------------------------------------------------------
    # Find all fire pixels in these ranges, grouped by day
    # ---------------------------------------------------------------------

    day_arr = read_fire_pixels(coords, dates[0], dates[-1])

    # ---------------------------------------------------------------------
    # Create a flat list of Fire objects with boundaries and overpasses
    # calculated
    # ---------------------------------------------------------------------

    fire_list = []
    for pixel_list in day_arr:
        fire_list.extend(find_fires(
            pixel_list=pixel_list,
            dist=args.radius,
            threshold=args.frp,
            delta_t=args.deltat,
            pad=args.padding
        ))

    # ---------------------------------------------------------------------
    # Save all fires to a pickle file
    # ---------------------------------------------------------------------

    pkl = './intercepts_s{sd:%Y%m%d}_e{ed:%Y%m%d}.p'.format(
        sd=dates[0], ed=dates[-1])

    if fire_list:
        with open(pkl, 'wb') as file:
            pickle.dump(fire_list, file, protocol=2)
