
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2018-02-23 11:21:21
# Contact:    Nolan.Dickson@Canada.ca

# TODO TEST EVERYTHING, error checking, documentation like fou, increase consistency
# TODO find a way to be able to read option list for all possible model attrs, not just height
# TODO make sure all results have same shape as attributes, as I claim
# TODO decide whether calipso methods should be in calipso or in plume methods, like model

import os
import sys
import glob
import math
import datetime
import itertools
import matplotlib.pyplot as plt

import warnings
warnings.simplefilter('ignore', FutureWarning)

stdout_loc = sys.stdout


def hide_output():
    '''Hide all standard output. Will hide *all* output, not just the output from this module'''
    sys.stdout = open(os.devnull, 'w')


def show_output():
    '''Show standard output again, for use after hide_output'''
    sys.stdout = stdout_loc


def correct_digits(num, dig):
    '''*Given an integer number (str or int) and amount of digits, return zero-padded string of number*'''
    try:
        num = int(num)
        return ('0' * (dig - len(str(num)))) + str(num) if len(str(num)) < dig else str(num)
    except ValueError as e:
        raise e


def str_or_flt(x):
    '''*Given an unknown value (single value, not iterable), return either float or string of value, depending on value*'''
    try:
        return float(x)
    except ValueError:
        return x


def arg_to_list(x):
    '''*Given arg value or iterable of args, return list of args*'''
    if isinstance(x, list):
        return x
    else:
        return [x]


def arg_to_gen(x):
    '''*Given arg value or iterable of args, return generator of args*'''
    if isinstance(x, list):
        return (n for n in x)
    else:
        return (n for n in [x])  # comprehensions are faster than everything, including list(x)


def haver_dist(long1, lat1, long2, lat2):
    '''Calculate distance in metres between two geographic coordinates using the Haversine formula'''

    R = 6371000  # Radius of Earth, in metres
    phi_c = lat2 * math.pi / 180.
    phi_p = lat1 * math.pi / 180.
    del_phi = (lat1 - lat2) * math.pi / 180.
    del_lam = (long1 - long2) * math.pi / 180.

    a = math.sin(del_phi / 2.)**2 + math.cos(phi_p) \
        * math.cos(phi_c) * math.sin(del_lam / 2.)**2
    dist = R * 2 * math.atan2(a**0.5, (1 - a)**0.5)

    return dist


def drawmap(minlon=-180, minlat=-90, maxlon=180, maxlat=90):
    '''TODO'''
    from mpl_toolkits.basemap import Basemap

    minlon, minlat, maxlon, maxlat = map(int, [minlon, minlat, maxlon, maxlat])

    m = Basemap(projection='cyl', resolution='l', llcrnrlat=minlat, urcrnrlat=maxlat, llcrnrlon=minlon, urcrnrlon=maxlon)
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels([x for x in range(minlat, maxlat + 1, round(((maxlat - minlat) / 4) + 0.5))], labels=[True, False, False, True])
    m.drawmeridians([x for x in range(minlon, maxlon + 1, round(((maxlon - minlon) / 4) + 0.5))], labels=[True, False, False, True])
    m.drawcountries(linewidth=0.5)

    return m


def width_to_degrees(distance, latitude):
    '''Given a width, in metres, and a latitude, convert to degrees of horizontal longitude'''
    import math
    return distance / (math.cos(latitude * math.pi / 180.) * 111320)


def drawswath(path, swath_width, swath_function=width_to_degrees):
    '''Given an orbit path, return a polygon of a determined width in metres around it. Swath function must accept swath width and a latitude'''

    poly_left, poly_right = [], []
    for pind, point in enumerate(path[:-1]):
        x1, y1 = point
        x2, y2 = path[pind + 1]

        dx = x2 - x1
        dy = y2 - y1
        d = math.sqrt((dx * dx) + (dy * dy)) / 2.

        # if being exact would use WGS84 ellipsoid equation but whatever
        try:
            r = swath_function(swath_width, y1 + dy)
        except TypeError:
            r = swath_width

        xo = x1 + math.sqrt((r * r) + (d * d))
        theta = math.atan(r / d)
        phi = math.atan(dy / dx)

        alpha_left = phi + theta if phi > 0. else phi - theta
        alpha_right = phi - theta if phi > 0. else phi + theta

        poly_left.append((
            (((xo - x1) * math.cos(alpha_left)) - ((y1 - y1) * math.sin(alpha_left)) + x1),
            (((xo - x1) * math.sin(alpha_left)) + ((y1 - y1) * math.cos(alpha_left)) + y1),
        ))
        poly_right.append((
            (((xo - x1) * math.cos(alpha_right)) - ((y1 - y1) * math.sin(alpha_right)) + x1),
            (((xo - x1) * math.sin(alpha_right)) + ((y1 - y1) * math.cos(alpha_right)) + y1),
        ))

    poly_right.reverse()
    swath_poly = poly_left + poly_right
    return swath_poly


def clean_inv(*listN, **kwargs):
    '''
    Given N number of lists (as *listN), remove indices of all requested invalid values from all.
    kwargs must be the same format as the defaults, if you want to change them.
    {'nans': True, 'infs': True, 'negatives': False, 'fills': (False, [-9999,])}
    I know it would make more sense to just have them as normal arguments, but I like passing the
    lists in separately and arbitrarily, so i made this hacky, sue me. #TODO this may be doable now that I switched to py3
    "Although practicality beats purity." - Zen of Python

    *Does not return lists in same shape as rest of attribute, strips them of invalid values indices completely*
    '''

    # TODO optimize, error checking if all vlues invalid or if only some of listN sublists are actually lists, better way to loop?
    # TODO this is hacky man, it hurts, atleast make the docstring prettier please

    if not kwargs:
        # set kwargs to defaults
        kwargs = {'nans': True, 'infs': True, 'negatives': False, 'fills': (False, [-9999])}

    listN = [list_x for list_x in listN if isinstance(list_x, list)]

    if not listN:
        raise TypeError("arguments supplied must be lists, with an optional **kwargs dict.")

    try:
        if kwargs['nans']:
            listN = [[x for i, x in enumerate(sublistN) if not any([math.isnan(chk_sublistN[i]) for chk_sublistN in listN])] for sublistN in listN]
        if kwargs['infs']:
            listN = [[x for i, x in enumerate(sublistN) if not any([math.isinf(chk_sublistN[i]) for chk_sublistN in listN])] for sublistN in listN]
        if kwargs['negatives']:
            listN = [[x for i, x in enumerate(sublistN) if not any([chk_sublistN[i] < 0. for chk_sublistN in listN])] for sublistN in listN]
        if kwargs['fills'][0]:
            for fill_val in kwargs['fills'][1]:
                listN = [[x for i, x in enumerate(sublistN) if not any([chk_sublistN[i] == fill_val for chk_sublistN in listN])] for sublistN in listN]
    except KeyError:
        raise KeyError("kwargs must be in same format as defaults, see `help(clean_inv)`.")
    except IndexError:
        raise IndexError("all supplied lists must be of same size.")

    return listN if len(listN) > 1 else listN[0]


class Plume(object):
    """*Plume object for the visualization and manipulation of MINX results and comparison with GEMMACH model output and CALIPSO CALIOP Lidar data*

    Wildfire plume object. Utilizes the results of NASA's `MISR INteractive eXplorer` software, given in a `Plumes` text file.
    You must have utilized MINX, and preferably the preMINX tool, to digitize what plumes you desire to inspect.

    If you wish to compare the MINX results to a corresponding GEMMACH output file, you must ensure that said file completely covers the plume
    and has a filename reflecting the date. To attach the model data, either utilize the get_model method, or include the model directory during init.

    This object contains methods to plot, in 3D and 2D, the MINX output and model output, methods to receive basic stats on the output, technical
    methods to move, cut and find output data and methods for comparing the MINX output to the model output.
    Most methods allow you to specify which dataset you wish to utilize, the options for which are given under `attributes` below.
    Alongside these attributes, this object will contain the majority of the MetaData which is present in the MINX output file.

    Unless explicitly stated, all results lists which may be returned by methods will match the length and sorting of all attributes, which are
    organized based on MINX coordinate points.

    Note
    ----
        All documentation for the module, as well as further explanation, guidance and examples, is available on the wiki page:
        lolnotyetsry.wiki

        All examples in the documentation of this module can be tested using the doctest module. The examples will utilize a test
        plume text file, which should be found in the /tests/ folder alongside the postMINX file.

        Plume fill values are NaN.

    **Initialization**
        *Plume.from_coord((long, lat))*

        *Plume.from_id(orb_num, blk_num, plm_num)*

        *Plume.from_filename(filepath)*

    **Attributes**
        Coords: *longitude, latitude*

        Geographic data: *distance, direction, terrain*

        Height data: *nowind_height, wind_height, filtered_height*

        Wind Speed data: *across_wind, along_wind, total_wind*

        Model data: *model_height, model_value (if given model)*

        points: *list of points containing the above values for each coordinate*

    **Methods**
        3D plots: *contour, colour, scatter_3D, surface_plot*

        2D plots: *direction_plot, scatter_2D, poly_2D, hist*

        Stats: *max, min, median, mean, std*

        Technicals: *get_model, set_reference, find_coord, coord_section, cross_section*

        Comparison: *table, abs_diff, RMSD, corr_coeff, linear_comp*

    """

    clr = itertools.cycle(['r', 'b', 'c', 'm', 'y'])

    def __getitem__(self, key):
        return self.__dict__[key]

    def __repr__(self):
        return 'Plume.from_filename("{out}/{orb}/{file}"{mod})'.format(
            out=self.out_dir,
            orb=self.orbit_num,
            file=self.filepath,
            mod='' if not getattr(self, 'model_dir', '') else ', model_dir="{mdir}"'.format(mdir=self.model_dir)
        )

    @classmethod
    def from_coord(cls, coord, date=False, out_dir=os.getcwd(), nearest=False, exact=False, model_dir="", cal_int=""):
        '''*Initialization method for reading plume based on a coordinate point*

        Initialization method which will search through a given directory for plume text files, and search said text files for the requested coordinate
        point. Utilizes the polygon checking functions from either Cpoly or preMINX, depending on the version of python, and passes necessary
        file on to from_filename for initialization.
        Used as a classmethod in order to namespace this function as a factory method, for consistency purposes.

        Args
        ----
        coord : tuple of float or iter of tuple of float
            Longitude-latitude coordinate point to search for. If list, must be of a N-sided convex polygon. Concave bbox will not work.
        date : str, optional
            Date of plume to search for. If given, must have format `YYYY-MM-DD`.
        out_dir : str, optional
            Path of either the desired MINX project directory, or anywhere within the `Digitize_Output` directory
            of said project.
        nearest : bool, optional
            If specified, will search for the nearest plume to the coordinate point, rather than searching for a plume
            which contains it. (Not implemented yet)
        exact : bool, optional
            If specified, will only accept a plume if the exact desired coordinate point contains results in the MINX output.
        model_dir : str, optional
            Path of the directory containing the model output files. Including this parameter will automatically attempt
            to open and attach the model data for this plume. See get_model for more information.

        Returns
        -------
        out : Plume
            A Plume object satisfying the specific requirements.

        Raises
        ------
        IOError: If no file satisfying the specific requirements is found.

        '''
        try:
            from Cpoly import polygon_check_single, polygon_check
        except ImportError:
            # running on python2, use python implementation of polygon checking from preMINX
            from preMINX import polygon_check_single, polygon_check

        initdict = {}

        initdict['out_dir'] = out_dir[:out_dir.rfind('Digitize_Output') + 15] if 'Digitize_Output' in out_dir else out_dir + '/Digitize_Output'
        txt_files = glob.glob('{0}/*/Plumes_*.txt'.format(initdict['out_dir']))
        for txt in txt_files:
            plm_flag = False
            date_acquired = None

            with open(txt, 'r') as file:
                cont = file.read().split('\n')

                polypoints = []
                poly_flag = False
                for ind, line in enumerate(cont):
                    if 'Date acquired' in line:
                        date_acquired = line.split(':', 1)[1].strip()  # date must be given in YYYY-MM-DD format

                    if 'DIRECTION:' in line:
                        break
                    if 'POLYGON:' in line:  # after to skip 'POLYGON' line itself
                        poly_flag = ind + 4

                    if poly_flag and ind > poly_flag:
                        try:
                            values = [x for x in map(float, line.split()[:3])]
                            polypoints.append((values[1], values[2]))
                        except IndexError:
                            continue

            if exact:
                # Coord must be a point in the file (TODO not polypoints, filepoints!)
                if coord in polypoints:
                    plm_flag = True
            elif nearest:
                # TODO find_coord code, but over all files?
                if False:
                    plm_flag = True
            else:
                try:
                    if polygon_check_single(polypoints, coord, ccwsorted=True):  # Not actually ccw sorted necessarily, but they're sorted so it's fine
                        plm_flag = True
                except TypeError:
                    try:
                        centBox = (sum([b[0] for b in coord]) / len(coord), sum([b[1] for b in coord]) / len(coord))
                        coord.sort(key=lambda b: math.atan2(b[1] - centBox[1], b[0] - centBox[0]))

                        if polygon_check(polypoints, coord, ccwsorted=True):
                            plm_flag = True
                    except TypeError:
                        mssg = "Both polygon_checks errored. Please ensure `coord` argument is entered correctly as either 2-tuple or list of 2-tuples."
                        raise TypeError(mssg)

            if plm_flag:
                if not date or (date and date_acquired == date):
                    initdict['filepath'] = txt
                    initdict['identifier'] = os.path.basename(txt)[7:26]
                    initdict['orbit_num'] = initdict['identifier'][1:7]
                    initdict['block_num'] = initdict['identifier'][9:12]
                    initdict['plume_num'] = initdict['identifier'][17:19]
                    initdict['band'] = initdict['identifier'][16:17]
                    break
                else:
                    mssg = 'No such file or directory: No plume results containing this coordinate point on {0}. Please check date'.format(date)
                    raise IOError(mssg)

        try:
            return cls.from_filename(initdict['filepath'], out_dir=out_dir, _initdict=initdict, model_dir=model_dir, cal_int=cal_int)
        except (KeyError, NameError):
            mssg = 'No such file or directory: No plume results containing this coordinate point in {0}.'.format(out_dir)
            raise IOError(mssg)

    @classmethod
    def from_id(cls, orbit_num, block_num, plume_num=1, band=None, out_dir=os.getcwd(), model_dir="", cal_int=""):
        '''*Initialization method for reading plume based on a plume identifier*

        Initialization method which will search through a given directory for a plume text file corresponding to to the plume identifier requested.
        A plume identifier is given, by MINX, as the orbit number, block number, plume number and band of a digitized plume. The first two come
        from the actual MISR satellite overpass, while the latter two are assigned by MINX. Every plume will have a unique identifier.
        Used as a class method in order to namespace this function as a factory method, for consistency purposes.

        Args
        ----
        orbit_num : str or int
            MISR orbit number of the plume.
        block_num : str or int
            MISR block number containing the plume. Blocks range from 0 to 180, decreasing from the north pole to the south.
        plume_num : str or int, optional
            The number of the plume. Plumes are numbered based on the order of their digitization, for each block.
        band : {'B', 'R'}, optional
            The band the plume was digitized in. The majority of the time, it will be done in both. Therefore, unless you desire a specific band,
            the default option is sufficient.
        out_dir : str, optional
            Path of either the desired MINX project directory, or anywhere within the `Digitize_Output` directory
            of said project.
        model_dir : str, optional
            Path of the directory containing the model output files. Including this parameter will automatically attempt
            to open and attach the model data for this plume. See get_model for more information.

        Returns
        -------
        out : Plume
            A Plume object satisfying the specific requirements.

        Raises
        ------
        IOError: If no file satisfying the specific requirements is found.

        '''
        initdict = {}

        # --------------------------------------------------------------------------
        # Read and format all arguments, and find corresponding text file
        # --------------------------------------------------------------------------

        initdict['orbit_num'] = correct_digits(orbit_num, 6)
        initdict['block_num'] = correct_digits(block_num, 3)
        initdict['plume_num'] = correct_digits(plume_num, 2)

        # TODO messes up if user uses relative directories
        initdict['out_dir'] = out_dir[:out_dir.rfind('Digitize_Output') + 15] if 'Digitize_Output' in out_dir else out_dir + '/Digitize_Output'
        identifier = "O{orb}-B{blk}-SPW[{BR}]{plm}".format(
            orb=initdict['orbit_num'],
            blk=initdict['block_num'],
            BR=band if band else 'BR',  # TODO check if this actually works
            plm=initdict['plume_num']
        )

        plume_search = glob.glob("{0}/{1}/Plumes_{2}.txt".format(initdict['out_dir'], initdict['orbit_num'], identifier))

        if plume_search and len(plume_search) <= 2:
            initdict['filepath'] = plume_search[0]
            initdict['identifier'] = os.path.basename(initdict['filepath'])[7:26]
            initdict['band'] = initdict['identifier'][16:17]

            return cls.from_filename(initdict['filepath'], out_dir=out_dir, _initdict=initdict, model_dir=model_dir, cal_int=cal_int)

        elif len(plume_search) > 2:
            mssg = "Too many Plumes text files found with identifier < {0} > in < {1}/{2}/ >. Please inspect your output directory, \
            this error should not occur with normal MINX results".format(identifier, initdict['out_dir'], initdict['orbit_num'])
        else:
            mssg = "No Plumes text file found with identifier < {0} > in < {1}/{2}/ >".format(identifier, initdict['out_dir'], initdict['orbit_num'])
        raise IOError(mssg)

    @classmethod
    def from_filename(cls, filepath, out_dir=os.getcwd(), _initdict={}, model_dir="", cal_int=""):
        '''*Initialization method for reading plume based on a plume text file name*

        Factory initialization method which directly reads a given file and initializes a plume based off of a given file path. This is the actual
        init method, and all other factory methods simply supply a file path to this method. This is also the method which is returned by __repr__.


        Args
        ----
        filepath : str
            Path to the MINX results file of the desired plume.
        out_dir : str, optional
            Path of either the desired MINX project directory, or anywhere within the `Digitize_Output` directory
            of said project.
        model_dir : str, optional
            Path of the directory containing the model output files. Including this parameter will automatically attempt
            to open and attach the model data for this plume. See get_model for more information.
        cal_int : str, optional
            Path to the pickle file containing the relevant MISR/CALIOP intercept information. Including this
            parameter will automatically attempt to donwload and attach the CALIOP data near this plume. See get_caliop
            for more information.

        Returns
        -------
        out : Plume
            A Plume object satisfying the specific requirements.

        Raises
        ------
        IOError: If no file satisfying the specific requirements is found.

        '''
        if not _initdict:
            # if using this method directly, grab identifier
            _initdict['filepath'] = filepath
            _initdict['orbit_num'] = filepath[-22:-16]
            _initdict['block_num'] = filepath[-14:-11]
            _initdict['band'] = filepath[-7:-6]
            _initdict['plume_num'] = filepath[-6:-4]

            _initdict['out_dir'] = out_dir[:out_dir.rfind('Digitize_Output') + 15] if 'Digitize_Output' in out_dir else out_dir + '/Digitize_Output'
            _initdict['identifier'] = "O{orb}-B{blk}-SPW[{BR}]{plm}".format(
                orb=_initdict['orbit_num'],
                blk=_initdict['block_num'],
                BR=_initdict['band'],
                plm=_initdict['plume_num']
            )

        # --------------------------------------------------------------------------
        # Open file and read data
        # --------------------------------------------------------------------------

        try:
            with open(filepath, 'r') as file:
                file_attrs = {}
                _initdict['points'] = []
                _initdict['polygon'] = []
                _initdict['units'] = {}

                cont = file.read().split('\n')

                # --------------------------------------------------------------------------
                # Collect MetaData (up until line containing 'level 1 radiance file')
                # --------------------------------------------------------------------------

                # TODO I made this index thing better in a method above i think
                for line in cont[:[ind for ind, substr in enumerate(cont) if 'Level 1 radiance file' in substr][0]]:
                    try:
                        file_attrs[line.split(':', 1)[0].strip()] = str_or_flt(line.split(':', 1)[1].strip())
                    except IndexError:
                        continue

                _initdict['datetime'] = datetime.datetime.strptime(file_attrs['Date acquired'] + file_attrs['UTC time'], '%Y-%m-%d%H:%M:%S')
                _initdict['origin'] = (file_attrs['First point longitude'], file_attrs['First point latitude'])
                _initdict['biome'] = tuple(str_or_flt(i) for i in file_attrs['Biome IGBP name, class'].split(','))
                _initdict['region'] = tuple(str_or_flt(i) for i in file_attrs['Geographic region'].split(','))
                _initdict['comments'] = file_attrs['Comments by digitizer']   # TODO collect comments from plumeprojorbitlist.txt
                _initdict['stddev_ht'] = file_attrs['Ht std. deviation (m)']  # TODO error bars or something?
                _initdict['frp'] = file_attrs['Total fire power (MW)']
                _initdict['max_ht'] = file_attrs['Max ht (m > fire)']
                _initdict['area'] = file_attrs['Area (sq km)']

                # --------------------------------------------------------------------------
                # Collect Polygon Data (all after: line containing 'POLYGON' + 4
                #   until 'DIRECTION')
                # --------------------------------------------------------------------------

                for line in cont[[ind for ind, substr in enumerate(cont) if 'POLYGON:' in substr][0] + 4:]:
                    try:
                        if 'DIRECTION' in line:
                            break
                        else:
                            values = [x for x in map(float, line.split()[:3])]
                            _initdict['polygon'].append((values[1], values[2]))
                    except IndexError:
                        continue

                # --------------------------------------------------------------------------
                # Collect Point Data (all after: line containing 'RESULTS' + 4)
                # --------------------------------------------------------------------------

                data = [
                    ('longitude', 1, '\u00b0'),            # degrees
                    ('latitude', 2, '\u00b0'),             # degrees

                    ('distance', 6, 'km'),                 # kilometres, distance to initial point
                    ('direction', 7, '\u00b0'),            # degrees, azimuthal orientation measured cw from north
                    ('terrain', 8, 'm'),                   # metres, terrain height above sea level

                    ('nowind_height', 9, 'm'),             # metres, zero-wind height
                    ('wind_height', 10, 'm'),              # metres, wind-corrected height
                    ('filtered_height', 11, 'm'),          # metres, filtered and smoothed height

                    ('across_wind', 12, 'm/s'),            # metres/second, windspeed across swath (E/W)
                    ('along_wind', 13, 'm/s'),             # metres/second, windspeed along swath (N/S)
                    ('total_wind', 14, 'm/s')              # metres/second, windspeed along direction of plume
                ]

                for datapoint in data:
                    _initdict[datapoint[0]] = []
                    _initdict['units'][datapoint[0]] = datapoint[2]

                for line in cont[[ind for ind, substr in enumerate(cont) if 'RESULTS:' in substr][0] + 4:]:
                    try:
                        # TODO ensure no other negative values are possible (other than lat long)
                        values = [x if (x >= 0 or i <= 2) else float('nan') for i, x in enumerate(map(float, line.split()))]
                        point_vals = {}

                        for datapoint in data:
                            point_vals[datapoint[0]] = values[datapoint[1]]
                            _initdict[datapoint[0]].append(values[datapoint[1]])

                        _initdict['points'].append(point_vals)

                    except IndexError:
                        continue
        except IOError:
            mssg = "No such Plume file: {0}".format(filepath)
            raise IOError(mssg)

        # --------------------------------------------------------------------------
        # Initialize Plume object with _initdict data
        # --------------------------------------------------------------------------

        return cls(model_dir=model_dir, cal_int=cal_int, **_initdict)

    def __init__(self, model_dir="", cal_int="", *args, **kwargs):
        super(Plume, self).__init__()

        if not kwargs:
            mssg = 'to create Plume object must use Plume.from_file(orb, blk, plm, band) or Plume.from_coord((lon, lat)).'
            raise TypeError(mssg)
        else:
            for option in kwargs:
                setattr(self, option, kwargs[option])

        self.reference = 'sea'

        if model_dir:
            self.model_dir = model_dir
            self.get_model(model_dir)

        if cal_int:
            self.get_caliop(cal_int)

    # --------------------------------------------------------------------------
    # Plotting methods
    # --------------------------------------------------------------------------
    # TODO error checking and argument parsing on all methods, prettify all plots, add terrain option to 3d?
    # makse sure to include reference and units in axis labels
    # TODO subplots into squares?

    def contour(self, option='filtered_height', fill=True, save=None, filename=None, *args, **kwargs):
        '''*3D Plotting method for creating a contour plot of a given option*

        3D plotting method for creating, on a latitude-longitude grid, a contour plot, for a given attribute.
        Uses natural neighbor interpolation based on Delaunay triangulation to grid the desired attributes.
        A normal contour, or a filled contour plot can both be created, although the non-fill contour will usually provide
        non-sensical results.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to be plotted. For a list of valid attributes see Plume class docstring.
        fill : bool, optional
            If specified, the plot will be a filled contour plot, rather than a line contour plot. It is recommmended to leave this option
            to the default value.
        save : ternary, optional
            If True, will save plot to a given `filename`, if False, will display plot, if None will return matplotlib object.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib contourf function

        Returns
        -------
        out : Nonetype or [array of] matplotlib.AxesSubplot
            Returns None or a [numpy array of] matplolib subplot axes objects
        (Will display or save plot if returning None)

        Examples
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.contour()  #doctest: +ELLIPSIS
        <matplotlib.axes.AxesSubplot object at 0x...>

        '''
        import numpy as np

        option = arg_to_list(option)

        fig, axarr = plt.subplots(1, len(option))
        if len(option) == 1:
            axarr = [axarr]

        x, y = self.longitude, self.latitude
        xgrid, ygrid = np.meshgrid(x, y)

        for ind, opt in enumerate(option):
            ax = axarr[ind]
            ax.set_title(opt)

            if not ind:
                ax.set_ylabel('Latitude ({0})'.format(self.units['latitude']))
            ax.set_xlabel('Longitude ({0})'.format(self.units['longitude']))

            if 'model' in opt.lower():
                try:
                    data = self.Model.height
                except AttributeError:
                    mssg = "This plume has no corresponding model data. Gather said data by using the 'get_model' method"
                    raise AttributeError(mssg)
            else:
                data = self[opt]

            zgrid = plt.mlab.griddata(x, y, data, xgrid, ygrid)

            if fill:
                CS = ax.contourf(xgrid, ygrid, zgrid, *args, **kwargs)
            else:
                CS = ax.contour(xgrid, ygrid, zgrid, *args, **kwargs)

            cbar = plt.colorbar(CS)
            if 'height' in opt:
                cbar.ax.set_ylabel("Metres above {lvl} level".format(lvl=self.reference))

        if save:
            if not filename:
                filename = raw_input('Enter desired filename: ')
            fig.savefig(filename)
        else:
            if save is None:
                return axarr if len(option) > 1 else axarr[0]

            plt.show()

        fig.clf()
        plt.close(fig)

        return

    def colour(self, option='filtered_height', save=None, filename=None, *args, **kwargs):
        '''*3D Plotting method for creating a colour plot of a given option*

        3D plotting method for creating, on a latitude-longitude grid, a colour plot, for a given attribute.
        Each point will be displayed as a coloured square, based on the desired attribute.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to be plotted. For a list of valid attributes see Plume class docstring.
        save : ternary, optional
            If True, will save plot to a given `filename`, if False, will display plot, if None will return matplotlib object.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib scatter function

        Returns
        -------
        out : Nonetype or [array of] matplotlib.AxesSubplot
            Returns None or a [numpy array of] matplolib subplot axes objects
        (Will display or save plot if returning None)

        Examples
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.colour()  #doctest: +ELLIPSIS
        <matplotlib.axes.AxesSubplot object at 0x...>

        '''
        option = arg_to_list(option)

        # TODO caliop check better
        try:
            del option[option.index('caliop')]
            caliop_comp = True
        except ValueError:
            caliop_comp = False

        # TODO share axes
        fig, axarr = plt.subplots(1, len(option))
        if len(option) == 1:
            axarr = [axarr]

        for ind, opt in enumerate(option):
            ax = axarr[ind]
            ax.set_title(opt, fontsize=14)

            if not ind:
                ax.set_ylabel('Latitude ({0})'.format(self.units['latitude']), fontsize=14)
            ax.set_xlabel('Longitude ({0})'.format(self.units['longitude']), fontsize=14)

            if 'model' in opt.lower():
                try:
                    data = self.Model.height
                except AttributeError:
                    mssg = "This plume has no corresponding model data. Gather said data by using the 'get_model' method"
                    raise AttributeError(mssg)
            else:
                data = self[opt]

            ax.scatter(self.longitude, self.latitude, c=data, s=90, marker='s', edgecolors='none', *args, **kwargs)

            if caliop_comp:
                try:
                    ax.scatter(self.Caliop.longitude, self.Caliop.latitude, c=self.Caliop.height, edgecolors='none', *args, **kwargs)
                except AttributeError:
                    mssg = "This plume has no corresponding Caliop data. Gather said data by using the 'get_caliop' method"
                    raise AttributeError(mssg)

        cur_axes = plt.gca()
        PCM = ax.get_children()[2]  # all labels should probs be dynamic
        plt.colorbar(PCM, ax=cur_axes, label="Metres above {lvl} level".format(lvl=self.reference))

        if save:
            if not filename:
                filename = raw_input('Enter desired filename: ')
            fig.savefig(filename)
        else:
            if save is None:
                return axarr if len(option) > 1 else axarr[0]

            plt.show()

        fig.clf()
        plt.close(fig)

        return

    def scatter_3D(self, option='filtered_height', origin=True, save=None, filename=None, *args, **kwargs):
        '''*3D Plotting method for creating a scatter plot of a given option*

        3D plotting method for creating, with a longitude-latitude x-y grid, a 3 dimensional scatter plot, for a given attribute.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to be plotted. For a list of valid attributes see Plume class docstring
        origin : bool, optional
            If specified, the location of the plume origin will be displayed as well, by a large `X`.
        save : ternary, optional
            If True, will save plot to a given `filename`, if False, will display plot, if None will return matplotlib object.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib plot function

        Returns
        -------
        out : Nonetype or matplotlib.AxesSubplot
            Returns None or a matplolib subplot axes objects
        (Will display or save plot if returning None)

        Examples
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.scatter_3D()  #doctest: +ELLIPSIS
        <matplotlib.axes.AxesSubplot object at 0x...>

        '''
        from mpl_toolkits.mplot3d import Axes3D

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.autoscale()
        # ax.set_zlim(zmin=0)

        # TODO caliop check better
        try:
            del option[option.index('caliop')]
            caliop_comp = True
        except ValueError:
            caliop_comp = False

        option = arg_to_gen(option)

        for opt in option:

            if 'model' in opt.lower():
                try:
                    data = self.Model.height
                except AttributeError:
                    mssg = "This plume has no corresponding model data. Gather said data by using the 'get_model' method"
                    raise AttributeError(mssg)
            else:
                data = self[opt]

            ax.plot(self.longitude, self.latitude, data, 'o', color=next(self.clr), label=opt, *args, **kwargs)

        if origin:
            ax.plot([self.origin[0]], [self.origin[1]], [self.terrain[0]], 'x', mew=4, ms=8, color=next(self.clr), label='Plume origin')

        if caliop_comp:
            # TODO axes when caliop
            try:
                ax.plot(self.Caliop.longitude, self.Caliop.latitude, self.Caliop.height, 'o', color=next(self.clr), label='caliop', *args, **kwargs)
            except AttributeError:
                mssg = "This plume has no corresponding Caliop data. Gather said data by using the 'get_caliop' method"
                raise AttributeError(mssg)

        ax.legend(numpoints=1)

        if save:
            if not filename:
                filename = raw_input('Enter desired filename: ')
            fig.savefig(filename)
        else:
            if save is None:
                return ax

            plt.show()

        fig.clf()
        plt.close(fig)

        return

    def surface_plot(self, option='filtered_height', save=None, filename=None, *args, **kwargs):
        '''*3D Plotting method for creating a surface plot of a given option*

        3D plotting method for creating, with a longitude-latitude x-y grid, a 3 dimensional surface plot, for a given attribute.
        Uses Delauney triangulation to connect all points of desired attribute. This method will not give very sensical results for most
        MINX plumes, unless they are extremely uniform.


        Args
        ----
        option : str or list of str
            Plume attribute(s) to be plotted. For a list of valid attributes see Plume class docstring
        save : ternary, optional
            If True, will save plot to a given `filename`, if False, will display plot, if None will return matplotlib object.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib plot_trisurf function

        Returns
        -------
        out : Nonetype or matplotlib.AxesSubplot
            Returns None or a matplolib subplot axes objects
        (Will display or save plot if returning None)

        Examples
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.surface_plot()  #doctest: +ELLIPSIS
        <matplotlib.axes.AxesSubplot object at 0x...>

        '''
        from mpl_toolkits.mplot3d import Axes3D
        # TODO this looks terrible

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.autoscale()
        # ax.set_zlim(zmin=0)

        option = arg_to_gen(option)

        for opt in option:

            if 'model' in opt.lower():
                try:
                    data = self.Model.height
                except AttributeError:
                    mssg = "This plume has no corresponding model data. Gather said data by using the 'get_model' method"
                    raise AttributeError(mssg)
            else:
                data = self[opt]

            x, y, z = clean_inv(self.longitude, self.latitude, data)
            ax.plot_trisurf(x, y, z, color=next(self.clr), label=opt, *args, **kwargs)

        ax.legend()

        if save:
            if not filename:
                filename = raw_input('Enter desired filename: ')
            fig.savefig(filename)
        else:
            if save is None:
                return ax

            plt.show()

        fig.clf()
        plt.close(fig)

        return

    def direction_plot(self, origin=True, boundaries=True, save=None, filename=None, *args, **kwargs):
        '''*2D Plotting method for creating a direction field of MINX wind data*

        2D plotting method for creating, with a longitude-latitude x-y grid, a direction field plot displaying the speed and
        direction of the MINX wind products.

        Args
        ----
        origin : bool, optional
            If specified, the location of the plume origin will be displayed as well, by a large `X`.
        save : ternary, optional
            If True, will save plot to a given `filename`, if False, will display plot, if None will return matplotlib object.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib quiver function

        Returns
        -------
        out : Nonetype or matplotlib.AxesSubplot
            Returns None or a matplolib subplot axes objects
        (Will display or save plot if returning None)

        Examples
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.direction_plot()  #doctest: +ELLIPSIS
        <matplotlib.axes.AxesSubplot object at 0x...>

        '''
        # not sure why I need all these negatives, probably don't
        u = [math.sin((90 - p['direction']) * math.pi / 180.) * p['total_wind'] for p in self.points]
        v = [math.cos((90 - p['direction']) * math.pi / 180.) * p['total_wind'] for p in self.points]

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.quiver(self.longitude, self.latitude, u, v, *args, **kwargs)

        if boundaries:
            plm_x, plm_y = zip(*self.polygon)
            ax.plot(plm_x, plm_y, 'o-', color=next(self.clr), label='Plume Boundaries', *args, **kwargs)

        if origin:
            ax.plot([self.origin[0]], [self.origin[1]], 'x', mew=4, ms=8, color=next(self.clr), label='Plume origin')

        if save:
            if not filename:
                filename = raw_input('Enter desired filename: ')
            fig.savefig(filename)
        else:
            if save is None:
                return ax

            plt.show()

        fig.clf()
        plt.close(fig)

        return

    def scatter_2D(self, option='filtered_height', axis='distance', show_terrain=True, save=None, filename=None, *args, **kwargs):
        '''*2D Plotting method for creating a scatter plot of a given option*

        2D plotting method for creating a 2 dimensional scatter plot, for a given attribute. The x-axis can be specified with the `axis` arg.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to be plotted. For a list of valid attributes see Plume class docstring.
        axis : {'distance', 'longitude', 'latitude'}, optional
            Specifies which attribute will be set as the x-axis the heights will be plotted versus.
        show_terrain : bool, optional
            If specified, will plot a green line indicating the terrain height for each point of the desired attribute.
        save : ternary, optional
            If True, will save plot to a given `filename`, if False, will display plot, if None will return matplotlib object.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib plot function

        Returns
        -------
        out : Nonetype or matplotlib.AxesSubplot
            Returns None or a matplolib subplot axes objects
        (Will display or save plot if returning None)

        Examples
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.scatter_2D()  #doctest: +ELLIPSIS
        <matplotlib.axes.AxesSubplot object at 0x...>

        '''

        fig = plt.figure()
        ax = fig.add_subplot(111)

        if show_terrain:
            terr_axis, terr_vals = zip(*sorted(zip(self[axis], self['terrain'])))
            ax.plot(terr_axis, terr_vals, 'g-', *args, **kwargs)

        option = arg_to_gen(option)

        for opt in option:

            if 'model' in opt.lower():
                try:
                    data = self.Model.height
                except AttributeError:
                    mssg = "This plume has no corresponding model data. Gather said data by using the 'get_model' method"
                    raise AttributeError(mssg)
            else:
                data = self[opt]

            if opt != 'terrain' or not show_terrain:
                ax.plot(self[axis], data, 'o', color=next(self.clr), label=opt, *args, **kwargs)
                #ax.set_ylim(ymin=0)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax.legend(numpoints=1, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=4)

        if save:
            if not filename:
                filename = raw_input('Enter desired filename: ')
            fig.savefig(filename)
        else:
            if save is None:
                return ax

            plt.show()

        fig.clf()
        plt.close(fig)

        return

    def poly_2D(self, option='filtered_height', axis='distance', deg=8, show_terrain=True, save=None, filename=None, *args, **kwargs):
        '''*2D Plotting method for creating a linear and polynomial plot of a given option*

        2D plotting method for creating a 2 dimensional linear plot, for a given attribute, as well as a N-th degree polynomial fit of said plot.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to be plotted. For a list of valid attributes see Plume class docstring
        axis : {'distance', 'longitude', 'latitude'}, optional
            Specifies which attribute will be set as the x-axis the heights will be plotted versus.
        deg : int, optional
            Specifies the degree of the fitting polynomial.
        show_terrain : bool, optional
            If specified, will plot a green line indicating the terrain height for each point of the desired attribute.
        save : ternary, optional
            If True, will save plot to a given `filename`, if False, will display plot, if None will return matplotlib object.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib plot function

        Returns
        -------
        out : Nonetype or matplotlib.AxesSubplot
            Returns None or a matplolib subplot axes objects
        (Will display or save plot if returning None)

        Examples
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.poly_2D()  #doctest: +ELLIPSIS
        <matplotlib.axes.AxesSubplot object at 0x...>

        '''
        import numpy as np

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.autoscale()
        # ax.set_ylim(ymin=0)

        if show_terrain:
            terr_axis, terr_vals = zip(*sorted(zip(self[axis], self['terrain'])))
            ax.plot(terr_axis, terr_vals, 'g-', *args, **kwargs)

        option = arg_to_gen(option)

        for opt in option:

            if opt != 'terrain' or not show_terrain:

                if 'model' in opt.lower():
                    try:
                        data = self.Model.height
                    except AttributeError:
                        mssg = "This plume has no corresponding model data. Gather said data by using the 'get_model' method"
                        raise AttributeError(mssg)
                else:
                    data = self[opt]

                opt_axis, opt_vals = zip(*sorted(zip(*clean_inv(self[axis], data))))

                ax.plot(opt_axis, opt_vals, color=next(self.clr), label=opt, *args, **kwargs)

                z = np.polyfit(opt_axis, opt_vals, deg)
                f = np.poly1d(z)

                axis_smooth = np.linspace(min(self[axis]), max(self[axis]), 300)
                option_smooth = f(axis_smooth)

                ax.plot(axis_smooth, option_smooth, '-', color=next(self.clr), label='smoothed ' + opt, *args, **kwargs)

        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
        ax.legend(numpoints=1, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=4)

        if save:
            if not filename:
                filename = raw_input('Enter desired filename: ')
            fig.savefig(filename)
        else:
            if save is None:
                return ax

            plt.show()

        fig.clf()
        plt.close(fig)

        return

    def shape(self, fill=True, save=None, filename=None, *args, **kwargs):
        '''*2D Plotting method for creating a filled in contour of the geographic shape of the digitized plume*

        2D plotting method for creating a 2 dimensional contour plot with no vertical information. The plot is only for showcasing the
        shape of the digitized plume. Plumes will be rotated so that a geographic line drawn directly between the origin and the
        furthest point will sit on the x-axis.

        Args
        ----
        fill : boolean, optional
            Determine wether the shape of the plume should be filled in with a semi-transparent colour or not.
        save : ternary, optional
            If True, will save plot to a given `filename`, if False, will display plot, if None will return matplotlib object.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib fill function

        Returns
        -------
        out : Nonetype or matplotlib.AxesSubplot
            Returns None or a matplolib subplot axes objects
        (Will display or save plot if returning None)

        Examples
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.shape()  #doctest: +ELLIPSIS
        <matplotlib.axes.AxesSubplot object at 0x...>

        '''
        # TODO return some value of 'linearity'? using end point as x-axis may be incorrect, maybe should use most linear part, or initial direction

        # plume stats methods only work for plume points, not polygon points
        # TODO is haverdistance necessary? dont think so, probs a distance in degrees is fine
        dist_arr = [(ind, haver_dist(pnt[0], pnt[1], self.origin[0], self.origin[1])) for ind, pnt in enumerate(self.polygon)]
        min_index, _ = min(dist_arr, key=lambda dist: dist[1])
        max_index, _ = max(dist_arr, key=lambda dist: dist[1])

        orig_x, orig_y = self.polygon[min_index]

        poly_x, poly_y = zip(*self.polygon)
        trans_x, trans_y = [x - orig_x for x in poly_x], [y - orig_y for y in poly_y]

        end_x, end_y = trans_x[max_index], trans_y[max_index]
        theta = math.atan2(end_y, end_x)

        rot_x = [(x * math.cos(-theta)) - (trans_y[ind] * math.sin(-theta)) for ind, x in enumerate(trans_x)]
        rot_y = [(trans_x[ind] * math.sin(-theta)) + (y * math.cos(-theta)) for ind, y in enumerate(trans_y)]

        fig = plt.figure()
        ax = fig.add_subplot(111)

        plt.plot(0, 0, 'ro')
        if fill:
            ax.fill(rot_x, rot_y, next(self.clr), alpha=0.5)
        else:
            ax.plot(rot_x, rot_y, next(self.clr))

        if save:
            if not filename:
                filename = raw_input('Enter desired filename: ')
            fig.savefig(filename)
        else:
            if save is None:
                return ax

            plt.show()

        fig.clf()
        plt.close(fig)

        return

    def hist(self, option='filtered_height', save=None, filename=None, *args, **kwargs):
        '''*2D Plotting method for creating a histogram plot of a given option*

        2D plotting method for creating a histogram plot, for a given attribute. If multiple attributes are requested, will be placed in
        different subplots.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to be plotted. For a list of valid attributes see Plume class docstring
        save : ternary, optional
            If True, will save plot to a given `filename`, if False, will display plot, if None will return matplotlib object.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib hist function

        Returns
        -------
        out : Nonetype or [array of] matplotlib.AxesSubplot
            Returns None or a [numpy array of] matplolib subplot axes objects
        (Will display or save plot if returning None)

        Examples
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.contour()  #doctest: +ELLIPSIS
        <matplotlib.axes.AxesSubplot object at 0x...>

        '''
        option = arg_to_list(option)

        fig, axarr = plt.subplots(1, len(option))
        if len(option) == 1:
            axarr = [axarr]

        for ind, opt in enumerate(option):

            if 'model' in opt.lower():
                try:
                    data = self.Model.height
                except AttributeError:
                    mssg = "This plume has no corresponding model data. Gather said data by using the 'get_model' method"
                    raise AttributeError(mssg)
            else:
                data = self[opt]

            ax = axarr[ind]
            ax.set_title(opt)
            ax.hist(clean_inv(data), *args, **kwargs)

        if save:
            if not filename:
                filename = raw_input('Enter desired filename: ')
            fig.savefig(filename)
        else:
            if save is None:
                return axarr if len(option) > 1 else axarr[0]

            plt.show()

        fig.clf()
        plt.close(fig)

        return

    # --------------------------------------------------------------------------
    # Statistical methods
    # --------------------------------------------------------------------------

    # TODO make work with model stuff, could make work with multiple options and weightings

    def max(self, option='filtered_height'):
        '''*Return point data for which a given option is the max for this plume*'''
        from numpy import nanmax
        ind = self[option].index(nanmax(self[option]))
        return self.points[ind], ind

    def min(self, option='filtered_height'):
        '''*Return point data for which a given option is the min for this plume*'''
        from numpy import nanmin
        ind = self[option].index(nanmin(self[option]))
        return self.points[ind], ind

    def median(self, option='filtered_height'):
        '''*Return the median value of a given option for all of this plume*'''
        from numpy import median
        return median(self[option])

    def mean(self, option='filtered_height'):
        '''*Return the mean value of a given option for all of this plume*'''
        from numpy import nanmean
        return nanmean(self[option])

    def std(self, option='filtered_height', ddof=0):
        '''*Return the population (0) or sample (1) standard deviation for a given option for all of this plume*'''
        from numpy import nanstd
        return nanstd(self[option], ddof=ddof)

    # --------------------------------------------------------------------------
    # Technical methods
    # --------------------------------------------------------------------------

    def get_model(self, fst_dir, var='AF', threshold=4.0):
        '''*Initializaton method for the Model Plume objects*

        Initializes a ModelPlume object corresponding to this MINX generated plume.
        Function does not return said ModelPlume, and it should only be accessed through the corresponding Plume object.

        Can either be called from an existing Plume, or during Plume initialization, if model_dir is provided to
        the initialization methods. ModelPlume's *must* be constructed using this method, as they should not be accessed
        apart from their corresponding MINx plumes.

        ModelPlume's are generated based on fst output files. In order to construct the plumes, you must ensure that the
        correct path to these files is provided and that they cover the timespan and location of the plume.

        See the documentation of the _ModelPlume class for a better description of the results of this function.

        Args
        ----
        fst_dir : str
            Path to the directory containing the model fst output files. Files must be named according to the time of the simulation.
        var : str, optional
            The name of the product variable which will be used to determine plume heights. It is recommended to leave this as PM2.5
            concentration ('AF').
        threshold : float, optional
            The minimum lower threshold value of the product, used to determine plume heights. It is recommended to leave this as the
            default 4.0, while being used with PM2.5 concentration.

        Returns
        -------
        out : Nonetype
            Returns None.
        (A reference to the model object will be stored in Plume().Model)

        '''
        # TODO return products if they already existed (if this is from a slice method)
        self.Model = _ModelPlume(plm=self, fst_dir=fst_dir, var=var, threshold=threshold)
        return

    def get_caliop(self, intersect, cal_path=os.getcwd(), verbose_error=False):
        '''
        get caliop data
        '''
        self.Caliop = _Caliop(plm=self, intersect=intersect, cal_path=cal_path, verbose_error=verbose_error)
        return

    def set_reference(self, reference=None):
        '''*Change reference level of all heights between `above-terrain` and `ASL`*

        Convert all heights stored within this plume and the corresponding model plume to the reference level desired.
        This reference level can be either `sea` (above sea level) or `terrain` (above terrain level). If no reference level is
        specified, will simply switch to the opposite of what is currently in use. To see what level is currently being used,
        see the Plume().reference attribute.

        Args
        ----
        reference : str, optional
            Reference level to convert all heights to. If not specified, will switch to level not currenlty in use.

        Returns
        -------
        out : Nonetype
            Returns None.
        (Will change all stored height values)

        '''
        # TODO make for all units?, make sure only doing height, do for other model products

        if not reference:
            reference = 'terrain' if self.reference == 'sea' else 'sea'
        elif reference == self.reference:
            UserWarning('Requested reference level is already in use, changing nothing.')
            return

        self.reference = reference

        # set the change to be either -terrain or +terrain
        change = [(-1) * terr for terr in self.terrain] if reference == 'terrain' else self.terrain

        for option, unit in self.units.items():
            if unit == 'm' and option != 'terrain':

                if 'model' in option.lower():
                    obj = self.Model
                    opt = 'height'
                elif 'caliop' in option.lower():
                    obj = self.Caliop
                    opt = 'height'
                else:
                    obj = self
                    opt = option

                new_vals = [val + change[ind] for ind, val in enumerate(obj[opt])]

                setattr(obj, opt, new_vals)

                try:
                    for pind, point in enumerate(obj['points']):
                        point[opt] = point[opt] + change[pind]
                except KeyError:
                    # TODO caliop, if possible (maybe you flags to find height of terrain starts)?
                    continue

    def find_coord(self, coord, exact=False):
        '''*Search for a given coordinate point within a plume*

        Find the nearest digitized point to a given longitude-latitude coordinate within a plume.
        Uses a linear exact nearest neighbour search, using the haversine formula to determine
        distance between points on a globe.

        Args
        ----
        coord : tuple
            Coordinate point to search plume for. Must be given as (longitude, latitude).
        exact : bool, optional
            If specified, will return None if `coord` does not match a digitized point in the plume exactly.

        Returns
        -------
        res : dict
            The point dictionary corresponding to the nearest point to `coord`. Contains all attributes pertaining to said point.
        model : float
            The model height corresponding to the nearest point to `coord`.
        dist : float
            The distance, in metres, between the nearest point and `coord`.
        ind : int
            The index corresponding to the nearest point to `coord`. Can be used in all plume attributes to return values
            pertaining to this point.

        Examples TODO
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.find_coord((100., 50))
        ({TODOdict}, 1000.3, 34.4, 6)

        '''
        dist = float('inf')
        res = None
        ind = None

        R = 6371000  # Radius of Earth, in metres
        phi_c = coord[1] * math.pi / 180.

        for i, p in enumerate(self.points):

            if exact and coord == (p['longitude'], p['latitude']):
                return (p, self.Model.height[i], 0.0, i)
            elif not exact:

                phi_p = p['latitude'] * math.pi / 180.
                del_phi = (p['latitude'] - coord[1]) * math.pi / 180.
                del_lam = (p['longitude'] - coord[0]) * math.pi / 180.

                a = math.sin(del_phi / 2.)**2 + math.cos(phi_p) * math.cos(phi_c) * math.sin(del_lam / 2.)**2
                c = R * 2 * math.atan2(a**0.5, (1 - a)**0.5)

                res, dist, ind = (p, c, i) if c < dist else (res, dist, ind)

        try:
            model = self.Model.height[ind]
        except AttributeError:
            model = None

        return (res, model, dist, ind) if not exact else None

    def coord_section(self, bbox, concave=False):
        '''*Given a bounding polygon, create a new plume using the new coordinate slice of points*

        Given an N-dimensional bounding polygon of longitude-latitude coordinate points, create a new plume object, utilizing all
        the same data and metadata from the old plume, and cutting out all points outside the geographic bounding box specified.
        Polygon must either be a convex polygon, or, if it is to be concave, must be sorted when passed to this method. This is
        because concave polygon sorting is impossible without source information.

        If the old plume had a corresponding model plume, said model plume will also be initialized.
        Be aware that this function may misbehave for coordinates near the ante-meridian and the poles, and that the metadata of
        the new plume may be slightly off in areas.

        Note
        ----
        Sorting of polygons is currently not supported, please sort your bounding boxes before passing to this function.

        Args
        ----
        bbox : list of tuples
            N-dimensional list of lat-long coordinates used to create a new plume from this geographic bounding polygon.
        concave : bool, optional
            If given bounding polygon is to be concave, must specify this flag and ensure the the polygon is already sorted.

        Returns
        -------
        out : Plume
            A Plume object satisfying the specific geographic bounding box requirements.

        Examples  TODO
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.coord_section([(1, 0), (0, 0), (0, 1), (1, 1)])  #doctest: +ELLIPSIS
        <postMINX.Plume object at 0x...>

        '''
        # TODO keep all model products too
        try:
            from Cpoly import polygon_check_single
        except ImportError:
            # running on python2, use python implementation of polygon checking from preMINX
            from preMINX import polygon_check_single
        from copy import deepcopy

        newattrs = deepcopy(self.__dict__)
        newattrs['polygon'] = bbox
        try:
            newattrs.pop('Model')
        except KeyError:
            pass

        inds = []

        for ind, p in enumerate(self.points):
            if not polygon_check_single(bbox, point=(p['longitude'], p['latitude']), ccwsorted=concave):
                inds.append(ind)

        for item in self.__dict__:
            if isinstance(self.__dict__[item], list) and not (item == 'polygon' or item == 'Model'):
                for ind in reversed(inds):
                    newattrs[item].pop(ind)

        try:
            model_dir = self.Model.dir
        except AttributeError:
            model_dir = ""

        return Plume(model_dir=model_dir, **newattrs)

    def cross_section(self, option='filtered_height', min_val=0., max_val=float('inf'), med_dev=False):
        '''*Given a bounding polygon, create a new plume using the new coordinate slice of points*

        Given minimum and maximum values, create a new plume object, utilizing all the same data and metadata from the old plume,
        and cutting out all points with an `option` value outside the cross section of values specified.
        If the old plume had a corresponding model plume, said model plume will also be initialized.

        Args
        ----
        option : str
            Plume attribute to be filtered on and the cross section cut from. For a list of valid attributes see Plume class docstring
        min_val : float, optional
            Minimum value of `option` for the cross section to be cut based on.
        max_val : float, optional
            Maximum value of `option` for the cross section to be cut based on.
        med_dev : bool, optional
            NotImplementedError

        Returns
        -------
        out : Plume
            A Plume object satisfying the specific `option` value requirements.

        Examples  TODO
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.cross_section(min_val=1000., max_val=2000.)  #doctest: +ELLIPSIS
        <postMINX.Plume object at 0x...>

        '''
        # TODO add sigma-esque method to cut outliers, using this method and std?, cut on multiple options?

        from copy import deepcopy

        if 'model' in option.lower():
            opt = self.Model.value if 'value' in option else self.Model.height
        else:
            opt = self[option]

        newattrs = deepcopy(self.__dict__)
        newattrs.pop('Model')
        newattrs.pop('model_dir')

        inds = []
        for ind, optval in enumerate(opt):
            if not (min_val < optval < max_val):
                inds.append(ind)

        for item in self.__dict__:
            if isinstance(self.__dict__[item], list) and not item == 'polygon':
                for ind in reversed(inds):
                    newattrs[item].pop(ind)

        try:
            model_dir = self.Model.dir
        except AttributeError:
            model_dir = ""

        return Plume(model_dir=model_dir, **newattrs)

    def set_release_time(self, intervalsize=0.1):
        """
        determine, using the distance, wind speed and direction, the *approximate* time that the smoke at each point may have been produced

        works by creating a 'simulation' line from the point to the origin and following it in <intervalsize> degrees intervals along the longitude axis.
        each time, finding the nearest digitized point, determines the windspeed velocity along the direction of the simulated line
        by rotating the nearest vector to have a parallel x-axis to the sim line, and taking its x-component.

        This is a huge approximation, for many reasons, and only currently functions for near-linear plumes.
        """
        import numpy as np
        import pyloader
        # TODO only for linearish plumes right now, need to find way to draw a polynomial line
        # TODO plot to test if getting Vwind correctly (think I'm doing this backwards)
        # TODO make sure the angles all work, no matter what they are (ie test different plumes)
        orgn = self.origin
        self.release_time = []
        load = pyloader.Load_bar(len(self.points))
        for pind, pnt in enumerate(self.points):
            load.update(pind)
            # draw line from point to origin

            if pnt['longitude'] - orgn[0]:
                m = (pnt['latitude'] - orgn[1]) / (pnt['longitude'] - orgn[0])
                b = pnt['latitude'] - ((pnt['latitude'] - orgn[1]) / (pnt['longitude'] - orgn[0])) * pnt['longitude']

            elif (pnt['longitude'], pnt['latitude']) == orgn:
                self.release_time.append(self['datetime'])
                self.points[pind]['release_time'] = self['datetime']
                continue
            else:
                m = b = None  # TODO handle the vert/hor lines

            # throughtout line, keep finding nearest point, on simulation from origin to point
            sim_range = np.arange(orgn[0], pnt['longitude'], intervalsize if orgn[0] < pnt['longitude'] else -intervalsize)

            # vector, magnitude from first sim_pnt to orgn  # to find angle of sim line, to find correct velocity components
            Vsx, Vsy = (orgn[0] - pnt['longitude'], orgn[1] - pnt['latitude'])
            dist = math.sqrt(((Vsx) * (Vsx)) + ((Vsy) * (Vsy)))

            # angle of sim line, to rotate nearest vector about, and also to find distance traveled, which is uniform about sim_range
            theta = math.acos((Vsx) / dist)

            # distance travelled each simulation *in metres*, does not change throughout sim line, *as long as this is linear*
            dist_per_sim = pnt['distance'] * 1000 / len(sim_range)

            last_nearest = None
            velocity_bins = []  # (velocity, distance) in each bin (changes each time a new nearest point, and thus speed, is found)
            vind = -1

            for sim_x in sim_range:
                sim_pnt = (sim_x, m * sim_x + b)
                nearest_pnt = self.find_coord(sim_pnt)[0]

                if nearest_pnt != last_nearest:
                    # if nearest point changes, find new velocity
                    # take speed at each sim_point to be along line component of speed of total wind of nearest point

                    # vector of nearest point wind velocity  # along sim component of this is the sim_velocity
                    phi = (90 - nearest_pnt['direction']) * math.pi / 180.  # TODO is this wrong?
                    Vnx, Vny = (nearest_pnt['total_wind'] * math.cos(phi), nearest_pnt['total_wind'] * math.sin(phi))

                    # rotated nearest vector
                    Vrx = Vnx * math.cos(-theta) - Vny * math.sin(-theta)
                    Vry = Vnx * math.sin(-theta) + Vny * math.cos(-theta)

                    # angle of rotated nearest vector
                    alpha = phi - theta

                    # x-component of rotated nearest vector (ie the component of nearest vector along sim line direction)
                    # for some reason the velocities are always negative, there's a flipped angle somewhere but idk where (doesnt matter tho)
                    Vwind = math.sqrt(Vrx * Vrx + Vry * Vry) * math.cos(alpha)
                    vind += 1
                    # TODO handle nans better (currently discounts the entire bin, should maybe find another point)
                    velocity_bins.append([Vwind, 0.])

                velocity_bins[vind][1] += dist_per_sim
                last_nearest = nearest_pnt

            # delta_t = sum(distance intervals / speed during intervals) (in seconds)
            delta_t = datetime.timedelta(seconds=abs(sum([vbin[1] / vbin[0] for vbin in velocity_bins if vbin[0] and not math.isnan(vbin[0])])))

            # initial release time = MISR time - delta_t
            self.release_time.append(self['datetime'] - delta_t)
            self.points[pind]['release_time'] = self['datetime'] - delta_t

        load.stop()

    def drawkml(self, filepath, *args):
        '''
        *must have simplekml installed*
        draw a kml with your choice of:
            - ('MISR_{cam}') -> MISR camera image (any camera)
            - ('CALIOP') -> calipso path (and height?)
            - ('FirePixel') -> MODIS fire pixels
            - ('dict({method}: [*args])') -> any of the above plotting methods which can be 2d overlaid (maybe even 3d, we'll see)
            -
        all options must be given as string *args
        '''
        try:
            import simplekml
        except ImportError:
            mssg = 'simplekml module must be installed to use drawkml (pip install simplekml)'

        kml = simplekml.Kml(name='Plume_{id}'.format(id=self.identifier))

        filepath = filepath if filepath else '{dir}/Plume_{id}.kml'.format(dir=os.getcwd(), id=self.identifier)
        kml = simplekml.Kml()

        for arg in args:
            if 'misr_' in arg.lower():
                import MISRutils as mu
                # draw camera image

                # callthe drawmisr fnuction and save the png then project the image on the kml

            elif 'caliop' in arg.lower():
                # draw caliop path

                cal_line = kml.newlinestring(
                    name='Calipso path',
                    description='Approximate orbital path of the Calipso satellite (+- 2km), projected on the earth',
                    coords=list(zip(self.Caliop.longitude, self.Caliop.latitude))
                )

            elif 'firepixel' in arg.lower():
                import netCDF4
                # draw fire pixels

                # grab granule dir from ?(text file)? and load corresponding granule to extract fire pixels

            elif isinstance(arg, dict):
                # draw the requested plotting method
                pass
                # call method with save arg, and then project the png on the kml

            else:
                mssg = "drawkml got an unexpected request '{arg}', see documentation for available options".format(arg=arg)
                raise TypeError(mssg)

        kml.save(filepath)

    # --------------------------------------------------------------------------
    # Comparison methods
    # --------------------------------------------------------------------------
    # TODO model comparisons should be default in these mehtods, not an option, ensure that
    # need to make these consistent, but still figure out a way to be able to do with other value and products than just height

    def table(self, option='filtered_height', model=True, output='json', full=False, custom=[('', []), ], verbose=True, save='', filename=None):
        '''*Function for printing a formatted table of values and returning said data in a desired output format*

        Prints a formatted table of values, for each lat/long point, comparing the plume attributes to the model directly.
        Can display multiple attributes and table can be saved in multiple different file formats, if requested.
        Data will be returned from function in the format specified in the `output` argument.
        Custom fields can be provided, by giving a list of tuples containing the name of the field and the data. Custom data
        should be of the same length as all other data. The custom names should be larger than four characters, for best results.
        As well, custom data must be numerical, lists of string will not function.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to be displayed. For a list of valid attributes see Plume class docstring
        model : bool, optional
            If specified, model plume heights will be displayed on left of table.
        output : str, optional
            Specifies the format the table data will be output as. [['json', 'dataframe', ]]
        full : bool, optional
            Specifies wether the entirety of the table will be printed to stdout, or just the first and last few lines.
        custom : tuple of str and iterable
            Custom dataset to be formatted and displayed in table. Must be formatted as (name, [data]).
        verbose : bool, optional
            Specifies wether a table will be printed to stdout.
        save : str, optional
            If desired, specifies a format to save the table data as. [['json', 'csv', 'text']]
        filename : str, optional
            Name or path to which you wish to save the JSON file. `save` arg must be specified in order to use `filename`.

        Returns
        -------
        tablejson : list of dict
            Flat JSON object (list of dictionaries representing each row) containing all requested attributes as well as their latitudes and longitudes.

        Raises
        ------
        ValueError: If custom fields are entered incorrectly.

        Examples  TODO
        --------
        >>> plm = Plume.from_filename('./tests/Plumes_O051339-B062-SPWB01.txt')
        >>> plm.table()

        '''
        # TODO include units, some options formatting still broken (terrain)
        # TODO maybe change from E to just spacing
        # TODO perhaps more output options

        # --------------------------------------------------------------------------
        # Argument Checking
        # --------------------------------------------------------------------------

        if model and not hasattr(self, 'Model'):
            mssg = "This plume has no corresponding model data. Gather said data by using the 'get_model' method"
            raise AttributeError(mssg)
        if custom[0][1] and not all([len(cus[1]) == len(self.points) for cus in custom]):
            mssg = 'custom fields must have same number of points as plume ({0})'.format(len(self.points))
            raise ValueError(mssg)

        option = arg_to_list(option)

        # --------------------------------------------------------------------------
        # Create header and break line.
        # Header has four spaces between each option name and line extends past each
        # option name by one dash on either side, with two spaces seperating each
        # --------------------------------------------------------------------------

        header = '\n {lat}    {lon}    {mod}    '.format(
            lat='latitude',
            lon='longitude',
            mod='model height' if model else '',
        )

        header += '    '.join(option + [cus[0] for cus in custom if cus[0]])

        brk = '{lat}  {lon}  {mod}  '.format(
            lat='-' * (len('latitude') + 2),
            lon='-' * (len('longitude') + 2),
            mod='-' * (len('model height') + 2),
        )

        brk += '  '.join(['-' * (len(x) + 2) for x in (option + [cus[0] for cus in custom if cus[0]])])

        # --------------------------------------------------------------------------
        # Loop through each point in the plume and create a line, saving the data
        # for printing and for the output
        # --------------------------------------------------------------------------

        lines = []
        tablejson = []
        for ind, point in enumerate(self.points):

            # --------------------------------------------------------------------------
            # Check if the table should be printed at this line. If `full` is not
            # specified, table should only be printed at beginning and ending rows.
            # --------------------------------------------------------------------------

            if (((len(self.points) - 6 < ind) or (ind < 5)) or full) and verbose:

                # --------------------------------------------------------------------------
                # Print all data, one column at a time, for this row, with two spaces in
                # between each section, as with the brk's.
                # The formatting of the data goes as follows:
                #   ^ -> center the data in the area
                #   {area} -> set the amount of space available to be the same as the brk's
                #   {size} -> the amount of significant figures the data will have
                # All fill values are replaced with `--------`
                # --------------------------------------------------------------------------

                # --------------------------------------------------------------------------
                # Latitude and Longitude, always present
                # --------------------------------------------------------------------------

                line = '{lat:^10.6}  {lon:^11.6}  '.format(lat=point['latitude'], lon=point['longitude'])

                # --------------------------------------------------------------------------
                # Model height
                # --------------------------------------------------------------------------

                if model:
                    line += '{mod:^14.8}  '.format(mod=self.Model.height[ind] if not math.isinf(self.Model.height[ind]) else '-' * 8)

                # --------------------------------------------------------------------------
                # All plume options
                # --------------------------------------------------------------------------

                line += '  '.join(["{data:^{area}.{size}}".format(data=point[op] if not math.isnan(point[op]) else '-' * 8,
                                                                  area=len(op) + 2, size=len(op) - 2) for op in option])

                # --------------------------------------------------------------------------
                # All custom options
                # --------------------------------------------------------------------------

                try:
                    if custom[0][1]:
                        line += '  '

                        line += '  '.join(["{data:^{area}.{size}}".format(data=cuval[ind] if not math.isnan(cuval[ind]) else '-' * 8,
                                                                          area=len(cuname) + 2, size=len(cuname) - 2) for cuname, cuval in custom])

                except (ValueError, TypeError):
                    mssg = 'invalid custom name or data. Please ensure all data is in float format and the custom names are valid strings'
                    raise ValueError(mssg)

                lines.append(line)

            elif ind == 5 and not full and verbose:
                lines.append('\n\t:\n\t:\n')

            # --------------------------------------------------------------------------
            # Save all data points for this row to a dictionary for the output
            # --------------------------------------------------------------------------

            tablejson.append({op: point[op] for op in (['latitude', 'longitude'] + option)})
            if custom[0][1]:
                for cuname, cuval in custom:
                    tablejson[ind][cuname] = cuval[ind]
            if model:
                tablejson[ind]['model'] = self.Model.height[ind]

        # --------------------------------------------------------------------------
        # If requested, print the table to the console
        # --------------------------------------------------------------------------

        if verbose:
            sys.stdout.write(header + '\n')
            sys.stdout.write(brk + '\n')
            for line in lines:
                sys.stdout.write(line + '\n')

        # --------------------------------------------------------------------------
        # If requested, save the output as the desired filetype
        # --------------------------------------------------------------------------

        if save:
            if not filename:
                filename = raw_input('Enter desired JSON filename: ')

            if save.lower() in ('json', 'dict'):
                import json
                with open(filename, 'w') as savefile:
                    json.dump(tablejson, savefile)

            elif save.lower() == 'csv':
                import csv
                with open(filename, 'w') as savefile:
                    dw = csv.DictWriter(savefile, tablejson[0].keys())
                    dw.writeheader()
                    dw.writerows(tablejson)

            elif save.lower() in ('text', 'txt', 'ascii'):
                with open(filename, 'w') as savefile:
                    savefile.write(header + '\n')
                    savefile.write(brk + '\n')
                    for line in lines:
                        savefile.write(line + '\n')

            else:
                mssg = 'incorrect save option "{0}", see documentation for options'.format(save)
                raise ValueError(mssg)

        # --------------------------------------------------------------------------
        # Return the output in the desired format
        # --------------------------------------------------------------------------

        if output.lower() == 'json':
            return tablejson
        elif output.lower() == 'dataframe':
            import pandas as pd
            return pd.DataFrame(tablejson)
        elif output is None:
            return
        else:
            mssg = 'incorrect output option "{0}", see documentation for options'.format(output)
            raise ValueError(mssg)

    def abs_diff(self, option='filtered_height', plot=True, axis='distance', style='-', save=False, filename=None, *args, **kwargs):
        '''*Calculate the absolute differences between given options and the model plume height*

        Calculate the absolute difference between given options and the model plume height, at every point.
        While any option is accepted, this method only makes sense for the height attributes of the plume.
        If either the option or model has an invalid value at a point, the abs. difference will be given as NaN. If desired,
        and by default, the absolute differences will be plotted, on a given axis, with a simple 2D plot.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to calculate absolute difference from plume height on. For a list of valid attributes see Plume class docstring.
        plot : bool, optional
            if specified, absolute difference values will be plotted on a given axis. All the following args depend on this flag.
        axis : {'distance', 'longitude', 'latitude'}, optional
            Specifies which attribute will be set as the x-axis the heights will be plotted versus.
        style : str
            Style of the absolute difference plot. See matplotlib `marker` documentation for list of valid options.
        save : bool, optional
            If specified, will save plot rather than simply displaying it.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib plot function

        Returns
        -------
        abs_list : list or list of list
            Returns list(s) containing all calulated absolute differences.
        (Will display or save plot if specified)

        '''
        if not hasattr(self, 'Model'):
            mssg = "This plume has no corresponding model data. Gather said data by using the 'get_model' method"
            raise AttributeError(mssg)

        option = arg_to_list(option)

        # dont use clean inv cause need same len as opt
        abs_list = [[abs(x - self.Model.height[i]) if not math.isinf(self.Model.height[i]) else float('nan')
                    for i, x in enumerate(self[opt])] for opt in option]

        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.autoscale()

            for ind, opt in enumerate(abs_list):
                opt_axis, opt_vals = zip(*sorted(zip(self[axis], opt)))
                plt.plot(opt_axis, opt_vals, style, color=next(self.clr), label=option[ind])

            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
            ax.legend(numpoints=1, loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=4)

            if save:
                if not filename:
                    filename = raw_input('Enter desired filename: ')
                fig.savefig(filename)
            else:
                plt.show()

            fig.clf()
            plt.close(fig)

        return abs_list[0] if len(abs_list) == 1 else abs_list

    def RMSD(self, option='filtered_height'):
        '''*Calculate the root-mean-squared deviation between given options and the model plume height*

        Calculate the root-mean-squared deviation between given options and the model plume height, at every point.

        RMSD formula -> sqrt(sum((ph - mh)^2) / N)

        While any option is accepted, this method only makes sense for the height attributes of the plume.
        Be aware that the returned lists do NOT match all other attribute shapes.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to calculate RMSD from plume height on. For a list of valid attributes see Plume class docstring.

        Returns
        -------
        RMS : list or list of list
            Returns list(s) containing all calculated RMSD.

        '''
        # TODO make this have same shape as attributes again
        if not hasattr(self, 'Model'):
            mssg = "This plume has no corresponding model data. Gather said data by using the 'get_model' method"
            raise AttributeError(mssg)

        option = arg_to_gen(option)

        RMS = []
        for opt in option:
            optval, modelval = clean_inv(self[opt], self.Model.height)

            numer = sum([(oval - modelval[i])**2 for i, oval in enumerate(optval)])
            RMS.append((numer / len(optval))**0.5)

        return RMS[0] if len(RMS) == 1 else RMS

    def corr_coeff(self, option=['filtered_height', 'model_height']):
        '''*Calculate a Pearson correlation coefficient between given options and the model plume height*

        Calculates a Pearson correlation coefficient and the corresponding p-value for testing of non-correlation, between given options and
        the model plume height. Utilizes the Scipy.stats pearsonr function. See Scipy documentation for a more in-depth explanation of
        their methodology.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to calculate RMSD from plume height on. For a list of valid attributes see Plume class docstring.

        Returns
        -------
        out : Tuple
            Returns Tuple containing the correlation coefficient, and the p-value for testing of non-correlation

        '''
        # TODO multiple options at once?, do for all possible model attrs (values and products not just height)
        # TODO think about taking model out of opt arg again, for consistency (although how will they specify if they want to compare two other opts then?)
        from scipy import stats

        option = arg_to_list(option)

        if len(option) == 1 and 'model' not in option[0]:
            option.append('model_height')

        if not len(option) == 2:
            raise NotImplementedError('must specify a list of either one or two options at this time.')

        arg1, arg2 = clean_inv(*[self[opt] if 'model' not in opt else self.Model.height for opt in option])

        return stats.pearsonr(arg1, arg2)

    def linear_comp(self, option='filtered_height', plot=True, save=False, filename=None, *args, **kwargs):
        '''*Plot model plume height versus given options, to determine correlation*

        Plot model plume height versus given options, to investigate linearity, using a linear least-squares regression,
        and correlation. The Scipy.stats linregress function is used to determine the slope & intercept of a regression line,
        the correlation coefficient, the p-value of a zero-slope null hypothesis, and the standard error.
        While any option is accepted, this method only makes sense for the height attributes of the plume.

        Args
        ----
        option : str or list of str
            Plume attribute(s) to plot and compare with model plume heights. For a list of valid attributes see Plume class docstring.
        plot : bool, optional
            if specified, the comparison plot will be displayed, along with the regression line. All the following args depend on this flag.
        save : bool, optional
            If specified, will save plot rather than simply displaying it.
        filename : str, optional
            Name or path to which you wish to save the plot. `save` arg must be specified in order to use `filename`.
        args, kwargs:
            All extra parameters passed to the matplotlib plot function

        Returns
        -------
        abs_list : list or list of list
            Returns list(s) containing all calulated absolute differences.
        (Will display or save plot if specified)

        '''
        import scipy.stats

        option = arg_to_list(option)

        if plot:
            fig = plt.figure()
            plt.suptitle('Linear Correlation Comparison')

        results = []
        for ind, opt in enumerate(option, 1):

            optval, modelval = clean_inv(self[opt], self.Model.height)

            m, b, r, p, stderr = scipy.stats.linregress(optval, modelval)
            results.append((m, b, r, p, stderr))  # better way?

            if plot:
                ax = fig.add_subplot(1, len(option), ind)

                if ind == 1:
                    ax.set_ylabel('Model Heights')
                ax.set_xlabel(opt)

                ax.plot(optval, modelval, 'o', color=next(self.clr), label=opt, *args, **kwargs)
                ax.plot(optval, [m * x + b for x in optval], '-', color=next(self.clr), label=opt)
                ax.text(min(optval), m * min(optval) + b, s='{m:.3f}x+{b:.1f}'.format(m=m, b=b), verticalalignment='top' if m > 0 else 'bottom')

        if plot and save:
            if not filename:
                filename = raw_input('Enter desired filename: ')
            fig.savefig(filename)

            fig.clf()
            plt.close(fig)

        elif plot:
            plt.show()

            fig.clf()
            plt.close(fig)

        return results[0] if len(results) == 1 else results

    def mean_bias(self, option='filtered_height', error=False, normalized=False):
        '''*Calculate the mean bias, or mean bias error (just adds absolutes inside sum)*
        ...
        '''
        # TODO make this have same shape as attributes again
        option = arg_to_gen(option)

        bias = []
        for opt in option:
            optval, modelval = clean_inv(self[opt], self.Model.height)

            if not error:
                numer = sum([(oval - modelval[i])**2 for i, oval in enumerate(optval)])
            else:
                numer = sum([abs(oval - modelval[i])**2 for i, oval in enumerate(optval)])

            denom = len(self[optval]) if not normalized else sum(modelval)

            bias.append(numer / denom)

        return bias[0] if len(bias) == 1 else bias


class _ModelPlume(object):
    """*ModelPlume object for the visualization and manipulation of GEMMACH model output and comparison with MINX results*

    Model plume object. Utilizes the model output of the GEMMACH FireWork or CFFEPS software, given in an output `fst` file.
    You must ensure that said file completely covers the plume and has a filename reflecting the date.

    *This class should only be created and accessed through a corresponding Plume object.
    All attributes and methods of a model plume should be accessed as follows: `Plume().Model.{attr/method}`.*

    The plume heights are calculated using local maximums of PM2.5 concentration, and the height and concentration values
    at this determined plume top are saved as `height` and `value` respectively. All arrays will correspond exactly with their
    Plume object counterparts, size and sorting wise.

    All other model products can also be loaded, as profiles at each MINX coordinate.


    Note
    ----
        Greater explanation of the plume height retrieval algorithm, as well as greater explanation of this class and its values
        is given on the wiki: lolnotyetsry.wiki

        Model fill values are -inf.

    **Initialization**
        *Plume().get_model(model_dir)*

        *Plume.from_coord((long, lat), model_dir=model_dir)*

        *Plume.from_id(orb_num, blk_num, plm_num, model_dir=model_dir)*

        *Plume.from_filename(filepath, model_dir=model_dir)*

    **Attributes**
        Plume top: *height, value*
        Extra Products: *product['height'], product['value'] (if product retrieved)*

    **Methods**
        Extra Products: *get_product*

    """

    def __getitem__(self, key):
        return self.__dict__[key]

    def __init__(self, plm, fst_dir, var='AF', threshold=4.0):
        super(_ModelPlume, self).__init__()

        import pyloader
        import rpnpy.librmn.all as rmn

        # --------------------------------------------------------------------------
        # Function for returning all relevant data given a key
        # --------------------------------------------------------------------------

        def data_from_key(key):
            '''given a key, return the metadata, grid info, xypos info and values at each lat/long of this plume, for this file'''

            # --------------------------------------------------------------------------
            # Determine the metadata corresponding to this key and add the file id
            # --------------------------------------------------------------------------

            meta = rmn.fstprm(key)
            meta['iunit'] = fid

            # --------------------------------------------------------------------------
            # Grab the actual data from this layer, as well as the relevant grid data
            # --------------------------------------------------------------------------

            data = rmn.fstluk(key)['d']
            grid = rmn.ezqkdef(meta)

            # --------------------------------------------------------------------------
            # Interpolate to, and calculate the values at, each MINX lat/lon point
            # --------------------------------------------------------------------------

            xypos = rmn.gdxyfll(grid, latarr, lonarr)
            val = rmn.gdxysval(grid, xypos['x'], xypos['y'], data)

            return meta, grid, xypos, val

        # --------------------------------------------------------------------------
        # Ensure that the ModelPlume is being instantiated from a MINX plume
        # --------------------------------------------------------------------------

        if not isinstance(plm, Plume):
            raise TypeError('ModelPlume object must be instantiated from a Plume object as Plume().get_model()')

        # --------------------------------------------------------------------------
        # Get time used in fst filename by rounding MINX plume time to the nearest
        # hour and then subtracting 12 hours, as fst file times are based off of
        # the time elapsed since 12 o'clock
        # --------------------------------------------------------------------------

        file_time = (plm.datetime if plm.datetime.minute < 30 else plm.datetime + datetime.timedelta(hours=1)) - datetime.timedelta(hours=12)
        file_path = "{0}/{1}".format(os.path.abspath(fst_dir), file_time.strftime('%Y%m%d12_0%H'))

        # --------------------------------------------------------------------------
        # Initialize interpolation options and open the fst file
        # --------------------------------------------------------------------------

        rmn.fstopt(rmn.FSTOP_MSGLVL, rmn.FSTOPI_MSG_CATAST)
        rmn.ezsetopt(rmn.EZ_OPT_INTERP_DEGREE, rmn.EZ_INTERP_NEAREST)

        fid = rmn.fstopenall(file_path, rmn.FST_RO)

        # --------------------------------------------------------------------------
        # Gather the list of points to check from the MINX plume and the list of
        # keys to check from this product entry in the fst file. each key
        # corresponds to a single layer of model output
        # --------------------------------------------------------------------------

        lonarr = plm.longitude
        latarr = plm.latitude
        keylist = rmn.fstinl(fid, nomvar=var)

        # --------------------------------------------------------------------------
        # Store some metadata
        # --------------------------------------------------------------------------

        self.Plume = plm
        self.Plume.units['model_height'] = 'm'
        self.Plume.units['model_value'] = 'ug/m^3'  # TODO fix unicode error, need dynamic units (as in get_product)

        ref_lvl = self.Plume.reference

        self.name = file_time.strftime('%Y%m%d12_0%H')
        self.dir = os.path.abspath(fst_dir)

        # --------------------------------------------------------------------------
        # Initialize the containers for the heights and values. Model fill values
        # are -inf and the containers are the same size as the lat/lon arrays
        # --------------------------------------------------------------------------

        self.height = [float('-inf')] * len(lonarr)
        self.value = [float('-inf')] * len(lonarr)
        self.points = [{} for _ in lonarr]

        # --------------------------------------------------------------------------
        # Decode the ip values for each layer and sort the keylist based on
        # ascending height. ip1 values are used to determine height, and they are
        # inversely porportional to it.
        # --------------------------------------------------------------------------

        iplist = []
        for key in keylist:
            meta = rmn.fstprm(key)
            iplist.append(rmn.DecodeIp(meta['ip1'], meta['ip2'], meta['ip3'])[0].v1)

        sorted_keylist = (x for _, x in sorted(zip(iplist, keylist), reverse=True))

        # --------------------------------------------------------------------------
        # Calculate the necessary data for the current layer to be checked
        # --------------------------------------------------------------------------
        for keyInd, keyRow in enumerate(sorted_keylist):
            cur_meta, cur_grid, cur_xypos, cur_val = data_from_key(keyRow)

            for ind, val in enumerate(cur_val):

                # --------------------------------------------------------------------------
                # Save the height as that of this layer if the var concentration is an
                # absolute maximum. Check the current concentration versus the
                # the concentration of both the preceeding and proceeding layers, and
                # check that it is both above the minimum threshold and greater than any
                # previously stored value.
                # As the keys are sorted by ascending height, this will return the height
                # at each point to be the highest (in altitude) absolute maximum of
                # PM2.5 concentration
                # --------------------------------------------------------------------------

                if (val >= self.value[ind]):

                    # --------------------------------------------------------------------------
                    # Calculate the actual height of the point in this accepted layer. Heights
                    # are stored within the GZ product.
                    # It is slightly wasteful to regenerate the heights for each point in the
                    # layer every time, but not enough to merit redoing this
                    # --------------------------------------------------------------------------

                    try:
                        gzkey = rmn.fstinf(fid, nomvar='GZ', ip1=cur_meta['ip1'])['key']
                        gzdata = rmn.fstluk(gzkey)['d']

                        height = rmn.gdxysval(cur_grid, cur_xypos['x'], cur_xypos['y'], gzdata) * 10
                        self.height[ind] = height[ind].item() if ref_lvl == 'sea' else (height[ind] - self.Plume.terrain[ind]).item()
                        self.value[ind] = val.item()

                        self.points[ind]['height'] = height[ind].item() if ref_lvl == 'sea' else (height[ind] - self.Plume.terrain[ind]).item()
                        self.points[ind]['value'] = val.item()

                    except TypeError:
                        print('GZ missing for ip1: {}'.format(cur_meta['ip1']))
                        # A GZ value is not present on certain layers, skipping this layer should be fine
                        continue


        rmn.fstcloseall(fid)
