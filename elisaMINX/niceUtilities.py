# functions that are used for plume analysis, but not directly related
import os
import errno
import math

def makeDir(path):
    ''' makeDir(path) makes the path a directory unless it already exists. '''

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def toFloats(row):
    floats = map(float, row.split())
    return floats

def readFile(path):
    f = open(path, 'r')
    lines = f.readlines()
    f.close()
    return lines

def clean_inv(*listN, **kwargs):
    if not kwargs:
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

def plotFuelMap(fuelDir, plmOrigins=None, pltFuels = False, pltPoints = True, modifier = ''):
    ''' Plot the canfuels map. '''

    import gdal
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    a = 6378137
    f = 1/298.257222101
    b = a*(1-f)
    w = 5351000
    h = 4558250
    cellWidth = 250

    fuelData = gdal.Open(fuelFile, gdal.GA_ReadOnly)
    fuelArray = fuelData.ReadAsArray()
    yList = []
    count = 0
    for y in range(len(fuelArray)):
        yRow = [count]*len(fuelArray[1])
        yList += [yRow]
        count += cellWidth
    yArray = np.array(yList)
    yArray = np.flipud(yArray)
    xList = []
    for x in range(len(fuelArray)):
        xRow = range(0,w,250)
        xList += [xRow]
    xArray = np.array(xList)

    fig = plt.figure(figsize = (6,6))

    m = Basemap(projection='lcc', rsphere = (a,b), lat_1=49, lat_2=77, lat_0=49, lon_0=-95, llcrnrlon = -121.112, llcrnrlat = 38.35766, urcrnrlon =-12.02121, urcrnrlat = 63.05306, resolution = 'i')
    if pltPoints:
        xoList = []
        yoList = []
        for o in plmOrigins:
            xoList += [o[0]]
            yoList += [o[1]]
        m.scatter(xoList, yoList, 1, marker = '^',color = 'k') #,zorder = 10)

    if pltFuels:
        m.pcolormesh(xArray,yArray,fuelArray, vmin = 100, vmax = 123, cmap=plt.cm.get_cmap('RdBu',22), alpha = 0.7)
    #m.drawrivers()
    m.drawcountries(linewidth =0.5, color = 'gray')
    m.drawcoastlines(linewidth = 0.5, color = 'gray')
    m.drawstates(linewidth = 0.5, color = 'gray')
    #m.drawmeridians(np.arange(-170, -20,10), linewidth = 0.2, labels=[0,0,0,1], fontsize=3)
    #m.drawparallels(np.arange(40,60,10), linewidth = 0.2, labels=[1,0,0,0],fontsize=3)
    #plt.colorbar(shrink=0.5)
    #assume passed in the xy origins

    fig.savefig(os.path.join(saveDir, 'fuel{}_distribution.png'.format(modifier)), dpi = 300)
    plt.close('all')
