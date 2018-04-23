## determine the mapType details for plotting

import matplotlib.pyplot as plt

cmaps = ['viridis', 'plasma', 'inferno', 'magma',
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
            'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c',
            'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg', 'hsv',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']

defaultCmap = 'viridis'
defaultBins = None
defaultExt = 'neither'

def reverseName(mapType, reverse):
    if reverse  == True:
        mapType = mapType + '_r'
        return mapType
    else:
        return mapType

def isCmap(mapType):
    if mapType in cmaps:
        return True
    else:
        print('Error: The mapType is not in the list of cmaps, please try again.')

def cmapType(plt,mapType=defaultCmap, reverse=False,totalBins=defaultBins):
    #note that extend should be one of 'neither','both','min' or 'max'
    # determine the color scheme of the map
    mapType = reverseName(mapType,reverse)
    # check if map is valid, can continue if is
    isCmap(mapType)
    if totaBins != None:
        return cmaps = plt.cm.get_cmap(mapType, totalBins)
    else:
        return cmaps = plt.cm.get_cmap(mapType)
    # if only one of min/max value are provided, clim should automatically compensate for the scaling

def colorbarInfo(plt,minVal=None,maxVal=None, extension=defaultExt):
    # returns details on the colorbar
    if (minVal==None) and (maxVal==None):
        plt.colorbar(extend=extension)
    else:
        plt.colorbar(extend=extension)
        plt.clim(minVal,maxVal)
