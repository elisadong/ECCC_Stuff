#! /usr/bin/env python

# script to run most postMINX features viable
# returns plume details and associated model details
# uses the postMINX program
# NOTE: function expects that there is a corresponding model file, and possibly a second one

import os
import postMINX as pm
import pandas as pd
import matplotlib.pyplot as plt
import errno

# example of use
# import speedyAltMINXdemo as sam
# sam.pltPlume(path = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/plumesP50C60', plume ='plume.txt', mod_fire='/fs/site2/dev/eccc/aq/r1/eld001/MINX/fireworkLowRes/mach', mod_nofire='/fs/site2/dev/eccc/aq/r1/eld001/MINX/fireworkLowRes_nofire', threshVal = 0, startTime = 12)

imageMethods = ['contour','colour', 'scatter_3D', 'surface_plot', 'direction_plot', 'scatter_2D', 'poly_2D', 'shape', 'hist']
numMethods = ['max', 'min','median','mean','std']
compMethods = ['table', 'abs_diff', 'linear_comp']
compNumMethods = ['RMSD', 'corr_coeff', 'mean_bias']
currentDir = os.getcwd()
picklePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/FW-GM_start00_pklPlmsOps/"
#picklePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/CF-GM_start00_pklPlmsOps/"
savePlumePlts = False

try:
    import cPickle as pickle
except ImportError:
    import pickle

# TODO: add option to use get products, then add to list of methods after when plotting

def makeDir(path):
    ''' makeDir(path) makes the path a directory unless it already exists. '''

    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def pltMethods(path, plume, plm, modelStats=False, saveMethod =True, modType=''):
    saveDir = os.path.join(currentDir, 'postAltOutputs',plume[7:-4] + modType)
    makeDir(saveDir)
    os.chdir(saveDir)
    for im in imageMethods:
        saveName = plume[7:-4] + im + '.png'
        try:
            eval('plm.{}(save={},filename="{}")'.format(im, saveMethod, saveName))
            plt.close('all')
        except:
            plt.close('all')
            continue
        print('Save {} success!'.format(saveName))

    # number methods for minx plume only
    numSave = plume[7:-4] + '.txt'
    option = 'filtered_height'
    numList = []
    for nu in numMethods:
        numVal = eval('plm.{}(option="{}")'.format(nu, option))
        numRow = [nu, numVal]
        numList += [numRow]

    numDF = pd.DataFrame(numList, columns=('method', option))
    numDF.to_csv(numSave, sep = ',', index=False)

    if modelStats:
        options = ['filtered_height', 'model_height']
        for im in imageMethods:
            try:
                saveName = 'Mod_{}{}.png'.format(plume[7:-4], im)
                eval('plm.{}(save={}, option={}, filename="{}")'.format(im,saveMethod, options, saveName))
                plt.close('all')
                print('Save {} success!'.format(saveName))
            except:
                plt.close('all')
                continue
        numSave='Mod_{}.txt'.format(plume[7:-4])
        numList = []
        option1 = 'filtered_height'
        option2 = 'model_height' # note that this won't actually get the model
        opt2 = plm.Model.height

    os.chdir(currentDir)

def pltPlume(path, plume, minFRP, mod_fire, mod_nofire=None, threshVal=4, startTime = 12, toSave=savePlumePlts, modelStats=True, modType = ''):
    '''
    pltPlume() returns a plume object with FRP, filtered_height, model height, model values(AF at height), and distance from plume origin.

    Additional plotting methods and stats can be turned on by setting toSave=True and modelStats=True.
    See postMINX.py documentation for other parameters.
    '''
    makeDir(picklePath)
    # try to unpickle first

    try:
        plm=pickle.load(open(os.path.join(picklePath, plume[0:5] + plume[6:-3]+'pkl'), 'rb'))
        if plm['frp'] >=minFRP:
            # Some model eval stuff, move back in when want to use, for now should be skipped
            if toSave  == True:
                try:
                    pltMethods(path, plume, plm, modelStats=True, modType = '')
                except:
                    pass
            return plm
        else:
            #print('{} FRP too low, plume will not be used.'.format(plume))
            return None
    except:
        plm = pm.Plume.from_filename(os.path.join(path,plume))
        if plm['frp'] >= minFRP:
            try:
                plm.get_model(fst_dir=mod_fire, ctrl_dir=mod_nofire, threshold=threshVal, filestart = startTime)
                plm.save_plume(os.path.join(picklePath, 'Plume_{}.pkl'.format(plume[7:26])))
                print('{} has been pickled.'.format(plume))
                if toSave == True:
                    try:
                        pltMethods(path, plume, plm, mod_fire, minFRP=0, mod_nofire=None, threshVal=4, startTime = 12, modelStats=True, modType = '')
                    except:
                        pass
                try:
                    return plm
                except:
                    #print('Plume and Model height were not returned.')
                    return None
            except:
                return None
        else:
            #print('{} FRP too low, plume will not be used.'.format(plume))
            return None

def plotFuelMap(fuelDir, plmOrigins=None, pltFuels = False, pltPoints = True, modifier = ''):
    import gdal
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    a = 6378137
    f = 1/298.257222101
    b = a*(1-f)
    w = 5351000
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
        m.pcolormesh(xArray,yArray,fuelArray, vmin = 100, vmax = 122, cmap=plt.cm.get_cmap('jet',21), alpha = 0.7)
    #m.drawrivers()
    m.drawcountries(linewidth =0.5, color = 'gray')
    m.drawcoastlines(linewidth = 0.5, color = 'gray')
    m.drawstates(linewidth = 0.5, color = 'gray')
    #m.drawmeridians(np.arange(-170, -20,10), linewidth = 0.2, labels=[0,0,0,1], fontsize=3)
    #m.drawparallels(np.arange(40,60,10), linewidth = 0.2, labels=[1,0,0,0],fontsize=3)
    #plt.colorbar(shrink=0.5)
    #assume passed in the xy origins

    fig.savefig(os.path.join(fuelDir, 'fuel{}_distribution.png'.format(modifier)), dpi = 300)
    plt.close('all')
