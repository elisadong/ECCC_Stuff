# script to run every single postMINX feature viable
# gets plume from filename

import sys
import os
import postMINX as pm
import pandas as pd
import matplotlib.pyplot as plt
import errno

plume = sys.argv[1]
# default model path /fs/site1/dev/eccc/aq/r5/jac001/PlumeRiseEval/ModelInput/HighRes2p5km/highResDupes
try:
    model = sys.argv[2]
    plm = pm.Plume.from_filename(plume, model_dir=model)
except:
    model = False
    plm = pm.Plume.from_filename(plume)

saveMethod = True
currentDir = os.getcwd()
saveDir = os.path.join(currentDir, 'postOutputs',plume[7:-4] + '_HighRes')

def makeDir(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


imageMethods = ['contour','colour', 'scatter_3D', 'surface_plot', 'direction_plot', 'scatter_2D', 'poly_2D', 'shape', 'hist']
numMethods = ['max', 'min','median','mean','std']
compMethods = ['table', 'abs_diff', 'linear_comp']
compNumMethods = ['RMSD', 'corr_coeff', 'mean_bias']

# add option to use get products, then add to list of methods after when plotting

makeDir(saveDir)
os.chdir(saveDir)

# image methods for minx plume only
for im in imageMethods:
    saveName = plume[7:-4] + im + '.png'
    eval('plm.{}(save={},filename="{}")'.format(im, saveMethod, saveName))
    plt.close()
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

if model:
    options = ['filtered_height', 'model_height']
    for im in imageMethods:
        try:
            saveName = 'Mod_{}{}.png'.format(plume[7:-4], im)
            eval('plm.{}(save={}, option={}, filename="{}")'.format(im,saveMethod, options, saveName))
            plt.close()
            print('Save {} success!'.format(saveName))
        except:
            continue
    numSave='Mod_{}.txt'.format(plume[7:-4])
    numList = []
    option1 = 'filtered_height'
    option2 = 'model_height' # note that this won't actually get the model
    opt2 = plm.Model.height
    # she's dead, Jim! Won't work in nice loops, this is due to be fixed. Best write own stats methods anyways
    #for nu in numMethods:
    #    numVal1 = eval('plm.{}(option="{}")'.format(nu, option1))
    #    numVal2 = eval('opt2.{}'.format(nu))
    #    numRow = [nu, numVal1, numVal2]
    #    numList += [numRow]
    #numDF = pd.DataFrame(numList, columns=('method', option1, option2))
    #numDF.to_csv(numSave,sep=',',index=False)

os.chdir(currentDir)
