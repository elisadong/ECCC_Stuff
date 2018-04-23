# check what the plumes over 20k are

import os
import postMINX as pm
import matplotlib.pyplot as plt
import numpy as np 
try:
    import cPickle as pickle
except ImportError:
    import pickle


CFpklPath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/CF-GM_start00_pklPlmsOps/"
FWpklPath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/FW-GM_start00_pklPlms/"
plumeList = ['Plume_O093120-B037-SPWB07', 'Plume_O093120-B037-SPWB05', 'Plume_O093120-B037-SPWB03', 'Plume_O093120-B037-SPWB06', 'Plume_O093120-B037-SPWB04']
plmCFs = []
plmFWs = []

for p in plumeList:
    pklFile = p + '.pkl'
    CFplm = pickle.load(open(os.path.join(CFpklPath, pklFile), 'rb'))
    FWplm = pickle.load(open(os.path.join(FWpklPath, pklFile), 'rb'))
    plmCFs += [CFplm]
    plmFWs += [FWplm]

# plot scatter plots

for pInd, (cPlm, fPlm, plm) in enumerate(zip(plmCFs, plmFWs, plumeList)):
    cName = 'CF-' + plm
    fName = 'FW-' + plm
    pm.Plume.scatter_3D(cPlm, option=['filtered_height', 'model_height'], save=True, filename='scatter_' + cName)
    plt.close('all')
    pm.Plume.scatter_3D(fPlm, option=['filtered_height', 'model_height'], save=True, filename='scatter_' + fName)
    plt.close('all')
    pm.Plume.hist(cPlm, option=['filtered_height', 'model_height'], save=True, filename='hist_' + cName)
    plt.close('all')
    pm.Plume.hist(fPlm, option=['filtered_height', 'model_height'], save=True, filename='hist_' + fName)
    plt.close('all')

cModVals = []
>>> for c in plmCFs:
...     cModVals += c.Model.value

# get average value
