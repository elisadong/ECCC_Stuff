# plot and save a time series for your choice of attributes
# should only be used after you've pickled the appropriate plumes until I fix how the pickle path is passed around

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import os
import datetime
import numpy as np
from postMINX import clean_inv

try:
    import cPickle as pickle
except ImportError:
    import pickle

FWpklPath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/FW-GM_start00_pklPlmsOps/"
CFpklPath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/CF-GM_start00_pklPlmsOps/"
plumePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/plumesP50C60/plumesCan/"

# for looking at the origin points
# assume everything in folder is plume
plmList = os.listdir(plumePath)

# make dictionary storing values
plumes = {}
for p in plmList:
    if p not in plumes.keys():
        try:
            cwPlm = pickle.load(open(os.path.join(CFpklPath, p[0:5] + p[6:-3]+'pkl'), 'rb'))
            fwPlm = pickle.load(open(os.path.join(FWpklPath, p[0:5] + p[6:-3]+'pkl'), 'rb'))
            plumes[p] = {'CFplm':cwPlm, 'FWplm': fwPlm}
        except:
            continue

# for plotting origin
dates = {}
for p in plumes:
    datePlm = plumes[p]['FWplm']['datetime']
    if datePlm not in dates.keys():
        dates[datePlm] = {'cfHeights':[], 'fwHeights': [], 'minxHeights': []}
    cfOrg = plumes[p]['CFplm'].Model.height[0]
    fwOrg = plumes[p]['FWplm'].Model.height[0]
    minxOrg = plumes[p]['FWplm']['filtered_height'][0]
    dates[datePlm]['cfHeights'] += [cfOrg]
    dates[datePlm]['fwHeights'] += [fwOrg]
    dates[datePlm]['minxHeights'] += [minxOrg]
    #
    if ((type(cfOrg) == int or type(cfOrg) == float) and (type(fwOrg) == int or type(fwOrg) == float) and (type(minxOrg) == int or type(minxOrg) == float)):
        dates[datePlm]['cfHeights'] += [cfOrg]
        dates[datePlm]['fwHeights'] += [fwOrg]
        dates[datePlm]['minxHeights'] += [minxOrg]

datekeys = dates.keys()
datekeys.sort()
dateList = []
cfHeightList = []
fwHeightList = []
minxHeightList = []

# origins done
for dk in datekeys:
    #dateList += [mdates.date2num(dk)]
    cfHeight = dates[dk]['cfHeights']
    fwHeight = dates[dk]['fwHeights']
    minxHeight = dates[dk]['minxHeights']
    for hInd, h in enumerate(cfHeight):
        if h < 20000:
            dateList += [mdates.date2num(dk)]
            cfHeightList += [cfHeight[hInd]]
            fwHeightList += [fwHeight[hInd]]
            minxHeightList += [minxHeight[hInd]]

dList, cList, fList, mList = clean_inv(dateList, cfHeightList, fwHeightList, minxHeightList)
fig = plt.figure()
plt.plot_date(dList, cList, 'bo', label = 'CFFEPS')
plt.plot_date(dList, fList, 'ro', label = 'FireWork')
plt.plot_date(dList, mList, 'ko', label = 'MINX')
plt.legend()

# for plotting last plume height listed.
# NOTE: There is a possibility that the last plume height is NOT the furthest point from the origin. 
dates = {}
for p in plumes:
    datePlm = plumes[p]['FWplm']['datetime']
    if datePlm not in dates.keys():
        dates[datePlm] = {'cfHeights':[], 'fwHeights': [], 'minxHeights': []}
    cfOrg = plumes[p]['CFplm'].Model.height[-1]
    fwOrg = plumes[p]['FWplm'].Model.height[-1]
    minxOrg = plumes[p]['FWplm']['filtered_height'][-1]
    dates[datePlm]['cfHeights'] += [cfOrg]
    dates[datePlm]['fwHeights'] += [fwOrg]
    dates[datePlm]['minxHeights'] += [minxOrg]

datekeys = dates.keys()
datekeys.sort()
dateList = []
cfHeightList = []
fwHeightList = []
minxHeightList = []

# tail points
for dk in datekeys:
    #dateList += [mdates.date2num(dk)]
    cfHeight = dates[dk]['cfHeights']
    fwHeight = dates[dk]['fwHeights']
    minxHeight = dates[dk]['minxHeights']
    for hInd, h in enumerate(cfHeight):
        if h < 20000:
            dateList += [mdates.date2num(dk)]
            cfHeightList += [cfHeight[hInd]]
            fwHeightList += [fwHeight[hInd]]
            minxHeightList += [minxHeight[hInd]]

dList, cList, fList, mList = clean_inv(dateList, cfHeightList, fwHeightList, minxHeightList)
fig = plt.figure()
plt.plot_date(dList, cList, 'bo', label = 'CFFEPS')
plt.plot_date(dList, fList, 'ro', label = 'FireWork')
plt.plot_date(dList, mList, 'ko', label = 'MINX')
plt.legend()

# median height
dates = {}
for p in plumes:
    datePlm = plumes[p]['FWplm']['datetime']
    if datePlm not in dates.keys():
        dates[datePlm] = {'cfHeights':[], 'fwHeights': [], 'minxHeights': []}
        cfs = plumes[p]['CFplm'].Model.height
        fws = plumes[p]['FWplm'].Model.height
        minxs = plumes[p]['FWplm']['filtered_height']
        try:
            cfList, fwList, minxList = clean_inv(cfs, fws, minxs)
            cfMed =np.median(cfList)
            fwMed= np.median(fwList)
            minxMed =np.median(minxList)
            dates[datePlm]['cfHeights'] += [cfMed]
            dates[datePlm]['fwHeights'] += [fwMed]
            dates[datePlm]['minxHeights'] += [minxMed]
        except:
            continue

datekeys = dates.keys()
datekeys.sort()
dateList = []
cfHeightList = []
fwHeightList = []
minxHeightList = []
for dk in datekeys:
    cfHeight = dates[dk]['cfHeights']
    fwHeight = dates[dk]['fwHeights']
    minxHeight = dates[dk]['minxHeights']
    for hInd, h in enumerate(cfHeight):
        if h < 20000:
            dateList += [mdates.date2num(dk)]
            cfHeightList += [cfHeight[hInd]]
            fwHeightList += [fwHeight[hInd]]
            minxHeightList += [minxHeight[hInd]]


fig = plt.figure()
plt.plot_date(dateList, cfHeightList, 'bo-', label = 'CFFEPS')
plt.plot_date(dateList, fwHeightList, 'ro-', label = 'FireWork')
plt.plot_date(dateList, minxHeightList, 'ko-', label = 'MINX')
plt.legend()

# plot average height
dates = {}
for p in plumes:
    datePlm = plumes[p]['FWplm']['datetime']
    if datePlm not in dates.keys():
        dates[datePlm] = {'cfHeights':[], 'fwHeights': [], 'minxHeights': []}
        cfs = plumes[p]['CFplm'].Model.height
        fws = plumes[p]['FWplm'].Model.height
        minxs = plumes[p]['FWplm']['filtered_height']
        try:
            cfList, fwList, minxList = clean_inv(cfs, fws, minxs)
            cfMed =np.average(cfList)
            fwMed= np.average(fwList)
            minxMed =np.average(minxList)
            dates[datePlm]['cfHeights'] += [cfMed]
            dates[datePlm]['fwHeights'] += [fwMed]
            dates[datePlm]['minxHeights'] += [minxMed]
        except:
            continue

datekeys = dates.keys()
datekeys.sort()
dateList = []
cfHeightList = []
fwHeightList = []
minxHeightList = []
for dk in datekeys:
    cfHeight = dates[dk]['cfHeights']
    fwHeight = dates[dk]['fwHeights']
    minxHeight = dates[dk]['minxHeights']
    for hInd, h in enumerate(cfHeight):
        if h < 20000:
            dateList += [mdates.date2num(dk)]
            cfHeightList += [cfHeight[hInd]]
            fwHeightList += [fwHeight[hInd]]
            minxHeightList += [minxHeight[hInd]]

fig = plt.figure()
plt.plot_date(dateList, cfHeightList, 'bo', label = 'CFFEPS')
plt.plot_date(dateList, fwHeightList, 'ro', label = 'FireWork')
plt.plot_date(dateList, minxHeightList, 'ko', label = 'MINX')
plt.legend()





dateList = []
cwHeightList = []
fwHeightList = []
minxHeightList = []

datekeys = dates.keys()
datekeys.sort()
for d in datekeys:
    dateList += [mdates.date2num(d)]
    cwHeightList += [np.average(dates[d]['cwHeights'])]
    fwHeightList += [np.average(dates[d]['fwHeights'])]
    minxHeightList += [np.average(dates[d]['minxHeights'])]

fig = plt.figure()
plt.plot_date(dateList, cwHeightList, 'bo-', label = 'CFFEPS')
plt.plot_date(dateList, fwHeightList, 'ro-', label = 'FireWork')
plt.plot_date(dateList, minxHeightList, 'ko-', label = 'MINX')
plt.legend()

# CW and MINX comparison
datesCW = {}
for p in plumes:
    datePlm = plumes[p]['CFplm']['datetime']
    if datePlm not in datesCW.keys():
        datesCW[datePlm] = {'cfHeights':[], 'fwHeights':[], 'minxHeights': []}
    cfOrg = plumes[p]['CFplm'].Model.height[0]
    minxOrg = plumes[p]['CFplm']['filtered_height'][0]
    if ((type(cfOrg) == int or type(cfOrg) == float) and (type(minxOrg) == int or type(minxOrg) == float)):
        datesCW[datePlm]['cwHeights'] += [cwOrg]
        datesCW[datePlm]['minxHeights'] += [minxOrg]

dateList = []
cwHeightList = []
minxHeightList = []
datekeys = datesCW.keys()
datekeys.sort()
for d in datekeys:
    dateList += [mdates.date2num(d)]
    cwHeightList += [np.average(dates[d]['cwHeights'])]
    minxHeightList += [np.average(dates[d]['minxHeights'])]

fig = plt.figure()
plt.plot_date(dateList, cwHeightList, 'bo-', label = 'CFFEPS')
plt.plot_date(dateList, minxHeightList, 'ko-', label = 'MINX')
plt.legend()

# compare FW and MINX
datesFW = {}
for p in plumes:
    datePlm = plumes[p]['FWplm']['datetime']
    if datePlm not in datesFW.keys():
        datesFW[datePlm] = {'fwHeights': [], 'minxHeights': []}
    fwOrg = plumes[p]['FWplm'].Model.height[0]
    minxOrg = plumes[p]['FWplm']['filtered_height'][0]
    if ((type(fwOrg) == int or type(fwOrg) == float) and (type(minxOrg) == int or type(minxOrg) == float)):
        datesFW[datePlm]['fwHeights'] += [fwOrg]
        datesFW[datePlm]['minxHeights'] += [minxOrg]

dateList = []
fwHeightList = []
minxHeightList = []

datekeys = datesFW.keys()
datekeys.sort()
for d in datekeys:
    dateList += [mdates.date2num(d)]
    cwHeightList += [np.average(dates[d]['cwHeights'])]
    fwHeightList += [np.average(dates[d]['fwHeights'])]
    minxHeightList += [np.average(dates[d]['minxHeights'])]

fig = plt.figure()
plt.plot_date(dateList, fwHeightList, 'ro-', label = 'FireWork')
plt.plot_date(dateList, minxHeightList, 'ko-', label = 'MINX')
plt.legend()
