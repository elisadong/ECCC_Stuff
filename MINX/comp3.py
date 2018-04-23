#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2018-01-23 11:27:51
# Contact:    Nolan.Dickson@Canada.ca

import math
import glob
import numpy as np
import postMINX as pm
import matplotlib.pyplot as plt

txts = glob.glob('/fs/site2/dev/eccc/aq/r1/nod001/MINX/Working/AB_SK_*2017/Digitize_Output/*/Plumes_O*-B*-SPWB*.txt')

plumes = []
for ind, file in enumerate(txts, 1):
    print '{0}: ({1}/{2})'.format(file, ind, len(txts))
    p = pm.Plume.from_filename(file, model_dir='/fs/site2/dev/eccc/aq/r2/jac001/FireWork/EvalCFFEPS_OS/model_output/OS2p5_FireWork')

    if all([0 < h < 60000 or math.isinf(h) for h in p.Model.height]):
        p.Model.get_product(['GZ', 'AF', 'TNO2', 'TNO'])
        plumes.append(p)

max_dist = max([plm.max('distance')['distance'] for plm in plumes])
num_bins = int(round(max_dist)) / 5
binned_rat = [[] for _ in range(num_bins + 1)]

for plm in plumes:
    for cind, coordinate in enumerate(plm.Model.GZ['value']):
        for vind, GZval in enumerate(coordinate):
            try:
                AFval = plm.Model.AF['value'][cind][vind]
            except IndexError:
                continue

            if AFval >= 5. and GZval < plm.Model.height[cind]:
                try:
                    NO2val = plm.Model.TNO2['value'][cind][vind]
                    NOval = plm.Model.TNO['value'][cind][vind]

                    binned_rat[int(round(plm.distance[cind])) / 5].append(NO2val / (NOval + NO2val))

                except IndexError:
                    continue

binned_avgs = [np.nanmean(section) for section in binned_rat]
binned_stddev = [np.nanstd(section) for section in binned_rat]
x = [(i * 5) + 2.5 for i, _ in enumerate(binned_avgs)]

plt.errorbar(x, binned_avgs, yerr=binned_stddev, fmt='o', color='k', label='mean_ratio= ' + str(np.nanmean(binned_avgs)))
plt.legend(numpoints=1, loc='lower right', fontsize=14)
plt.title('Ratio over all plumes and layers using FireWork output', fontsize=14)
plt.xlim(0, 160)
plt.ylim(0, 1)
plt.xticks(range(0, 160, 20))
plt.xlabel('Distance from fire (km)', fontsize=14)
plt.ylabel('NO2 / (NO2 + NO)', fontsize=14)

plt.show()

plt.close('all')
