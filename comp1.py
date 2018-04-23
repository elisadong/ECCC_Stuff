#! /usr/bin/env python2.7
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2018-03-01 15:50:27
# Contact:    Nolan.Dickson@Canada.ca

import postMINX as pm
import matplotlib.pyplot as plt
import glob
import math
import collections
import scipy.stats as stats
import numpy as np
import rpnpy


'''
Plots first a binned average graph of all MINX plume heights,
then FRP sorted model/MINX comparisons,
then biome sorted model/MINX comparisons
'''


def dist_binned(plm_list, model=False):
    # returns all points from all plm_list plumes grouped in 5km distance bins
    max_dist = max([plm.max('distance')[0]['distance'] for plm in plm_list])
    num_bins = int(round(max_dist)) / 5
    binned_list = [[] for _ in range(num_bins + 1)]

    for plm in plm_list:
        for pind, pnt in enumerate(plm.points):
            binned_list[int(round(pnt['distance'])) / 5].append(pnt if not model else plm.Model.points[pind])

    return binned_list


txts = glob.glob('/fs/site2/dev/eccc/aq/r1/nod001/MINX/Working/AB_SK_Sep2017/Digitize_Output/*/Plumes_O*-B*-SPWB*.txt')

plumes = []
for ind, file in enumerate(txts):
    print('{0}: ({1}/{2})'.format(file, ind, len(txts)))

    try:
        # CFFEPS #
        p = pm.Plume.from_filename(file, model_dir='/fs/site2/dev/eccc/aq/r2/jac001/FireWork/EvalCFFEPS_OS/model_output/OS2p5_CFFEPS')

        # FireWork #
        # p = pm.Plume.from_filename(file, model_dir='/fs/site2/dev/eccc/aq/r2/jac001/FireWork/EvalCFFEPS_OS/model_output/OS2p5_FireWork')

        # Newest CFFEPS Run #
        # p = pm.Plume.from_filename(file)
        # p.get_model('/space/hall1/sitestore/eccc/aq/r1/pic001/os_2p5km_dailyruns_ppp1hare_CFFEPS/outp_fst/OS2p5km/{d:%Y%m%d}120000/model'.format(d=p.datetime))

    except rpnpy.librmn.fstd98.FSTDError:
        continue

    if all([0 < h < 60000 or math.isinf(h) for h in p.Model.height]):
        # exclude plumes with model outliers (*used to* occur occasionally)

        p.set_release_time()
        plumes.append(p)

frp_grouped = collections.OrderedDict([
    ('all', plumes),  # not going to be altering any data so reference copy is aight here
    ('high ( (2000, inf) MW )', [plm for plm in plumes if plm.frp > 2000]),
    ('med ( [1000, 2000] MW )', [plm for plm in plumes if 1000 <= plm.frp <= 2000]),
    ('low ( [0, 1000) MW )', [plm for plm in plumes if plm.frp < 1000])
])

biome_grouped = {}
for plm in plumes:
    if int(plm.biome[1]) == 1 or int(plm.biome[1]) == 9:
        try:
            biome_grouped[plm.biome[0]].append(plm)
        except KeyError:
            biome_grouped[plm.biome[0]] = [plm]

# --------------------------------
# plots
# --------------------------------

binned_avgs = [np.nanmean([pnt['filtered_height'] - pnt['terrain'] for pnt in section]) for section in dist_binned(plumes)]

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_title('Average Filtered Height over 5km bins\nAB-SK Summer 2017')
x = [i * 5 for i, _ in enumerate(binned_avgs)]

ax.bar(x, binned_avgs, width=5, align='edge')
plt.xticks(x)
ax.set_ylabel('Height ASL (m)')
ax.set_xlabel('Distance from fire (km)')

plt.show()

fig = plt.figure()
for plm in plumes:
    fig = plm.shape(fig=fig)
fig.suptitle('Plume Shapes\nAB-SK Summer 2017')

plt.show()

# --------------------

for grouped in [(frp_grouped, 'FRP'), (biome_grouped, '')]:

    # fig = plt.figure()
    fig, axarr = plt.subplots(1, len(grouped[0]), sharey=True, sharex=True)
    if len(grouped[0]) == 1:
        axarr = [axarr]

    fig.suptitle('Average Height over 5km bins, Sorted by {0} (AB-SK Summer 2017)'.format(grouped[1]))

    for ind, level in enumerate(grouped[0], 1):
        ax = axarr[ind - 1]

        # MINX #

        # if ind == 1:  # if only want to do a single plot
        # ax = fig.add_subplot(1, 1, 1)

        # ax = fig.add_subplot(2, len(grouped[0]), ind)
        ax.set_title('{0} {1}\n{2} plumes'.format(level, grouped[1], len(grouped[0][level])))
        ax.set_ylim(0, 5)
        ax.set_xlim(0, 150)
        if ind == 1:
            ax.set_ylabel('Plume Height (km above terrain) (CFFEPS)')

        X_dbins = dist_binned(grouped[0][level])
        X_binned_avgs = [np.nanmean([pnt['filtered_height'] - pnt['terrain'] for pnt in section if pnt]) / 1000. for section in X_dbins]
        X_binned_stddev = [np.nanstd([pnt['filtered_height'] - pnt['terrain'] for pnt in section if pnt]) / 1000. for section in X_dbins]
        x = [i * 5 for i, _ in enumerate(X_binned_avgs)]

        # ax.bar(x, X_binned_avgs, width=5, align='edge')
        # ax.plot([xi + 2.5 for xi in x], X_binned_avgs, '-', color='k', label='Satellite')
        ax.errorbar([xi + 2.5 for xi in x], X_binned_avgs, yerr=X_binned_stddev, fmt='o', color='k', label='Satellite', linewidth=2.0)

        # plot number of points in each bin above the line  # TODO test and tweak with error bars, maybe near points instead
        for i, xi in enumerate(x):
            ax.text(xi, -1, str(len(X_dbins[i])), rotation=90, fontsize=10)

        # Model #

        # ax = fig.add_subplot(2, len(grouped[0]), ind)  # + len(grouped[0]))
        # ax.set_ylim(0, 5)
        # ax.set_xlim(0, 150)
        ax.set_xlabel('Distance from fire (km)')
        # if ind == 1:
        #     ax.set_ylabel('Model Height above terrain (km)')

        M_dbins = dist_binned(grouped[0][level], model=True)
        M_binned_avgs = [np.nanmean([pnt['height'] - X_dbins[si][pi]['terrain'] for pi, pnt in enumerate(section) if pnt]) / 1000.
                         for si, section in enumerate(M_dbins)]
        M_binned_stddev = [np.nanstd([pnt['height'] - X_dbins[si][pi]['terrain'] for pi, pnt in enumerate(section) if pnt]) / 1000.
                           for si, section in enumerate(M_dbins)]

        x = [i * 5 for i, _ in enumerate(M_binned_avgs)]

        # ax.bar(x, M_binned_avgs, width=5, align='edge', color='r')
        # ax.plot([xi + 2.5 for xi in x], M_binned_avgs, '-', color='r', label='Model')
        ax.errorbar([xi + 2.5 for xi in x], M_binned_avgs, yerr=M_binned_stddev, fmt='o', color='r', label='Model', linewidth=2.0)

        X_binned_avgs, M_binned_avgs = pm.clean_inv(X_binned_avgs, M_binned_avgs)

        mean_bias = (1. / len(M_binned_avgs)) * sum([MISR - M_binned_avgs[i] for i, MISR in enumerate(X_binned_avgs)])
        corr = stats.pearsonr(X_binned_avgs, M_binned_avgs)[0]
        ax.text(10, 4.5, 'mean bias = {0:0.2f}'.format(mean_bias))
        ax.text(10, 4.2, 'corr. coeff. = {0:0.2f}'.format(corr))
        print('{0}: {1}: bias: {2:0.2f}'.format(grouped[1] if grouped[1] else 'biome', level, mean_bias))
        print('{0}: {1}: corr: {2:0.2f}'.format(grouped[1] if grouped[1] else 'biome', level, corr))

        # plt.setp(ax.get_xticklabels())
        # plt.setp(ax.get_yticklabels())
        # plt.subplots_adjust(bottom=0.2)

    fig.subplots_adjust(left=0.04, right=0.94, bottom=0.2, top=0.85)
    ax.legend(numpoints=1, loc='upper center', bbox_to_anchor=(1.05, 0.5), ncol=1)

    plt.show()
