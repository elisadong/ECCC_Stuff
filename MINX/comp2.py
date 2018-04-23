#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2017-12-18 17:13:10
# Contact:    Nolan.Dickson@Canada.ca

import os
import math
import glob
import numpy as np
import postMINX as pm
import rpnpy.librmn.all as rmn
import matplotlib.pyplot as plt

txts = glob.glob('/fs/site2/dev/eccc/aq/r1/nod001/MINX/Working/AB_SK_*2017/Digitize_Output/*/Plumes_O*-B*-SPWB*.txt')

plumes = []
for ind, file in enumerate(txts, 1):
    print '{0}: ({1}/{2})'.format(file, ind, len(txts))
    p = pm.Plume.from_filename(file, model_dir='/fs/site2/dev/eccc/aq/r2/jac001/FireWork/EvalCFFEPS_OS/model_output/OS2p5_CFFEPS')

    if all([0 < h < 60000 or math.isinf(h) for h in p.Model.height]):
        plumes.append(p)

res = {
    # # key for every layer
    # ip1_vals: [ # dict for every coordinate
    #    {'NO2': val, 'NOx': val, 'ratio': val, 'height': val, 'distance': val},
    # ]
}

for plm in plumes:
    lon_list, lat_list = zip(*[(lon, plm.latitude[i]) for i, lon in enumerate(plm.longitude)])

    lon_list = np.array(lon_list, dtype=np.float32)
    lat_list = np.array(lat_list, dtype=np.float32)

    model_file = os.path.join(plm.Model.dir, plm.Model.name)

    funit = rmn.fstopenall(model_file, rmn.FST_RO)
    NO2_keylist = rmn.fstinl(funit, nomvar='TNO2')

    for ind, key in enumerate(NO2_keylist):    # a key for each layer of this var
        meta = rmn.fstprm(key)
        meta['iunit'] = funit
        NO2data = rmn.fstluk(key)['d']
        grid = rmn.ezqkdef(meta)

        res[meta['ip1']] = [] if not res.get(meta['ip1'], False) else res[meta['ip1']]

        try:
            gzkey = rmn.fstinf(funit, nomvar='GZ', ip1=meta['ip1'])['key']
            gzdata = rmn.fstluk(gzkey)['d']

            NOkey = rmn.fstinf(funit, nomvar='TNO', ip1=meta['ip1'])['key']
            NOdata = rmn.fstluk(NOkey)['d']

            inttype = rmn.EZ_INTERP_CUBIC
            rmn.ezsetopt(rmn.EZ_OPT_INTERP_DEGREE, inttype)
            xypos = rmn.gdxyfll(grid, lat_list, lon_list)

            NO2val = rmn.gdxysval(grid, xypos['x'], xypos['y'], NO2data)          # arrays of values and heights for every desired latlong in this layer
            height = rmn.gdxysval(grid, xypos['x'], xypos['y'], gzdata) * 10
            NOval = rmn.gdxysval(grid, xypos['x'], xypos['y'], NOdata)

            for vind, val in enumerate(NO2val):
                res[meta['ip1']].append({
                    'NO2': val,
                    'NOx': val + NOval[vind],
                    'height': height[vind] - plm.terrain[vind],
                    'ratio': val / (val + NOval[vind]),

                    'distance': plm.distance[vind],  # fst functions retain exact order, so can use index to access plume data as well
                    'plume_ht': plm.Model.height[vind] - plm.terrain[vind]
                })

        except TypeError, e:
            print meta['ip1'], 'is not in GZ or TNO ip1s'
            continue

# fig = plt.figure()
# fig.suptitle('NO2 vs NOx in smoke plume')

print 'mean height of layer    mean ratio of layer'
print '--------------------    -------------------'

totx, toty = [], []
for key in res:
    # key for every layer
    try:
        x, y, lay = zip(*[(coord['distance'], coord['ratio'], coord['height']) for coord in res[key] if coord['plume_ht'] - 1000 < coord['height'] < coord['plume_ht'] + 500])

        mean_ht = np.mean(lay)
        mean_ratio = np.mean(y)

        totx += x
        toty += y

        # ax = fig.add_subplot(2, 3, math.ceil(mean_ht / 7500.))
        # ax = fig.add_subplot(1, 1, 1)

        plt.plot(x, y, 'o', label='mean_ratio= ' + str(mean_ratio))
        plt.title('ratio at layer with mean height: ' + str(mean_ht))
        plt.xlim(0, 160)
        plt.ylim(0, 1)
        plt.xlabel('Distance from fire')
        plt.ylabel('NO2 / (NO2 + NO)')
        plt.legend(numpoints=1, loc='lower right')
        plt.show()

        print mean_ht, '          ', mean_ratio

    except (ValueError, KeyError), e:
        # print '-------'
        # print key, "doesn't have either NO2 or NOx or height or something who knows."
        # print "If its 76696048 then just ignore it, I honestly don't know what that layer is."
        # print e
        # print 'addendum, its actually cause one of the heights isnt with 1 km of plume height, for now'
        # print '-------'
        continue

tot_mean_ratio = np.mean(toty)
print '---------------------------------------------'
print 'total mean ratio: ', tot_mean_ratio

plt.plot(totx, toty, 'o', label='total_mean_ratio= ' + str(tot_mean_ratio))
plt.legend(numpoints=1, loc='lower right')
plt.title('Ratio at all layers')
plt.xlim(0, 160)
plt.ylim(0, 1)
plt.xlabel('Distance from fire (km)')
plt.ylabel('NO2 / (NO2 + NO)')

plt.show()

plt.close('all')
