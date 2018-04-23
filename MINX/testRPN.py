#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author:     Nolan Dickson
# Last Edit:  2017-12-13 10:07:45
# Contact:    Nolan.Dickson@Canada.ca

import os
import numpy as np
import rpnpy.librmn.all as rmn
import datetime


# try:
#     file_id = rmn.fstopenall(file_name, rmn.FST_RO)
# except:
#     raise rmn.FSTDError("file not ofund/readable: {0}".format(file_name))

# try:
#     k1 = rmn.fstinf(file_id)['key']
#     k1_oth = rmn.fstinf(file_id)
#     # print k1_oth

#     TA3_key = rmn.fstinf(file_id, nomvar='TA3')['key']
#     TA3_info = rmn.fstinf(file_id, nomvar='TA3')
#     TA3_inl = rmn.fstinl(file_id, nomvar='TA3')

#     print TA3_info
#     # print TA3_inl
#     print '-'*8
# except Exception, e:
#     raise e
#     # raise rmn.FSTDError("go away we cant find those things")

# try:
#     # for key in TA3_inl:
#     #     print rmn.fstprm(key)['ip2']
#     TA3_meta = rmn.fstprm(TA3_key)
#     print TA3_meta
# except:
#     raise rmn.FSTDError("Error: Problem getting record metadata")

# try:
#     TA3_rec = rmn.fstluk(TA3_key)
#     # print TA3_rec
# except:
#     raise rmn.FSTDError("Error: Problem getting record data")

# rmn.fstcloseall(file_id)

inttype = rmn.EZ_INTERP_CUBIC

file_name = "/fs/site2/dev/eccc/aq/r2/jac001/FireWork/EvalCFFEPS_OS/model_output/OS2p5_CFFEPS/2017090312_017"
varname = 'AF'
lon = np.array([-115.3], dtype=np.float32)
lat = np.array([47.9], dtype=np.float32)

# try:
rmn.fstopt(rmn.FSTOP_MSGLVL, rmn.FSTOPI_MSG_CATAST)
funit = rmn.fstopenall(file_name, rmn.FST_RO)
keylist = rmn.fstinl(funit, nomvar=varname)

results = [{}]
v = []
h = []
ip = []
hip = []

# --------
for ind, key in enumerate(keylist):    # a key for each layer of this var
    meta = rmn.fstprm(key)
    meta['iunit'] = funit
    data = rmn.fstluk(key)['d']
    grid = rmn.ezqkdef(meta)

    try:
        gzkey = rmn.fstinf(funit, nomvar='GZ', ip1=meta['ip1'])['key']
        gzdata = rmn.fstluk(gzkey)['d']

        inttype = rmn.EZ_INTERP_CUBIC
        rmn.ezsetopt(rmn.EZ_OPT_INTERP_DEGREE, inttype)

        xypos = rmn.gdxyfll(grid, lat, lon)
        dataval = rmn.gdxysval(grid, xypos['x'], xypos['y'], data)          # arrays of values and heights for every desired latlong in this layer
        height = rmn.gdxysval(grid, xypos['x'], xypos['y'], gzdata) * 10

        v.append(dataval)
        h.append(height)
        ip.append(meta['ip1'])

        # print '--float--', rmn.DecodeIp(meta['ip1'], meta['ip2'], meta['ip3'])[0].v1, type(rmn.DecodeIp(meta['ip1'], meta['ip2'], meta['ip3'])[0])

        hip.append((height[0], meta['ip1'], rmn.DecodeIp(meta['ip1'], meta['ip2'], meta['ip3'])[0].v1))

    except TypeError, e:
        print meta['ip1'], 'is not in GZ ip1s'
        continue

# print lat[0]
# print lon[0]
# print np.array(v)[:, 0]
# print np.array(h)[:, 0]

# print ip
# print h
# print v

for x in sorted(hip): print x

for i, x in enumerate(results):
    x['datetime'] = datetime.datetime(2017, 8, 11, 19, 00).strftime('%Y-%m-%d, %H:%M')
    x['latitude'] = lat[i]
    x['longitude'] = lon[i]
    x['values'] = np.array(v)[:, i]
    x['heights'] = np.array(h)[:, i]
    # x['ip1'] = np.array(h)[:, i].fill(ip[i])
    # x['ip1'] = np.array(ip)[:, i]

    # print x['ip1']

    x['heights'], x['values'] = zip(*sorted(zip(x['heights'], x['values'])))


for x in results:

    print '\n', '-' * 20
    print 'latitude, longitude, date, time'
    print '{0}, {1}, {2}'.format(x['latitude'], x['longitude'], x['datetime'])
    print '\nvalue (µg/m³), height (m)'
    for i, v in enumerate(x['values']):
        print '{0}, {1}'.format(v, x['heights'][i])


# --------

    # gzk = rmn.fstinf(funit, nomvar='GZ')['key']

    # data = rmn.fstluk(k)['d']
    # gzdata = rmn.fstluk(gzk)['d']

    # meta = rmn.fstprm(k)
    # # print meta, '\n'
    # # print data, '\n'
    # # print rmn.DecodeIp(meta['ip1'], meta['ip2'], meta['ip3'])[0], '\n'
# except Exception, e:
#     print 'first part broke yo'
#     raise e
    # raise rmn.RMNError('Problem opening/reading var=%s in File=%s' % (varname,file_name))

# try:
#     meta['iunit'] = funit
#     grid = rmn.ezqkdef(meta)
# except Exception, e:
#     print 'couldnt define grid'
#     raise e

# try:
#     rmn.ezsetopt(rmn.EZ_OPT_INTERP_DEGREE, inttype)
# except Exception, e:
#     print 'couldnt set ezsetopt'
#     raise e

# # i
# xypos = rmn.gdxyfll(grid, lat, lon)

# # ii
# dataval = rmn.gdxysval(grid, xypos['x'], xypos['y'], data)

# # iii
# height = rmn.gdxysval(grid, xypos['x'], xypos['y'], gzdata) * 10


# lldata = rmn.gdllsval(grid, lat, lon, data)
# latlonpos = rmn.gdll(grid)

    # print '------RESULTS------'
    # print '| x,y: ', xypos['x'], xypos['y'], '|'
    # print '| val: ', dataval, '|'
    # print '| height: ', height, '|'
    # print '| level: ', rmn.DecodeIp(meta['ip1'], meta['ip2'], meta['ip3'])[0], '|'

# print '| grid: {0} -------|'.format(grid)
# print '| data: {0}'.format(data)
# print '|----- lldata: {0} ------|'.format(lldata)
# print '|------ llind: {0} -------|'.format(llind)
# print '|---- lats: {0}, \nlons: {1} ------|'.format(latlonpos['lat'], latlonpos['lon'])
# print '|------ posind: {0} -------|'.format(posind)

rmn.fstcloseall(funit)

# Old Method --------------

# for ind, val in enumerate(dataval):     # loop over all desired latlong values and heights for this layer

#     if val >= threshold and height[ind] > self.height[ind]:

#         self.height[ind] = height[ind]
#         self.value[ind] = val

#         self.points[ind]['height'] = height[ind]
#         self.points[ind]['value'] = val
