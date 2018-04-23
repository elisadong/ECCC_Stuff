#this calculation assumes that TH and QQ are in the file
#appends file to include QPOT
import rpnpy.librmn.all as rmn
import rpnpy.vgd.all as vgd
import numpy as np
import pygeode as pyg
import sys

#These are the coefficients for the vertical grid from MACC. I don't know
#how to get them out of the netcdf file (not even in ncdump!)
Ak = np.array([ 0.0,20.000000, 38.425343, 63.647804, 95.636963, 134.483307, 180.584351, 234.779053, \
       298.495789, 373.971924, 464.618134, 575.651001, 713.218079, 883.660522, 1094.834717, \
       1356.474609,1680.640259,2082.273926,2579.888672,3196.421631,3960.291504,4906.708496, \
       6018.019531, 7306.631348,8765.053711,10376.126953,12077.446289,13775.325195,15379.805664, \
       16819.474609,18045.183594,19027.695313,19755.109375,20222.205078,20429.863281,20384.480469, \
       20097.402344,19584.330078,18864.750000,17961.357422,16899.468750,15706.447266,14411.124023, \
       13043.218750,11632.758789,10209.500977,8802.356445,7438.803223,6144.314941,4941.778320, \
       3850.913330,2887.696533,2063.779785,1385.912598,855.361755,467.333588,210.393890,65.889244, \
       7.367743,0.000000,0.000000 ])
Bk = np.array([0.0,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000, \
      0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000, \
      0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,0.00007582,0.00046139,0.00181516, \
      0.00508112,0.01114291,0.02067788,0.03412116,0.05169041,0.07353383,0.09967469,0.13002251,0.16438432, \
      0.20247594,0.24393314,0.28832296,0.33515489,0.38389215,0.43396294,0.48477158,0.53570992,0.58616841, \
      0.63554746,0.68326861,0.72878581,0.77159661,0.81125343,0.84737492,0.87965691,0.90788388,0.93194032, \
      0.95182151,0.96764523,0.97966272,0.98827010,0.99401945,0.99763012,1.00000000])
#PS = 101325.
#now define the dictionary for the variables to get RPN tags
varnames = { ('aermr01') : {'nomvar' : 'TSS1', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 1.19689950054894e-08, 'scale_factor' : 3.65287035509047e-13 }, \
             ('aermr02') : {'nomvar' : 'TSS2', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 4.37358414179047e-07 , 'scale_factor' : 1.33479342665887e-11 }, \
             ('aermr03') : {'nomvar' : 'TSS3', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 1.6015011439134e-06, 'scale_factor' : 4.88769194870719e-11 }, \
             ('aermr04') : {'nomvar' : 'TCM1', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 4.52156804428452e-11 , 'scale_factor' : 1.48153698539026e-06 }, \
             ('aermr05') : {'nomvar' : 'TCM2', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 3.48243845183938e-07, 'scale_factor' : 1.06282074462534e-11 }, \
             ('aermr06') : {'nomvar' : 'TCM3', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 2.47611611041928e-07, 'scale_factor' : 7.55696792534724e-12 }, \
             ('aermr11') : {'nomvar' : 'TSU1', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 5.09227407801796e-08, 'scale_factor' : 1.55413357688395e-12 }, \
             ('go3')     : {'nomvar' : 'TO3',  'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' :  9.85818869816391e-06, 'scale_factor' : 3.00866406890188e-10 }, \
             ('co')      : {'nomvar' : 'TCO',  'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 6.44896118390205e-06, 'scale_factor' :  1.96818689309103e-10 }, \
             ('t')       : {'nomvar' : 'TT',   'conv_fac' : 1., 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 259.446889881307, 'scale_factor' : 0.00252150691776352 }, \
             ('lnsp')    : {'nomvar' : 'P0',   'conv_fac' : 1., 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 11.2055747073771, 'scale_factor' : 1.10543009042722e-05 }, \
             ('ch4')     : {'nomvar' : 'TCH4', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 1, 'scale_factor' : 1}, \
             ('nox')     : {'nomvar' : 'TNOX', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 3.99218436371045e-07, 'scale_factor' : 1.21839231633721e-11 }, \
             ('hcho')    : {'nomvar' : 'TNOX', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 1, 'scale_factor' : 1}, \
             ('so2')     : {'nomvar' : 'TSO2', 'conv_fac' : 1.e9, 'tv' : 'C', 'grtyp': 'Z', 'DTYPE' : 'f16', \
                                'add_offset' : 3.10082039836176e-07, 'scale_factor' : 9.46353017872734e-12 }, \
}


NK=1
DEET=0
NPAS=0
hfore=0
ETIKET= 'MACCREAN'
#define GEM grid so it is the same as in the MOZART files
#First read in a file from Verica that contains thegrid used for regional GEM_MACH
#This is a link in the directory you are in to Verica's file (wherever that is)
filename = "/users/dor/arqj/dip/mozart/moz4_JJ/mozart4geos5_201001_4p6.rpn"   ***********************<- change
try:
    fid = rmn.fstopenall(filename,rmn.FST_RO)  #read only
except:
    raise rmn.FSTDError("ERROR: File not found/readable: "+filename)

#get the grid that we want to use for the output
try:
    grec = rmn.fstlir(fid,nomvar='GZ')
    grid = rmn.readGrid(fid,grec)
except:
    raise rmn.FSTDError('ERROR: Could not read grid.')
latlon = rmn.gdll(grid)
ni = grid['ni'] ; nj = grid['nj']
rmn.fstfrm(fid)


#get the lnsp to convert vertical coordinates
ncidp = pyg.open('netcdf-atls06-a562cefde8a29a7288fa0b8b7f9413f7-BaXcgu.nc')
lnsp = ncidp.vars[0][:]                 #this is an array of (124,81,320)
mnsp = np.mean(lnsp,axis=0)
phalf = np.zeros([61,81,320])
for i in np.arange(0,61,1) :
    phalf[i,:,:] = Ak[i]+Bk[i]*np.exp(mnsp[:,:])
pfull = np.zeros([60,81,320])
pfull[0,:,:]=phalf[1,:,:]/2.
for i in np.arange(1,60,1) :
    pfull[i,:,:] = (phalf[i,:,:]+phalf[i+1,:,:])/2.
#PS = 101325.
#hyb_lvl = phalf[:]/PS

#now go through all the files and convert each variable
#for filename in filenames :
ncid = pyg.openall(filenames)

#get the axes
lon = ncid.axes[3].values
lat = ncid.axes[2].values
level = ncid.axes[1].values
times = ncid.axes[0].values
NI = len(lon)
NJ = len(lat)
#define grid for this file - note that the MACC grid in the file is
#defined for lons -180 to 180, but the python defGrid_L can't deal
#with that and defines the grid from 0 to 360 so will have to reorder
#the MACC fields a bit, or they end up 180 deg out of phase
# Also necessary to add one more longitude to wrap around
dlatlon=float(1.125)    ******************<- check
MACCgrid = rmn.defGrid_L(NI,NJ,lat0=0,lon0=0.,dlat=dlatlon,dlon=dlatlon)


#also need to put pfull onto the same horizontal grid as gem
outprs = np.zeros([gem_ni,gem_nj,60])
tmp = np.zeros([NJ,NI])
for k in np.arange(0,60,1) :
    tmp[:,0:160]=pfull[k,:,160:320]
    tmp[:,160:320]=pfull[k,:,0:160]
    outprs[:,:,k] = rmn.ezsint(gemgrid['id'],MACCgrid['id'],np.transpose(tmp))

#get the variables
nvars = len(ncid.vars)
outdat = np.zeros([gem_ni,gem_nj,60])
for ivar in np.arange(0,nvars,1) :
    varname = ncid.vars[ivar].name
    if varname == 'lnsp' :
        #skip this one
        print 'skipping lnsp....'
    else :

        TV = varnames[varname]['tv']
        NOMVAR = varnames[varname]['nomvar']
        DTYPE = varnames[varname]['DTYPE']
        GRTYP = varnames[varname]['grtyp']
        #do monthly mean level by level to avoid memory blowout
        for k in np.arange(0,60,1) :
            data = ncid.vars[ivar][:,k,:,:]*varnames[varname]['conv_fac']    #reads field (t,z,y,x)
            data_mn = np.mean(data,axis=0)   #takes monthly average which is what we want
            #interpolate to gem grid
            tmp[:,0:160]=data_mn[:,160:320]
            tmp[:,160:320]=data_mn[:,0:160]
            outdat[:,:,k] = rmn.ezsint(gemgrid['id'],MACCgrid['id'],np.transpose(tmp))

        #Now have the whole "cube" in on the correct horizontal grid
        #interpolate data to the correct hybrid level
        #NOTE: all these quantities should be on the thermo grid so use the ph levels from gem
        for ihyb in np.arange(0,82,1) :
            for i in np.arange(0,gem_ni,1) :
                for j in np.arange(0,gem_nj,1) :
                    gem_p0[i,j] = np.interp(ph[i,j,ihyb],outprs[i,j,:],outdat[i,j,:])

            rp1 = rmn.FLOAT_IP(hyh_mod[ihyb],hyh_mod[ihyb],rmn.LEVEL_KIND_HYB)
            rp2 = rmn.FLOAT_IP(0,0,rmn.TIME_KIND_HR)
            rp3 = rmn.FLOAT_IP(0,0,rmn.KIND_ARBITRARY)
            (IP1, IP2, IP3) = rmn.EncodeIp(rp1,rp2,rp3)
            FST_METADATA = {'dateo' : mydate,'deet' : 0, 'npas' : 0, 'ni' : gem_ni, 'nj' : gem_nj, 'nk' : 1, \
                            'ip1' : IP1, 'ip2' : 0, 'ip3' : 0, 'typvar' : TV, 'nomvar' : NOMVAR, 'etiket' : ETIKET, \
                            'grtyp' : GRTYP, 'ig1' : gem_lond['ip1'], 'ig2' : gem_lond['ip2'], \
                            'ig3' : gem_lond['ip3'], 'ig4' : 0 }
            rmn.fstecr(fid,gem_p0.astype('float32'),FST_METADATA,rewrite=True)

rmn.fstfrm(fid)
rmn.fclos(fid)

#*********************Calculation of pv using u, v, t and pressure***************************
**************change reading from files ************************

#set some stuff up for calculating the pv
a     = 6.371229e6        #earth radius in m
g0    = 9.81
kn2ms = 0.51444444     #conversion from knots to m/s
kappa = 0.287514
#actual lat-lon coords (not rotated)
cosphi = np.cos(lat*np.pi/180.)
f0_p   = 2*np.pi/86400.*np.sin(lat*np.pi/180.)
dx     = np.zeros([ni,nj])
dy     = np.zeros([ni,nj])
for i in np.arange(1,ni-1,1) :
    dx[i,:] = a*cosphi[i,:]*(lon[i+1,:]-lon[i-1,:])*np.pi/360.
dx[0,:] = a*cosphi[0,:]*(lon[1,:]-lon[0,:])*np.pi/180.
dx[-1,:] = a*cosphi[-1,:]*(lon[-1,:]-lon[-2,:])*np.pi/180.
for j in np.arange(1,nj-1,1) :
    dy[:,j] = a*(lat[:,j+1]-lat[:,j-1])*np.pi/360.
dy[:,0] = a*(lat[:,1]-lat[:,0])*np.pi/180.
dy[:,-1] = a*(lat[:,-1]-lat[:,-2])*np.pi/180.

#stuff on level k (uu and vv will be on this level)
dvdx=np.zeros([ni,nj])
dudy=np.zeros([ni,nj])
term1 = np.zeros([ni,nj,nk_mom])
#stuff on level k-1
uulast = np.zeros([ni,nj])
vvlast = np.zeros([ni,nj])
qqlast = np.zeros([ni_p,nj_p])
#stuff on level k+1/2 (tt will be on this level)
dthdx=np.zeros([ni,nj])
dthdy=np.zeros([ni,nj])
term2 = np.zeros([ni,nj,nk_thm])
#stuff on level k-1/2
thlast = np.zeros([ni,nj])
#Potential vorticity field
PV=np.zeros([ni,nj,nk_mom])


#now go level by level to calculate term1 and term2
for k in np.arange(0,nk_mom-1,1) :

    #get the necessary fields for the calcuation
    #uu and vv in m/s, qq in s^-1 and th in K
    if k > 0 :
        uulast = uu      #mom level k-1
        vvlast = vv      #mom level k-1
        qqlast = qq      #mom level k-1
        thlast = th      #thm level k-1/2
    try:
        uu = rmn.fstlir(fid,ip1=ip1_mom[k],nomvar='UU')['d']*kn2ms
        vv = rmn.fstlir(fid,ip1=ip1_mom[k],nomvar='VV')['d']*kn2ms
        qq = rmn.fstlir(fid,ip1=ip1_mom[k],nomvar='QQ')['d']
        #Note:  QQ is listed in the dictionary as absolute vorticity but it is sometimes relative vorticity
        th = rmn.fstlir(fid,ip1=ip1_thm[k],nomvar='TH')['d']   #thm level k+1/2
    except :
        raise rmn.FSTDError('ERROR: Missing field for PV calculation')

    #dthdy, dthdx on thm level k-1/2 and on the diagnostic grid because of the way I set up
    #the dx and dy (centred on i,j)
    dthdy[:,0]  = (thlast[:,1]-thlast[:,0])/dy[:,0]
    dthdy[:,-1] = (thlast[:,-1]-thlast[:,-2])/dy[:,-1]
    for j in np.arange(1,nj-1,1) :
        dthdy[:,j] = (thlast[:,j+1]-thlast[:,j-1])/dy[:,j]/2.
    dthdx[0,:]  = (thlast[1,:]-thlast[0,:])/dx[0,:]
    dthdx[-1,:] = (thlast[-1,:]-thlast[-2,:])/dx[-1,:]
    for i in np.arange(1,ni-1,1) :
        dthdx[i,:] = (thlast[i+1,:]-thlast[i-1,:])/dx[i,:]/2.


    #term1 is on mom level k
    #term2 is on htm level k-1/2
    #We want these on the diagnostic grid, not the physics grid so interpolate where necessary
    if k == 0 :
        term2[:,:,k] = np.zeros([ni,nj])
        term1[:,:,k] = np.zeros([ni,nj])
    else :
        term2[:,:,k] = (uu-uulast)/(hyb_mom[k]-hyb_mom[k-1])*dthdy - (vv-vvlast)/(hyb_mom[k]-hyb_mom[k-1])*dthdx
        
        dudy[:,0]  = (uu[:,1]*cosphi[:,1]-uu[:,0]*cosphi[:,1])/dy[:,0]
        dudy[:,-1] = (uu[:,-1]*cosphi[:,-1]-uu[:,-2]*cosphi[:,-2])/dy[:,-1]
        for j in np.arange(1,nj-1,1) :
            dudy[:,j] = (uu[:,j+1]*cosphi[:,j+1]-uu[:,j-1]*cosphi[:,j-1])/dy[:,j]/2.
        dvdx[0,:]  = (vv[1,:]-vv[0,:])/dx[0,:]
        dvdx[-1,:] = (vv[-1,:]-vv[-2,:])/dx[-1,:]
        for i in np.arange(1,ni-1,1) :
            dvdx[i,:] = (vv[i+1,:]-vv[i-1,:])/dx[i,:]/2.
        qq  = f0_p + dvdx - dudy/cosphi

        term1[:,:,k] = (th-thlast)/(hyb_thm[k]-hyb_thm[k-1]) * tmp

#now the terms for k=nk
uulast = uu      #mom level N
vvlast = vv      #mom level N
qqlast = qq      #mom level N
thlast = th      #thm level N-1/2
#get the necessary fields for the calculation
try:
    #get the diagnostic level to stand in as hybrid 1 (5002 is okay, 5005 is approximation)
    th = rmn.fstlir(fid,ip1=dg_thm,nomvar='TH')['d']
    uu = rmn.fstlir(fid,ip1=dg_mom,nomvar='UU')['d']*kn2ms
    vv = rmn.fstlir(fid,ip1=dg_mom,nomvar='VV')['d']*kn2ms
    #note that QQ is not on this surface for 5005
except :
    raise rmn.FSTDError('ERROR: Missing field for PV calculation')

#dtdy, dtdx on diagnostic level
dthdy[:,0] = (th[:,1]-th[:,0])/dy[:,0]
dthdy[:,-1] = (th[:,-1]-th[:,-2])/dy[:,-1]
for j in np.arange(1,nj-1,1) :
    dthdy[:,j] = (th[:,j+1]-th[:,j-1])/dy[:,j]/2.

dthdx[0,:] = (th[1,:]-th[0,:])/dx[0,:]
dthdx[-1,:] = (th[-1,:]-th[-2,:])/dx[-1,:]
for i in np.arange(1,ni-1,1) :
    dthdx[i,:] = (th[i+1,:]-th[i-1,:])/dx[i,:]/2.

#term2 is on full level N+1/2 (note one more level for term2
term2[:,:,-1] = (uu-uulast)/(hyb_mom[-1]-hyb_mom[-2])*dthdy - (vv-vvlast)/(hyb_mom[-1]-hyb_mom[-2])*dthdx

PV=np.zeros([ni,nj,nk_mom])
for k in np.arange(1,nk_mom,1) :
    PV[:,:,k] = -g0/dpdnu_mom[:,:,k] * (term1[:,:,k] + (term2[:,:,k]+term2[:,:,k-1])/2. )

