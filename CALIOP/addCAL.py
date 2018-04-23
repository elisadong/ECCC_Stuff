# Add CALIOP data

import os
import postMINX as pm
import datetime
import niceUtilities as nu

try:
    import cPickle as pickle
except ImportError:
    import pickle
# for now do the september cases

picklePath = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/Working/apr_oct/FW-GM_start00_pklPlmsOps/"
savePklPath = './FW-GM_CAL_start00_pklPlmsOps'
plumeList = os.listdir(picklePath)

septInts = "/fs/site2/dev/eccc/aq/r1/eld001/MINX/CALIOP/intercepts_s20170901_e20170930.p"
nu.makeDir(savePklPath)
septArray = []
for p in  plumeList:
    plm=pickle.load(open(os.path.join(picklePath, p), 'rb'))
    print(plm)
    plmDate =plm['datetime']
    mnth = plmDate.month
    if mnth == 9:
        septArray +=[plm]

firePath = '/fs/site2/dev/eccc/aq/r1/eld001/MINX/CALIOP/intercepts_s20170901_e20170930.p'
fires = pickle.load(open(firePath,'rb'))
firePasses = []
fireBounds = []
keptFires = []
fDays = []

for f in fires:
    if f.overpasses != []:
        firePasses += [f.overpasses]
        fireBounds += [f.boundaries]
        keptFires += [f]

for o in firePasses:
    day = o[0].day
    fDays += [day]

matchArray = []
for s in septArray:
    sdate = s.datetime
    if sdate.day in fDays:
        matchArray += [s]









def get_caliop(self, fire, cal_path=os.getcwd()):
    self.Caliop = _Caliop(plm=self, fire = fire, cal_path=cal_path)
    return

matchArray = []
for s in septArray:
    sdate = s.datetime
    if sdate.day in fDays:
        matchArray += [s]

#check if origin inside of caliop boundary

class _Caliop(object):
    def __getitem__(self, key):
        return self.__dict__[key]
    def __init__(self, plm, fire, cal_path=os.getcwd(), verbose_error=False, auth=('elisadong', 'elisaPass1')):
        import copy
        import ftplib
        import MISRutils as mu
        import pyhdf.SD as SD
        from Cpoly.Cpoly import polygon_check_single, polygon_check

        self.Plume = plm
        ftp_url = 'ftp.icare.univ-lille1.fr'
        ftp_dir = 'SPACEBORNE/CALIOP/VFM.v4.10/{yr}/{date}/'.format(yr=fire.overpasses[0].year, date=fire.overpasses[0].strftime('%Y_%m_%d'))

        icare = ftplib.FTP(ftp_url, user=auth[0], passwd=auth[1])
        icare.cwd(ftp_dir)
        file_list = icare.nlst()

        vind, vfile = min(enumerate(file_list),
                          key=lambda filename: abs(datetime.datetime.strptime(filename[1][30:49], '%Y-%m-%dT%H-%M-%S') - fire.overpasses[0]
                                                   if datetime.datetime.strptime(filename[1][30:49], '%Y-%m-%dT%H-%M-%S') < fire.overpasses[0]
                                                   else datetime.timedelta(999)))
        delta_t = datetime.datetime.strptime(vfile[30:49], '%Y-%m-%dT%H-%M-%S') - fire.overpasses[0]
        if 0 < delta_t.seconds < 300:
            # file is within 5 minutes greater than intersect, grab later one
            inter_files = [vfile, file_list[vind + 1]]
        elif -300 < delta_t.seconds < 0:
            # file is within 5 minutes less than intersect, grab earlier one
            inter_files = [file_list[vind - 1], vfile]
        else:
            inter_files = [vfile]

        lonarr, latarr = np.empty((2, 0))
        fcf_arr = np.empty((0, 5515))

        for find, ftp_VFM_file in enumerate(inter_files):
            local_VFM_path = '{dir}/{file}'.format(dir=cal_path, file=ftp_VFM_file)
            if not os.path.exists(local_VFM_path):  # if dont have locally, download
                with open(local_VFM_path, 'wb') as f:
                    icare.retrbinary('RETR ' + ftp_VFM_file, f.write)
            VFM_file = SD.SD(local_VFM_path, SD.SDC.READ)
            # TODO ss stands for subsatellite, must ensure its what I think it is
            lonarr = np.concatenate((lonarr, VFM_file.select('ssLongitude')[:, 0]))
            latarr = np.concatenate((latarr, VFM_file.select('ssLatitude')[:, 0]))
            fcf_arr = np.concatenate((fcf_arr, VFM_file.select('Feature_Classification_Flags').get()))

        icare.close()

        plm_inter_flag = 0
        polypoints = fire.boundaries
        for pind, (lon, lat) in enumerate(zip(lonarr, latarr)):
            # TODO meridian problem again
            try:
                if polygon_check(polypoints, [(lon, lat), (lonarr[pind + 1], latarr[pind + 1])]):
                    plm_inter_flag = pind + 1
                    break
            except IndexError:
                break

        if plm_inter_flag:
            # intersects, grab a bunch of nearby points (+- 1500 covers about 2 degrees of lat/long usually)
            # (need to modulize to nearest multiple of 15 for fcf array indices, which are only 3744 long)
            cal_cut_ind = (
                plm_inter_flag - 1500 - (plm_inter_flag % 15) if plm_inter_flag - 1500 - (plm_inter_flag % 15) > 0 else 0,
                plm_inter_flag + 1500 - (plm_inter_flag % 15)
            )

            lonarr = lonarr[cal_cut_ind[0]: cal_cut_ind[1]]
            latarr = latarr[cal_cut_ind[0]: cal_cut_ind[1]]
            fcf_arr = np.uint(fcf_arr[int(cal_cut_ind[0] / 15.): int(cal_cut_ind[1] / 15.), 1165:])  # cut off the profile at 8.2km alt, for simplicity sake

            fcf_arr = fcf_arr & 7
            # fcf_arr[(fcf_arr != 3) & (fcf_arr != 4)] = 0  # only focus on aerosols (not sure if faster to get rid of others here or check better down there)

            # TODO docs explaining this weird data format
            # cut the fcf array (which is only till 8.2km) into the 290 high cols, find the index of the highest flag, find its alt from range
            alt_range = [-0.5 + i * 0.03 for i in range(289, -1, -1)]
            fcf_arr = fcf_arr.reshape(fcf_arr.shape[0] * 15, 290)

            # ref_lvl = self.Plume.reference

            self.longitude = lonarr
            self.latitude = latarr
            self.flags = fcf_arr
            self.aerosols = copy.deepcopy(fcf_arr)
            self.aerosols[(self.aerosols != 3) & (self.aerosols != 4)] = 0
            self.height = [float('Nan')] * len(fcf_arr)
            for rind, r in enumerate(fcf_arr):
                for zind, z in enumerate(fcf_arr[rind]):
                    if z == 3 or z == 4:  # search for aerosols
                        self.height[rind] = alt_range[zind] * 1000.
                        break
            self.Plume.units['caliop_height'] = 'm'
