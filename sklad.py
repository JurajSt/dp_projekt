import datetime
import os
import utility

def fdayofyear(y, m, d):
    d0 = datetime.date(y, 1, 1)
    d1 = datetime.date(y, m, d)
    delta = (d1 - d0)
    return delta.days+1

obs_path = 'data/kame_2/'
eph_path = 'data/eph/test/eph/'#igr18644.sp3'
#path_obs = 'data/kame_1/KAME2740.15o'#'data/eph/1864/KAME2740.15o'

#obsfiles = next(os.walk(obs_path))[2]
ephfiles = next(os.walk(eph_path))[2]
#obsfiles = ['KAME2740.15o']
doy_list = []
vyska = []

for f in ephfiles:

    if 'sp3' != f.split('.')[1]:
        continue

    e = eph_path+f
    print e
    satellites = ['G05']#['G05', 'G07', 'G08', 'G15', 'G21', 'G22', 'G26', 'G28', 'G30']
    #satellites = ['G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G11', 'G12', 'G13', 'G14', 'G15', 'G16',
    #              'G17', 'G18', 'G19', 'G20', 'G21', 'G22', 'G23', 'G24', 'G25', 'G26', 'G27', 'G28', 'G29', 'G30', 'G31',
    #              'G32']

    sp3, date_e = utility.fReadSP3(e)
    doy_eph = fdayofyear(int(date_e[0]), int(date_e[1]), int(date_e[2]))
    if doy_eph > 0 and doy_eph < 10:
        namefile = '00' + str(doy_eph)+ '0'
    elif doy_eph >= 10 and doy_eph < 100:
        namefile = '0' + str(doy_eph) + '0'
    else: namefile = str(doy_eph) + '0'
    obsfile = obs_path + 'KAME' + namefile + '.17o'
    if not os.path.exists(obsfile):
        obsfile = obs_path + 'KAME' + namefile + '.16o'
    #else: continue

    print doy_eph, obsfile