# -*- coding: utf-8 -*-

import os, sys
import datetime, time
import numpy as np
import xlsxwriter
import utility
#from input_data import *
import math as m
import model
import scipy.signal as signal
import matplotlib.pyplot as plt
import funkcie


obs_path = 'data/kame_2/'  # 'D:/diplomka/telg_2/'
eph_path = 'data/eph/test/eph/'  # igr18644.sp3'
# path_obs = 'data/kame_1/KAME2740.15o'#'data/eph/1864/KAME2740.15o'

# obsfiles = next(os.walk(obs_path))[2]
ephfiles = next(os.walk(eph_path))[2]
# obsfiles = ['KAME2740.15o']
doy_list = []
vyska = []
kame = [3892532.358, 1572220.333, 4785952.565]
for f in ephfiles:

    if 'sp3' != f.split('.')[1]:
        continue

    e = eph_path + f
    print e
    satellites = ['G21']  # ['G05', 'G07', 'G08', 'G15', 'G21', 'G22', 'G26', 'G28', 'G30']
    # satellites = ['G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G11', 'G12', 'G13', 'G14', 'G15', 'G16',
    #              'G17', 'G18', 'G19', 'G20', 'G21', 'G22', 'G23', 'G24', 'G25', 'G26', 'G27', 'G28', 'G29', 'G30', 'G31',
    #              'G32']

    sp3, date_e = utility.fReadSP3(e)

    doy_eph = funkcie.fdayofyear(int(date_e[0]), int(date_e[1]), int(date_e[2]))

    if doy_eph > 0 and doy_eph < 10:
        namefile = '00' + str(doy_eph) + '0'
    elif doy_eph >= 10 and doy_eph < 100:
        namefile = '0' + str(doy_eph) + '0'
    else:
        namefile = str(doy_eph) + '0'
    obsfile = obs_path + 'KAME' + namefile + '.17o'
    if not os.path.exists(obsfile):
        obsfile = obs_path + 'KAME' + namefile + '.16o'
    # else: continue
    r_obs = open(obsfile, "r")
    lines = r_obs.readlines()
    r_obs.close()

    header, version, headlines, headlength, obstimes, sats, numallsvs, numsvs, date_o = funkcie.scan(lines)
    data, obstypes = funkcie.processBlocks(lines, header, obstimes, headlines, headlength, sats, numallsvs, satellites)

    del lines
    deg = 6
    indexS1 = obstypes.index('S1') + 2
    indexS2 = obstypes.index('S2') + 2

    # svs = []
    den_hodiny = []
    for sat in satellites:
        print ("pocitam pre: " + sat)
        t = []
        signals1 = []
        signals2 = []
        sinele = []

        for k, d in enumerate(data):

            if sat not in sats[k]:
                continue
            # print "pocitam pre: " + sat, k

            for i in d:  # jedna sekunda pre x druzic
                # print i
                sv = i[0]
                if sat not in sv:
                    continue
                # print sv
                # signals = (i[indexS1], i[indexS2])
                lg = utility.fInterpolateSatelliteXYZ(sp3, sv, i[1])
                azi, ele = utility.fComputeAziEle(kame, [lg[0], lg[1], lg[2]])

                if ele >= 5 and ele <= 30 and azi >= 60 and azi <= 100:

                    t.append(i[1])
                    sinele.append(m.sin(m.radians(ele)))
                    # signals1.append(m.pow(10,float(i[indexS1])/20))
                    signals2.append(m.pow(10, float(i[indexS2]) / 20))
                    # print ele, m.sin(ele), i[1]
                    break

        # print max(sinele), len(signals1), len(signals2)
        ind_max = sinele.index(max(sinele))             # rozdelenie primaneho signalu na vzostupny  a zostupny
        for l in range(0, 2):
            if l == 0:
                s1 = signals1[:ind_max + 1]
                s2 = signals2[:ind_max + 1]
                sinele2 = sinele[:ind_max + 1]
                t2 = t[:ind_max + 1]
                # print "pre " + sat + " počet observacií: " + str(len(sinele2))
                # print len(sinele2), sinele2[-1]
            else:
                s1 = signals1[ind_max + 1:]
                s2 = signals2[ind_max + 1:]
                sinele2 = sinele[ind_max + 1:]
                t2 = t[ind_max + 1:]
                # print len(sinele2), sinele2[0]
                # print "pre " + sat + " počet observacií: " + str(len(sinele2))
            if len(sinele2) < 2000:
                print "pre ", sat, " málo observácií (" + str(len(sinele2)) + ") - vypocet neprebehne"
                continue

            print "pre " + sat + " počet observacií: " + str(len(sinele2))
            doy_list.append(doy_eph)
            SNR1 = np.array(s1)  # zoznam na np.array
            SNR2 = np.array(s2)
            x = np.array(sinele2)
            # den_sek = np.array(cas[j])  # cas dna merania v sekundach
            timenum = np.array(t2)
            f = []

            alt_prob = np.linspace(1.0, 2.0, np.size(timenum))

            # print alt_prob
            for i in range(len(alt_prob)):
                f.append(((2 / 0.24) * alt_prob[i]))

            modelSNR2, residuals2, SNR2_vector = model.modelSNR(timenum, SNR2, deg)
            modelSNR2_1, residuals2_1, SNR2_vector_1 = model.modelSNR(x, SNR2_vector, deg)

            pgram = signal.lombscargle(timenum, SNR2_vector_1, f, )  # vypocet periodogramu

            ind = np.ndarray.argmax(pgram)

            vyska.append(alt_prob[ind])
            print alt_prob[ind]

            plt.subplot(3, 1, 1)
            plt.plot(x, SNR2_vector_1, x, modelSNR2_1, 'r')
            plt.xlabel('sin ')
            plt.ylabel('residuals (volt)')

            plt.subplot(3, 1, 2)
            plt.plot(x, SNR2_vector_1)
            plt.xlabel('sin ')
            plt.ylabel(' SNR2 (db)')
            # plt.show()

            plt.subplot(3, 1, 3)
            plt.plot(alt_prob, pgram)
            plt.xlabel('Reflector height (m)')
            plt.ylabel('Spectral Ampl. ')
            plt.show()

print len(doy_list), len(vyska)
plt.plot(doy_list, vyska, 'o')
plt.show()



