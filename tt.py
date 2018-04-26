# -*- coding: utf-8 -*-

import os, sys
import datetime, time
import numpy as np
import psycopg2
from numpy import *
import xlsxwriter
import utility
from input_data import *
import math as m
import model
from azi_ele_coor import *
import scipy.signal as signal
import matplotlib.pyplot as plt
import ogr
import t3
import matplotlib.patches as mpatches
import lomb
#from astropy.stats import LombScargle

obs_path = 'data/ganp/obs_jan/'#'D:/diplomka/telg_2/'
eph_path = 'data/ganp/eph_jan/'#igr18644.sp3'
#eph_path = 'D:/diplomka/dp_projekt/data/eph/1921/'#'data/eph/1864/KAME2740.15o'

#obsfiles = next(os.walk(obs_path))[2]
ephfiles = next(os.walk(eph_path))[2]
#ephfiles = ['igr19326.sp3'] # 'igr19181.sp3'
doy_list = []
vyska = []
a = 0.0
func= 1000

connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()

cursor.execute("DROP TABLE IF EXISTS test")
cursor.execute("CREATE TABLE test (stname CHAR (4), azimuthM REAL, elevangleM REAL, azimuthJ REAL, elevangleJ REAL, azi REAL, geom GEOMETRY)")


for f in ephfiles:

    if 'sp3' != f.split('.')[1]:
        continue

    e = eph_path+f
    print e
    #telg = sat 5,6,9,16, 21,23,29,30,31
    satellites = ['G19']#, 'G07', 'G08', 'G15', 'G21', 'G22', 'G26', 'G28', 'G30']
    #satellites = ['G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G11', 'G12', 'G13', 'G14', 'G15', 'G16',
    #              'G17', ']G18', 'G19', 'G20', 'G21', 'G22', 'G23', 'G24', 'G25', 'G26', 'G27', 'G28', 'G29', 'G30', 'G31',
    #              'G32']

    sp3, date_e = utility.fReadSP3(e)
    doy_eph = fdayofyear(int(date_e[0]), int(date_e[1]), int(date_e[2]))
    dat = datetime.date(int(date_e[0]), int(date_e[1]), int(date_e[2]))

    print dat, doy_eph
    if doy_eph > 0 and doy_eph < 10:
        namefile = '00' + str(doy_eph)+ '0'
    elif doy_eph >= 10 and doy_eph < 100:
        namefile = '0' + str(doy_eph) + '0'
    else: namefile = str(doy_eph) + '0'
    obsfile = obs_path + 'ganp' + namefile + '.12o'
    if not os.path.exists(obsfile):
        obsfile = obs_path + 'ganp' + namefile + '.12o'
        if not os.path.exists(obsfile):
            continue
    #else: continue
    print obsfile
    r_obs = open(obsfile, "r")
    lines = r_obs.readlines()
    r_obs.close()

    """
    row = 1
    col = 0
    workbook = xlsxwriter.Workbook('D:\diplomka\matlabCode\hofn'+namefile+'.xlsx')
    worksheet = workbook.add_worksheet()
    worksheet.write(0, 0, 's2')
    worksheet.write(0, 1, 'sin')
    worksheet.write(0, 2, 't2')"""

    header, version, headlines, headlength, obstimes, sats, numallsvs, numsvs, date_o = scan(lines)
    data, obstypes = processBlocks(lines, header, obstimes, headlines, headlength, sats, numallsvs, satellites)

    del lines
    deg = 6
    indexS1 = obstypes.index('S1')+2
    indexS2 = obstypes.index('S2')+2

    #svs = []
    den_hodiny = []
    for sat in satellites:
        print ("pocitam pre: " + sat)
        t = []
        signals1 = []
        signals2 = []
        #ss = []
        #eles = []
        sinele = []
        #tims = []

        for k, d in enumerate(data):
            #print d
            #print k, sats[k]

            if sat not in sats[k]:
                continue
            #print "pocitam pre: " + sat, k

            for i in d:             # jedna sekunda pre x druzic
                #print i
                sv = i[0]
                if sat not in sv:
                   continue
                #print sv
                #signals = (i[indexS1], i[indexS2])
                lg = utility.fInterpolateSatelliteXYZ(sp3, sv, i[1])
                azi, ele = utility.fComputeAziEle(ganp[0], [lg[0], lg[1], lg[2]])
                azi2, ele2 = t3.aziEle(ganp[0], [lg[0], lg[1], lg[2]])
                Lat, Lon, hel = fXYZ_to_LatLonH(ganp[0][0], ganp[0][1], ganp[0][2])
                Latsv, Lonsv, helsv = fXYZ_to_LatLonH(lg[0], lg[1], lg[2])
                WGSmercatorX, WGSmercatorY = fwgs4326towgs3857(Lat, Lon)
                WGSmercatorXsv, WGSmercatorYsv = fwgs4326towgs3857(Latsv, Lonsv)
                azi3 = fCalculateAzimuth([WGSmercatorX, WGSmercatorY], [WGSmercatorXsv, WGSmercatorYsv])
                azi4 = fCalculateAzimuth(hofn[0], [lg[0], lg[1], lg[2]])
                #print azi, azi2, azi3
                line = ogr.Geometry(ogr.wkbLineString)
                line.AddPoint(WGSmercatorX, WGSmercatorY)
                line.AddPoint(WGSmercatorXsv, WGSmercatorYsv)
                wkt_line = line.ExportToWkt()

                cursor.execute("INSERT INTO test (stname, azimuthM, elevangleM, azimuthJ, elevangleJ, azi, geom)"
                               " VALUES ('%s',  %s, %s, %s, %s, %s, ST_GeometryFromText('%s'))" % (sv, azi, ele, azi2, ele2, azi3, wkt_line))
                del line, wkt_line
                continue
                #ss.append(float(i[indexS2]))
                #eles.append(ele)
                #tims.append(i[1])
                if 5 <= ele <= 30 and 180 <= azi3 <= 220:
                    #print azi, ele, i[1]
                    t.append(i[1])
                    sinele.append(m.sin(m.radians(ele)))
                    #signals1.append(m.pow(10,float(i[indexS1])/20))
                    signals2.append(m.pow(10,float(i[indexS2])/20))
                    break

        if len(sinele) < 70:
            print "pre ", sat, " málo observácií - vypocet neprebehne", len(sinele)
            continue
        """
        maxim = max(sinele)
        print maxim
        n=1

        sinele4 = []
        a = 0
        for i in range(len(sinele)):
            if i < maxim and maxim < 1:
                c = sinele[i]*(-1)
            elif i > maxim and maxim < 1:
                c = sinele[i]
            if i == a:
                sinele4.append(c)


        print sinele[0], sinele[-1]
        print t[0], t[-1]
        
        fig = plt.figure()
        ax1 = fig.add_axes((0.1, 0.2, 0.8, 0.6))  # create an Axes with some room below
        ax1.set_title('SNR GPS L2')
        ax1.plot(tims, ss)
        ax1.set_xlabel('time (h)')
        ax1.set_ylabel('SNR (dB)')
        ax1_ticks = ax1.get_xticks()
        diff_ax1_ticks = ax1_ticks[-1] - ax1_ticks[0]
        a = 0
        b = int(len(tims) / 4)
        tt = []
        sinele4 = []
        for i in range(len(tims)):
            if i == a:
                tt.append((tims[i] / diff_ax1_ticks) - (ax1_ticks[0] / diff_ax1_ticks))
                sinele4.append(round(eles[i],2))
                a = a+b
        tt.append((t[-1] / diff_ax1_ticks) - (ax1_ticks[0] / diff_ax1_ticks))
        sinele4.append(round(eles[-1], 2))
        # create second Axes. Note the 0.0 height
        ax2 = fig.add_axes((0.1, 0.1, 0.8, 0.0))
        ax2.yaxis.set_visible(False)  #
        ax2.set_xlabel('elevation angle')
        ax2.set_xticks(tt)
        ax2.set_xticklabels(sinele4)
        plt.show()"""


        #print max(sinele), len(signals1), len(signals2)
        ind_max = sinele.index(max(sinele))
        for l in range(0,2):
            if l == 0:
                s1 = signals1[:ind_max+1]
                s2 = signals2[:ind_max+1]
                sinele2 = sinele[:ind_max+1]
                t2 = t[:ind_max+1]
                #print "pre " + sat + " počet observacií: " + str(len(sinele2))
                #print len(sinele2), sinele2[-1]
            else:
                s1 = signals1[ind_max+1:]
                s2 = signals2[ind_max+1:]
                sinele2 = sinele[ind_max+1:]
                t2 = t[ind_max+1:]
                #print len(sinele2), sinele2[0]
                #print "pre " + sat + " počet observacií: " + str(len(sinele2))
            if len(sinele2) < 70:
                print "pre ", sat, " málo observácií (" + str(len(sinele2)) + ") - vypocet neprebehne"
                continue

            print "pre " + sat + " počet observacií: " + str(len(sinele2))

            SNR1 = np.array(s1)  # zoznam na np.array
            SNR2 = np.array(s2)
            x = np.array(sinele2)
            #den_sek = np.array(cas[j])  # cas dna merania v sekundach
            timenum = np.array(t2)
            vyska_odrazec = []

            #print np.size(timenum)

            alt_prob = np.linspace(1, 2, np.size(timenum))   # odhadovaná výška odrážača
            funkcia = []

            #print alt_prob
            for i in range(len(alt_prob)):
                funkcia.append((2*alt_prob[i])/0.244)

            #print vyska_odrazec
            #print funkcia
            modelSNR2, residuals2, SNR2_vector = model.modelSNR(timenum, SNR2, deg)
            modelSNR2_1, residuals2_1, SNR2_vector_1 = model.modelSNR(x, SNR2_vector, deg)
            """
            p, hist = np.histogram(modelSNR2_1)
            a = plt.hist(SNR2_vector_1)
            plt.show()
            sdata = []
            tim = []
            SNR2L = SNR2_vector_1.tolist()
            for d in SNR2_vector_1:
                if  max(SNR2L) >= d >= 5-max(SNR2L) or 5+min(SNR2L) >= d >= min(SNR2L):
                    dind = SNR2L.index(d)
                    tim.append(timenum[dind])
                    sdata.append(d)
            
            for i in range (np.size(timenum)):
                worksheet.write(row, col, SNR2[i])
                worksheet.write(row, col+1, x[i])
                worksheet.write(row, col+2, timenum[i])
                row = row +1
            workbook.close()
            """
            pgram = signal.lombscargle(timenum, SNR2_vector_1, funkcia)  # vypocet periodogramu

            #freq, welch = signal.welch(SNR2_vector_1)
            ind = np.ndarray.argmax(pgram)
            print alt_prob[ind], funkcia[ind]
            plt.plot(alt_prob, pgram)
            plt.show()
            #doy_list.append(dat)
            #vyska.append(alt_prob[ind])
            """
            # that's the line, you need:
            #a = diff(sign(diff(pgram))).nonzero()[0] + 1  # local min+max
            #b = (diff(sign(diff(pgram))) > 0).nonzero()[0] + 1  # local min
            c = (diff(sign(diff(pgram))) < 0).nonzero()[0] + 1  # local max

            #print pgram[ind], (np.mean(pgram) * konstanta)

            if pgram[ind] > (np.mean(pgram)*2):
                if pgram[ind] > np.mean(pgram[c])*1.5:
                    if 3.9 <= alt_prob[ind] <= 4.6:
                        doy_list.append(dat)
                        vyska.append(alt_prob[ind])
                        print alt_prob[ind]

            plt.plot(alt_prob, pgram)
            plt.plot([2, 6], [(np.mean(pgram)) * 2, (np.mean(pgram) * 2)], 'r')
            plt.plot([2, 6], [np.mean(pgram[c])*1.5, np.mean(pgram[c])*1.5], 'g')
            plt.show()
            
            doy_list.append(doy_eph)
            vyska.append(alt_prob[ind])
            
            plt.subplot(3, 1, 1)
            plt.plot(x, SNR2)
            plt.xlabel('sin ')
            plt.ylabel('residuals (volt)')
            
            plt.xlabel('sinus elevation angle')
            plt.ylabel('residuals (volts)')
            #red_patch = mpatches.Patch(color='red', label='model')
            #blue_patch = mpatches.Patch(color='blue', label='residuals')
            plt.legend(handles=[blue_patch, red_patch]) 
            
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
            
            vyska.append(x)
            vyska.append(SNR2_vector_1)"""
           


plt.plot(doy_list, vyska, 'o')
#plt.plot([min(doy_list),4.4], [max(doy_list),4.4])
plt.show()
"""
plt.plot(vyska[0],vyska[1],vyska[2], vyska[3],'r')
plt.xlabel('sinus elevation angle')
plt.ylabel('residuals (volts)')
red_patch = mpatches.Patch(color='red', label='SNR2 obdobie zo snehom')
blue_patch = mpatches.Patch(color='blue', label='SNR2 obdobie bez snehu')
plt.legend(handles=[blue_patch, red_patch])
plt.show()"""



