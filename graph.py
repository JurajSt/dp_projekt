# -*- coding: utf-8 -*-
import psycopg2
import xlsxwriter as xls
import scipy.signal as signal
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from input_data import *
import model
import copy_srnestimate
import numpy as np
import math as m
import time

#import pyAstronomy
pi = m.pi
connection = psycopg2.connect(DB)  # pripoenie k db
connection.autocommit = True
cursor = connection.cursor()

tablename = "ganp3060"

cursor.execute('SELECT svname, COUNT(svname) FROM %s WHERE azimuth > 170 and azimuth < 200 and '
               'elevangle > 5 and elevangle < 30 GROUP BY svname' % tablename)
k = 1
count_rows = cursor.fetchall()
workbook = xls.Workbook('D:\diplomka\excel\_'+tablename+'.xlsx')
for row in count_rows:            # kontrolujem pocet meranich zaznamov
    if row[1] < 70:
        continue
    sv = row[0]
    print sv
    cursor.execute('SELECT s1, s2, elevangle, azimuth, datetime FROM '+ tablename +' WHERE svname = '+ str(9))#+' and '
                   #'azimuth > 170 and azimuth < 200 and elevangle > 5 and elevangle < 30')
    rows = cursor.fetchall()

    N = len(rows)
    SNR1 = []
    SNR2 = []
    angle = []
    timenum = []
    timetext = []
    azimuth = []
    sin = []

    deg = 5

    for row in rows:
        #SNR1.append(m.pow(10, row[0] / 20))
        #SNR2.append(m.pow(10, row[1] / 20))
        if row[1] == 0:
            continue
        SNR1.append(row[0])
        SNR2.append(row[1])
        angle.append(m.radians(row[2]))
        sin.append(row[1])
        azimuth.append(row[3])
        timenum.append(row[4])
        t = time.asctime(time.localtime(row[4]))
        txt = t.split()
        cas = txt[3].split(':')
        s = int(cas[0]) * 60 * 60 + int(cas[1]) * 60 + float(cas[2])
        timetext.append(s)
    SNR1 = np.array(SNR1)
    SNR2 = np.array(SNR2)
    x = np.array(sin)

    #print cas, s, timetext

    modelSNR1, residuals1 = model.modelSNR(timenum, SNR1, deg)
    modelSNR2, residuals2 = model.modelSNR(timetext, SNR2, deg)
    #print residuals2
    cas = []
    for i in range(len(timetext)):
        cas.append(timetext[i]/60/60)
    funkcia = []
    f = 1
    while f <= 21:
        funkcia.append(f)
        f = f + 1
    vyska_odrazec2 = []
    for i in range(0, len(funkcia)):
        vyska_odrazec2.append(funkcia[i] * 0.244)

    per = signal.periodogram(timetext, 30)
    # print per

    pgramN = signal.lombscargle(timetext, residuals2, funkcia, normalize=True)
    print pgramN
    #pgram = pgram.tolist()
    """pgramN = pgramN.tolist()
    plt.subplot(2, 1, 1)
    plt.plot(timenum, SNR2)
    plt.subplot(2, 1, 2)
    plt.plot(vyska_odrazec2, pgramN)
    #plt.show()"""

    #print cas
    #print SNR2

    plt.plot(cas, SNR2)
    plt.xlabel('cas(h)')
    plt.ylabel('SNR2 (db)')
    blue_patch = mpatches.Patch(color='blue', label='SNR2')
    plt.legend(handles=[blue_patch])
    plt.show()

    plt.plot(cas, SNR2, cas, modelSNR2, "r")
    plt.xlabel('cas (h)')
    plt.ylabel('SNR2 (db)')
    red_patch = mpatches.Patch(color='red', label='model')
    blue_patch = mpatches.Patch(color='blue', label='SNR2')
    plt.legend(handles=[red_patch, blue_patch])
    #plt.legend(handles=)
    plt.show()

    worksheet = workbook.add_worksheet(str(sv))
    worksheet.write(0, 0, "SNR1")
    worksheet.write(0, 1, "SNR1 MODEL")
    worksheet.write(0, 2, "RESIDUAL SRN1")
    worksheet.write(0, 3, "SNR2")
    worksheet.write(0, 4, "SNR2 MODEL")
    worksheet.write(0, 5, "RESIDUAL SRN2")
    worksheet.write(0, 6, "AZIMUTH")
    worksheet.write(0, 7, "SIN_ELEV")
    worksheet.write(0, 8, "TIME SEC")
    worksheet.write(0, 9, "TIME")

    for i in range(len(SNR1)):
        worksheet.write(i + 1, 0, SNR1[i])
        worksheet.write(i + 1, 1, modelSNR1[i])
        worksheet.write(i + 1, 2, residuals1[i])
        worksheet.write(i + 1, 3, SNR2[i])
        worksheet.write(i + 1, 4, modelSNR2[i])
        worksheet.write(i + 1, 5, residuals2[i])
        worksheet.write(i + 1, 6, azimuth[i])
        worksheet.write(i + 1, 7, angle[i])
        worksheet.write(i + 1, 8, timenum[i])
        worksheet.write(i + 1, 9, timetext[i])

workbook.close()
for i in range(len(azimuth)):
    if round(angle[i], 2) == 0.5:
        print angle[i]
        h = 2.9
        lambda1 = 19.0  # vlnova dlzka L1 v cm
        lambda2 = 0.244  # vlnova dlzka L2
        n=1
        d =n*lambda2/2
        R = h/m.tan(angle[i]) + (d/m.sin(angle[i]))/m.tan(angle[i])
        b = m.sqrt(((2*d*h)/m.sin(angle[i])) + (m.pow(d/m.sin(angle[i]),2)))
        a = b/m.sin(angle[i])
        k = 0
        st = m.radians(1)      # stupen na radiany
        ze = []
        fzx = []
        fzy = []
        while k <= 2*pi:
            xi = a * m.cos(k) + R
            yi = b * m.sin(k)
            xe = m.sin(azimuth[i])*xi - m.cos(azimuth[i])*yi
            ye = m.sin(azimuth[i])*yi + m.cos(azimuth[i])*xi
            k = k+st
            zze = [xe, ye]
            fzx.append(xe)
            fzy.append(ye)
            ze.append(zze)

        plt.plot(fzx, fzy)
        plt.show()