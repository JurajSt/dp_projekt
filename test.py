import psycopg2
import xlsxwriter as xls
import scipy.signal as signal
import matplotlib.pyplot as plt
from input_data import *
import model
import copy_srnestimate
import numpy as np
import math as m
import time
#import pyAstronomy

connection = psycopg2.connect(DB)  # pripoenie k db
connection.autocommit = True
cursor = connection.cursor()
SNR1 = []
SNR2 = []
sin = []
timenum = []
timetext = []
azimuth = []
deg = 5
cursor.execute('SELECT s1, s2, elevangle, azimuth, datetime FROM ganp3060 WHERE svname = 23 and azimuth > 80 and azimuth < 110 and '
               'elevangle > 5 and elevangle < 30')
rows = cursor.fetchall()
N = len(rows)
print rows
for row in rows:
    SNR1.append(m.pow(10, row[0] / 20))
    SNR2.append(m.pow(10, row[1] / 20))
    sin.append(m.sin(m.radians(row[2])))
    #sin.append(row[1])
    azimuth.append(row[3])
    timenum.append(row[4])
    timetext.append(time.asctime(time.localtime(row[4])))
SNR1 = np.array(SNR1)
SNR2 = np.array(SNR2)
x = np.array(sin)

normval = x.shape[0] # For normalization of the periodogram
modelSNR1, residuals1 = model.modelSNR(x, SNR1, deg)
modelSNR2, residuals2 = model.modelSNR(x, SNR2, deg)
vyska_odrazec = []
funkcia = []
f = 1
while f <= 20:
    funkcia.append(f)
    f = f+1

for i in range(0,len(funkcia)):
    vyska_odrazec.append(funkcia[i]*0.19)

pgram = signal.lombscargle(timenum, residuals1, funkcia, normalize=True)

plt.subplot(2, 1, 1)
plt.plot(x, residuals1, 'b+')
plt.subplot(2, 1, 2)
plt.plot(vyska_odrazec, pgram)
plt.show()
workbook = xls.Workbook("D:\diplomka\modelSNR_PG23All.xlsx")
worksheet = workbook.add_worksheet()
worksheet.write(0, 0, "SNR1")
worksheet.write(0, 1, "SNR1 MODEL")
worksheet.write(0, 2, "RESIDUAL SRN1")
worksheet.write(0, 3, "SNR2")
worksheet.write(0, 4, "SNR2 MODEL")
worksheet.write(0, 5, "RESIDUAL SRN2")
worksheet.write(0, 6, "AZIMUTH")
worksheet.write(0, 7, "SIN_ELEV")
worksheet.write(0, 8, "TIME")
for i in range(len(SNR1)):
    worksheet.write(i+1, 0, SNR1[i])
    worksheet.write(i + 1, 1, modelSNR1[i])
    worksheet.write(i + 1, 2, residuals1[i])
    worksheet.write(i + 1, 3, SNR2[i])
    worksheet.write(i + 1, 4, modelSNR2[i])
    worksheet.write(i + 1, 5, residuals2[i])
    worksheet.write(i + 1, 6, azimuth[i])
    worksheet.write(i + 1, 7, sin[i])
    worksheet.write(i + 1, 8, timetext[i])

workbook.close()
