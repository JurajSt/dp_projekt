# -*- coding: utf-8 -*-
import scipy.signal as signal
import matplotlib.pyplot as plt
from input_data import *
import psycopg2
import model
import numpy as np
import PreparingData

connection = psycopg2.connect(DB)  # pripoenie k db
connection.autocommit = True
cursor = connection.cursor()

def fPgram(tablename, azimin, azimax, elevmin, elevmax):
    print tablename
    # vyber druzice podla poctu zaznamov (zistujem len cislo druzice)
    cursor.execute('SELECT svname, COUNT(svname) FROM ' + tablename + ' WHERE azimuth > ' + str(azimin) + ''
                   'and azimuth < ' + str(azimax) + ' and elevangle > ' + str(elevmin) + ' and elevangle < '
                   + str(elevmax) + ' GROUP BY svname')

    count_rows = cursor.fetchall()

    for row in count_rows:  # kontrolujem pocet meranich zaznamov
        if row[1] < 70:  # pri 30s rinex datach minimalny pocet observacii
            continue
        sv = row[0]
        print sv

        # priprava dat - odtranenie opakujucich, urcenie klesajuci resp stupajucich druzic, kontrola minimalneho poctu
        # observacii
        SNR1, SNR2, sin, timenum, timetext, azimuth, cas = \
            PreparingData.fPreparingData(tablename, azimin, azimax, elevmin, elevmax, sv)

        deg = 5

        for j in range(len(timenum)):
            SNR1 = np.array(SNR1[j])   # zoznam na np.array
            SNR2 = np.array(SNR2[j])
            x = np.array(sin[j])
            den_sek= np.array(cas[j])  # cas dna merania v sekundach
            timenum = np.array(timenum[j])

            #modelSNR1, residuals1 = model.modelSNR(x, SNR1, deg)    # vypocet modelu
            modelSNR2, residuals2 = model.modelSNR(x, SNR2, deg)

            vyska_odrazec = []
            funkcia = []
            den_hodiny = []

            for i in range(len(cas)):
                den_hodiny.append(den_sek[i] / 60 / 60)

            f = 1
            while f <= 21:
                funkcia.append(f)
                f = f+1

            for i in range(0,len(funkcia)):
                vyska_odrazec.append(funkcia[i]*0.244)

            pgram = signal.lombscargle(timenum, residuals2, funkcia)  # vypocet periodogramu

            plt.subplot(2, 1, 1)
            plt.plot(x, residuals2)
            #plt.show()
            plt.subplot(2, 1, 2)
            plt.plot(vyska_odrazec, pgram)
            plt.show()

azimin = 170
azimax = 200
elevmin = 5
elevmax = 30
i = 2740

#t = fPgram('hofn2880', azimin, azimax, elevmin, elevmax)
while i <= 2880:
    tablename = 'hofn' + str(i)
    t = fPgram(tablename, azimin, azimax, elevmin, elevmax)
    i = i+10