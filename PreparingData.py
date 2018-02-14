import psycopg2
from input_data import *
import math as m
import time


connection = psycopg2.connect(DB)  # pripoenie k db
connection.autocommit = True
cursor = connection.cursor()

def fPreparingData(tablename, azimin, azimax, elevmin, elevmax, sv):
    #print 'SELECT s1, s2, elevangle, azimuth, datetime FROM '+ tablename +' WHERE svname = '+"'" + str(sv) +"'" + ' and azimuth > '+ str(azimin) +' and azimuth < '+ str(azimax) +' and elevangle > '+ str(elevmin) +' and elevangle < '+ str(elevmax)
    cursor.execute('SELECT s1, s2, elevangle, azimuth, datetime FROM '+ tablename +' WHERE svname = '+"'" + str(sv) +"'"
                   ' and azimuth > '+ str(azimin) +' and azimuth < '+ str(azimax) +' and elevangle > '+ str(elevmin) +
                   ' and elevangle < '+ str(elevmax))
    rows = cursor.fetchall()

    SNR1u = []
    SNR2u = []
    sinu = []
    timenumu = []
    timetextu = []
    azimuthu = []
    casu = []

    for row in rows:
        if row[1] == 0:
            continue
        try:
            timenumu.index(row[4])           # na odstranenie zdovjenych dat
        except ValueError as e:
            timetextu.append(time.asctime(time.localtime(row[4])))
            SNR1u.append(m.pow(10, row[0] / 20))  # snr data prepocitane  do linearnej skali
            SNR2u.append(m.pow(10, row[1] / 20))
            sinu.append(m.sin(m.radians(row[2])))  # sinus elevacneho uhla
            azimuthu.append(row[3])
            timenumu.append(row[4])
            t = time.asctime(time.localtime(row[4]))
            txt = t.split()
            tt = txt[3].split(':')
            s = int(tt[0]) * 60 * 60 + int(tt[1]) * 60 + float(tt[2])
            casu.append(s)      # cas dna v sekundach

    print SNR1u
    SNR1 = [[],[]]
    SNR2 = [[],[]]
    sin = [[],[]]
    timenum = [[],[]]
    timetext = [[],[]]
    azimuth = [[],[]]
    cas = [[],[]]
    c = 0

    for i in range(len(timenumu)):  # rozdelenie dat na stupania a klesanie druzice

        j = 1 + i

        while j <= (len(timenumu) - 1):

            if timenumu[j] - timenumu[i] == 30:
                c = 0
                SNR1[c].append(SNR1u[i])
                SNR2[c].append(SNR2u[i])
                sin[c].append(sinu[i])
                timenum[c].append(timenumu[i])
                timetext[c].append(timetextu[i])
                azimuth[c].append(azimuthu[i])
                cas[c].append(casu[i])
                break

            else:
                c = 1
                #print sinu[i] , sinu[a]
                SNR1[c].append(SNR1u[i])
                SNR2[c].append(SNR2u[i])
                sin[c].append(sinu[i])
                timenum[c].append(timenumu[i])
                timetext[c].append(timetextu[i])
                azimuth[c].append(azimuthu[i])
                cas[c].append(casu[i])
                #c = c + 1
                break

    SNR1[c].append(SNR1u[i])
    SNR2[c].append(SNR2u[i])
    sin[c].append(sinu[i])
    timenum[c].append(timenumu[i])
    timetext[c].append(timetextu[i])
    azimuth[c].append(azimuthu[i])
    cas[c].append(casu[i])

    i = 0
    while i < (len(timenum)):

        if len(timenum[i]) < 70:            # pri 30s rinex datach minimalny pocet observacii
            SNR1.remove(SNR1[i])
            SNR2.remove(SNR2[i])
            sin.remove(sin[i])
            timenum.remove(timenum[i])
            timetext.remove(timetext[i])
            azimuth.remove(azimuth[i])
            cas.remove(cas[i])
            i = 0
        else:
            i = i+1

    return SNR1, SNR2, sin, timenum, timetext, azimuth, cas