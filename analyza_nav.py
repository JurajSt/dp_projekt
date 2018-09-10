import os, sys
import datetime, time
import numpy as np
from numpy import *
from input_data import *
import math as m
import model
from azi_ele_coor import *
import scipy.signal as signal
import matplotlib.pyplot as plt
from dataSHMU import *
import reader
import ogr
import t3
import matplotlib.patches as mpatches
import poloha_druzice

try:
    import psycopg2
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"

try:
    import matplotlib.pyplot as plt
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"

try:
    import xlsxwriter
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"

connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()

def fAnalyzaUlozenie():

    cursor.execute("DROP TABLE IF EXISTS VyberDruzice")
    cursor.execute("CREATE TABLE VyberDruzice (sat_cislo CHAR (4), st_nazov CHAR (25), azimut REAL, ele_uhol REAL, cas FLOAT, geom GEOMETRY)")
    print "DB pre vyber druzice vytvorena."

def fAnalyzaZapis(sat, sur_druzice, sur_stanice, st_nazov, ele_uhol, azimut, cas):

    line = ogr.Geometry(ogr.wkbLineString)
    if sursystem == 'xyz' or sursystem == 'blh':
        line.AddPoint(round(sur_stanice[0], 3), round(sur_stanice[1], 3), round(sur_stanice[2], 3))
        line.AddPoint(round(sur_druzice[0], 3), round(sur_druzice[1], 3), round(sur_druzice[2], 3))

    else:
        line.AddPoint(round(sur_stanice[0], 3), round(sur_stanice[1], 3))
        line.AddPoint(round(sur_druzice[0], 3), round(sur_druzice[1], 3))

    wkt_line = line.ExportToWkt()

    cursor.execute("INSERT INTO VyberDruzice (sat_cislo, st_nazov, azimut, ele_uhol, cas, geom)"
                   " VALUES ('%s', '%s', %s, %s, %s, ST_GeometryFromText('%s'))" % (sat, st_nazov, azimut, ele_uhol, cas,  wkt_line))

def fAnalyza():

    row = 1
    col = 0
    workbook = xlsxwriter.Workbook('D:/diplomka/vystupcsv/vystup.' + ulozenie)
    worksheet = workbook.add_worksheet()
    worksheet.write(0, 0, 'cislo druzice')
    worksheet.write(0, 1, 'pocet_vstkych_observacii')
    worksheet.write(0, 2, 'pocet_vybranych_observacii')

    obsfiles = np.sort(next(os.walk(obs_path))[2])

    sat_all = ['G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G11', 'G12', 'G13', 'G14', 'G15',
               'G16', 'G17', 'G18', 'G19', 'G20', 'G21', 'G22', 'G23', 'G24', 'G25', 'G26', 'G27', 'G28', 'G29',
               'G30', 'G31', 'G32']

    # datum zaciatku a konca vypoctu
    dat_start = datetime.date(rok_start, mes_start, den_start)
    dat_end = datetime.date(rok_konec, mes_konec, den_konec)

    c = 0
    while c < len(stanice):

        for f in obsfiles:

            if 'o' != f.split('.')[1][-1]:
                continue

            o = obs_path + f # cesta k SP3 suboru
            print
            print o

            r_obs = open(o, "r")
            lines = r_obs.readlines()
            r_obs.close()

            #sp3, date_e = utility.fReadSP3(e)   # nacitanie presnych efemerid
            header, version, headlines, headlength, obstimes, sats, numallsvs, numsvs, date_o = reader.scan(lines)
            prva_obs = header['TIME OF FIRST OBS'].split()
            dat = datetime.date(int(prva_obs[0]), int(prva_obs[1]), int(prva_obs[2]))  # datum prvej observacie

            # print int(date_e[0]), int(date_e[1]), int(date_e[2])

            if dat > dat_end: # kontrola ci sa datum zaznamu sa nachadza v zadanom intervale datumu
                break
            if dat_start > dat:
                continue

            rok = prva_obs[0][2:]
            navfile = nav_path + f.split('.')[0] + '.' + rok + 'n'
            if not os.path.exists(navfile):
                print 'ziadne navigacne data'
                continue
            # else: continue

            print navfile
            print dat

            #parsrovanie observacnych dat
            data, obstypes = reader.processBlocks(lines, header, obstimes, headlines, headlength, sats, numallsvs, sat_all)
            navdata = reader.rinexnav(navfile)
            navvalues = navdata.values  # navigacne data
            navtime = navdata['t'].to_pandas()  # vyber casu
            del lines

            max_obs = []

            for sat in sat_all:

                pocitadlo = 0
                ele_pocitadlo = 0


                for k, d in enumerate(data):
                    # print d
                    # print k, sats[k]

                    if sat not in sats[k]:
                        continue
                    # print "pocitam pre: " + sat, k

                    for i in d:  # jedna sekunda pre x druzic
                        # print i
                        sv = i[0]
                        obssek = i[1]
                        #print obssek
                        if sat not in sv:
                            continue

                        for j in range(len(navtime)):
                            navsv = int(navvalues[j][0])

                            if int(sat[1:]) != navsv:
                                continue

                            navdatetime = time.strptime(str(navtime[j]), "%Y-%m-%d %H:%M:%S")  # delenie casu
                            navtimesecond = datetime.datetime(navdatetime[0], navdatetime[1], navdatetime[2],
                                                              navdatetime[3], navdatetime[4], navdatetime[5])
                            dat_nav = datetime.date(int(navdatetime[0]), int( navdatetime[1]), int( navdatetime[2]))  # datum prvej observacie
                            navsek = (navdatetime[3]*60*60) + (navdatetime[4]*60) + navdatetime[5]
                            if dat != dat_nav:
                                continue
                            Dt = obssek-navsek
                            #print navsv, navsek, obssek
                            if Dt < 0:
                                break
                            if Dt >= platnos_spravy:    # platonst spravy 7200 sekund = 2 hodiny uveden v input_data
                                continue

                            sattelite = poloha_druzice.fvypocet_poloha(navvalues[j], Dt)  # vypocet polohy druzice

                            azi, ele = fComputeAziEle(stanice[c][0], [sattelite[0], sattelite[1], sattelite[2]])   # vypocet elevacneho uhla, azimutu

                            # zobrazenie podkladovych map je v WGS 84 Web Mercator (EPSG:3857) tj. suradnice su prepocitane tiez a
                            # tym musi byt aj azimut prepocitany
                            latst, lonst, hel = fXYZ_to_LatLonH(stanice[c][0][0], stanice[c][0][1], stanice[c][0][2])
                            surstx, sursty = fwgs4326towgs3857(latst, lonst)
                            latsv, lonsv, hel2 = fXYZ_to_LatLonH(sattelite[0], sattelite[1], sattelite[2])
                            sursvx, sursvy = fwgs4326towgs3857(latsv, lonsv)

                            azi2 = fCalculateAzimuth([surstx, sursty], [sursvx, sursvy])
                            pocitadlo = pocitadlo + 1

                            ele_min = stanice[c][1][0]
                            ele_max = stanice[c][1][1]
                            azi_min = stanice[c][2][0]
                            azi_max = stanice[c][2][1]

                            if ele_min <= ele <= ele_max and azi_min <= azi2 <= azi_max:
                                ele_pocitadlo = ele_pocitadlo+1

                            if sursystem == 'xy':
                                # cislo druzice, sur druzice, sur stanice,nazov stanice, ele uhol, azimut
                                fAnalyzaZapis(sat, [sursvx, sursvy], [surstx, sursty], stanice[c][5], ele, azi2, i[1])

                            if sursystem == 'xyz':
                                fAnalyzaZapis(sat, sattelite, stanice[c][0], stanice[c][5], ele, azi, i[1])

                            if sursystem == 'llh':
                                fAnalyzaZapis(sat, [latsv, lonsv, hel2], [latst, lonst, hel], stanice[c][5], ele, azi, i[1])

                #print 'pre druzicu cislo ' + sat + ' na '  + stanice[c][4] + 'je ' + str(pocitadlo) + ' observacii z toho ' + str(ele_pocitadlo) + ' v zadanom rozsahu azimutu a elev. uhla'
                worksheet.write(row, col, sat)
                worksheet.write(row, col + 1, pocitadlo)
                worksheet.write(row, col + 2, ele_pocitadlo)
                max_obs.append(ele_pocitadlo)
                row = row + 1
            ind = max_obs.index(max(max_obs))
            doporucena_druzica = sat_all[ind]
            print ' doporucena druzica je ' + doporucena_druzica + ' pre ' + stanice[c][4] + ' s ' + str(max(max_obs)) + ' observaciami'
        c = c+1
