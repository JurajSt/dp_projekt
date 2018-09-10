import os,sys
from datetime import date
import reader
from input_data import *
import ogr
from azi_ele_coor import *
import numpy as np

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
    if sursystem == 'xyz' or sursystem == 'llh':
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

    ephfiles = np.sort(next(os.walk(eph_path))[2])
    #obsfiles = np.sort(next(os.walk(obs_path))[2])

    dat_start = date(rok_start, mes_start, den_start)
    dat_end = date(rok_konec, mes_konec, den_konec)

    sat_all = ['G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G11', 'G12', 'G13', 'G14', 'G15',
               'G16', 'G17', 'G18', 'G19', 'G20', 'G21', 'G22', 'G23', 'G24', 'G25', 'G26', 'G27', 'G28', 'G29',
               'G30', 'G31', 'G32']

    fAnalyzaUlozenie()
    c = 0
    while c < len(stanice):
        print "ukladam data do db"
        for f in ephfiles:

            if 'sp3' != f.split('.')[1]:
                continue

            e = eph_path + f # cesta k SP3 suboru
            print
            print e

            sp3, date_e = reader.fReadSP3(e)   # nacitanie presnych efemerid
            doy_eph = reader.fdayofyear(int(date_e[0]), int(date_e[1]), int(date_e[2]))  # prepocet datumu na den v roku
            dat = date(int(date_e[0]), int(date_e[1]), int(date_e[2])) # datum prveho zaznamu sp3 suboru
            # print int(date_e[0]), int(date_e[1]), int(date_e[2])

            if dat > dat_end: # kontrola ci sa datum zaznamu sa nachadza v zadanom intervale datumu
                break
            if dat_start > dat:
                continue

            if 0 < doy_eph < 10:            # zistenie nazvu pre observacne data pomocou doy
                namefile = '00' + str(doy_eph) + '0'
            elif 10 <= doy_eph < 100:
                namefile = '0' + str(doy_eph) + '0'
            else:
                namefile = str(doy_eph) + '0'

            rok = date_e[0][2:]
            obsfile = obs_path + stanice[c][4] + namefile + '.' + rok + 'o'
            if not os.path.exists(obsfile):
                print 'ziadny observacny subor'
                continue
            # else: continue
            r_obs = open(obsfile, "r")
            lines = r_obs.readlines()
            r_obs.close()

            print obsfile
            print dat

            #parsrovanie observacnych dat
            header, version, headlines, headlength, obstimes, sats, numallsvs, numsvs, date_o = reader.scan(lines)
            data, obstypes = reader.processBlocks(lines, header, obstimes, headlines, headlength, sats, numallsvs, sat_all)

            del lines

            max_obs = []
            for sat in sat_all:
                pocitadlo = 0
                ele_pocitadlo = 0

                for k, d in enumerate(data):
                    #print d
                    # print k, sats[k]

                    if sat not in sats[k]:
                        continue
                    #print "pocitam pre: " + sat, k

                    for i in d:  # jedna sekunda pre x druzic
                        # print i
                        sv = i[0]
                        if sat not in sv:
                            continue

                        lg = reader.fInterpolateSatelliteXYZ(sp3, sv, i[1])  # vypocet lagrangeovej interpolacie
                        azi, ele = fComputeAziEle(stanice[c][0], [lg[0], lg[1], lg[2]])   # vypocet elevacneho uhla, azimutu

                        # zobrazenie podkladovych map je v WGS 84 Web Mercator (EPSG:3857) tj. suradnice su prepocitane tiez a
                        # tym musi byt aj azimut prepocitany
                        latst, lonst, hel = fXYZ_to_LatLonH(stanice[c][0][0], stanice[c][0][1], stanice[c][0][2])
                        surstx, sursty = fwgs4326towgs3857(latst, lonst)
                        latsv, lonsv, hel2 = fXYZ_to_LatLonH(lg[0], lg[1], lg[2])
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
                            fAnalyzaZapis(sat, lg, stanice[c][0], stanice[c][5], ele, azi, i[1])

                        if sursystem == 'llh':
                            fAnalyzaZapis(sat, [latsv, lonsv, hel2], [latst, lonst, hel], stanice[c][5], ele, azi, i[1])

                #print 'pre druzicu ' + sat + ' ' + str(pocitadlo) + ' observacii z toho ' + str(ele_pocitadlo) + ' v zadanom rozsahu azimutu a elev. uhla'
                worksheet.write(row, col, sat)
                worksheet.write(row, col + 1, pocitadlo)
                worksheet.write(row, col + 2, ele_pocitadlo)
                max_obs.append(ele_pocitadlo)
                row = row + 1
            ind = max_obs.index(max(max_obs))
            doporucena_druzica = sat_all[ind]
            print ' doporucena druzica ' + doporucena_druzica + ' s ' + str(max(max_obs)) + ' observaciami'
            workbook.close()
        c = c+1
