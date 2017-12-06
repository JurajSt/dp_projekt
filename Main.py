# -*- coding: utf-8 -*-

import psycopg2
import rinex_reader as reader
import poloha_druzice
import numpy as np
import time, datetime
from input_data import *
import ogr
import azi_ele_coor as aec
import os
#import test


tic = time.time()

connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()

id = 1
navigfiles = next(os.walk(rinexNavpath))[2]   #dir is your directory path as string
obsfiles =next(os.walk(rinexObspath))[2]
#print  len(navigfiles), len(obsfiles)


for nf in navigfiles:
    for of in obsfiles:
        if nf[:8] != of[:8]:
            continue

        navdata = reader.rinexnav(rinexNavpath+nf)             # citanie navigacnych sprav
        obsdata = reader.rinexobs(rinexObspath+of)
        navvalues = navdata.values                              # navigacne data
        obsvalues = obsdata[0].values                           # observacne data

        navtime = navdata['t'].to_pandas()                      # vyber casu
        obstime = obsdata[0]['t'].to_pandas()                   # [0] index aby som dostal z tuple > arrayDataset, t > pre cas
        obstype = obsdata[0]['type'].values.tolist()            # list s obs. signalmi
        sv = obsdata[0]['sv'].values.tolist()                   # zoznam indexov maxialneho poctu druzic

        indexL1 = obstype.index('L1')                           # index signalu
        indexL2 = obstype.index('L2')
        indexS1 = obstype.index('S1')
        indexS2 = obstype.index('S2')


        #print 'nav time',len(navtime),'nav data', len(navvalues), 'sv', len(sv)
        #print 'obs time',len(obstime),'obs data', len(obsvalues)
        cursor.execute("DROP TABLE IF EXISTS hofn_line")
        cursor.execute("CREATE TABLE hofn_line (id INT NOT NULL, svname CHAR (4), azimuth REAL, azimuth1 REAL, elevangle REAL, geom GEOMETRY, geom_xyz GEOMETRY, PRIMARY KEY (ID))")

        cursor.execute("DROP TABLE IF EXISTS " + nf[:8])
        cursor.execute("CREATE TABLE "+ nf[:8] +"(id INT NOT NULL, svname INT, datetime BIGINT, azimuth1 REAL, azimuth REAL, elevangle REAL,l1 "
                       "DOUBLE PRECISION,l2 DOUBLE PRECISION, s1 DOUBLE PRECISION, s2 DOUBLE PRECISION, geom_xyz GEOMETRY, "
                       "geom_webmercator GEOMETRY, geom_latlonh GEOMETRY, PRIMARY KEY (ID))")

        for i in range(len(navtime)):
            Dt = 0          # vynulovanie Dt
            navdatetime = time.strptime(str(navtime[i]), "%Y-%m-%d %H:%M:%S")  # delenie casu
            # print date[0], date[1], date[2], date[3], date[4], date[5]  # datum: Y, M, D cas: H, M, S
            # navtimesecond = (navdatetime[3]*60*60)+(navdatetime[4]*60)+navdatetime[5]  # cas navig.sprav v sekundach
            navtimesecond = datetime.datetime(navdatetime[0], navdatetime[1], navdatetime[2], navdatetime[3],
                                                navdatetime[4], navdatetime[5])
            navsecond = time.mktime(navtimesecond.timetuple())  # cas a datum v navig.sprav v sekundach
            navsv = int(navvalues[i][0])  # cislo druzice v navig. sprave
            #print 'nav sv', navsv, navdatetime

            for j in range(len(obstime)):
                #if Dt < 0:
                 #   break


                obsdatetime = time.strptime(str(obstime[j]), "%Y-%m-%d %H:%M:%S")
                #obstimesecond2 = (obsdatetime[3] * 60 * 60) + (obsdatetime[4] * 60) + obsdatetime[5]  # cas obs. dat v sekundach
                obstimesecond = datetime.datetime(obsdatetime[0], obsdatetime[1], obsdatetime[2], obsdatetime[3],
                                                  obsdatetime[4], obsdatetime[5])
                obssecond = time.mktime(obstimesecond.timetuple())  # cas a den v  obs.datach v sekundach

                for k in range(len(sv)):
                    # print str(obsvalues[i][k][0])
                    if str(obsvalues[j][k][0]) == 'nan':  # kontrola prazdnych hodnot
                        break

                    obssv = int(obsvalues[j][k][0])  # cislo druzice v obs. datach
                    #print "cislo druzice v obs. datach je: ", obssv
                    if int(str(obssv)[0]) == 1 and int(str(obssv)[1:]) == int(navsv):
                        obssv = int(str(obssv)[1:])
                        Dt = obssecond - navsecond                 # cas na vypocet druzice
                        #print "cislo druzice v obs. datach je: ", obssv, Dt
                        if Dt < 0:
                            break

                        if Dt < platnos_spravy:
                            #if id == 18266:  # :18746
                             #   print id
                            point = ogr.Geometry(ogr.wkbPoint)
                            #print navvalues[i]
                            sattelite = poloha_druzice.fvypocet_poloha(navvalues[i], Dt)       # vypocet polohy druzice
                            point.AddPoint(sattelite[0], sattelite[1], sattelite[2])
                            wkt_xyz = point.ExportToWkt()           # wgs84 xyz

                            latlonh = aec.fXYZ_to_LatLonH(sattelite[0], sattelite[1], sattelite[2])
                            point.AddPoint(latlonh[1], latlonh[0], latlonh[2])
                            wkt_latlonh = point.ExportToWkt()       # wgs84 lat lon h

                            coor2 = aec.fwgs4326towgs3857(latlonh[0], latlonh[1])
                            point.AddPoint(coor2[0], coor2[1])
                            wkt_webmercator = point.ExportToWkt()   # wgs84 web mercator xyz

                            angles = aec.fComputeAziEle(hofn, [sattelite[0],sattelite[1], sattelite[2]])
                            azimuth = aec.fCalculateAzimuth(hofn,[coor2[0], coor2[1]])
                            cursor.execute('INSERT INTO '+ nf[:8] +' (id, svname, datetime, azimuth1, azimuth, elevangle,' 
                                           'l1, l2, s1, s2, geom_xyz, geom_webmercator, geom_latlonh)' 
                                           'VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, ST_GeometryFromText(%s),' 
                                           'ST_GeometryFromText(%s), ST_GeometryFromText(%s))',
                                            (id, obssv, obssecond, azimuth, angles[0], angles[1], obsvalues[j][k][indexL1],
                                             obsvalues[j][k][indexL2],obsvalues[j][k][indexS1], obsvalues[j][k][indexS2],
                                             wkt_xyz, wkt_webmercator, wkt_latlonh))
                            line = ogr.Geometry(ogr.wkbLineString)
                            line.AddPoint(-1691824.83965334, 9417966.19245271)
                            line.AddPoint(coor2[0], coor2[1])
                            wkt_line = line.ExportToWkt()
                            line2 = ogr.Geometry(ogr.wkbLineString)
                            line2.AddPoint(2679690.298, -727951.336, 5722789.244)
                            line2.AddPoint(sattelite[0], sattelite[1])
                            wkt_line_xyz = line2.ExportToWkt()
                            cursor.execute('INSERT INTO hofn_line (id, svname, azimuth1, azimuth, elevangle, geom, geom_xyz) '
                                           'VALUES (%s, %s, %s, %s, %s, ST_GeometryFromText(%s), ST_GeometryFromText(%s))',
                                (id, obssv,  azimuth, angles[0], angles[1], wkt_line, wkt_line_xyz))
                            id = id+1
                            break
                        break
        print("finished in {:.2f} seconds".format(time.time() - tic))
        #t = test.fTest(nf[:8])
        break




print("finished in {:.2f} seconds".format(time.time() - tic))








