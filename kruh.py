import psycopg2
from input_data import *
import math as m
import ogr
import azi_ele_coor as aec
import t3


def fwktpoint2xyz(wkt):
    point = wkt.replace(')','').split('(')
    coor = point[1].split()
    coor = map(float, coor)
    return coor

def kruh (sur):
    krok = 10 # krok v stupnoch po kruznici
    uhol = 0
    while uhol < 360:

        if uhol in range(0,91):
            a = (m.sin(m.radians(uhol))*r)/m.sin(m.radians(90))
            b = (m.sin(m.radians(90-uhol))*r)/m.sin(m.radians(90))
        elif uhol in range(92,181):
            a = (m.sin(m.radians(180-uhol)) * r) / m.sin(m.radians(90))
            b = -1*((m.sin(m.radians(uhol-90)) * r) / m.sin(m.radians(90)))
        elif uhol in range(182,271):
            a = -1*((m.sin(m.radians(uhol-180)) * r) / m.sin(m.radians(90)))
            b = -1*((m.sin(m.radians(90-(uhol-180))) * r) / m.sin(m.radians(90)))
        elif uhol in range(272,359):
            a = -1*((m.sin(m.radians(uhol-270)) * r) / m.sin(m.radians(90)))
            b = (m.sin(m.radians(90-(uhol-270))) * r) / m.sin(m.radians(90))
        x = sur[0] + b    #stationcoor[0] + b      #
        y = sur[1] + a    #stationcoor[1] + a      #
        xx = stationcoor[0] + b      #
        yy = stationcoor[1] + a      #
        azimuth = aec.fCalculateAzimuth([sur[0], sur[1]], [x, y])
        #blh = aec.fXYZ_to_LatLonH(stationcoor[0]+b, stationcoor[1]+a, stationcoor[2])
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x, y)
        wkt_kruh = point.ExportToWkt()
        cursor.execute("INSERT INTO kruh (stname, azimuth, geom) VALUES ('%s', %s, ST_GeometryFromText('%s'))" %
                       (stname, azimuth, wkt_kruh))

        uhol = uhol + krok


connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()

r = 100 # polomer kruznice v metroch

# vytvorenie tabulky pre body azimutu v epsg:3857 wgs84 web mercator
cursor.execute('DROP TABLE IF EXISTS kruh')
cursor.execute('CREATE TABLE  kruh (id INT NOT NULL, stname CHAR (4), azimuth REAL, azimuth1 REAL, geom GEOMETRY, geom_xyz GEOMETRY, PRIMARY KEY (id))')

if sursystem == 'xy':
    # vytvorenie tabulky pre stanice v epsg:3857 wgs84 web mercator
    cursor.execute('DROP TABLE IF EXISTS stanice_webmercator')
    cursor.execute('CREATE TABLE  stanice_webmercator (stname CHAR (4) NOT NULL, geom GEOMETRY, PRIMARY KEY (stname))')

if  0 < abs(stanice[c][0]) < 180 and 0 < abs(stanice[c][1]) < 180 or sursystem == 'xyz':
    # vytvorenie tabulky pre stanice vo wgs84
    cursor.execute('DROP TABLE IF EXISTS stanice_xyz')
    cursor.execute('CREATE TABLE  stanice_xyz (stname CHAR (4) NOT NULL, geom GEOMETRY, PRIMARY KEY (stname))')
else:
    # vytvorenie tabulky pre stanice vo wgs84
    cursor.execute('DROP TABLE IF EXISTS stanice_blh')
    cursor.execute('CREATE TABLE  stanice_blh (stname CHAR (4) NOT NULL, geom GEOMETRY, PRIMARY KEY (stname))')

cursor.execute('SELECT stname, ST_AsText(geom) FROM station')
stationrows = cursor.fetchall()

id = 0

for s in stanice:

    stname = telg[4]
    stationcoor = s[0]
    point = ogr.Geometry(ogr.wkbPoint)
    if sursystem == 'xy':
        if 0 < abs(stanice[c][0]) < 180 and 0 < abs(stanice[c][1]) < 180:
            blh = aec.fXYZ_to_LatLonH(stationcoor[0], stationcoor[1], stationcoor[2])
            coor2 = aec.fwgs4326towgs3857(blh[0], blh[1])
        else:
            coor2 = aec.fwgs4326towgs3857(s[0], s[1])
        point.AddPoint(coor2[0], coor2[1])
        wkt = point.ExportToWkt()
        cursor.execute("INSERT INTO stanice_webmercator (stname, geom) VALUES ('%s', ST_GeometryFromText('%s'))"% (stname, wkt))
        kruh(coor2)

    elif sursystem == 'xyz':
        if 0 < abs(stanice[c][0]) < 180 and 0 < abs(stanice[c][1]) < 90:
            point.AddPoint(s[0], s[1], s[2])
            kruh(s)
        else:
            xyz = aec.fLatLonH_to_XYZ(s[0], s[1], s[2])
            point.AddPoint(xyz[0], xyz[1], xyz[2])
        wkt = point.ExportToWkt()
        cursor.execute("INSERT INTO stanice_xyz (stname, geom) VALUES ('%s', ST_GeometryFromText('%s'))" % (stname, wkt))
        kruh(xyz)

    else:
        if 0 < abs(stanice[c][0]) < 180 and 0 < abs(stanice[c][1]) < 90:
            blh = aec.fXYZ_to_LatLonH(stationcoor[0], stationcoor[1], stationcoor[2])
            point.AddPoint(blh[0], blh[1], blh[2])
            kruh(blh)
        else:
            point.AddPoint(blh[0], blh[1], blh[2])
            kruh(s)
        wkt = point.ExportToWkt()
        cursor.execute("INSERT INTO stanice_blh (stname, geom) VALUES ('%s', ST_GeometryFromText('%s'))" % (stname, wkt))


