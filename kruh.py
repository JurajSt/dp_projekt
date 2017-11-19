import psycopg2
from input_data import *
import math as m
import ogr
import azi_ele_coor as aec


def fwktpoint2xyz(wkt):
    point = wkt.replace(')','').split('(')
    coor = point[1].split()
    coor = map(float, coor)
    return coor

connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()

r = 100 # polomer kruznice v metroch

# vytvorenie tabulky pre body azimutu v epsg:3857 wgs84 web mercator
cursor.execute('DROP TABLE IF EXISTS kruh')
cursor.execute('CREATE TABLE  kruh (id INT NOT NULL, stname CHAR (4), azimuth REAL, azimuth1 REAL, geom GEOMETRY, PRIMARY KEY (id))')
# vytvorenie tabulky pre liniu kruhu v epsg:3857 wgs84 web mercator
cursor.execute('DROP TABLE IF EXISTS kruh_line')
cursor.execute('CREATE TABLE  kruh_line (id INT NOT NULL, stname CHAR (4), geom GEOMETRY, PRIMARY KEY (id))')
# vytvorenie tabulky pre stanice v epsg:3857 wgs84 web mercator
cursor.execute('DROP TABLE IF EXISTS station_webmercator')
cursor.execute('CREATE TABLE  station_webmercator (stname CHAR (4) NOT NULL, geom GEOMETRY, PRIMARY KEY (stname))')

cursor.execute('SELECT stname, ST_AsText(geom) FROM station')
stationrows = cursor.fetchall()

id = 0

for row in stationrows:

    stname = row[0]
    stationcoor = fwktpoint2xyz(row[1])
    blh = aec.fXYZ_to_LatLonH(stationcoor[0], stationcoor[1], stationcoor[2])
    coor2 = aec.fwgs4326towgs3857(blh[0],blh[1])
    point = ogr.Geometry(ogr.wkbPoint)
    point.AddPoint(coor2[0], coor2[1])
    wkt = point.ExportToWkt()
    cursor.execute('INSERT INTO station_webmercator (stname, geom) VALUES (%s, ST_GeometryFromText(%s))', (stname, wkt))
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
        x = coor2[0] + b    #stationcoor[0] + b      #
        y = coor2[1] + a    #stationcoor[1] + a      #
        azimuth = aec.fCalculateAzimuth([coor2[0], coor2[1]], [x, y])
        angle = aec.fComputeAziEle(stationcoor, [stationcoor[0]+b, stationcoor[1]+a, stationcoor[2]])
        #blh = aec.fXYZ_to_LatLonH(stationcoor[0]+b, stationcoor[1]+a, stationcoor[2])
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x, y)
        wkt = point.ExportToWkt()
        cursor.execute('INSERT INTO kruh (id, stname, azimuth, azimuth1, geom) VALUES (%s, %s, %s, %s, ST_GeometryFromText(%s))',
                       (id, stname, azimuth, angle[0], wkt))
        id = id +1
        uhol = uhol + krok
id = 0
"""
for strow in stationrows:
    stname = strow[0]
    cursor.execute("SELECT ST_AsText(geom) FROM kruh WHERE stname = '%s' ORDER BY azimuth ASC" % stname)
    rows = cursor.fetchall()
    line = ogr.Geometry(ogr.wkbLineString)
    for krow in rows:
        coor = fwktpoint2xyz(krow[0])
        line.AddPoint(coor[0], coor[1])
        wkt_line = line.ExportToWkt()
    coor = fwktpoint2xyz(rows[0][0])
    line.AddPoint(coor[0], coor[1])
    wkt_line = line.ExportToWkt()
    cursor.execute('INSERT INTO kruh_line (id, stname, geom) VALUES (%s, %s, ST_GeometryFromText(%s))', (id, stname, wkt_line))
    id = id + 1

cursor.execute('SELECT ST_AsText(geom) FROM ganp3060_16')
satteliterows = cursor.fetchall()
cursor.execute('ALTER TABLE ganp3060_16 DROP COLUMN IF EXISTS geom_webmercator')
cursor.execute('ALTER TABLE ganp3060_16 ADD COLUMN geom_webmercator GEOMETRY')
print len(satteliterows)
i = 0
for row in satteliterows:
    print i
    point = ogr.Geometry(ogr.wkbPoint)
    coor = fwktpoint2xyz(row[0])
    point.AddPoint(coor[0], coor[1])
    wkt= point.ExportToWkt()
    cursor.execute("UPDATE ganp3060_16 SET geom_webmercator = ST_GeometryFromText('%s')" % wkt)
    i = i+1"""