#-*- coding: utf-8 -*-
from input_data import *
import psycopg2
import numpy as np
import datetime
import time
import os
import ogr

connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()

def fEphRead(path):

    cursor.execute("DROP TABLE IF EXISTS igr18644")
    cursor.execute("CREATE TABLE igr18644 (id INT NOT NULL, svname CHAR (4), datetime FLOAT (30), "
                   "coorX FLOAT (25), coorY FLOAT (25), coorZ FLOAT (25), geom_xyz GEOMETRY, PRIMARY KEY (ID))")

    c_eph = os.path.join(path)
    r_eph = open(c_eph, "r")
    eph = r_eph.readlines()

    sv = []
    epoch = []
    raws = ''
    #kml = simplekml.Kml()
    #skip header, which has non-constant number of rows
    c=0
    k=1
    #m = np.zeros([32, 86400/900, 4])
    for data in eph:
        if "*" in data[0]:
            eph = eph[c:]
            break
        c=c+1
    #now read data
    points = []
    t = []
    id = 1

    for data in eph:
        if "EOF" in data:
            break
        if "*" in data[0]:
            #date_row = data#.split()
            date = data.split()
            year = int(date[1])
            month = int(date[2])
            day = int(date[3])
            hour = int(date[4])
            minute = int(date[5])
            second = float(date[6])
            timesec = hour*60*60 + minute*60 +second
            t.append(timesec)
            date_row = str(hour)+":"+str(minute)+":"+str(second)
            dat = datetime.datetime(year, month, day, hour, minute, int(second))
            dat = time.mktime(dat.timetuple())
            t.append(dat)
            continue

        data_row = data.split()
        sat = data_row[0]
        if sat[:2] == 'PR':
            continue
        x = float(data_row[1])*1000
        y = float(data_row[2])*1000
        z = float(data_row[3])*1000
        clk = float(data_row[4])/1e6
        point = [sat, timesec, x, y, z, clk]
        points.append(point)
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x,y,z)
        wkt = point.ExportToWkt()

        cursor.execute('INSERT INTO igr18644 (id, svname, datetime, coorX, coorY, coorZ, geom_xyz) '
                       'VALUES (%s, %s, %s, %s, %s, %s, ST_GeometryFromText(%s))', (id, sat, dat, x, y, z, wkt))
        id = id+1

        #return [points, t]


