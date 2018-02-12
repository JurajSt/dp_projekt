# -*- coding: utf-8 -*-
import os
import datetime, time
import numpy as np
import ogr
import utility
from input_data import *
import psycopg2

#connection = psycopg2.connect(DB)                   # pripojenie k DB
#connection.autocommit = True
#cursor = connection.cursor()

#cursor.execute("DROP TABLE IF EXISTS interp1864")
#cursor.execute("DROP TABLE IF EXISTS interp1864_1894")
#cursor.execute("CREATE TABLE interp1864_1894 (id INT NOT NULL, svname CHAR(3), datetime BIGINT, geom_xyz GEOMETRY, PRIMARY KEY (ID))")

id = 0
path_dat = 'data/eph/'#1864/ #igr18644.sp3'
dirs = next(os.walk(path_dat))[1]

for dir in dirs:
    path = path_dat+dir
    files = next(os.walk(path))[2]

    for file in files:
        if file.split('.')[1] == 'sp3':
            path_file = path+'/'+file
        else:
            continue
        sp3 = utility.fReadSP3(path_file)
        print path_file

        for d in sp3[0]:
            sv = d[0]
            i = 0
            print sv, d[1][0]
            point = ogr.Geometry(ogr.wkbPoint)
            while i < 86400:
                t = d[1][0]+i
                lg = utility.fInterpolateSatelliteXYZ(sp3, sv, t)
                i = i+1
                point.AddPoint(lg[0], lg[1], lg[2])
                wkt = point.ExportToWkt()
                #cursor.execute('INSERT INTO interp1864_1894 (id, svname, datetime, geom_xyz) '
                #               'VALUES (%s, %s, %s, ST_GeometryFromText(%s))',
                #               (id, sv, d[1][0], wkt))
                id = id+1











