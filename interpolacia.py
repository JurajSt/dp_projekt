import rinex_reader as reader
import poloha_druzice
import CubicSplineInterpolation
from LagrangeInterpolation import *
from CubicSplineInterpolation import *
import numpy as np
import matplotlib.pyplot as plt
import read_eph
import psycopg2
from input_data import *

connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()

c_eph = 'data/eph/1864/igr18644.sp3'
rEph = read_eph.fEphRead(c_eph)

#cursor.execute('SELECT svname, COUNT(svname) FROM igr18644 GROUP BY svname')
sv = 'PG01'
cursor.execute("SELECT coorx, coory,coorz, datetime FROM igr18644 WHERE svname = '%s' " % sv)
rows = cursor.fetchall()
x, y, z, t = [],[],[],[]
for row in rows:
    x.append(row[0])
    y.append(row[1])
    z.append(row[2])
    t.append(row[3])
    #print row

cs = fCubicSplineInterpolation(x,y,z,t)

plt.plot(x, y, 'o', label='data')
plt.plot(cs[0], cs[1], 'r', label='data')
plt.show()
