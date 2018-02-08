import rinex_reader as reader
import sys, os
import numpy as np
import time, datetime
import psycopg2
from input_data import *
from CubicSplineInterpolation import *
import ogr
import azi_ele_coor as aec

connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()

cursor.execute("DROP TABLE IF EXISTS interp18644")
cursor.execute("CREATE TABLE interp18644 (id INT NOT NULL, svname INT, datetime BIGINT, azimuth REAL, elevangle REAL, "
               "s1 DOUBLE PRECISION, s2 DOUBLE PRECISION, geom_xyz GEOMETRY, PRIMARY KEY (ID))")

id = 0
#A = np.nan* np.empty([3600,3600,11])
c_eph = os.path.join('data/eph/1864/KAME2740.15o')
r_eph = open(c_eph, "r")
eph = r_eph.readlines()
#"""
header = []
# Capture header info
for i, l in enumerate(eph):
    header.append(l)
    if "END OF HEADER" in l:
        i += 1  # skip to data
        break
body = eph[i:]
value = []
j = 0
while j < len(body):

    value.append(body[j])
    if "MARKER NAME" in body[j]:
        value.remove(value[-1])
        value.remove(value[-1])
        print value[0]
        value = []
        for k, h in enumerate(header):
            value.insert(k, h)
        break
        obsdata = reader.rinexobs(value)
        if "MARKER NUMBER" in body[j+1]:
            j = j+6
            continue
        j = j+5
        continue
    j = j+1

#for k, h in enumerate(header):
#    value.insert(k, h)
obsdata = reader.rinexobs(value)
obsvalues = obsdata[0].values                           # observacne data
obstime = obsdata[0]['t'].to_pandas()                   # [0] index aby som dostal z tuple > arrayDataset, t > pre cas
obstype = obsdata[0]['type'].values.tolist()            # list s obs. signalmi
osv = obsdata[0]['sv'].values.tolist()                   # zoznam indexov maxialneho poctu druzic
indexS1 = obstype.index('S1')
indexS2 = obstype.index('S2')


path_eph = 'data/eph/1864/igr18644.sp3'
path_obs ='data/eph/1864/KAME2740.15o'

sp3 = utility.fReadSP3(path_eph)
#obs = utility.fReadObservationRinex(path_obs,['G'], 1)

obst = []
for date in obstime:
    obsdatetime = time.strptime(str(date), "%Y-%m-%d %H:%M:%S")
    obstimesecond2 = (obsdatetime[3] * 60 * 60) + (obsdatetime[4] * 60) + obsdatetime[5]  # cas obs. dat v sekundach
    obst.append(obssecond2)

for j in range(len(obst)):
    for k in range(len(osv)):
        # print str(obsvalues[i][k][0])
        if str(obsvalues[j][k][0]) == 'nan':  # kontrola prazdnych hodnot
            break
        obssv = int(obsvalues[j][k][0])
        if int(str(obssv)[0]) == 1 and int(str(obssv)[1:]) == int(sv[2:]):