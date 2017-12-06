#-*- coding: utf-8 -*-

from pathlib2 import Path
import numpy as np
from datetime import datetime
from pandas import read_hdf
import xarray
from io import BytesIO
from time import time
import re
import os
import azi_ele_coor as aec
import simplekml
import xlsxwriter as xls
# %% eph file

c_eph = os.path.join('data/bordel/cod18221.eph')
r_eph = open(c_eph, "r")
eph = r_eph.readlines()

sv = [];
epoch = [];
raws = ''
kml = simplekml.Kml()
#skip header, which has non-constant number of rows
c=0
k=1
for data in eph:
    if "*" in data[0]:
        eph = eph[c:]
        break
    c=c+1
#now read data
workbook = xls.Workbook("D:\diplomka\eph.xls")
worksheet = workbook.add_worksheet()
worksheet.write(0, 0, "svname")
worksheet.write(0, 1, "lat")
worksheet.write(0, 2, "lon")
for data in eph:
    if "EOF" in data:
        break
    if "*" in data[0]:
        #date_row = data#.split()
        date = data.split()
        year = int(date[1])
        month = int(date[2])
        day = int(data[3])
        hour = int(date[4])
        minute = int(date[5])
        second = float(date[6])
        timesec = hour*60*60 + minute*60 +second
        date_row = str(hour)+":"+str(minute)+":"+str(second)
        continue
    data_row = data.split()[:4]
    sat = data_row[0]
    x = float(data_row[1])*100
    y = float(data_row[2])*100
    z = float(data_row[3])*100

    if sat[:2] == 'PG' or sat[:2] == 'PR':
        timest = str(year) + "-" + str(month) + "-" + str(day) + "T" + str(hour) + ":" + str(minute) + ":" + \
                   str(second) + "Z"
        latlonh = aec.fXYZ_to_LatLonH(x,y,z)
        xyz = aec.fwgs4326towgs3857(latlonh[0], latlonh[1])
        worksheet.write(k, 0, data_row[0])
        worksheet.write(k, 1,latlonh[0])
        worksheet.write(k, 2, latlonh[1])
        k=k+1
        pnt = kml.newpoint(name=sat[2:], description=sat[:2]+","+str(timesec), coords=[(latlonh[1], latlonh[0])])  # lon, lat optional height
        #pnt.timestamp.when = times
    else:
        continue

    a = (latlonh[1],latlonh[0])
    sv.append(a)


#lin = kml.newlinestring(name="Pathway", description="A pathway in Kirstenbosch",coords=sv)


#print kml.kml()

kml.save('testkml.kml')
workbook.close()