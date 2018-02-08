# -*- coding: utf-8 -*-
import os
import datetime, time
import numpy as np
import ogr
import utility
from input_data import *
import psycopg2

import Tkinter

connection = psycopg2.connect(DB)                   # pripojenie k DB
connection.autocommit = True
cursor = connection.cursor()

cursor.execute("DROP TABLE IF EXISTS interp18644")
cursor.execute("CREATE TABLE interp18644 (id INT NOT NULL, svname CHAR(3), datetime BIGINT, azimuth REAL, elevangle REAL, "
               "s1 DOUBLE PRECISION, s2 DOUBLE PRECISION, geom_xyz GEOMETRY, PRIMARY KEY (ID))")

path_eph = 'data/eph/1864/igr18644.sp3'
c_obs = os.path.join('data/eph/1864/KAME2740.15o')
r_eph = open(c_obs, "r")
lines = r_eph.readlines()

def _obstime(fol):
    year = int(fol[0])
    if 80 <= year <= 99:
        year += 1900
    elif year < 80:  # because we might pass in four-digit year
        year += 2000
    dat = datetime.datetime(year=year, month=int(fol[1]), day=int(fol[2]),
                    hour=int(fol[3]), minute=int(fol[4]),
                    second=int(float(fol[5])),
                    microsecond=int(float(fol[5]) % 1 * 100000))
    sek = time.mktime(dat.timetuple())
    sod = (int(fol[3]) * 60 * 60) + (int(fol[4]) * 60) + float(fol[5])
    return sod

def scan(L):
    header = {}
    # Capture header info
    for i, l in enumerate(L):
        if "END OF HEADER" in l:
            i += 1  # skip to data
            break
        if l[60:80].strip() not in header:  # Header label
            header[l[60:80].strip()] = l[:60]  # don't strip for fixed-width parsers
            # string with info
        else:
            header[l[60:80].strip()] += " " + l[:60]
            # concatenate to the existing string

    verRinex = float(header['RINEX VERSION / TYPE'][:9])  # %9.2f
    # list with x,y,z cartesian
    header['APPROX POSITION XYZ'] = [float(j) for j in header['APPROX POSITION XYZ'].split()]
    # observation types
    header['# / TYPES OF OBSERV'] = header['# / TYPES OF OBSERV'].split()
    #print header['# / TYPES OF OBSERV']
    # turn into int number of observations
    header['# / TYPES OF OBSERV'][0] = int(header['# / TYPES OF OBSERV'][0])
    header['INTERVAL'] = float(header['INTERVAL'][:10])

    headlines = []
    headlength = []
    rowlength = []
    obstimes = []

    numallsvs = []
    sats = []
    svset = set()
    # %%

    row = 1 + (header['# / TYPES OF OBSERV'][0] - 1) // 5  # number of rows observed values
    #print row, len(L), len(L)/row
    while i < len(L):
        if len(L[i][0:26].split()) == 6:  # then its headerline
            if int(L[i][28]) in (0, 1, 5, 6):  # CHECK EPOCH FLAG STATUS
                headlines.append(i)
                year, month, day, hour = L[i][1:3], L[i][4:6], L[i][7:9], L[i][10:12]
                minute, second = L[i][13:15], L[i][16:26]
                obstimes.append(_obstime([year, month,
                                          day, hour,
                                          minute, second]))
                # ONLY GPS SATELLITES
                numsvs = int(L[i][29:32])  # Number of visible satellites %i3
                numallsvs.append(numsvs)
                headlength.append(1 + (numsvs - 1) // 12)  # number of lines in header, depends on how many svs on view
                row = 1 + (header['# / TYPES OF OBSERV'][0] - 1) // 5   # number of rows observed values
                rowlength.append(row)
                if numsvs > 12:
                    sv = []
                    for s in range(numsvs):
                        if s > 0 and s % 12 == 0:
                            i += 1  # every 12th sat  will add new headline row ex >12 2 rows
                        if L[i][33 + (s % 12) * 3 - 1] == 'G' or 'R' or 'E' or 'S' or 'J' or 'C':
                            sv.append(L[i][32 + (s % 12) * 3:35 + (s % 12) * 3].replace(' ', '0'))
                    sats.append(sv)
                    i += numsvs*row + 1

                else:
                    sats.append([L[i][32 + s * 3:35 + s * 3].replace(' ', '0') for s in range(numsvs) if
                                 L[i][33 + s * 3 - 1] == 'G' or 'R' or 'E' or 'S' or 'J' or 'C'])  # lista de satelites (numeros prn)

                    i += numsvs*row + 1

            else:  # there was a comment or some header info
                flag = int(L[i][28])
                if (flag != 4):
                    print(flag)
                skip = int(L[i][30:32])
                i += skip + 1
                # %% get every SV that appears at any time in the file, for master index


    return header, verRinex, headlines, headlength, obstimes, sats, numallsvs, numsvs

def processBlocks(lines, header, obstimes, ihead, headlength, sats, numallsvs, numsvs):
    obstypes = header['# / TYPES OF OBSERV'][1:]
    data = []
    for i in range(len(ihead)):
        linesinblock = numallsvs[i] * int(np.ceil(header['# / TYPES OF OBSERV'][0] / 5.))
        block = ''.join((lines[ihead[i] + headlength[i]:ihead[i] + linesinblock + headlength[i]]))# nsats x observations
        bdf = _block2df(block, obstypes, sats[i], len(sats[i]))
        for t in bdf:
            t.insert(1, obstimes[i])
        #data.append(bdf)
        #print data
        #break
    return bdf, obstypes

def _block2df(block, obstypes, svnames, svnum):
    """
    input: block of text corresponding to one time increment INTERVAL of RINEX file
    output: 2-D array of float64 data from block. Future: consider whether best to use Numpy, Pandas, or Xray.
    """
    assert isinstance(svnum, int)
    N = len(obstypes)
    e = 0
    text = ""
    data2 = []
    if N%5>0:
        row = 1+(N//5)                  #pocet riadkov
    else:
        row = N/5

    block = block.split("\n")           # odstranujem zlomy stran
    for i in range(len(block)):         # vytvorenie zo zaznamu jeden riadok
        while len(block[i]) < 80:
            block[i] = block[i] + " "
        if i == (row + e)-1:
            text = text + block[i] + "\n"
            e = e + row
        else:
            text = text + block[i]

    text = text.split("\n")

    i,j = 0,0
    while i < (len(text)-1): # cyklus na vytvorenie chybajucich hodnôt (nahradi nulov)
        txt = ""
        txt = txt + svnames[i] + " "  #.replace("G","1").replace("R","3").replace("E","5").replace("S","7").replace("C","9").replace("J","21")+ " "
        for j in range(len(text[i])/16):
            val = text[i][(16)*j:16*(j+1)].split()    # odsranuje prazdne miesta a hodnoty sily signalu
            if len(val) == 0:
                txt = txt + "0.0" + " "
            elif len(val) > 1:
                txt = txt+val[0]+" "
            else:
                txt = txt + val[0] + " "
            if j == N-1:
                txt = txt+"\n"
                break
        i = i+1
        data = txt.split()
        data2.append(data)
        data = ""
    return data2

header, version, headlines, headlength, obstimes, sats, numallsvs, numsvs = scan(lines)
data, obstypes = processBlocks(lines, header, obstimes, headlines, headlength, sats, numallsvs, numsvs)
sp3 = utility.fReadSP3(path_eph)
id = 0

indexS1 = obstypes.index('S1')
indexS2 = obstypes.index('S2')
point = ogr.Geometry(ogr.wkbPoint)
for d in range(len(data)):
    if 'G' != data[d][0][0]:
        continue
    for i in range(len(obstimes)):
        if data[d][0] in sats[i]:
            lg = utility.fInterpolateSatelliteXYZ(sp3, data[d][0], i)
            print data[d][0], i, id
            azi, ele = utility.fComputeAziEle(kame, [lg[0], lg[1], lg[2]])
            point.AddPoint(lg[0], lg[1], lg[2])
            wkt = point.ExportToWkt()
            cursor.execute('INSERT INTO interp18644 (id, svname, datetime, azimuth, elevangle, s1, s2, geom_xyz) '
                           'VALUES (%s, %s, %s, %s, %s, %s, %s, ST_GeometryFromText(%s))',
                           (id, data[d][0], obstimes[i], azi, ele, indexS1, indexS2, wkt))
            id = id+1
        else:
            continue















