# -*- coding: utf-8 -*-
#citanie excelu zo snehovou pokrivkou od SHMU

from input_data import *
from datetime import date
import matplotlib.pyplot as plt
import math as m
import numpy as np
#import csv
import os, sys

# funkcia nacita csv s vyskou snehu.
# data musia byt strukturovane do troch stlpcov kde v prvom je nazov meteo stanice,v druhom  datum,v tretom vyska snehu
def fSHMU(c):

    meteo = [[],[]]

    vstup = open(csvpath, "r")
    subor = vstup.readlines()
    vstup.close()

    hlavicka = subor[0]
    #print hlavicka
    for data in subor[1:]:
        data1 = data.strip('\n')
        datas = data1.split(';')
        #print datas[2]
        if datas[2] == '\r':
            datas[2] = 0
        elif len(datas[2]) == 0:
            datas[2] = 0
        elif int(datas[2]) == 999 or int(datas[2]) == 995:
            datas[2] = 0
        if datas[0] == stanice[c][5]:
            dat = datas[1].split('.')
            d = date(int(dat[2]), int(dat[1]), int(dat[0]))
            meteo[0].append(d)
            meteo[1].append(int(datas[2]))

    data = np.array(meteo[1])
    suma = 0
    pocet = 0
    for d in data:
        if float(d) == 0.0:
            continue
        else:
            suma = suma+d
            pocet = pocet+1

    primer = suma / pocet
    med = np.median(data)

    return meteo
"""
a = fSHMU()

plt.plot(a[0],a[1])
plt.show()
"""