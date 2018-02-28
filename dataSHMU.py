# -*- coding: utf-8 -*-
#citanie excelu zo snehovou pokrivkou od SHMU

from input_data import *
from datetime import date
import matplotlib.pyplot as plt
import math as m
import numpy as np
import csv
import os, sys

Lies = [[],[]]
Telg = [[],[]]
Kame = [[],[]]

csvpath = 'data/SHMU/KAME_1617.csv'
vstup = open(csvpath, "r")    # r - read, w - write, a - updaate
subor = vstup.readlines()
vstup.close()

hlavicka = subor[0]
print hlavicka
for data in subor[1:]:
    data1 = data.strip('\n')
    datas = data1.split(';')
    #print datas[2], len(datas[2])
    if len(datas[2]) == 0:
        datas[2] = 0
    if int(datas[2]) == 999 or int(datas[2]) == 995:
        datas[2] = 0
    if datas[0] == 'Liesek':
        dat = datas[1].split('.')
        d = date(int(dat[2]), int(dat[1]), int(dat[0]))
        Lies[0].append(d)
        Lies[1].append(int(datas[2]))
    elif datas[0] == 'Telgart':
        dat = datas[1].split('.')
        d = date(int(dat[2]), int(dat[1]), int(dat[0]))
        Telg[0].append(d)
        Telg[1].append(int(datas[2]))
    elif datas[0] == 'Kamenica nad Cirochou':
        dat = datas[1].split('.')
        d = date(int(dat[2]), int(dat[1]), int(dat[0]))
        Kame[0].append(d)
        Kame[1].append(int(datas[2]))

data = np.array(Kame[1])
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
histo = np.histogram(data)
plt.plot(Kame[0],Kame[1])
plt.show()

plt.hist(data)
plt.show()