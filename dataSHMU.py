#citanie excelu zo snehovou pokrivkou od SHMU

from input_data import *
import csv
import os, sys

vstup = open(csvpath, "r")    # r - read, w - write, a - updaate
subor = vstup.readlines()
vstup.close()

hlavicka = subor[0]
print hlavicka
for data in subor[1:]:
    data = data.strip('\n')
    data = data.split('\t')
