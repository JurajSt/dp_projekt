import rinex_reader as reader
import poloha_druzice
import numpy as np
import time
from input_data import *

'''
navdata = reader.rinexnav(rinexNavfilename)     # citanie navigacnych sprav
data = navdata.values                           # vyber hodnot
print len(data)

#poloha = poloha_druzice.fvypocet_poloha(data, 0)
#print poloha

navcas = navdata['t'].to_pandas()               # vyber casu
for i in range(len(navcas)):
    date = time.strptime(str(navcas[i]), "%Y-%m-%d %H:%M:%S")       # delenie casu
    #print date[0], date[1], date[2], date[3], date[4], date[5]  # datum: Y, M, D cas: H, M, S
    timesecond = (date[3]*60*60)+(date[4]*60)+date[5]
    if date <> 0.0:
        dodatok = 60-date[5]
'''
obsdata = reader.rinexobs(rinexObsfilename)
cas = obsdata[0]['t'].to_pandas()# [0] index aby som dostal z tuple > arrayDataset, t > pre cas
obs = obsdata[0]['type'].values.tolist()
sv = obsdata[0]['sv'].values.tolist()
data = obsdata[0].values
#print obsdata
ss = []
for i in range(len(cas)):
    date = time.strptime(str(cas[i]), "%Y-%m-%d %H:%M:%S")
   #print date[0], date[1], date[2], date[3], date[4], date[5]  #datum: Y, M, D cas: H, M, S
    c = [date[0], date[1], date[2], date[3], date[4], date[5]]
    ss.append(c)
print 'obs', obs
print 'sv', sv
print 'cas', ss
signal = 'L1'
obs_i = obs.index(signal)
for cas in range(len(ss)):
    for sat in range(len(sv)):
        print data[cas][sat][0],ss[cas], data[cas][sat][obs_i]



