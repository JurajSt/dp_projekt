# -*- coding: utf-8 -*-

import os, sys
import datetime, time
import numpy as np
import xlsxwriter
import utility
from input_data import *
import math as m
import model
import scipy.signal as signal
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
    #sek = time.mktime(dat.timetuple())
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
    #svset = set()
    # %%

    row = 1 + (header['# / TYPES OF OBSERV'][0] - 1) // 5  # number of rows observed values
    #print row, len(L), len(L)/row
    while i < len(L):
        if len(L[i][0:26].split()) == 6:  # then its headerline
            if int(L[i][28]) in (0, 1, 5, 6):  # CHECK EPOCH FLAG STATUS
                headlines.append(i)
                year, month, day, hour = L[i][1:3], L[i][4:6], L[i][7:9], L[i][10:12]
                minute, second = L[i][13:15], L[i][16:26]
                date = [year, month, day, hour, minute, second]
                obstimes.append(_obstime([year, month,day, hour, minute, second]))
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


    return header, verRinex, headlines, headlength, obstimes, sats, numallsvs, numsvs, date

def processBlocks(lines, header, obstimes, ihead, headlength, sats, numallsvs, numsvs):
    obstypes = header['# / TYPES OF OBSERV'][1:]
    data = []
    for i in range(len(ihead)):
        linesinblock = numallsvs[i] * int(np.ceil(header['# / TYPES OF OBSERV'][0] / 5.))
        block = ''.join((lines[ihead[i] + headlength[i]:ihead[i] + linesinblock + headlength[i]]))# nsats x observations
        bdf = _block2df(block, obstypes, sats[i], len(sats[i]), satellites)
        for t in bdf:
            t.insert(1, obstimes[i])
        data.append(bdf)
        #print data
        #break
    #print("read rinex in {:.2f} seconds".format(time.time() - tic))
    return data, obstypes

def _block2df(block, obstypes, svnames, svnum, satellites):
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
    for sat in satellites:
        i,j = 0,0
        while i < (len(text)-1): # cyklus na vytvorenie chybajucich hodnôt (nahradi nulov)
            txt = ""
            if sat not in svnames[i]:
                i = i + 1
                continue
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

def fdayofyear(y, m, d):
    d0 = datetime.date(y, 1, 1)
    d1 = datetime.date(y, m, d)
    delta = (d1 - d0)
    return delta.days+1

obs_path = 'data/kame_2/'#'D:/diplomka/telg_2/'
eph_path = 'data/eph/test/eph/'#igr18644.sp3'
#path_obs = 'data/kame_1/KAME2740.15o'#'data/eph/1864/KAME2740.15o'

#obsfiles = next(os.walk(obs_path))[2]
ephfiles = next(os.walk(eph_path))[2]
#obsfiles = ['KAME2740.15o']
doy_list = []
vyska = []
a = 0.0
func= 1000
row = 1
col = 0
workbook = xlsxwriter.Workbook('Expenses01T.xlsx')
worksheet = workbook.add_worksheet()
worksheet.write(0, 0, 's2')
worksheet.write(0, 1, 'sin')
worksheet.write(0, 2, 't2')
for f in ephfiles:

    if 'sp3' != f.split('.')[1]:
        continue

    e = eph_path+f
    print e
    satellites = ['G21']#['G05', 'G07', 'G08', 'G15', 'G21', 'G22', 'G26', 'G28', 'G30']
    #satellites = ['G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G11', 'G12', 'G13', 'G14', 'G15', 'G16',
    #              'G17', 'G18', 'G19', 'G20', 'G21', 'G22', 'G23', 'G24', 'G25', 'G26', 'G27', 'G28', 'G29', 'G30', 'G31',
    #              'G32']

    sp3, date_e = utility.fReadSP3(e)
    doy_eph = fdayofyear(int(date_e[0]), int(date_e[1]), int(date_e[2]))
    if doy_eph > 0 and doy_eph < 10:
        namefile = '00' + str(doy_eph)+ '0'
    elif doy_eph >= 10 and doy_eph < 100:
        namefile = '0' + str(doy_eph) + '0'
    else: namefile = str(doy_eph) + '0'
    obsfile = obs_path + 'KAME' + namefile + '.17o'
    if not os.path.exists(obsfile):
        obsfile = obs_path + 'KAME' + namefile + '.16o'
    #else: continue
    r_obs = open(obsfile, "r")
    lines = r_obs.readlines()
    r_obs.close()

    header, version, headlines, headlength, obstimes, sats, numallsvs, numsvs, date_o = scan(lines)
    data, obstypes = processBlocks(lines, header, obstimes, headlines, headlength, sats, numallsvs, numsvs)

    del lines
    deg = 6
    indexS1 = obstypes.index('S1')+2
    indexS2 = obstypes.index('S2')+2

    #svs = []
    den_hodiny = []
    for sat in satellites:
        print ("pocitam pre: " + sat)
        t = []
        signals1 = []
        signals2 = []
        sinele = []

        for k, d in enumerate(data):
            #print d
            #print k, sats[k]

            if sat not in sats[k]:
                continue
            #print "pocitam pre: " + sat, k

            for i in d:             # jedna sekunda pre x druzic
                #print i
                sv = i[0]
                if sat not in sv:
                   continue
                #print sv
                #signals = (i[indexS1], i[indexS2])
                lg = utility.fInterpolateSatelliteXYZ(sp3, sv, i[1])
                azi, ele = utility.fComputeAziEle(kame, [lg[0], lg[1], lg[2]])

                if ele >= 5 and ele <= 30 and azi >= 60 and azi <= 100:
                    if float(i[indexS2]) == 0.0:
                        break
                        #if i[1] == a:
                    t.append(i[1]/60)
                    #dataD.append((sv, t))
                    #dataD.append(lg)
                    #dataD.append(signals)
                    sinele.append(m.sin(m.radians(ele)))
                    #signals1.append(m.pow(10,float(i[indexS1])/20))
                    signals2.append(m.pow(10,float(i[indexS2])/20))
                    #print ele, m.sin(ele), i[1]
                    break
        """
        maxim = max(sinele)
        print maxim
        n=1

        sinele4 = []
        a = 0
        for i in range(len(sinele)):
            if i < maxim and maxim < 1:
                c = sinele[i]*(-1)
            elif i > maxim and maxim < 1:
                c = sinele[i]
            if i == a:
                sinele4.append(c)


        print sinele[0], sinele[-1]
        print t[0], t[-1]
        #if len(sinele) < 2000:
        #    print "pre ", sat, " málo observácií - vypocet neprebehne"
        #    continue
        fig = plt.figure()
        ax1 = fig.add_axes((0.1, 0.2, 0.8, 0.6))  # create an Axes with some room below
        ax1.set_title('SNR GPS L2')
        ax1.plot(t, signals2)
        ax1.set_xlabel('time (h)')
        ax1.set_ylabel('SNR (dB)')
        ax1_ticks = ax1.get_xticks()
        diff_ax1_ticks = ax1_ticks[-1] - ax1_ticks[0]
        a = 0
        b = int(len(t) / 4)
        tt = []
        sinele4 = []
        for i in range(len(t)):
            if i == a:
                tt.append((t[i] / diff_ax1_ticks) - (ax1_ticks[0] / diff_ax1_ticks))
                sinele4.append(round(sinele[i],2))
                a = a+b
        tt.append((t[-1] / diff_ax1_ticks) - (ax1_ticks[0] / diff_ax1_ticks))
        sinele4.append(round(sinele[-1], 2))
        # create second Axes. Note the 0.0 height
        ax2 = fig.add_axes((0.1, 0.1, 0.8, 0.0))
        ax2.yaxis.set_visible(False)  #
        ax2.set_xlabel('elevation angle')
        ax2.set_xticks(tt)
        ax2.set_xticklabels(sinele4)
        #plt.show()"""


        #print max(sinele), len(signals1), len(signals2)
        ind_max = sinele.index(max(sinele))
        for l in range(0,2):
            if l == 0:
                s1 = signals1[:ind_max+1]
                s2 = signals2[:ind_max+1]
                sinele2 = sinele[:ind_max+1]
                t2 = t[:ind_max+1]
                print "pre " + sat + " počet observacií: " + str(len(sinele2))
                #print len(sinele2), sinele2[-1]
            else:
                s1 = signals1[ind_max+1:]
                s2 = signals2[ind_max+1:]
                sinele2 = sinele[ind_max+1:]
                t2 = t[ind_max+1:]
                #print len(sinele2), sinele2[0]
                print "pre " + sat + " počet observacií: " + str(len(sinele2))
            if len(sinele2) < 2000:
                print "pre ", sat, " málo observácií - vypocet neprebehne"
                continue

            doy_list.append(doy_eph)
            SNR1 = np.array(s1)  # zoznam na np.array
            SNR2 = np.array(s2)
            x = np.array(sinele2)
            #den_sek = np.array(cas[j])  # cas dna merania v sekundach
            timenum = np.array(t2)
            vyska_odrazec = []


            #print np.size(timenum)

            for i in range (np.size(timenum)):
                worksheet.write(row, col, SNR2[i])
                worksheet.write(row, col+1, x[i])
                worksheet.write(row, col+2, timenum[i])
                row = row +1
            workbook.close()
            alt_prob = np.linspace(1.0, 2.0, 100)
            funkcia = []
            """f = 1
            while f < 2:
                funkcia.append(f)
                f = f+0.24

            funkcia.append(2)"""
            #print alt_prob
            for i in range(len(alt_prob)):
                vyska_odrazec.append(((2/0.24)*alt_prob[i]))
            print vyska_odrazec
            #print funkcia
            modelSNR2, residuals2, SNR2_vector = model.modelSNR(timenum, SNR2, deg)
            modelSNR2_1, residuals2_1, SNR2_vector_1 = model.modelSNR(x, SNR2_vector, deg)

            pgram = signal.lombscargle(timenum, SNR2_vector_1, vyska_odrazec)  # vypocet periodogramu

            ind = np.ndarray.argmax(pgram)
            #print ind
            #svs.append(sat)
            #vyska.append(alt_prob[ind])
            print alt_prob[ind]
            #funkcia.append(x)
            #funkcia.append(SNR2_vector_1)

            plt.subplot(2, 1, 1)
            plt.plot(x, SNR2_vector_1, x, modelSNR2_1, 'r')
            plt.xlabel('sin ')
            plt.ylabel('residuals (volt)')
            plt.xlabel('sinus elevation angle')
            plt.ylabel('residuals (volts)')
            red_patch = mpatches.Patch(color='red', label='model')
            blue_patch = mpatches.Patch(color='blue', label='residuals')
            plt.legend(handles=[blue_patch, red_patch])

            """
            plt.subplot(2, 1, 2)
            plt.plot(x, SNR2_vector_1)
            plt.xlabel('sin ')
            plt.ylabel(' SNR2 (db)')
            # plt.show()
            
            plt.subplot(3, 1, 3)
            plt.plot(alt_prob, pgram)
            plt.xlabel('Reflector height (m)')
            plt.ylabel('Spectral Ampl. ')"""
            plt.show()
            vyska.append(x)
            vyska.append(SNR2_vector_1)
print len(vyska[0]), len(vyska[1]), len(vyska[2]), len(vyska[3])
print len(doy_list), len(vyska)
print len(funkcia)
plt.plot(vyska[0],vyska[1],vyska[2], vyska[3],'r')
plt.xlabel('sinus elevation angle')
plt.ylabel('residuals (volts)')
red_patch = mpatches.Patch(color='red', label='SNR2 obdobie zo snehom')
blue_patch = mpatches.Patch(color='blue', label='SNR2 obdobie bez snehu')
plt.legend(handles=[blue_patch, red_patch])
plt.show()



