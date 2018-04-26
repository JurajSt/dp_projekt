# -*- coding: utf-8 -*-
"""
Read navigation and observation rinex files
@author: pyrinex (modified)

"""
from pathlib2 import Path
import numpy as np
from datetime import datetime
from pandas import read_hdf
import xarray
from io import BytesIO
from time import time
import re
import sys


# %% Navigation file

def rinexnav(fn, ofn=None):
    """
    Reads RINEX 2.11 NAV files
    Michael Hirsch
    It may actually be faster to read the entire file via f.read() and then .split()
    and asarray().reshape() to the final result, but I did it frame by frame.
    http://gage14.upc.es/gLAB/HTML/GPS_Navigation_Rinex_v2.11.html
    """
    fn = Path(fn).expanduser()

    startcol = 3  # column where numerical data starts
    N = 7  # number of lines per record

    sv = [];
    epoch = [];
    raws = ''

    with fn.open('r') as f:
        """
        skip header, which has non-constant number of rows
        """
        while True:
            if 'END OF HEADER' in f.readline():
                break
        """
        now read data
        """
        for l in f:
            # format I2 http://gage.upc.edu/sites/default/files/gLAB/HTML/GPS_Navigation_Rinex_v2.11.html
            sv.append(int(l[:2]))
            # format I2
            year = int(l[3:5])
            if 80 <= year <= 99:
                year += 1900
            elif year < 80:  # good till year 2180
                year += 2000
            epoch.append(datetime(year=year,
                                  month=int(l[6:8]),
                                  day=int(l[9:11]),
                                  hour=int(l[12:14]),
                                  minute=int(l[15:17]),
                                  second=int(l[17:20]),  # python reads second and fraction in parts
                                  microsecond=int(l[21]) * 100000))


            """
            now get the data as one big long string per SV
            """
            raw = l[22:80]
            for _ in range(N):
                raw += f.readline()[startcol:80]
            # one line per SV
            # raws += raw + '\n'
            raws += raw + ' '

    raws = raws.replace('D-', 'Em')
    raws = raws.replace('D', 'E')
    raws = raws.replace('-', ' -')
    raws = raws.replace('m', '-')
    raws = re.sub(r'\n', r' ', raws)
    #print raws


    lista = [float(i) for i in raws.split(' ') if len(i) != 0]
    sat_info = np.array(lista)
    #print sat_info
    sat_info = sat_info.reshape(len(lista) / 29, 29)
    nav = xarray.DataArray(data=np.concatenate((np.atleast_2d(sv).T, sat_info), axis=1),
                           coords={'t': epoch,
                                   'data': ['sv', 'SVclockBias', 'SVclockDrift', 'SVclockDriftRate', 'IODE',
                                            'Crs', 'DeltaN', 'M0', 'Cuc', 'Eccentricity', 'Cus', 'sqrtA', 'TimeEph',
                                            'Cic', 'OMEGA', 'CIS', 'Io', 'Crc', 'omega', 'OMEGA DOT', 'IDOT',
                                            'CodesL2', 'GPSWeek', 'L2Pflag', 'SVacc', 'SVhealth', 'TGD', 'IODC',
                                            'TransTime', 'FitIntvl']},    #, ''
                           dims=['t', 'data'])

    if ofn:
        ofn = Path(ofn).expanduser()
        print('saving NAV data to', ofn)
        if ofn.is_file():
            wmode = 'a'
        else:
            wmode = 'w'
        nav.to_hdf(ofn, key='NAV', mode=wmode, complevel=6)

    return nav


# %% Observation File
def rinexobs(fn, ofn=None):
    """
    Program overviw:
    1) scan the whole file for the header and other information using scan(lines)
    2) each epoch is read and the information is put in a 4-D xarray.DataArray
    3)  rinexobs can also be sped up with if an h5 file is provided,
        also rinexobs can save the rinex file as an h5. The header will
        be returned only if specified.

    rinexobs() returns the data in a 4-D xarray.DataArray, [Parameter,Sat #,time,data/loss of lock/signal strength]?????
    """
    # open file, get header info, possibly speed up reading data with a premade h5 file
    if type(fn) is list:
        lines = fn
        header, version, headlines, headlength, obstimes, sats, svset, numallsvs, numsvs = scan(lines)
        #print('{} is a RINEX {} file, {} kB.'.format(fn, version, fn.stat().st_size // 1000))
        data = processBlocks(lines, header, obstimes, svset, headlines, headlength, sats, numallsvs, numsvs)
    else:
        fn = Path(fn).expanduser()
        with fn.open('r') as f:
            #tic = time()
            #lines = f.readlines()
            lines = f.read().splitlines(True)
            header, version, headlines, headlength, obstimes, sats, svset, numallsvs, numsvs = scan(lines)
            print('{} is a RINEX {} file, {} kB.'.format(fn, version, fn.stat().st_size // 1000))
            if fn.suffix == '.h5':
                data = read_hdf(fn, key='data')
            else:
                data = processBlocks(lines, header, obstimes, svset, headlines, headlength, sats, numallsvs, numsvs)

            #print("finished in {:.2f} seconds".format(time() - tic))

        # write an h5 file if specified
        if ofn:
            ofn = Path(ofn).expanduser()
            print('saving OBS data to', str(ofn))
            if ofn.is_file():
                wmode = 'a'
            else:
                wmode = 'w'
                # https://github.com/pandas-dev/pandas/issues/5444
            data.to_hdf(ofn, key='OBS', mode=wmode, complevel=6, format='table')

    return data, header


# this will scan the document for the header info and for the line on
# which each block starts
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
    a = 1               # v GANP sa nachadza po hodine komentar
    row = 1 + (header['# / TYPES OF OBSERV'][0] - 1) // 5  # number of rows observed values
    #print row, len(L), len(L)/row
    while i < len(L):
        if "COMMENT" in L[i+a]:
            a = 0
            i = i+1
            continue
        a = 1  # v GANP sa nachadza po hodine komentar
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
                        if L[i][33 + (s % 12) * 3 - 1] == 'G' or 'R' or 'E' or 'S':
                            sv.append(L[i][32 + (s % 12) * 3:35 + (s % 12) * 3].replace(' ', '0'))
                    sats.append(sv)
                    i += numsvs*row + 1

                else:
                    sats.append([L[i][32 + s * 3:35 + s * 3].replace(' ', '0') for s in range(numsvs) if
                                 L[i][33 + s * 3 - 1] == 'G' or 'R' or 'E' or 'S'])  # lista de satelites (numeros prn)

                    i += numsvs*row + 1

            else:  # there was a comment or some header info
                flag = int(L[i][28])
                if (flag != 4):
                    print(flag)
                skip = int(L[i][30:32])
                i += skip + 1
                # %% get every SV that appears at any time in the file, for master index

    for sv in sats:
        svset = svset.union(set(sv))

    return header, verRinex, headlines, headlength, obstimes, sats, svset, numallsvs, numsvs


def processBlocks(lines, header, obstimes, svset, ihead, headlength, sats, numallsvs, numsvs):
    # lines,header,obstimes,svset,ihead, headlength,sats
    obstypes = header['# / TYPES OF OBSERV'][1:]
    blocks = np.nan * np.empty((len(obstimes),      # tvorba matice
                                max(numallsvs),
                                len(obstypes)+1))  # por que max

    for i in range(len(ihead)):
        linesinblock = numallsvs[i] * int(np.ceil(header['# / TYPES OF OBSERV'][0] / 5.))  # nsats x observations
        # / 5 there is space for 5 observables per line
        block = ''.join(lines[ihead[i] + headlength[i]:ihead[i] + linesinblock + headlength[i]])
        bdf = _block2df(block, obstypes, sats[i], len(sats[i]))

        for x in range(len(bdf[0])):                # naplnenie matice blocks
            for y in range(len(bdf[0][0])):
                blocks[i][x][y] = bdf[0][x][y]

    obstypes.insert(0, 'sv')                        # number of satellites
    blocks = xarray.DataArray(data=blocks,          # tvorba matice, kde sv = falosne cislo druzice > len index koli poctu
                              coords={'t': obstimes,
                                      'sv': np.arange(max(numallsvs)),
                                      'type': obstypes},                 #['sv', 'L1', 'L2', 'P1', 'P2', 'C1', 'S1', 'S2']
                              dims=['t','sv', 'type'])

    return blocks


def _obstime(fol):
    year = int(fol[0])
    if 80 <= year <= 99:
        year += 1900
    elif year < 80:  # because we might pass in four-digit year
        year += 2000
    return datetime(year=year, month=int(fol[1]), day=int(fol[2]),
                    hour=int(fol[3]), minute=int(fol[4]),
                    second=int(float(fol[5])),
                    microsecond=int(float(fol[5]) % 1 * 100000)
                    )


def _block2df(block, obstypes, svnames, svnum):
    """
    input: block of text corresponding to one time increment INTERVAL of RINEX file
    output: 2-D array of float64 data from block. Future: consider whether best to use Numpy, Pandas, or Xray.
    """
    assert isinstance(svnum, int)
    N = len(obstypes)
    e = 0
    text = ""

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
    txt = ""
    while i < (len(text)-1):                           # cyklus na vytvorenie chybajucich hodnôt (nahradi nulov)
        txt = txt + svnames[i].replace("G","1").replace("R","3").replace("E","5").replace("S","7").replace("C","9").replace("J","21")+ " "
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

    sio = BytesIO(txt.encode('ascii'))
    barr = np.genfromtxt(sio, delimiter=" ",dtype=str, autostrip=True) # delenie riadka na hodnoty podla rinex   16,1,15,1,14,1,15,1,15,1,15,1,15
    data = np.zeros([1,len(barr),len(barr[0])],dtype=float)  # vytvorenie prazdnej matice

    for x in range(len(data[0])):                   # naplnenie matice
        for y in range(len(data[0][0])):
            data[0][x][y] = round(float(barr[x][y]), 5)


    return data

nav_path = 'D:/diplomka/dp_projekt/data/hofn/navig/hofn2740.15n'
navdata = rinexnav(nav_path)