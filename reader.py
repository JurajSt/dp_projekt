from datetime import datetime
from datetime import date
import numpy as np
import sys

try:
    import xarray
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"

try:
    from pathlib2 import Path
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"

try:
    import re
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"

try:
    import gzip
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def _obstime(fol):
    year = int(fol[0])
    if 80 <= year <= 99:
        year += 1900
    elif year < 80:  # because we might pass in four-digit year
        year += 2000
    dat = datetime(year=year, month=int(fol[1]), day=int(fol[2]),
                    hour=int(fol[3]), minute=int(fol[4]),
                    second=int(float(fol[5])),
                    microsecond=int(float(fol[5]) % 1 * 100000))
    #sek = time.mktime(dat.timetuple())
    sod = (int(fol[3]) * 60 * 60) + (int(fol[4]) * 60) + float(fol[5])
    return sod

#lagrange polynomial interpolation
#https://ilrs.cddis.eosdis.nasa.gov/data_and_products/dfpwg/pfsg/eph_interpolation.pdf
#https://gist.github.com/melpomene/2482930
#http://www.navipedia.net/index.php/Precise_GNSS_Satellite_Coordinates_Computation
def fLagrangeInterpolation(points):
    def P(x):
		total = 0
		n = len(points)
		for i in xrange(n):
			xi, yi = points[i]
			def g(i, n):
				tot_mul = 1
				for j in xrange(n):
					if i == j:
						continue
					xj, yj = points[j]
					tot_mul *= (x - xj) / float(xi - xj)
				return tot_mul
			total += yi * g(i, n)
		return total
    return P

#reads a list with SP3 file and returns interpolated XYZ coordinates of input
#satellite (given by its PRN code at given epoch (seconds of day)
def fInterpolateSatelliteXYZ(sp3_list, prn, epoch):
    sour_sat = []
##    print prn, epoch
    sour_x = None
    for satellite_x in sp3_list[0]:
        if satellite_x[0] == prn:
            for i in range (1, len(satellite_x)):
                if abs(satellite_x[i][0] - epoch) < 900:
                    index_epoch = i
                    break
            if index_epoch < 6:
                index_epoch = 6
            elif index_epoch > len(satellite_x)-1-4:
                index_epoch = len(satellite_x)-1-4
##            print index_epoch
##            print satellite_x[index_epoch-5:index_epoch+5]
            funkce_sour_x = fLagrangeInterpolation(satellite_x[index_epoch-5:index_epoch+5])
            sour_x = funkce_sour_x(epoch)
    if sour_x != None:
        for satellite_y in sp3_list[1]:
            if satellite_y[0] == prn:
                funkce_sour_y = fLagrangeInterpolation(satellite_y[index_epoch-5:index_epoch+5])
                sour_y = funkce_sour_y(epoch)
        for satellite_z in sp3_list[2]:
            if satellite_z[0] == prn:
                funkce_sour_z = fLagrangeInterpolation(satellite_z[index_epoch-5:index_epoch+5])
                sour_z = funkce_sour_z(epoch)
        sour_sat = [sour_x,sour_y,sour_z]
        return sour_sat
    else:
        return -1

#read an SP3 file and create a list of all coordinates in this format:
#[[X],[Y],[Z]] where [X] or [Y] or [Z] contains items of following structure:
#[PRN, (epoch, coordinate), (epoch, coordinate), ...]
def fReadSP3(sp3_input_path):
    if '.gz' in sp3_input_path or '.GZ' in sp3_input_path or '.z' in sp3_input_path or '.Z' in sp3_input_path or '.zip' in sp3_input_path or '.ZIP' in sp3_input_path:
        with gzip.open(sp3_input_path, 'r') as f:
            lines = f.readlines()
    else:
        cti = open(sp3_input_path, 'r')
        lines = cti.readlines()
        cti.close()
    date = lines[0][3:].split()[:6]
    #get the list of satellites
    satellite_PRN_list = []
    for i in range(0, len(lines)):
        if lines[i][0] == '*':
            j=i+1
            while '*' not in lines[j]:
                prn = lines[j].split()[0]
                satellite_PRN_list.append(prn)
                j=j+1
            break

    positions = [[],[],[]] #one list for each coordinate (x,y,z)
    for satellite in satellite_PRN_list:
        sat_x = []
        sat_y = []
        sat_z = []
        sat_x.append(satellite.replace('P',''))
        sat_y.append(satellite.replace('P',''))
        sat_z.append(satellite.replace('P',''))
        for i in range(0, len(lines)):
            if lines[i][0] == '*':
                j=i+1
                time_list = lines[i].split()
                time_seconds = int(time_list[4]) * 3600 + int(time_list[5]) * 60 + float(time_list[6])
                #obstimesecond = datetime.datetime(int(time_list[1]), int(time_list[2]), int(time_list[3]), int(time_list[4]),
                #                                  int(time_list[5]), int(float(time_list[6])))
                #time_seconds = time.mktime(obstimesecond.timetuple())
                while '*' not in lines[j]:
                    if satellite in lines[j]:
                        sat_x.append((int(time_seconds),float(lines[j].split()[1])*1000))
                        sat_y.append((int(time_seconds),float(lines[j].split()[2])*1000))
                        sat_z.append((int(time_seconds),float(lines[j].split()[3])*1000))
                        break
                    j=j+1
        positions[0].append(sat_x)
        positions[1].append(sat_y)
        positions[2].append(sat_z)
    return positions, date

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
            if 'COMMENT' in L[i]:
                i = i + 1
                continue
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

        elif len(L[i][0:26].split()) <> 6:
            i = i+1

    return header, verRinex, headlines, headlength, obstimes, sats, numallsvs, numsvs, date

def processBlocks(lines, header, obstimes, ihead, headlength, sats, numallsvs, satelittes):
    obstypes = header['# / TYPES OF OBSERV'][1:]
    data = []
    satt = satelittes
    #print satt
    for i in range(len(ihead)):

        linesinblock = numallsvs[i] * int(np.ceil(header['# / TYPES OF OBSERV'][0] / 5.))
        block = ''.join((lines[ihead[i] + headlength[i]:ihead[i] + linesinblock + headlength[i]]))# nsats x observations
        bdf = _block2df(block, obstypes, sats[i], len(sats[i]), satt)
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
        while i < (len(text)-1): # cyklus na vytvorenie chybajucich hodnot (nahradi nulov)
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

# %% Navigation file
def rinexnav(fn, ofn=None):
    """
    Reads RINEX 2.11 and 2.10 NAV files
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
    try:
        sat_info = sat_info.reshape(len(lista) / 29, 29)
    except ValueError as e:
        nav = rinexnav210(fn)
        return nav
    nav = xarray.DataArray(data=np.concatenate((np.atleast_2d(sv).T, sat_info), axis=1),
                           coords={'t': epoch,
                                   'data': ['sv', 'SVclockBias', 'SVclockDrift', 'SVclockDriftRate', 'IODE',
                                            'Crs', 'DeltaN', 'M0', 'Cuc', 'Eccentricity', 'Cus', 'sqrtA', 'TimeEph',
                                            'Cic', 'OMEGA', 'CIS', 'Io', 'Crc', 'omega', 'OMEGA DOT', 'IDOT',
                                            'CodesL2', 'GPSWeek', 'L2Pflag', 'SVacc', 'SVhealth', 'TGD', 'IODC',
                                            'TransTime', 'FitIntvl']},    #, ''
                           dims=['t', 'data'])

    return nav


def rinexnav210(fn, ofn=None):
    """
    Reads RINEX 2.11 and 2.10 NAV files
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
    # print raws


    lista = [float(i) for i in raws.split(' ') if len(i) != 0]
    sat_info = np.array(lista)
    # print sat_info
    sat_info = sat_info.reshape(len(lista) / 28, 28)
    nav = xarray.DataArray(data=np.concatenate((np.atleast_2d(sv).T, sat_info), axis=1),
                           coords={'t': epoch,
                                   'data': ['sv', 'SVclockBias', 'SVclockDrift', 'SVclockDriftRate', 'IODE',
                                            'Crs', 'DeltaN', 'M0', 'Cuc', 'Eccentricity', 'Cus', 'sqrtA', 'TimeEph',
                                            'Cic', 'OMEGA', 'CIS', 'Io', 'Crc', 'omega', 'OMEGA DOT', 'IDOT',
                                            'CodesL2', 'GPSWeek', 'L2Pflag', 'SVacc', 'SVhealth', 'TGD', 'IODC',
                                            'TransTime']},  # , ''
                           dims=['t', 'data'])

    return nav

def fdayofyear(y, m, d):
    d0 = date(y, 1, 1)
    d1 = date(y, m, d)
    delta = (d1 - d0)
    return delta.days+1

#nav_path = 'D:/diplomka/dp_projekt/data/hofn/navig/hofn2740.15n'
#navdata = rinexnav(nav_path)