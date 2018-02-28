# GPS Time is a uniformly counting time scale beginning at the 1/5/1980 to 1/6/1980 midnight. January 6, 1980
# is a Sunday. GPS Time counts in weeks and seconds of a week from this instant. The weeks begin at
# the Saturday/Sunday transition.

#Note that previous result can be used to calculate UTC offset of your current timezone.
# In this example this is +1h, i.e. UTC+0100.
import datetime, time
import math as m
import urllib
import os, sys
import subprocess
isfile = os.path.isfile
join = os.path.join

#subprocess.call('set PATH=%PATH%;C:\\Program Files\\7-Zip\\', shell=False)  # for unzip

def greg2gps(epochT, gregT):
    days = abs(epochT - gregT) / one_day
    weeks = m.floor(days / 7)
    seconds = m.floor((abs(gregT - epochT) - weeks * 7 * one_day) / 1000)
    return int(weeks)

def gps2greg(week, sow):
    gregDate = epochT + week * 7 * one_day + sow * 1000
    timeZone = 3600
    time = datetime.datetime.fromtimestamp((gregDate/1000)-timeZone)
    return time


import ftplib


def download(ftp, directory, file):
    ftp.cwd(directory)
    f = open(file, "wb")
    ftp.retrbinary("RETR " + file, f.write)
    f.close()


one_day = 1000 * 60 * 60 * 24
# epoch time
epochYear = 1980
epochMonth = 1
epochDay = 6
epochHours = 0
epochMinutes = 0
epochSeconds = 0
epochMilliseconds = 0

#start time for calculate gps week
gregYear = 2016
gregMonth = 10
gregDay = 1
gregHours = 0
gregMinutes = 0
gregSeconds = 0
gregMilliseconds = 0

#end time for calculate gps week
EndYear = 2017
EndMonth = 4
EndDay = 30
EndHours = 0
EndMinutes = 0
EndSeconds = 0
EndMilliseconds = 0

#sow = 0 # seconds of the week
#week = 1864 # gps week

# all times in milliseconds
epochT = datetime.datetime(epochYear, epochMonth, epochDay, epochHours  , epochMinutes, epochSeconds, epochMilliseconds)
epochT = time.mktime(epochT.timetuple()) * 1000

gregT = datetime.datetime(gregYear, gregMonth, gregDay, gregHours, gregMinutes, gregSeconds, gregMilliseconds)
gregT = time.mktime(gregT.timetuple()) * 1000

EndT = datetime.datetime(EndYear, EndMonth, EndDay, EndHours, EndMinutes, EndSeconds, EndMilliseconds)
EndT = time.mktime(EndT.timetuple()) * 1000

if gregT < epochT:
    print 'Date is before GPS epoch'
    sys.exit()

if EndT < gregT:
    print 'End before start'
    sys.exit()

#g = gps2greg(week, sow)
startw = greg2gps(epochT, gregT)
w = startw
endw = greg2gps(epochT, EndT)

path = 'data/eph/'

try:
    os.makedirs(path)
except OSError:
    if not os.path.isdir(path):
        raise

while w <= endw:
    dir = 'data\\eph\\' + str(w)

    try:
        os.makedirs(dir)
    except OSError:
        if not os.path.isdir(dir):
            raise
    d = 0

    while d <=6:  # day of the week
        file = 'igr' + str(w) + str(d) + '.sp3.Z'
        ftp = 'ftp://www.igs.org/pub/product/' + str(w) + '/' + file
        #print ftp
        urllib.urlretrieve(ftp, path + str(w) +'/' + file)
        urllib.urlcleanup()    # activate for python 2.7 (https://stackoverflow.com/questions/44733710/downloading-second-file-from-ftp-fails)
        args0 = '7z e -o' + os.getcwd() + '\\' + dir + ' ' + os.getcwd() + '\\' + dir + '\\'+ file
        #print args0
        subprocess.call(args0, shell=False)
        os.remove(dir + '/' + file)
        d = d+1
    print 'GPS week ' + str(w) + ' has been downloaded'
    w = w+1