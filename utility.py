#-------------------------------------------------------------------------------
# Name:        supportive functions for slants computations
# Purpose:     give support for main.py
#
# Author:      Michal Kacmarik, Institute of Geoinformatics, VSB-TUO
#
# Created:     2015-2017
# Copyright:   (c) Michal Kacmarik 2017
# Licence:     GNU
#-------------------------------------------------------------------------------

import math
import os
import subprocess #to call external exe files
import datetime, time
from datetime import date
import re
import sys
import numpy
import string
import gzip

#-------------------------------------------------------------------------------
#returns current date and time
def fCurrentDateTime():
    current_date = datetime.datetime.now()
    current_year = str(current_date).split()[0][2:4]
    current_doy = current_date.timetuple().tm_yday
    current_time = (str(current_date).split()[1]).split(':')
    current_time_seconds = int((float(current_time[0]) * 3600) + (float(current_time[1]) * 60) + float(current_time[2]))
    date_list = [current_year, current_doy, current_time_seconds]
    return date_list
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#prepocet casu v sekundach do podoby hhmm, pouziva se pri praci se ZTD ve formatu TRO
def seconds_to_hours_minutes_seconds(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    cas= "%02d%02d%02d" % (h, m, s)
    return cas

#-------------------------------------------------------------------------------
#prevod XYZ do geografickych souradnic, pokud se vezmou souradnice z TRO, presnost
#zpetne trasnformace do XYZ je dle testovani do 1 mm ve vsech slozkach souradnic
#postup prevzat z Hoffman-Wellenhof, provadi se dve iterace
def XYZ_to_LatLonH(X,Y,Z):
    X=float(X)
    Y=float(Y)
    Z=float(Z)

    a=6378137.0
    b=6356752.3142
    e_2=((a*a)-(b*b))/(a*a)
    e=math.sqrt(e_2)

    lon = math.atan(Y/X)
    lon_degrees = math.degrees(lon)
    p = math.sqrt(X*X+Y*Y)

    lat0=math.atan((Z/p)*(1/(1-e_2)))
    cos2lat0 = (1+math.cos(2*lat0))/2
    sin2lat0 = (1-math.cos(2*lat0))/2
    N0 = (a*a)/math.sqrt(a*a*cos2lat0 + b*b*sin2lat0)
    h0=(p/math.cos(lat0))-N0
    lat1=math.atan((Z/p)*(1/(1-(e_2)*(N0/(N0+h0)))))

    cos2lat1 = (1+math.cos(2*lat1))/2
    sin2lat1 = (1-math.cos(2*lat1))/2
    N1 = (a*a)/math.sqrt(a*a*cos2lat1 + b*b*sin2lat1)
    h1=(p/math.cos(lat1))-N1
    lat2=math.atan((Z/p)*(1/(1-(e_2)*(N1/(N1+h1)))))
    lat2_degrees = math.degrees(lat2)

    cos2lat2 = (1+math.cos(2*lat2))/2
    sin2lat2 = (1-math.cos(2*lat2))/2
    N2 = (a*a)/math.sqrt(a*a*cos2lat2 + b*b*sin2lat2)
    h2=(p/math.cos(lat2))-N2
    lat3=math.atan((Z/p)*(1/(1-(e_2)*(N2/(N2+h2)))))
    lat3_degrees = math.degrees(lat3)

    return lat3_degrees, lon_degrees, h2

#-------------------------------------------------------------------------------
#konverze geografickych souradnic v DDMMSS do DD.DDD
def dd_mm_ss_to_ddd(myList):
  D,M,S = myList.split(' ')
  DD = ((60/float(S)) + 60/float(M)) + float(D)
  return DD

#-------------------------------------------------------------------------------
#prevod geografickych souradnic do XYZ
def latlon_to_XYZ(latitude,longtitude,elevation):
    lat_rad = math.radians(latitude)
    lon_rad = math.radians(longtitude)
    a=6378137.0
##    b=6356752.3142 #wgs84
    b=6356752.314140 #etrs

    Njmenovatel1 = (a*a)*((1+math.cos(2*lat_rad))/2)
    Njmenovatel2 = (b*b)*((1-math.cos(2*lat_rad))/2)
    N=a*a/math.sqrt(Njmenovatel1+Njmenovatel2)

    X=(N+elevation)*math.cos(lat_rad)*math.cos(lon_rad)
    Y=(N+elevation)*math.cos(lat_rad)*math.sin(lon_rad)

    Z1 = (((b*b)/(a*a))*N) + elevation
    Z=Z1*(math.sin(lat_rad))
    return [X,Y,Z]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#function for meteo RINEX files conversion to own format, converts all files in given folder
#structure of own MET format = station_name, date, time, pressure, temp
def fMeteoRINEX_myMETEO(Met_RIN,Met_MET):
    print '--------------------------------------------------------------------------------------------------'
    print "convertion of meteo RINEX files to own MET format was started - function "+fMeteoRINEX_myMETEO.__name__
    print '--------------------------------------------------------------------------------------------------'

    flag_error = 0
    if len(os.listdir(Met_RIN)) > 0:
        for file in os.listdir(Met_RIN):
            # Nacteni obsahu souboru
            f = open(Met_RIN+'/'+file, 'r',)
            myLines = f.readlines()
            f.close()

            myMeteoParam=[]
            soubor=[]

            for myLine in myLines:
                # Zjisteni nazvu meteo stanice
                if not re.search('MARKER NAME', myLine) == None:
                  myLine = myLine.strip('\n')
                  myMarkerName = (re.split('\s+',myLine.strip(' '))[0][0:4]).lower()

                # Zjisteni typu observovanym meteo dat a jejich poradi
                if not re.search('TYPES OF OBSERV', myLine) == None:
                  myLine = myLine.strip('\n')
                  myHeaderLine = re.split('\s+',myLine.strip(' '))

                  parNumber = int(myHeaderLine[0])
                  for i in range(parNumber):
                    myMeteoParam.append(myHeaderLine[i+1])

                # Detekce konce hlavicky a zacatku meteo data
                if not re.search('END OF HEADER', myLine) == None:
                  myHeaderFlag = myLines.index(myLine)

            #priprava prace se souborem a otevirani vystupu
##            print file, myMarkerName, myMeteoParam
            #zamezeni nacteni prvniho zaznamu souboru, kde je udaj pro 23:xx = dela to pak bordel pri dalsim zpracovani
            if myLines[myHeaderFlag+1].split()[3] == '23':
                telo = myLines[myHeaderFlag+2:]
            else:
                telo = myLines[myHeaderFlag+1:]
            datum = (telo[0][1:3].replace(' ','0'))+(telo[0][4:6].replace(' ','0'))+(telo[0][7:9].replace(' ','0'))
            zapis = open(Met_MET+'/'+myMarkerName+datum+'.MET', 'w')

            #zjisteni indexu na ktereych se v meteo RINEXu nachazi hodnoty tlaku (PR)
            #a teploty (TD) - ruzni se to mezi soubory, tato cast skriptu zajisti to
            #ze do vystupu se vzdy tiskne nejprve tlak a pak teplota
            #je to reseno diky pevne strukture zaznamu v meteo RINEXu
            prvek_pr = myMeteoParam.index('PR')
            index_PR = [(18+(prvek_pr*7)),(25+(prvek_pr*7))]
            prvek_td = myMeteoParam.index('TD')
            index_TD = [(18+(prvek_td*7)),(25+(prvek_td*7))]
##            print index_PR, index_TD

            #prochazeni jednotlive radky vstupnich dat
            for myLine in telo:
                myLine = myLine.strip('\n')
                #parsovani datumu a textu do podoby yymmdd
                date = (myLine[1:3].replace(' ','0'))+(myLine[4:6].replace(' ','0'))+(myLine[7:9].replace(' ','0'))
                time = (myLine[10:12].replace(' ','0'))+(myLine[13:15].replace(' ','0'))
                pressure = float(myLine[index_PR[0]:index_PR[1]])
                temperature = float(myLine[index_TD[0]:index_TD[1]])
                if temperature > 100:
                    temperature = temperature - 273.1
                myLineNew = myMarkerName + ' ' + date + ' ' + time + ' ' + str(pressure) + ' ' +str(temperature) + '\n'

                #pokud je interval zaznamu v meteo RINEX kratsi nez 1 minuta, tak se
                #do vystupu ukladaji jen ty v celou minutu
                #nasledujici radky filtruji pridavani radku do pole soubor tak, ze
                #pokud se v nem jiz nachazi radek s hodnotou hhmm
                if len(soubor)>1:
                    posledni_radek = soubor[-1]
                    if not posledni_radek[0:16] == myLineNew[0:16]:
                        soubor.append(myLineNew)
                else:
                    soubor.append(myLineNew)

            #doplnek predchoziho filtrovani zaznamu s kratsim nez minutovym intervalem,
            #resi, ze se do vystupu nezapisuji dva zaznamy pro prvni minutu, jelikoz
            #predchozi filtr vzdy vyhodi dva zaznamy pro prvni minutu a az dalsi zaznamy
            #filtruje korektne
            if soubor[0][0:16] != soubor[1][0:16]:
                zapis.writelines(soubor)
            else:
                zapis.writelines(soubor[1:])
            zapis.close()
    else:
        print '------------------------------------------------------------------------------------------'
        print 'Error in function '+fMeteoRINEX_myMETEO.__name__ + ': no Meteo RINEX file with observations'
        print '------------------------------------------------------------------------------------------'
        flag_error = 1
        sys.exit()

    if flag_error == 0:
        print '--------------------------------------------------------------------------------------------------------------------'
        print "data import from Meteo RINEX files to own MET files was successful, result files are stored in: "+Met_MET
        print '--------------------------------------------------------------------------------------------------------------------'
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#get a list of stations for which the processing will be run = observation RINEX
#files are used
def fGet_Station_List(input_RINOBS):
    stations_list_GNSS = []
    if len(os.listdir(input_RINOBS)) > 0:
        for file in os.listdir(input_RINOBS):
            stations_list_GNSS.append(file[0:4].lower())
    else:
        print '--------------------------------------------------------------------'
        print 'Error in function '+fGet_Station_List.__name__ + ': no RINEX OBS file'
        print '--------------------------------------------------------------------'
        sys.exit()
    return stations_list_GNSS
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#read stations coordinates for processed stations from CRD file stored in STA directory
def fRead_CRD_File_Create_STA(STA, station_list):
    for file in os.listdir(STA):
        if 'CRD' in file:
            if '.gz' in file or '.GZ' in file:
                with gzip.open(STA+'/'+file, 'r') as f:
                    lines = f.readlines()
            else:
                read = open(STA+'/'+file, 'r')
                lines = read.readlines()
                read.close()

    station_coordinates = []
    for line in lines:
        if 'LOCAL GEODETIC DATUM' in line:
            geodetic_datum = line.split(':')[1].split()[0]
            station_coordinates.append(geodetic_datum)
            break
    for line in lines:
        if len(line.split())>2:
            if line.split()[1].upper() in station_list or line.split()[1].lower() in station_list:
                line_list = line.split()
                station_coordinates.append(line.split()[1] + ' ' + line.split()[3] + ' ' + line.split()[4] + ' ' + line.split()[5])

    #create gnss.STA file - necessary for compute_SIWV (last step of processing)
    d = open(STA+'/gnss.STA', 'w')
    for radek in station_coordinates[1:]:
        radek_coord = radek.split()
        geograf = XYZ_to_LatLonH(radek_coord[1],radek_coord[2],radek_coord[3])
        #for LON values interval 0-360 is used (not -180-0, 0-180)
        if geograf[1] < 0:
            lon = 360.0 + geograf[1]
        else:
            lon = geograf[1]
        d.write(radek_coord[0] + ' ' + radek_coord[1] + ' ' + radek_coord[2] + ' ' + radek_coord[3] + ' '
        + str(geograf[0]) + ' ' + str(lon) + ' ' + str(geograf[2]) +'\n')
    d.close()

    return station_coordinates
#-------------------------------------------------------------------------------

#read ZHD, ZWD, ZTD values from Bernese TRP file, only data for necessary
#stations are obtained
#-------------------------------------------------------------------------------
def fRead_TRP_File(input_ZTD, station_list, ZTD_input_format):
    TRP_list = []
    stations_processed = []
    stations_formal_error_dict = {}
    for file in os.listdir(input_ZTD):
        if ZTD_input_format in file:
            if '.gz' in file or '.GZ' in file:
                with gzip.open(input_ZTD+'/'+file, 'r') as f:
                    lines = f.readlines()
            else:
                read = open(input_ZTD+'/'+file, 'r')
                lines = read.readlines()
                read.close()

            if len(TRP_list) == 0:
                for line in lines:
                    if len(line.split())>3 and line.split()[3].isdigit():
                        year_day = line.split()[3][2:] + line.split()[4] + line.split()[5]
                        TRP_list.append(year_day)
                        break
            date_to = 19000101000000

            file_station_list = []
            stations_indexes = []
            for i in range(6, len(lines)):
                if len(lines[i].split()) > 1 and lines[i].split()[0].lower() not in file_station_list:
                    file_station_list.append(lines[i].split()[0].lower())
                    stations_indexes.append(i)
                    if len(lines[i].split()) == 23:
                        index_plus = 6
                    else:
                        index_plus = 0
            stations_indexes.append(len(lines))

            #if station occurs in more than one TRP file, the solution with the minimum
            #daily average formal error of ZTD is written to result files
            for j in range(0, len(stations_indexes)-1):
                if file_station_list[j] not in stations_processed:
                    if file_station_list[j].lower() in station_list:
                        formal_error_one_station_list = []
                        for line in lines[stations_indexes[j]:stations_indexes[j+1]]:
                            if len(line.split()) > 1:
                                line_list = line.split()
                                date_to_line = int(line_list[3]+line_list[4]+line_list[5]+line_list[6]+line_list[7]+line_list[8])
                                if date_to_line > date_to:
                                    date_to = date_to_line
                                ZTD_record = [line_list[0], line_list[3][2:]+line_list[4]+line_list[5],
                                line_list[6]+line_list[7]+line_list[8], line_list[9+index_plus], line_list[10+index_plus], line_list[12+index_plus], line_list[11+index_plus],
                                line_list[13+index_plus], line_list[14+index_plus], line_list[15+index_plus], line_list[16+index_plus]]
                                TRP_list.append(ZTD_record)
                                formal_error_one_station_list.append(float(line_list[11+index_plus]))
                        formal_error_station = sum(formal_error_one_station_list)/len(formal_error_one_station_list)
                        stations_formal_error_dict.update({file_station_list[j]:formal_error_station})
                        stations_processed.append(file_station_list[j])
                else:
                    if file_station_list[j].lower() in station_list:
                        formal_error_one_station_list = []
                        for line in lines[stations_indexes[j]:stations_indexes[j+1]]:
                            if len(line.split()) > 1:
                                formal_error_one_station_list.append(float(line.split()[11+index_plus]))
                        formal_error_station = sum(formal_error_one_station_list)/len(formal_error_one_station_list)
                        if formal_error_station < stations_formal_error_dict.get(file_station_list[j]):
                            stations_formal_error_dict[file_station_list[j]] = formal_error_station
                            TRP_list = [x for x in TRP_list if x[0] != file_station_list[j]]
                            for line in lines[stations_indexes[j]:stations_indexes[j+1]]:
                                if len(line.split()) > 1:
                                    line_list = line.split()
                                    date_to_line = int(line_list[3]+line_list[4]+line_list[5]+line_list[6]+line_list[7]+line_list[8])
                                    if date_to_line > date_to:
                                        date_to = date_to_line
                                    ZTD_record = [line_list[0], line_list[3][2:]+line_list[4]+line_list[5],
                                    line_list[6]+line_list[7]+line_list[8], line_list[9+index_plus], line_list[10+index_plus], line_list[12+index_plus], line_list[11+index_plus],
                                    line_list[13+index_plus], line_list[14+index_plus], line_list[15+index_plus], line_list[16+index_plus]]
                                    TRP_list.append(ZTD_record)
    TRP_list.insert(1,str(date_to))
    return TRP_list
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def fConvert_TRP_to_IWV_files(ZTD_TRP_list,station_coordinates,limit_sdev_ZTD,result_IWV,result_trosnx,trosinex_flag,SNX_basis_path,solution_name):
    print '-----------------------------------------------------------------------------------------------------------'
    print "TRP to own ZWD/IWV format conversion was started - function "+fConvert_TRP_to_IWV_files.__name__
    print '-----------------------------------------------------------------------------------------------------------'
    flag_error = 0

    if trosinex_flag == 'Y':
        read_snx_base = open(SNX_basis_path, 'r')
        snx_base = read_snx_base.readlines()
        read_snx_base.close()

    for prvek in station_coordinates[1:]:
        station = prvek.split()[0]
        year_day = ZTD_TRP_list[0]
        zapis = open(result_IWV+'/'+station.upper()+year_day+'.IWV', 'w')
        doy = datetime.datetime(int(year_day[0:2]),int(year_day[2:4]),int(year_day[4:6])).timetuple().tm_yday
        year_doy = year_day[0:2]+str(doy)

        if trosinex_flag == 'Y':
            list_trosinex_ztd_data = []

        for ztd_record in ZTD_TRP_list[2:]:
            if station in ztd_record and float(ztd_record[6])*1000 < limit_sdev_ZTD:
                zapis.write(station.upper() + ' ' + ztd_record[1] + ' ' + ztd_record[2][0:4] + str("{0:.1f}".format(float(ztd_record[5])*1000)).rjust(8, ' ') + str("{0:.1f}".format(float(ztd_record[6])*1000)).rjust(5, ' ')
                + str("{0:.3f}".format(float(ztd_record[7])*1000)).rjust(8, ' ') + str("{0:.3f}".format(float(ztd_record[8])*1000)).rjust(7, ' ') + str("{0:.3f}".format(float(ztd_record[9])*1000)).rjust(8, ' ') + str("{0:.3f}".format(float(ztd_record[10])*1000)).rjust(7, ' ')
                + ' -9999' + '  -9999' + '  -9999' + str("{0:.1f}".format(float(ztd_record[3])*1000)).rjust(8, ' ') + str("{0:.1f}".format(float(ztd_record[4])*1000)).rjust(8, ' ') + '    -9999' + '\n')

                if trosinex_flag == 'Y':
                    ztd_date_object = datetime.datetime(int(ztd_record[1][0:2])+2000,int(ztd_record[1][2:4]),int(ztd_record[1][4:6]))
                    seconds = int(ztd_record[2][0:2])*3600 + int(ztd_record[2][2:4])*60 + int(ztd_record[2][4:6])
                    date_print = ztd_record[1][0:2]+':'+str(ztd_date_object.timetuple().tm_yday)+':'+str(seconds).rjust(5,'0')
                    snx_line = (' ' + station.upper() + ' ' + date_print + str("{0:.1f}".format(float(ztd_record[5])*1000)).rjust(8, ' ') + str("{0:.1f}".format(float(ztd_record[6])*1000)).rjust(10, ' ')
                    + str("{0:.1f}".format(float(ztd_record[7])*1000)).rjust(10, ' ') + str("{0:.1f}".format(float(ztd_record[8])*1000)).rjust(10, ' ') + str("{0:.1f}".format(float(ztd_record[9])*1000)).rjust(10, ' ')
                    + str("{0:.1f}".format(float(ztd_record[10])*1000)).rjust(10, ' ') + str("{0:.1f}".format(float(ztd_record[3])*1000)).rjust(10, ' ') + str("{0:.1f}".format(float(ztd_record[4])*1000)).rjust(10, ' ') +'\n')
                    list_trosinex_ztd_data.append(snx_line)

        #TROSNX = cast resici tisk ZTD reseni + ZHD, ZWD do TROSNX souboru
        separator = '*-------------------------------------------------------------------------------\n'
        if trosinex_flag == 'Y':
            date_from = year_doy[0:2]+':'+year_doy[2:5]+':00000'
            date_to_input = ZTD_TRP_list[1]
            date_to_object = datetime.datetime(int(date_to_input[0:4]),int(date_to_input[4:6]),int(date_to_input[6:8]))
            date_to = date_to_input[2:4]+':'+str(date_to_object.timetuple().tm_yday)+':00000'

            write_snx = open(result_trosnx+'/'+station.upper()+'_'+year_doy+'.SNX', 'w')
            current_date_time = fCurrentDateTime()
            #write header
            #%=TRO 0.09 GOP 15:225:31188 GOP 10:360:00000 10:360:85500 P  MIX
            write_snx.write('%=TRO 0.09 '+solution_name+' '+current_date_time[0]+':'+str(current_date_time[1])+':'+str(current_date_time[2]).rjust(5,'0')+' '+solution_name+' '+ date_from + ' ' + date_to + ' P  MIX\n' )
            write_snx.write(separator)

            #write base from file
            write_snx.writelines(snx_base)

            #write coordinates
            write_snx.write('+TROP/STA_COORDINATES\n*SITE PT SOLN T __STA_X_____ __STA_Y_____ __STA_Z_____ SYSTEM REMRK\n')
            write_snx.write(' '+prvek.split()[0] + '  A    1 P' + str("{0:.3f}".format(float(prvek.split()[1]))).rjust(13, ' ') + str("{0:.3f}".format(float(prvek.split()[2]))).rjust(13, ' ') + str("{0:.3f}".format(float(prvek.split()[3]))).rjust(13, ' ') + station_coordinates[0].rjust(6, ' ') + solution_name.rjust(5, ' ')+'\n')
            write_snx.write('-TROP/STA_COORDINATES\n')
            write_snx.write(separator)

            write_snx.write('+TROP/SOLUTION\n*SITE ____EPOCH___  TROTOT    STDDEV    TGNTOT    STDDEV    TGETOT    STDDEV    TROHYD    TROWET\n')
            write_snx.writelines(list_trosinex_ztd_data)
            write_snx.write('-TROP/SOLUTION\n')
            write_snx.write(separator)
            write_snx.close()
        zapis.close()

    #delete files for which RINEX observations were available, but not ZTD estimates
    for file in os.listdir(result_IWV):
        if os.stat(result_IWV+'/'+file).st_size == 0:
            os.remove(result_IWV+'/'+file)
    return year_doy
    print '-----------------------------------------------------------------------------------------------------------'
    print "TRP to own ZWD/IWV format conversion was succesfull - function "+fConvert_TRP_to_IWV_files.__name__
    print '-----------------------------------------------------------------------------------------------------------'
#-------------------------------------------------------------------------------

#read ZTD+grad values from TROSNX file, only data for necessary
#stations are obtained
#-------------------------------------------------------------------------------
def fRead_TRO_File(input_ZTD,STA,meteo_STA_file,Met_MET,ZTD_input_format):
    station_list = []
    stations_coord_meteo = []
    #read list of METEO stations from meteo_STA file, create a list of stations
    #from it
    cti_MET_STA = open(STA+'/'+meteo_STA_file, 'r')
    metSTAlines = cti_MET_STA.readlines()
    cti_MET_STA.close()
    for radek in metSTAlines[1:]:
        meteo_sour = []
        radek_coord = radek.split()
        meteo_sour.append(radek_coord[0].lower())
        meteo_sour.append(radek_coord[1].lower())
        meteo_sour.append(float(radek_coord[2]))
        meteo_sour.append(float(radek_coord[3]))
        meteo_sour.append(float(radek_coord[4]))
        stations_coord_meteo.append(meteo_sour)
        station_list.append(radek_coord[1].lower())

    #check if for station is METEO file available, if not its ZTD values are not loaded
    #since the station is removed from station_list
    station_list_MET = []
    for file1 in os.listdir(Met_MET):
        station_list_MET.append(file1[0:4].lower())

    station_list_2 = station_list
    station_list = []
    for station in station_list_2:
        for prvek in stations_coord_meteo:
            if station == prvek[1]:
                if prvek[0] in station_list_MET:
                    station_list.append(station)

    zapis_GNSS_STA = open(STA+'/gnss.STA','w')
    TRO_list = []
    stations_coord_GNSS = []
    stations_processed = []
    stations_formal_error_dict = {}
    for file in os.listdir(input_ZTD):
        if ZTD_input_format in file:
            if '.gz' in file or '.GZ' in file:
                with gzip.open(input_ZTD+'/'+file, 'r') as f:
                    lines = f.readlines()
            else:
                read = open(input_ZTD+'/'+file, 'r')
                lines = read.readlines()
                read.close()
            if len(TRO_list) == 0:
                first_line = lines[0]
                if len(first_line.split())>3 and (first_line.split()[5]).split(':')[0].isdigit():
                    date_from = first_line.split()[5].split(':')[0] + first_line.split()[5].split(':')[1]
                    date_to = first_line.split()[6].split(':')[0] + first_line.split()[6].split(':')[1]
                    TRO_list.append(date_from)
                    TRO_list.append(date_to)

            #prochazeni TRO souboru a hledani indexu ohranicujicich jeho jednotlive casti
            for myLine in lines:
                #popis nastaveni urcovani parametru troposfery
                if '+TROP/DESCRIPTION' in myLine:
                    TropStartFlag = lines.index(myLine)
                if '-TROP/DESCRIPTION' in myLine:
                    TropStoFlag = lines.index(myLine)
                #souradnice stanic v XYZ
                if '+TROP/STA_COORDINATES' in myLine:
                    CoordStartFlag = lines.index(myLine)
                if '-TROP/STA_COORDINATES' in myLine:
                    CoordStopFlag = lines.index(myLine)
                #samotne reseni troposfery
                if '+TROP/SOLUTION' in myLine:
                    bodyStartFlag = lines.index(myLine)
                if '-TROP/SOLUTION' in myLine:
                    bodyStopFlag = lines.index(myLine)

            file_station_list = []
            stations_indexes = []
            for i in range(bodyStartFlag+1, bodyStopFlag):
                if len(lines[i].split()) > 1 and lines[i].split()[0].lower() not in file_station_list:
                    file_station_list.append(lines[i].split()[0].lower())
                    stations_indexes.append(i)
            stations_indexes.append(len(lines))

            #if station occurs in more than one TRP file, the solution with the minimum
            #daily average formal error of ZTD is written to result files
            for j in range(0, len(stations_indexes)-1):
                if file_station_list[j] not in stations_processed:
                    if file_station_list[j].lower() in station_list:
                        formal_error_one_station_list = []
                        for line in lines[stations_indexes[j]:stations_indexes[j+1]]:
                            if len(line.split()) > 1:
                                TRO_list.append(line)
                                formal_error_one_station_list.append(float(line.split()[3]))
                        formal_error_station = sum(formal_error_one_station_list)/len(formal_error_one_station_list)
                        stations_formal_error_dict.update({file_station_list[j]:formal_error_station})
                        stations_processed.append(file_station_list[j])
                        for radek in lines[CoordStartFlag+2:CoordStopFlag]:
                            stan_sour = []
                            radek_coord = radek.split()
                            if radek_coord[0].lower() == file_station_list[j]:
                                geograf = XYZ_to_LatLonH(radek_coord[4],radek_coord[5],radek_coord[6])
                                zapis_GNSS_STA.write(radek_coord[0] + ' ' + radek_coord[4] + ' ' + radek_coord[5] + ' ' + radek_coord[6] + ' '
                                + str(geograf[0]) + ' ' + str(geograf[1]) + ' ' + str(geograf[2]) +'\n')
                                if len(stations_coord_GNSS) == 0:
                                    stations_coord_GNSS.append(radek_coord[7])
                                stations_coord_GNSS.append([radek_coord[0].lower(), radek_coord[4], radek_coord[5], radek_coord[6], geograf[0], geograf[1], geograf[2]])
                                break
                else:
                    if file_station_list[j].lower() in station_list:
                        formal_error_one_station_list = []
                        for line in lines[stations_indexes[j]:stations_indexes[j+1]]:
                            if len(line.split()) > 1:
                                formal_error_one_station_list.append(float(line.split()[3]))
                        formal_error_station = sum(formal_error_one_station_list)/len(formal_error_one_station_list)
                        if formal_error_station < stations_formal_error_dict.get(file_station_list[j]):
                            stations_formal_error_dict[file_station_list[j]] = formal_error_station
                            TRO_list = [x for x in TRO_list if x[1:5].lower() != file_station_list[j]]
                            for line in lines[stations_indexes[j]:stations_indexes[j+1]]:
                                if len(line.split()) > 1:
                                    TRO_list.append(line)
                            stations_coord_GNSS = [x for x in stations_coord_GNSS if x[0] != file_station_list[j]]
                            for radek in lines[CoordStartFlag+2:CoordStopFlag]:
                                stan_sour = []
                                radek_coord = radek.split()
                                if radek_coord[0].lower() == file_station_list[j]:
                                    geograf = XYZ_to_LatLonH(radek_coord[4],radek_coord[5],radek_coord[6])
                                    zapis_GNSS_STA.write(radek_coord[0] + ' ' + radek_coord[4] + ' ' + radek_coord[5] + ' ' + radek_coord[6] + ' '
                                    + str(geograf[0]) + ' ' + str(geograf[1]) + ' ' + str(geograf[2]) +'\n')
                                    stations_coord_GNSS.append([radek_coord[0].lower(), radek_coord[4], radek_coord[5], radek_coord[6], geograf[0], geograf[1], geograf[2]])
                                    break
    zapis_GNSS_STA.close()
    return TRO_list, station_list, stations_coord_GNSS, stations_coord_meteo
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#if meteo RINEX files plus TRO files are used, they are firstly concatenated and then ZHD, ZWD, IWV computed
def fConcat_Meteo_ZTD_Get_IWV(TRO_list,stations_list_GNSS,stations_coord_GNSS,stations_coord_meteo,Met_MET,limit_sdev_ZTD,max_cas_rozdil_meteo_ZTD,STA,result_IWV,trosinex_flag,result_trosnx,SNX_basis_path,solution_name):
    print '-----------------------------------------------------------------------------------------------------------'
    print "ZWD/IWV computation was started - function "+fConcat_Meteo_ZTD_Get_IWV.__name__
    print '-----------------------------------------------------------------------------------------------------------'
    flag_error = 0

    stations_list_meteo = []
    station_meteo_dictionary = {}
    if len(os.listdir(Met_MET)) > 0:
        for file in os.listdir(Met_MET):
            stations_list_meteo.append(file[0:4].lower())
            station_meteo_dictionary.update({file[0:4].lower():file})
            date_meteo = file[4:10]
    else:
        print '--------------------------------------------------------------------'
        print 'Error in function '+fConcat_Meteo_ZTD_Get_IWV.__name__ + ': no Meteo MET file'
        print '--------------------------------------------------------------------'
        flag_error = 1
        sys.exit()

    #prochazeni jednotlivych meteo stanic = nejvyssi uroven cyklu funkce
    for station in stations_list_GNSS:
        print "current proccessed station is: ", station
        zapis = open(result_IWV+'/'+station.upper()+date_meteo+'.IWV', 'w')

        linesMeteo = []
        MeteoKey = {}

        year_doy = TRO_list[0]
        year_doy_plus1 = TRO_list[1]
        year_doy_print = year_doy[0:2]+':'+year_doy[2:5]+':00000'
        year_doy_plus1_print = year_doy_plus1[0:2]+':'+year_doy_plus1[2:5]+':00000'

        #TROSNX = cast resici tisk ZTD reseni + ZHD, ZWD do TROSNX souboru
        separator = '*-------------------------------------------------------------------------------\n'
        if trosinex_flag == 'Y':
            write_snx = open(result_trosnx+'/'+station.upper()+'_'+year_doy+'.SNX', 'w')

            read_snx_base = open(SNX_basis_path, 'r')
            snx_base = read_snx_base.readlines()
            read_snx_base.close()

            current_date_time = fCurrentDateTime()
            #write header
            #%=TRO 0.09 GOP 15:225:31188 GOP 10:360:00000 10:360:85500 P  MIX
            write_snx.write('%=TRO 0.09 '+solution_name+' '+current_date_time[0]+':'+str(current_date_time[1])+':'+str(current_date_time[2]).rjust(5,'0')+' '+solution_name+' '+ year_doy_print + ' ' + year_doy_plus1_print + ' P  MIX\n' )
            write_snx.write(separator)

            #write base from file
            write_snx.writelines(snx_base)

            write_snx.write('+TROP/STA_COORDINATES\n*SITE PT SOLN T __STA_X_____ __STA_Y_____ __STA_Z_____ SYSTEM REMRK\n')
            #write coordinates
            for prvek in stations_coord_GNSS[1:]:
                if station in prvek:
                    write_snx.write(' '+prvek[0].upper() + '  A    1 P' + prvek[1].rjust(13, ' ') + prvek[2].rjust(13, ' ') + prvek[3].rjust(13, ' ') + stations_coord_GNSS[0].rjust(6, ' ') + solution_name.rjust(5, ' ')+'\n')
            write_snx.write('-TROP/STA_COORDINATES\n')
            write_snx.write(separator)
            write_snx.write('+TROP/SOLUTION\n*SITE ____EPOCH___  TROTOT    STDDEV    TGNTOT    STDDEV    TGETOT    STDDEV    TROHYD    TROWET\n')

        #nacteni meteo souboru pro stanici, se kterou se aktualne pracuje
        for prvek in stations_coord_meteo:
            if prvek[1] == station:
                station_meteo = prvek[0]
        soubor = station_meteo_dictionary[station_meteo]
        g = open(Met_MET+'/'+soubor, 'r',)
        linesMeteo = g.readlines()
        g.close()
        #zjisteni vsech intervalu, pro ktere jsou v meteo souboru udaje, ulozeni
        #do slovniku (MeteoKey) do podoby 'index_radku_souboru':'cas zaznamu'
        #slouzi pri naslednem vybirani nejblizsiho meteo mereni pro ZTD
        for i in range (0,len(linesMeteo)):
            MeteoKey.update({i:linesMeteo[i][12:16]})
        #ziskani lat, lon pro GNSS a meteo pro aktualne zpracovavanou stanici
        for j in range(len(stations_coord_GNSS)):
            if stations_coord_GNSS[j][0].lower() == station or stations_coord_GNSS[j][0].upper() == station:
                coord_GNSS = stations_coord_GNSS[j]
        for k in range(len(stations_coord_meteo)):
            if stations_coord_meteo[k][1].lower() == station or stations_coord_meteo[k][1].upper() == station:
                coord_meteo = stations_coord_meteo[k]

        #rozdil vysky H pro GNSS a meteo stanici
##        print coord_GNSS[3], coord_meteo[3]
        height_difference = coord_meteo[4] - coord_GNSS[6]

        #porovnani lat, lon mezi GNSS a meteo
        if abs(coord_GNSS[4]-coord_meteo[2]) > 0.01 and abs(coord_GNSS[5]-coord_meteo[3]) > 0.01:
            print station + ': position of GNSS reference station and meteo stations differs in latitude of:' + str(coord_GNSS[4]-coord_meteo[2]) + ' , in longitude of: ' + str(coord_GNSS[5]-coord_meteo[3]) + ', altitude difference is: ' + str(height_difference)
        else:
            print station + ': position of GNSS reference station corresponds to position of meteo station, their altitude difference is: ' + str(height_difference)
        #prochazeni radku v tele TRO
        for lineZTD in TRO_list:
##            print lineZTD[1:5], station
            if lineZTD[1:5].lower() == station:
                radek_ztd_pole = lineZTD.split()

                #filtr na ZTD hodnoty prevysujici zadanou maximalni povolenou STDSDEV z TRO souboru (defaultne 3 mm)
                sdev_ZTD = float(radek_ztd_pole[3])
                if sdev_ZTD <= limit_sdev_ZTD:

                    year = int('20'+lineZTD[6:8])
                    doy = int(lineZTD[9:12])
                    #reseni situace s poslednim zaznamem v ZTD, ktery ma datum pro nasledujici den a hodinu 00:00
                    #nahradi se stavajicim dnem a casem 24:00 (meteo RINEX obsahuji posledni zaznam 23:59)
                    if doy == int(year_doy_plus1[2:5]):
                        date_ztd = (date.fromordinal(date(year, 1, 1).toordinal() + doy - 2)).strftime("%y%m%d")
                        time_ztd = 2400
                    else:
                        date_ztd = (date.fromordinal(date(year, 1, 1).toordinal() + doy - 1)).strftime("%y%m%d")
                        time_ztd_seconds = int(lineZTD[13:18])
                        time_ztd_60 = seconds_to_hours_minutes_seconds(time_ztd_seconds) #volani funkce pro prevod casu v s
                        #na tvar hhmm
                        hodiny = int(time_ztd_60[0:2])*100
                        minuty = int((float(time_ztd_60[2:4])/60)*100) #minuty se z sedesatkove
                        #soustavy transformuji na rozsah 0-100, aby nevznikal problem,
                        #ze se pro ZTD v case 2300 vybere zaznam 2314 misto 2259, jelikoz
                        #2314-2300 je mene nez 2259. proto se 2259 prevede na cca 2298
                        time_ztd = hodiny+minuty
    ##                print date_ztd, time_ztd

                    #vyber zaznamu s meteo udaji, ktery je nejblize ZTD hodnote
                    minDelta = 86400
    ##                print MeteoKey
                    for klic, hodnota in MeteoKey.iteritems():
                        hodiny_met = int(hodnota[0:2])*100
                        minuty_met = int((float(hodnota[2:4])/60)*100)
                        time_met = hodiny_met + minuty_met
                        time_delta = int(time_ztd) - time_met

                        if abs(time_delta) < minDelta:
                            minDelta = time_delta
                            memKey = klic

                    #priprava casu pro zapis - jelikoz posledni radek ZTD obsahuje udaj pro nasledujici den
                    #s casem 0000, musi se to resit takto (silene)
                    cas_tisk = str(time_ztd).rjust(4,'0')
                    radek_meteo_pole = re.split('\s+',linesMeteo[memKey].strip(' '))
                    #volani funkce pro vypocet ZHD, ZWD, IWV
                    IWV_pole = computeIWV(coord_GNSS[4],coord_GNSS[6],height_difference,float(radek_meteo_pole[3]),float(radek_meteo_pole[4]),float(radek_ztd_pole[2])/1000)
##                    print coord_GNSS[1],coord_GNSS[3],height_difference,float(radek_meteo_pole[3]),float(radek_meteo_pole[4]),float(radek_ztd_pole[2])/1000
##                    print IWV_pole

                    #pokud je rozdil mezi casem v ZTD a meteo zaznamem mensi nez 2 hodiny (1 hod = 100) a
                    #datum obou zaznamu je totozne, zaznam se zapise do souboru
                    if abs(int(time_ztd)-int(radek_meteo_pole[2])) < max_cas_rozdil_meteo_ZTD and date_ztd==radek_meteo_pole[1]:
                        zapis.write(station.upper() + ' ' + str(date_ztd) + ' ' + str(cas_tisk)[0:2] + str(int((float(cas_tisk[2:4])/100)*60)).rjust(2,'0') + '  '
                        + radek_ztd_pole[2] + '  ' + radek_ztd_pole[3] + '  ' + radek_ztd_pole[4].rjust(6, ' ') + '  ' + radek_ztd_pole[5] + '  ' + radek_ztd_pole[6].rjust(6, ' ') + '  '
                        + radek_ztd_pole[7] + '  ' + str("{0:.1f}".format(float(radek_meteo_pole[4]))).rjust(4,' ') + '  '
                        + str("{0:.1f}".format(float(radek_meteo_pole[3]))).rjust(4,' ') + '  '  + str("{0:.1f}".format(float(IWV_pole[0]))).rjust(4,' ') + '  '
                        + str("{0:.1f}".format(round(IWV_pole[1]*1000,2)).rjust(3,' ')) + '  ' + str("{0:.1f}".format(round(IWV_pole[2]*1000,2)).rjust(6,' ')) + '  ' + str("{0:.4f}".format(round(IWV_pole[3],4)).rjust(7,' ')+ '\n'))

                        if trosinex_flag == 'Y':
                            time_seconds = int(str(cas_tisk)[0:2])*3600 + (int((float(cas_tisk[2:4])/100)*60))*60
                            if time_seconds == 86400:
                                time_seconds = 0
                            write_snx.write(' ' + station.upper() + ' ' + str(date_ztd)[0:2]+ ':' + str(doy) + ':' + str(time_seconds).rjust(5,'0') + radek_ztd_pole[2].rjust(10, ' ')
                            + radek_ztd_pole[3].rjust(10, ' ') + radek_ztd_pole[4].rjust(10, ' ') + radek_ztd_pole[5].rjust(10, ' ') + radek_ztd_pole[6].rjust(10, ' ')
                            + radek_ztd_pole[7].rjust(10, ' ') + str("{0:.1f}".format(round(IWV_pole[1]*1000,1)).rjust(10,' ')) + str("{0:.1f}".format(round(IWV_pole[2]*1000,1)).rjust(10,' ')) +'\n')
##                Body2.remove(lineZTD)
        zapis.close()
        if trosinex_flag == 'Y':
            write_snx.write('-TROP/SOLUTION\n')
            write_snx.write(separator)
            write_snx.close()
    ##                print memKey, minDelta
    ##                print date_ztd, time_ztd, linesMeteo[memKey]

    return year_doy
    if flag_error == 0:
        print '-----------------------------------------------------------------------------------------------------------------------------'
        print "ZWD/IWV computation successful; results stored in v: "+result_IWV
        print '-----------------------------------------------------------------------------------------------------------------------------'
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def computeIWV(Lat,H_GNSS,H_delta,Pres,Temp,ZTD):
  Lat = float(Lat)
  H_GNSS = float(H_GNSS)
  H_delta = float(H_delta)
  Temp = float(Temp)
  Pres = float(Pres)
  ZTD = float(ZTD)
##  print Lat,H_GNSS,H_delta,Pres,Temp,ZTD

  # Lattitude from DD.DDDDDDD to RADIANS
  LatRad = math.radians(Lat)

  # Heigth from m to km
  H = H_GNSS / 1000.0

  #Pressure at Meteo Station recalibrated to GNSS station altitutde - Berg model
  Pres_GNSS = Pres * (math.pow((1-(0.0000226*(H_delta*-1))),5.225))

  # Temperature stupne to 273.15K
  Temp += 273.15
  #Temperature at Meteo Station recalibrated to GNSS station altitude
  Temp_GNSS = Temp-0.0065*(H_delta*-1)

  ZHD = (0.0022768 * Pres_GNSS) /((1-0.00266 * (math.cos(2*LatRad)))-0.00028*H)
  ZWD = ZTD - ZHD

  #k - verze pouzivana vsemi skupinami v GNSS4SWEC WG3 (vztah z Askne and Nordius, 1987,
  #konstanty z Bevis et al. 1994)
  c2apo = (70.4 - 77.6)*0.621977
  k = math.pow(10,6) / ((373900/((70.2 + 0.72 * Temp)+c2apo))*461.5)

  #k - moje verze pouzivana na DP, disertaci, rozdily ve vysledcich - nepouzivat
##  k = math.pow(10,6) / ((8.31434 / (18.0152 * math.pow(10,-3))) * ((3.776 * math.pow(10,5)) / (70.2 + 0.72 * Temp) + 17 - ((18.0152 * math.pow(10,-3)) / (28.9644 / 1000.0)) * 77.695))

  IWV = ZWD * k / 999.9720 * 100000
  return [Pres_GNSS,ZHD,ZWD,IWV]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#vypocet AZ, EL pro vsechny observace z RINEX souboru
#na vstupu jsou vyuzivany meteo RINEX soubory, produkty presnych efemerid druzic,
#souradnice prijimace pro vypocet AZ, EL
def fGenerateAZEL(input_RINOBS, input_EPH, STA, temp_AZEL, gnss_list, interval_slants, cut_off):
    sp3_list = []
    print '---------------------------------------------------------------------------------------------------'
    print "computation of AZI/ELE angles for observations from observation RINEX files was started - function "+fGenerateAZEL.__name__
    print '---------------------------------------------------------------------------------------------------'

    flag_error = 0
    gnss_coord = []
    #nacitani souradnic GNSS stanic z gnss.STA souboru (puvodne souradnice pochazeji z
    #TRO souboru), souradnice se pouzvaji v args1 pri vypoctu az a el druzic
    if len(os.listdir(STA)) >0:
        for file1 in os.listdir(STA):
##            if file1 != None:
            if file1 == 'gnss.STA':
                cti = open(STA+'/gnss.STA', 'r')
                radky = cti.readlines()
                cti.close()
                for radek in radky:
                    pole = radek.split()
                    gnss_coord.append(pole[0:4])
    else:
        print '------------------------------------------------------------------------------------------'
        print 'Error in function '+fGenerateAZEL.__name__ + ': no file with GNSS station coordinates found'
        print '------------------------------------------------------------------------------------------'
        flag_error = 1
        sys.exit()

    for prvek in gnss_coord:
        station_coord = prvek[0]
        coord_station = prvek[1:]
        coord_station = map(float, coord_station)
        if len(os.listdir(input_RINOBS)) >0:
            for soubor_rinex in os.listdir(input_RINOBS):
                if station_coord.lower() == soubor_rinex[0:4].lower():
                    rinex_list = fReadObservationRinex(input_RINOBS+'/'+soubor_rinex, gnss_list, interval_slants)
                    if len(rinex_list) > 1: #tests if there are observations for selected station
                        azel_file_name = soubor_rinex[0:4]+rinex_list[0][0]+rinex_list[0][1].rjust(2, '0')+rinex_list[0][2].rjust(2, '0')+'.AZEL'
                        zapis = open(temp_AZEL+'/'+azel_file_name, 'w')
                        if len(sp3_list) == 0:
                            for soubor_sp3 in os.listdir(input_EPH):
                                if '.gz' in soubor_sp3 or '.GZ' in soubor_sp3:
                                    with gzip.open(input_EPH+'/'+soubor_sp3, 'r') as f:
                                        for line in f:
                                            datum_sp3 = line.split()[0:3]
                                            break
                                else:
                                    cti_sp3 = open(input_EPH+'/'+soubor_sp3, 'r')
                                    first_line_sp3 = cti_sp3.readline()
                                    cti_sp3.close()
                                    datum_sp3 = first_line_sp3.split()[0:3]
                                if int(datum_sp3[0][-2:]) == int(rinex_list[0][0]) and int(datum_sp3[1]) == int(rinex_list[0][1]) and int(datum_sp3[2]) == int(rinex_list[0][2]):
            ##                        print soubor_sp3
                                    sp3_list = fReadSP3(input_EPH+'/'+soubor_sp3)
                ##                    print sp3_list
                                    break
                        if len(sp3_list) != 0:
                            for epoch_list in rinex_list[1:]:
                                epoch = epoch_list[0]
                                for satellite in epoch_list[1:]:
                                    coord_sat = fInterpolateSatelliteXYZ(sp3_list, satellite, epoch)
                                    if coord_sat != -1: #solve a problem when satellite ephemerides for needed epoch are not present in the SP3 file
                                        obs_azi_ele = fComputeAziEle(coord_station, coord_sat)
                                        if obs_azi_ele[1] >= cut_off:
                                           zapis.write(str(epoch).ljust(5, ' ') + satellite.rjust(4, ' ') + str("{0:.3f}".format(coord_sat[0])).rjust(15,' ') + str("{0:.3f}".format(coord_sat[1])).rjust(15,' ') +
                                           str("{0:.3f}".format(coord_sat[2])).rjust(15,' ') + str("{0:.5f}".format(obs_azi_ele[0])).rjust(12,' ') + str("{0:.5f}".format(obs_azi_ele[1])).rjust(12,' ') + '\n')
                        else:
                            print '-----------------------------------------------------------------------------------'
                            print 'Error in function '+compute_AZEL.__name__ + 'An SP3 file for processed date', str(rinex_list[0]) ,'is missing.'
                            print '-----------------------------------------------------------------------------------'
                            flag_error = 1
                            sys.exit()
                        zapis.close()
                        break
        else:
            print '-----------------------------------------------------------------------------------'
            print 'Error in function '+compute_AZEL.__name__ + ': no file with GNSS RINEX observations'
            print '-----------------------------------------------------------------------------------'
            flag_error = 1
            sys.exit()
    if flag_error == 0:
        print '--------------------------------------------------------------------------------------------------------------------------'
        print "computation of AZI/ELE angles for observations from observation RINEX files was successful - function  "+temp_AZEL
        print '--------------------------------------------------------------------------------------------------------------------------'
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#vypocet parametru a, b, c WET mapovaci funkce Niell
def compute_param_WET_niell(zem_sirka):
    pole_lat = {15:15, 30:30, 45:45, 60:60, 75:75}
    #matice hodnot parametru a, b, c pro WET Niell funkci prevzata z Niell, 1996
    matice_a_wet = {15:5.8021897e-4, 30:5.6794847e-4, 45:5.8118019e-4, 60:5.9727542e-4, 75:6.1641693e-4}
    matice_b_wet = {15:1.4275268e-3, 30:1.5138625e-3, 45:1.4572752e-3, 60:1.5007428e-3, 75:1.7599082e-3}
    matice_c_wet = {15:4.3472961e-2, 30:4.6729510e-2, 45:4.3908931e-2, 60:4.4626982e-2, 75:5.4736038e-2}

    #wet mf:
    lat_min = 15
    for klic, hodnota in pole_lat.iteritems():
        rozdil_lat = hodnota - zem_sirka
        if rozdil_lat <= 0 and abs(rozdil_lat) <lat_min:
            lat_min = rozdil_lat
            memKey = klic

    if zem_sirka >= 15.0 and zem_sirka <=75.0:
        a_wet = (zem_sirka - float(memKey))/15.0 * (matice_a_wet.get(memKey+15)-matice_a_wet.get(memKey)) + matice_a_wet.get(memKey)
        b_wet = (zem_sirka - float(memKey))/15.0 * (matice_b_wet.get(memKey+15)-matice_b_wet.get(memKey)) + matice_b_wet.get(memKey)
        c_wet = (zem_sirka - float(memKey))/15.0 * (matice_c_wet.get(memKey+15)-matice_c_wet.get(memKey)) + matice_c_wet.get(memKey)
    else:
        print "Parameters of Wet Niell mapping function can't be computed for given latitude of:" + str(zem_sirka)
        return

    return a_wet, b_wet, c_wet

#-------------------------------------------------------------------------------
#vypocet parametru a, b, c WET mapovaci funkce GMF (prepsano z fortran dle http://ggosatm.hg.tuwien.ac.at/DELAY/SOURCE/gmf.f)
def compute_param_WET_GMF(zem_sirka, zem_delka, doy):
    #pole parametru prevzaty z gmf.f z http://ggosatm.hg.tuwien.ac.at/DELAY/SOURCE/
    aw_mean = (56.4,1.555,-1.011,-3.975,.03171,.1065, .6175,.1376,.04229,.003028,1.688,-.1692,.05478,.02473,6.059e-4,
    2.278,.006614,-3.505e-4,-.006697,8.402e-4,7.033e-4,-3.236,.2184, -.04611,-.01613,-.001604,5.42e-5,7.922e-5,-.2711,+
    -.4406,-.03376,-.002801,-4.09e-4,-2.056e-5,6.894e-6,2.317e-6,1.941,-.2562,.01598,.005449,3.544e-4,1.148e-5,7.503e-6,
    -5.667e-7,-3.66e-8,.8683,-.05931,-.001864,-1.277e-4,2.029e-4,1.269e-5,1.629e-6,9.66e-8,-1.015e-7,-5e-10)
    bw_mean = (0.,0.,.2592,0.,.02974,-.5471,0., -.5926,-.103,-.01567,0.,.171,.09025,.02689,.002243,0.,.3439,
    .02402,.00541,.001601,9.669e-5,0.,.09502,-.03063,-.001055,-1.067e-4,-1.13e-4,2.124e-5,0.,-.3129,.008463,
    2.253e-4,7.413e-5,-9.376e-5,-1.606e-6,2.06e-6,0.,.2739,.001167,-2.246e-5,-1.287e-4,-2.438e-5,-7.561e-7,
    1.158e-6,4.95e-8,0.,-.1344,.005342,3.775e-4,-6.756e-5,-1.686e-6,-1.184e-6,2.768e-7,2.73e-8,5.7e-9)
    aw_amp = (.1023,-2.695,.3417,-.1405,.3175,.2116,3.536,-.1505,-.0166,.02967,.3819,-.1695,-.07444,.007409,
    -.006262,-1.836,-.01759,-.06256,-.002371,7.947e-4,1.501e-4,-.8603,-.136,-.03629,-.003706,-2.976e-4,
    1.857e-5,3.021e-5,2.248,-.1178,.01255,.001134,-2.161e-4,-5.817e-6,8.836e-7,-1.769e-7,.7313,-.1188,
    .01145,.001011,1.083e-4,2.57e-6,-2.14e-6,-5.71e-8,2e-8,-1.632,-.006948,-.003893,8.592e-4,7.577e-5,
    4.539e-6,-3.852e-7,-2.213e-7,-1.37e-8,5.8e-9)
    bw_amp = (0.,0.,-.08865,0.,-.4309,.0634,0.,.1162,.06176,-.004234,0.,.253,.04017,-.006204,.004977,0.,-.1737,
    -.005638,1.488e-4,4.857e-4,-1.809e-4,0.,-.1514,-.01685,.005333,-7.611e-5,2.394e-5,8.195e-6,0.,.09326,-.01275,
    -3.071e-4,5.374e-5,-3.391e-5,-7.436e-6,6.747e-7,0.,-.08637,-.003807,-6.833e-4,-3.861e-5,-2.268e-5,1.454e-6,
    3.86e-7,-1.068e-7,0.,-.02658,-.001947,7.131e-4,-3.506e-5,1.885e-7,5.792e-7,3.99e-8,2e-8,-5.7e-9)

    zem_sirka = math.radians(zem_sirka)
    zem_delka = math.radians(zem_delka)
    doy = float(doy) + 0.5 #pricitam 0.5, aby se hodnoty pocitaly pro poledne zpracovavaneho dne,
    #tedy pro stred zpracovavaneho obdobi (24 hodin)
    t = math.sin(zem_sirka)

    #parametry b a c se nemeni, dopocitava se pouze hodnota parametru a
    b_wet = 0.00146
    c_wet = 0.04391

    #pomocne promenne
    m = 9 #order m
    n = 9 #degree n
    dfac = []
    P = [None] * 100
    aP = []
    bP = []
    awm = 0.0
    awa = 0.0

    #determine n!  (faktorial)  moved by 1
    dfac.append(1.0)
    for i in range(0,(2*n)+1):
        dfac_i = dfac[i]*(i+1)
        dfac.append(dfac_i)

    #determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
    for i in range(0,n+1):
        for j in range(0,min(i,m)+1):
            ir = int((i - j)/2)
            sum = numpy.longdouble(0.0)
            for k in range(0,ir+1):
                sum = (sum + (-1.0)**k*dfac[2*i - 2*k]/dfac[k]/dfac[i - k]/dfac[i - j - 2*k]*pow(t,(i - j - 2*k)))
            p = ((1.0/2**i*math.sqrt((1 - t**2)**(j))*sum))
            P[i + 1 + (j+1) * 10 - 11] = (p)

    #spherical harmonics
    for j in range(0,10):
        for k in range(0,j+1):
            aP_i = ((P[j + 1 + (k + 1) * 10 - 11])*(math.cos(k*zem_delka)))
            aP.append(aP_i)

            bP_i = (P[j + 1 + (k + 1) * 10 - 11]*(math.sin(k*zem_delka)))
            bP.append(bP_i)

    for i in range(0,55):
        awm = (awm + (aw_mean[i]*aP[i] + bw_mean[i]*bP[i])*1e-5)
        awa = (awa + (aw_amp[i] *aP[i] + bw_amp[i]*bP[i])*1e-5)

    a_wet = awm + awa*(math.cos(2*math.pi*((doy-28)/365.25)))
    return a_wet, b_wet, c_wet

#-------------------------------------------------------------------------------
#vypocet WET mapovaci funkce (Niell ci GMF), navazuje na compute_param_WET_xxx
def compute_mf_WET(a_wet, b_wet, c_wet, elevace):
    sinus_el = math.sin(math.radians(float(elevace)))
    m_wet = (1 + (a_wet/(1 + (b_wet/(1 + c_wet)))))/(sinus_el + (a_wet/(sinus_el + (b_wet/(sinus_el + c_wet)))))
    return m_wet

#-------------------------------------------------------------------------------
#vypocet parametru a, b, c DRY mapovaci funkce Niell
def compute_param_DRY_Niell(zem_sirka,doy):
    pole_lat = {15:15, 30:30, 45:45, 60:60, 75:75}
    #matice hodnot parametru a, b, c pro DRY Niell mapovaci funkci, prevzato z Niell, 1996
    #average
    matice_a_dry_avg = {15:1.2769934e-3, 30:1.2683230e-3, 45:1.2465397e-3, 60:1.2196049e-3, 75:1.2045996e-3}
    matice_b_dry_avg = {15:2.9153695e-3, 30:2.9152299e-3, 45:2.9288445e-3, 60:2.9022565e-3, 75:2.9024912e-3}
    matice_c_dry_avg = {15:62.610505e-3, 30:62.837393e-3, 45:63.721774e-3, 60:63.824265e-3, 75:64.258455e-3}
    #amplitude
    matice_a_dry_amp = {15:0.0, 30:1.2709626e-5, 45:2.6523662e-5, 60:3.4000452e-5, 75:4.1202191e-5}
    matice_b_dry_amp = {15:0.0, 30:2.1414979e-5, 45:3.0160779e-5, 60:7.2562722e-5, 75:11.723375e-5}
    matice_c_dry_amp = {15:0.0, 30:9.0128400e-5, 45:4.3497037e-5, 60:84.795348e-5, 75:170.37206e-5}

    doy = float(doy) + 0.5 #pricitam 0.5, aby se hodnoty pocitaly pro poledne zpracovavaneho dne,
    #tedy pro stred zpracovavaneho obdobi (24 hodin)
    lat_min = 15

    #hledani odpovidajici zem_sirky z tabulky
    for klic, hodnota in pole_lat.iteritems():
        rozdil_lat = hodnota - zem_sirka
        if rozdil_lat <= 0 and abs(rozdil_lat) <lat_min:
            lat_min = rozdil_lat
            memKey = klic

    #interpolace pro jednotlive parametry pro konkretni zadanou zem_sirku
    #na zaklade udaju v maticich vyse (15, 30, 45, 60, 75 stupnu)
    if zem_sirka >= 15.0 and zem_sirka <=75.0:
        a_dry_avg = (zem_sirka - float(memKey))/15.0 * (matice_a_dry_avg.get(memKey+15)-matice_a_dry_avg.get(memKey)) + matice_a_dry_avg.get(memKey)
        b_dry_avg = (zem_sirka - float(memKey))/15.0 * (matice_b_dry_avg.get(memKey+15)-matice_b_dry_avg.get(memKey)) + matice_b_dry_avg.get(memKey)
        c_dry_avg = (zem_sirka - float(memKey))/15.0 * (matice_c_dry_avg.get(memKey+15)-matice_c_dry_avg.get(memKey)) + matice_c_dry_avg.get(memKey)

        a_dry_amp = (zem_sirka - float(memKey))/15.0 * (matice_a_dry_amp.get(memKey+15)-matice_a_dry_amp.get(memKey)) + matice_a_dry_amp.get(memKey)
        b_dry_amp = (zem_sirka - float(memKey))/15.0 * (matice_b_dry_amp.get(memKey+15)-matice_b_dry_amp.get(memKey)) + matice_b_dry_amp.get(memKey)
        c_dry_amp = (zem_sirka - float(memKey))/15.0 * (matice_c_dry_amp.get(memKey+15)-matice_c_dry_amp.get(memKey)) + matice_c_dry_amp.get(memKey)
    else:
        print "Parameters of Dry Niell mapping function can't be computed for given latitude of:"
        return

    a_dry = a_dry_avg - (a_dry_amp * (math.cos(2*math.pi*((doy-28)/365.25))))
    b_dry = b_dry_avg - (c_dry_amp * (math.cos(2*math.pi*((doy-28)/365.25))))
    c_dry = b_dry_avg - (c_dry_amp * (math.cos(2*math.pi*((doy-28)/365.25))))

    return a_dry, b_dry, c_dry

#-------------------------------------------------------------------------------
#vypocet parametru a, b, c DRY mapovaci funkce GMF (prepsano z fortran dle http://ggosatm.hg.tuwien.ac.at/DELAY/SOURCE/gmf.f)
def compute_param_DRY_GMF(zem_sirka, zem_delka, doy):
    #pole parametru prevzaty z http://ggosatm.hg.tuwien.ac.at/DELAY/SOURCE/gmf.f
    ah_mean = (125.17,.8503,.06936,-6.76,.1771,.0113, .5963,.01808,.002801,-.001414,-1.212,.093,.003683,.001095,
    4.671e-5,.3959,-.03867,.005413,-5.289e-4,3.229e-4,2.067e-5,.3, .02031,.0059,4.573e-4,-7.619e-5,2.327e-6,3.845e-6,.1182,.01158,
    .005445,6.219e-5,4.204e-6,-2.093e-6,1.54e-7,-4.28e-8,-.4751, -.0349,.001758,4.019e-4,-2.799e-6,-1.287e-6,5.468e-7,7.58e-8,
    -6.3e-9,-.116,.008301,8.771e-4,9.955e-5,-1.718e-6,-2.012e-6, 1.17e-8,1.79e-8,-1.3e-9,1e-10)

    bh_mean= (0.,0.,.03249,0.,.03324,.0185,0.,-.1115,.02519,.004923,0.,.02737,.01595,-7.332e-4,1.933e-4,0.,
    -.04796,.006381,-1.599e-4,-3.685e-4,1.815e-5,0.,.07033,.002426,-.001111,-1.357e-4,-7.828e-6,2.547e-6,0.,.005779,.003133,
    -5.312e-4,-2.028e-5,2.323e-7,-9.1e-8,-1.65e-8,0.,.03688,-8.638e-4,-8.514e-5,-2.828e-5,5.403e-7,4.39e-7,1.35e-8,1.8e-9,0.,
    -.02736,-2.977e-4,8.113e-5,2.329e-7,8.451e-7,4.49e-8,-8.1e-9,-1.5e-9,2e-10)

    ah_amp = (-.2738,-2.837,.01298,-.3588,.02413,.03427,-.7624,.07272,.0216,-.003385,.4424,.03722,.02195,-.001503,
    2.426e-4,.3013,.05762,.01019,-4.476e-4,6.79e-5,3.227e-5,.3123,-.03535,.00484,3.025e-6,-4.363e-5,2.854e-7,-1.286e-6,-.6725,
    -.0373,8.964e-4,1.399e-4,-3.99e-6,7.431e-6,-2.796e-7,-1.601e-7,.04068,-.01352,7.282e-4,9.594e-5,2.07e-6,-9.62e-8,-2.742e-7,
    -6.37e-8,-6.3e-9,.08625,-.005971,4.705e-4,2.335e-5,4.226e-6,2.475e-7,-8.85e-8,-3.6e-8,-2.9e-9,0.)

    bh_amp = (0.,0.,-.1136,0.,-.1868,-.01399,0.,-.1043,.01175,-.00224,0.,-.03222,.01333,-.002647,-2.316e-5,0.,
    .05339,.01107,-.003116,-1.079e-4,-1.299e-5,0.,.004861,.008891,-6.448e-4,-1.279e-5,6.358e-6,-1.417e-7,0.,.03041,.00115,-8.743e-4,
    -2.781e-5,6.367e-7,-1.14e-8,-4.2e-8,0.,-.02982,-.003,1.394e-5,-3.29e-5,-1.705e-7,7.44e-8,2.72e-8,-6.6e-9,0.,.01236,-9.981e-4,
    -3.792e-5,-1.355e-5,1.162e-6,-1.789e-7,1.47e-8,-2.4e-9,-4e-10)

    zem_sirka = math.radians(zem_sirka)
    zem_delka = math.radians(zem_delka)
    doy = float(doy) + 0.5 #pricitam 0.5, aby se hodnoty pocitaly pro poledne zpracovavaneho dne,
    #tedy pro stred zpracovavaneho obdobi (24 hodin)
    t = math.sin(zem_sirka)

    #pomocne promenne
    m = 9 #order m
    n = 9 #degree n
    dfac = []
    P = [None] * 100
    aP = []
    bP = []

    #determine n!  (faktorial)  moved by 1
    dfac.append(1.0)
    for i in range(0,(2*n)+1):
        dfac_i = dfac[i]*(i+1)
        dfac.append(dfac_i)

    #determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
    for i in range(0,n+1):
        for j in range(0,min(i,m)+1):
            ir = int((i - j)/2)
            sum = numpy.longdouble(0.0)
            for k in range(0,ir+1):
                sum = (sum + (-1.0)**k*dfac[2*i - 2*k]/dfac[k]/dfac[i - k]/dfac[i - j - 2*k]*pow(t,(i - j - 2*k)))
            p = ((1.0/2**i*math.sqrt((1 - t**2)**(j))*sum))
            P[i + 1 + (j+1) * 10 - 11] = (p)

    #spherical harmonics
    for j in range(0,10):
        for k in range(0,j+1):
            aP_i = ((P[j + 1 + (k + 1) * 10 - 11])*(math.cos(k*zem_delka)))
            aP.append(aP_i)

            bP_i = (P[j + 1 + (k + 1) * 10 - 11]*(math.sin(k*zem_delka)))
            bP.append(bP_i)

    b_dry = 0.0029
    c0dry = 0.062
    if zem_sirka < 0:   #southern hemisphere
        phh  = math.pi
        c11h = 0.007
        c10h = 0.002
    else:               #northern hemisphere
        phh  = 0
        c11h = 0.005
        c10h = 0.001
    c_dry = c0dry + ((math.cos(2*math.pi*((doy-28)/365.25) + phh)+1)*c11h/2 + c10h)*(1-math.cos(zem_sirka))

    ahm = 0.0
    aha = 0.0
    for i in range(0,55):
        ahm = (ahm + (ah_mean[i]*aP[i] + bh_mean[i]*bP[i])*1e-5)
        aha = (aha + (ah_amp[i] *aP[i] + bh_amp[i]*bP[i])*1e-5)

    a_dry = ahm + aha*(math.cos(2*math.pi*((doy-28)/365.25)))

    return a_dry, b_dry, c_dry

#-------------------------------------------------------------------------------
#vypocet DRY mapovaci funkce (Niell ci GMF), navazuje na funkci compute_param_DRY_xxx
def compute_mf_DRY(a_dry, b_dry, c_dry, elevace, H_vyska):
    #height correction
    a_height = 2.53e-5
    b_height = 5.49e-3
    c_height = 1.14e-3
    H_vyska = float(H_vyska)/1000 #prepocet na km

    sinus_el = math.sin(math.radians(float(elevace)))
    m_dry_apriori = (1 + (a_dry/(1 + (b_dry/(1 + c_dry)))))/(sinus_el + (a_dry/(sinus_el + (b_dry/(sinus_el + c_dry)))))
    delta_h = ((1 / sinus_el) - (1 + (a_height/(1 + (b_height/(1 + c_height)))))/(sinus_el + (a_height/(sinus_el + (b_height/(sinus_el + c_height)))))) * H_vyska
    m_dry = m_dry_apriori + delta_h

    return m_dry
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#funkce pro bilinearni interpolaci = pouziva se pri praci s VMF1 gridy
def bilinear_interpolation(x, y, points):
    '''Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)
#-------------------------------------------------------------------------------

#obtain LAT/LON grid size of used VMF1 grid files
#-------------------------------------------------------------------------------
def fVMF1_Get_Grid_Size(VMF1_grid):
    if len(os.listdir(VMF1_grid)) > 0:
        for file1 in os.listdir(VMF1_grid):
            if '.gz' in file1 or '.GZ' in file1:
                with gzip.open(VMF1_grid+'/'+file1, 'r') as f:
                    radky = f.readlines()
            else:
                cti = open(VMF1_grid+'/'+file1, 'r')
                radky = cti.readlines()
                cti.close()

            for radek in radky:
                if radek[0] == '!' and 'Range' in radek:
                    lat_size = float(radek.split(':')[1].split()[-2])
                    lon_size = float(radek.split(':')[1].split()[-1])
                    break
            break
        VMF1_grid_size = [lat_size, lon_size]
        return VMF1_grid_size
    else:
        print '------------------------------------------------------------------------------------------'
        print 'Error in function '+fVMF1_Get_Grid_Size.__name__ + ': no files with VMF1 grids found'
        print '------------------------------------------------------------------------------------------'
        sys.exit()
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#read content off VMF1 grid file(s) for processed day and return a_dry, a_wet
#parameters for a station position for all available epochs of VMF grid files
def fRead_vmf1_grid_file(VMF1_grid, VMF1_grid_size, zem_sirka, zem_delka, year_doy):
    list_time = []
    list_adry = []
    list_awet = []

    date_object = datetime.datetime.strptime(year_doy, '%y%j')
    year_day = str(date_object).split('-')[0][2:]+str(date_object).split('-')[1]+str(date_object).split('-')[2].split()[0]
    date_object_plus_day = date_object + datetime.timedelta(days=1)
    year_day_plus_day = str(date_object_plus_day).split('-')[0][2:]+str(date_object_plus_day).split('-')[1]+str(date_object_plus_day).split('-')[2].split()[0]

    #tvorba poli s rozsahy rohovych bodu bunek globalni site ve smeru lat a lon
    pole_lat = []
    lat_size = VMF1_grid_size[0]
    pole_lon = []
    lon_size = VMF1_grid_size[1]
    for i in range(0, 91):
        hodnota = -90.0 + i * lat_size
        pole_lat.append(hodnota)
    for i in range(0, 144):
        hodnota = 0 + i * lon_size
        pole_lon.append(hodnota)

    #hledani 4 rohovych bodu globalni site v okoli mista vypoctu (=souradnice stanice)
    index_lat_nizsi = min(range(len(pole_lat)), key=lambda i: abs(pole_lat[i]-zem_sirka))
    if pole_lat[index_lat_nizsi] > zem_sirka:
        index_lat_vyssi = index_lat_nizsi
        index_lat_nizsi = index_lat_nizsi - 1
    else:
        index_lat_vyssi = index_lat_nizsi + 1
    index_lon_nizsi = min(range(len(pole_lon)), key=lambda i: abs(pole_lon[i]-zem_delka))
    if pole_lon[index_lon_nizsi] > zem_delka:
        index_lon_vyssi = index_lon_nizsi
        index_lon_nizsi = index_lon_nizsi - 1
    else:
        index_lon_vyssi = index_lon_nizsi + 1
##    print pole_lat[index_lat_nizsi], pole_lat[index_lat_vyssi]
##    print pole_lon[index_lon_nizsi], pole_lon[index_lon_vyssi]

    for file1 in os.listdir(VMF1_grid):
        if year_day in file1 or (year_day_plus_day in file1 and '00' in file1.split('.')[1]):
            if '.gz' in file1 or '.GZ' in file1:
                with gzip.open(VMF1_grid+'/'+file1, 'r') as f:
                    radky = f.readlines()
            else:
                cti = open(VMF1_grid+'/'+file1, 'r')
                radky = cti.readlines()
                cti.close()

            indexy_epochy = []
            for i in range(0, len(radky)):
                if radky[i][0] == '!' and 'Epoch' in radky[i]:
                    indexy_epochy.append(i)
            indexy_epochy.append(len(radky))

            for i in range(0, len(indexy_epochy)-1):
                minLatminLon = [None] * 4
                minLatmaxLon = [None] * 4
                maxLatminLon = [None] * 4
                maxLatmaxLon = [None] * 4

                for radek in radky[indexy_epochy[i]:indexy_epochy[i+1]]:
                    if radek[0] == '!' and 'Epoch' in radek:
                        datum_cas_list = radek.split(':')[1].split()
                        datum_VMF1 = float(datum_cas_list[0][2:]+datum_cas_list[1]+datum_cas_list[2])
                        time_VMF1 = float(datum_cas_list[3])*100 + ((float(datum_cas_list[4])/60)*100)
                        if datum_VMF1 == float(year_day_plus_day) and time_VMF1 == 0.0:
                            datum_VMF1 = float(year_day)
                            time_VMF1 = 2400.0
                        list_time.append(time_VMF1)

                    elif '!' not in radek:
                        radek_pole = radek.split()
                        if float(radek[0:5]) == pole_lat[index_lat_nizsi] and float(radek[6:11]) == pole_lon[index_lon_nizsi]:
                            minLatminLon[0] = float(radek_pole[0])
                            minLatminLon[1] = float(radek_pole[1])
                            minLatminLon[2] = float(radek_pole[2])
                            minLatminLon[3] = float(radek_pole[3])
                        if float(radek[0:5]) == pole_lat[index_lat_nizsi] and float(radek[6:11]) == pole_lon[index_lon_vyssi]:
                            minLatmaxLon[0] = float(radek_pole[0])
                            minLatmaxLon[1] = float(radek_pole[1])
                            minLatmaxLon[2] = float(radek_pole[2])
                            minLatmaxLon[3] = float(radek_pole[3])
                        if float(radek[0:5]) == pole_lat[index_lat_vyssi] and float(radek[6:11]) == pole_lon[index_lon_nizsi]:
                            maxLatminLon[0] = float(radek_pole[0])
                            maxLatminLon[1] = float(radek_pole[1])
                            maxLatminLon[2] = float(radek_pole[2])
                            maxLatminLon[3] = float(radek_pole[3])
                        if float(radek[0:5]) == pole_lat[index_lat_vyssi] and float(radek[6:11]) == pole_lon[index_lon_vyssi]:
                            maxLatmaxLon[0] = float(radek_pole[0])
                            maxLatmaxLon[1] = float(radek_pole[1])
                            maxLatmaxLon[2] = float(radek_pole[2])
                            maxLatmaxLon[3] = float(radek_pole[3])

                input_dry = minLatminLon[0:3], minLatmaxLon[0:3], maxLatminLon[0:3], maxLatmaxLon[0:3]
                a_dry = bilinear_interpolation(zem_sirka,zem_delka,input_dry)
                minLatminLon.pop(2), minLatmaxLon.pop(2), maxLatminLon.pop(2), maxLatmaxLon.pop(2)
                input_wet = minLatminLon, minLatmaxLon, maxLatminLon, maxLatmaxLon
                a_wet = bilinear_interpolation(zem_sirka,zem_delka,input_wet)
                list_adry.append(a_dry)
                list_awet.append(a_wet)
    #sort values from starting epoch (0.0) to ending epoch (2400.0)
    xyz = zip(list_time, list_adry, list_awet)
    xyz_sorted = sorted(xyz)
    list_time = [item[0] for item in xyz_sorted]
    list_adry = [item[1] for item in xyz_sorted]
    list_awet = [item[2] for item in xyz_sorted]
    return [list_time, list_adry, list_awet]
#-------------------------------------------------------------------------------

#compute c_dry parameter of VMF1 function depending only on LAT and DOY
#-------------------------------------------------------------------------------
def fCompute_VMF1_c_dry_parameter(zem_sirka, doy):
    c0dry = 0.062
    if zem_sirka < 0:   #southern hemisphere
        phh  = math.pi
        c11h = 0.007
        c10h = 0.002
    else:               #northern hemisphere
        phh  = 0
        c11h = 0.005
        c10h = 0.001
    c_dry = c0dry + ((math.cos(2*math.pi*((doy-28)/365.25) + phh)+1)*c11h/2 + c10h)*(1-math.cos(zem_sirka))
    return c_dry

#-------------------------------------------------------------------------------
#compute parameters of VMF1 gridded function, a linear interpolation in time is performed
#uses results of previous functions fRead_vmf1_grid_file, fCompute_VMF1_c_dry_parameter
def compute_param_VMF1(cas, VMF1_a_parameters_all_epochs, c_dry):
    cas = float(cas)
    #osetreni vstupniho casu, aby nepresahl 23:59:59
    if cas == 2400.0:
        cas = 2399.999

    pole_cas_grid = VMF1_a_parameters_all_epochs[0]
    #hledani dvou nejblizsich epoch se zaznamy pro zadanou
    #observaci = index predchazejici epochy se ulozi do index_cas_nizsi a index
    #v case nasledujici epochy do promenne index_cas_vyssi
    #jinak receno - je potreba najit mezi ktere dve epochy hodnot GRIDU spada
    #cas observace, aby se dalo mezi temito dvema epochami v case interpolovat
    index_cas_nizsi = min(range(len(pole_cas_grid)), key=lambda i: abs(pole_cas_grid[i]-cas))
    if pole_cas_grid[index_cas_nizsi] > cas:
        index_cas_vyssi = index_cas_nizsi
        index_cas_nizsi =  index_cas_nizsi - 1
    else:
        index_cas_vyssi = index_cas_nizsi + 1
##    print index_cas_nizsi, index_cas_vyssi
##    print index_radky_datum[index_cas_nizsi], index_radky_datum[index_cas_vyssi]

    a_dry1 = VMF1_a_parameters_all_epochs[1][index_cas_nizsi]
    a_dry2 = VMF1_a_parameters_all_epochs[1][index_cas_vyssi]
    a_wet1 = VMF1_a_parameters_all_epochs[2][index_cas_nizsi]
    a_wet2 = VMF1_a_parameters_all_epochs[2][index_cas_vyssi]
    #linearni interpolace v case pro hodnoty a_dry/wet1 a a_dry/wet2 platne pro predchazejici a
    #nasledujici epochu zaznamu v GRID file - dela se pro cas observace
    a_dry = a_dry1 + (a_dry2 - a_dry1)*((cas-pole_cas_grid[index_cas_nizsi])/(pole_cas_grid[index_cas_vyssi]-pole_cas_grid[index_cas_nizsi]))
    a_wet = a_wet1 + (a_wet2 - a_wet1)*((cas-pole_cas_grid[index_cas_nizsi])/(pole_cas_grid[index_cas_vyssi]-pole_cas_grid[index_cas_nizsi]))
##    print a_dry, a_wet
    #dry - parametr b se nemeni, c zavisi na zem. sirce, a se dopocitava
    b_dry = 0.0029
    #wet - parametry b a c se nemeni, dopocitava se pouze hodnota parametru a
    b_wet = 0.00146
    c_wet = 0.04391
    return a_dry, b_dry, c_dry, a_wet, b_wet, c_wet
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#create a list of observations from observation RINEX file with this structure:
#[[epoch in seconds, PRN, PRN, PRN, ...],...]
#input parameters = path to rinex file, list of wanted GNSS ('G','R','E',...),
#wanted time interval of observations in seconds = use the interval in which
#you want to deliver slants (i.e. 150 seconds)
def fReadObservationRinex(rinex_input_path, gnss_list, interval):
    observations_list = []

    cti = open(rinex_input_path, 'r')
    lines = cti.readlines()
    cti.close()

    for i in range(0, len(lines)):
        if 'TIME OF FIRST OBS' in lines[i]:
            index_date = i
        if 'END OF HEADER' in lines[i]:
            index_end_header = i
            break
    date = [lines[index_date].split()[0][2:], lines[index_date].split()[1], lines[index_date].split()[2]]
    observations_list.append(date)

    for i in range(index_end_header+1, len(lines)):
        if len(lines[i].split()) > 3:
            line_list = lines[i].split()
            if line_list[0].rjust(2, '0') == date[0].rjust(2, '0') and line_list[1].rjust(2, '0') == date[1].rjust(2, '0') and line_list[2].rjust(2, '0') == date[2].rjust(2, '0'):
                epoch_list = []
                epoch = int(int(line_list[3])*3600 + int(line_list[4])*60 + float(line_list[5]))
                if float(epoch)% interval == 0:
                    epoch_list.append(epoch)
                    list_satellites = []

                    #it is necessary to find start and end index of list of satellites since in some RINEXes PRN codes G 1, G 2, G 3, ..., G11, .. are used and in others codes
                    #G01, G02, G03, ..., G11, G12 are used
                    radek = lines[i]
                    for j in range(0, len(radek)):
                        if radek[j].isalpha() == True:
                            index_radek_start = j
                            break
                    for j in range(index_radek_start, len(radek)):
                        if radek[j].isalpha() == True:
                            index_radek_end = j
                    prvek = radek[index_radek_start:index_radek_end+3].replace(' ','0')
                    list_satellites_input = prvek.replace('G','_G').replace('R','_R').replace('S','_S').replace('E','_E').replace('J','_J')
                    list_satellites = list_satellites_input.split('_')[1:]

                    #solves problem where GNSS abbreviation is not given in RINEX file, the observation is considered as from GPS satellite
                    if len(list_satellites) < 1:
                        if len(line_list)-8 > 1:
                            for prvek in line_list[8:]:
                                sat = 'G'+prvek
                                list_satellites.append(sat)
                    for prvek in list_satellites[:]:
                        if prvek[0] not in gnss_list:
                            list_satellites.remove(prvek)

                    list_satellites2 = []
                    for k in range(i+1, i+5):
                        if 'G' in lines[k] or 'R' in lines[k] or 'E' in lines[k]:
                            radek = lines[k]
                            for j in range(0, len(radek)):
                                if radek[j].isalpha() == True:
                                    index_radek_start = j
                                    break
                            for j in range(index_radek_start, len(radek)):
                                if radek[j].isalpha() == True:
                                    index_radek_end = j
                            prvek = radek[index_radek_start:index_radek_end+3].replace(' ','0')
                            list_satellites_input2 = prvek.replace('G','_G').replace('R','_R').replace('S','_S').replace('E','_E').replace('J','_J')
                            list_satellites2 = list_satellites2 + list_satellites_input2.split('_')[1:]
                        for prvek in list_satellites2[:]:
                            if prvek[0] not in gnss_list:
                                list_satellites2.remove(prvek)

                    epoch_list = epoch_list + sorted(list_satellites+list_satellites2)
                    observations_list.append(epoch_list)
    return observations_list
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def fComputeAziEle(receiver, satellite): #compute azimuth and elevation angles for given coordinats of receiver and satellite,
    #coordinates are given in XYZ format as a Python list = [X,Y,Z]
    #https://www.mathworks.com/matlabcentral/fileexchange/41364-gps-navigation-toolbox/content/GNT08.1.2/final-GPS/Calc_Azimuth_Elevation.m
    #http://www.navipedia.net/index.php/Transformations_between_ECEF_and_ENU_coordinates
    vector =[satellite[0]-receiver[0],satellite[1]-receiver[1],satellite[2]-receiver[2]]
    vector_size = math.sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2])
    unit_vector = numpy.matrix([vector[0]/vector_size, vector[1]/vector_size, vector[2]/vector_size])
    unit_vector_T = unit_vector.T

    rec_latlonH = XYZ_to_LatLonH(receiver[0],receiver[1],receiver[2])
    rec_lat = math.radians(rec_latlonH[0])
    rec_lon = math.radians(rec_latlonH[1])
    ENU = numpy.matrix([[-1*math.sin(rec_lon), math.cos(rec_lon), 0], [-1*math.sin(rec_lat)*math.cos(rec_lon), -1*math.sin(rec_lat)*math.sin(rec_lon), math.cos(rec_lat)],
    [math.cos(rec_lat)*math.cos(rec_lon), math.cos(rec_lat)*math.sin(rec_lon), math.sin(rec_lat)]])

    elevation = math.degrees(math.asin(float(ENU[2]*unit_vector_T)))
    azimuth_first = math.degrees(math.atan(float(ENU[0]*unit_vector_T)/float(ENU[1]*unit_vector_T)))
    #azimuth correction to place it in a correct quadrant
    if (float(ENU[0]*unit_vector_T) > 0 and float(ENU[1]*unit_vector_T) < 0) or (float(ENU[0]*unit_vector_T) < 0 and float(ENU[1]*unit_vector_T) < 0):
        azi_correction = 180
    elif float(ENU[0]*unit_vector_T) < 0 and float(ENU[1]*unit_vector_T) > 0:
        azi_correction = 360
    else:
        azi_correction = 0

    azimuth = azimuth_first + azi_correction
    return [azimuth, elevation]
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
def fAdd_Residuals_toTROSINEX(input_RES, result_trosnx, interval_slants):
    print '---------------------------------------------------------------------------------------------------------------------'
    print "Adding of post-fit residuals to TRO-SINEX files was started - function "+fAdd_Residuals_toTROSINEX.__name__
    print '---------------------------------------------------------------------------------------------------------------------'
    pocet_stanic = 0
    for soubor in os.listdir(result_trosnx):
        if 'SNX' in soubor:
            station_name_snx = soubor[0:4]
            doy_snx = soubor.split('.')[0][-3:]
##            print station_name_snx, doy_snx
            #find RESIDUAL file corresponding to TROSNX file
            for soubor_res in os.listdir(input_RES):
                doy_res = soubor_res[0:3]
                if doy_snx == doy_res and station_name_snx in soubor_res:
                    print 'processed TRO-SINEX file is: ',soubor,', processed RESIDUAL file is:', soubor_res
                    pocet_stanic = pocet_stanic + 1
                    #read TROSNX
                    cti_snx = open(result_trosnx+'/'+soubor, 'r')
                    radky_snx = cti_snx.readlines()
                    cti_snx.close()

                    #read only lines with reasiduals which correspond to selected interval for slants
                    radky_res = []
                    if '.gz' in soubor_res or '.GZ' in soubor_res or '.z' in soubor_res or '.Z' in soubor_res or '.zip' in soubor_res or '.ZIP' in soubor_res:
                        with gzip.open(input_RES+'/'+soubor_res, 'r') as f:
                            for radek_res in f:
                                if float(radek_res.split()[0]) % interval_slants == 0:
                                    radky_res.append(radek_res)
                    else:
                        with open(input_RES+'/'+soubor_res, 'r') as f:
                            for radek_res in f:
                                if float(radek_res.split()[0]) % interval_slants == 0:
                                    radky_res.append(radek_res)

                    #create output SNX file
                    zapis_snx = open(result_trosnx+'/'+soubor, 'w')
                    #detect start of slants part in input TROSNX
                    for radek in radky_snx:
                        if '+SLANT/SOLUTION' in radek:
                            myHeaderFlag = radky_snx.index(radek)
                            break
                    #write zenith part of TROSNX into results file
                    zapis_snx.writelines(radky_snx[0:myHeaderFlag+2])
                    #loop through slants in TROSNX
                    for i in range (myHeaderFlag+2,len(radky_snx)-2):
                        sat_snx_priprava = radky_snx[i].split()[10]
                        if sat_snx_priprava[0] == 'G':
                            sat_snx = str(int(sat_snx_priprava[1:]))
                        elif sat_snx_priprava[0] == 'R':
                            sat_snx = str(int(sat_snx_priprava[1:])+100)
                        seconds_snx = str(int(radky_snx[i].split()[1].split(':')[2]))
                        #loop through residuals in RES file
                        for radek_res in radky_res:
                            seconds_res = radek_res.split()[0]
                            sat_res = radek_res.split()[1]
                            #testing if satellite PRN and epoch are identical in slant and res line
                            if sat_snx == sat_res and seconds_snx == seconds_res:
                                residuum = float(radek_res.split()[4])
                                slant_with_res = float(radky_snx[i].split()[2]) + residuum
                                radky_snx[i] = radky_snx[i][0:18] + str("{0:.1f}".format(slant_with_res)).rjust(10,' ') + radky_snx[i][28:38] + str("{0:.1f}".format(residuum)).rjust(10,' ') + radky_snx[i][48:]
                                del radky_res[radky_res.index(radek_res)]
                                break
                    zapis_snx.writelines(radky_snx[myHeaderFlag+2:len(radky_snx)-2])
##                                        print len(radky_res)
                    zapis_snx.write('-SLANT/SOLUTION\n%=ENDTRO\n')
                    zapis_snx.close()
                    radky_snx = None
    print '---------------------------------------------------------------------------------------------------------------------'
    print "Adding of post-fit residuals to TRO-SINEX files ended - function "+fAdd_Residuals_toTROSINEX.__name__
    print "Number of processed stations is:"+str(pocet_stanic)
    print '---------------------------------------------------------------------------------------------------------------------'

#-------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
def read_config(ConfData,MyParameter):
  Values = []
  Flag = 0

  for line in ConfData:
    param = str.split(line,' ')
    name = str.strip(param[0])
    parValues = str.strip(line)

    if Flag == 1 and name != '[END]' and name[0:1] != '#':
      Values.append(parValues)
    elif Flag == 1 and name == '[END]':
      return Values
    if string.upper(name) == str.upper('[' + MyParameter + ']'):
      Flag = 1

  # Whole file was read but parameter with a given value wasn't found
  else:
    print(u'Non valid parameter name [' + MyParameter + u']!')
    sys.exit(u'Script run was untimely ended!')
# ------------------------------------------------------------------------------