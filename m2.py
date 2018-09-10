import os, sys
import datetime, time
import numpy as np
from numpy import *
from input_data import *
import math as m
import model
from azi_ele_coor import *

from dataSHMU import *
import reader
import ogr
import poloha_druzice

try:
    import matplotlib.pyplot as plt
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"

try:
    #ephfiles = np.sort(next(os.walk(eph_path))[2])
    obsfiles = np.sort(next(os.walk(obs_path))[2])

    doy_list1 = []
    doy_list2 = []
    SHMUdata = [[], []]
    vyska1 = []
    vyska2 = []
    vyska3 = []

    deg = 6   # stupen polynomu pri odstranovani trendu

    # datum zaciatku a konca vypoctu
    dat_start = datetime.date(rok_start, mes_start, den_start)
    dat_end = datetime.date(rok_konec, mes_konec, den_konec)
    delta = dat_end-dat_start
    pocet_dni = int((delta.total_seconds())/86400)

    c = 0
    while c < len(stanice):

        SHMU = fSHMU(c)

        for f in obsfiles:

            if 'o' != f.split('.')[1][-1]:
                continue

            o = obs_path + f # cesta k SP3 suboru
            print
            print o

            r_obs = open(o, "r")
            lines = r_obs.readlines()
            r_obs.close()

            #sp3, date_e = utility.fReadSP3(e)   # nacitanie presnych efemerid
            header, version, headlines, headlength, obstimes, sats, numallsvs, numsvs, date_o = reader.scan(lines)
            prva_obs = header['TIME OF FIRST OBS'].split()
            dat = datetime.date(int(prva_obs[0]), int(prva_obs[1]), int(prva_obs[2]))  # datum prvej observacie

            # print int(date_e[0]), int(date_e[1]), int(date_e[2])

            if dat > dat_end: # kontrola ci sa datum zaznamu sa nachadza v zadanom intervale datumu
                break
            if dat_start > dat:
                continue

            rok = prva_obs[0][2:]
            navfile = nav_path + f.split('.')[0] + '.' + rok + 'n'
            if not os.path.exists(navfile):
                print 'ziadne navigacne data'
                continue
            # else: continue

            print navfile
            print dat

            #parsrovanie observacnych dat
            data, obstypes = reader.processBlocks(lines, header, obstimes, headlines, headlength, sats, numallsvs,satellite)
            navdata = reader.rinexnav(navfile)
            navvalues = navdata.values  # navigacne data
            navtime = navdata['t'].to_pandas()  # vyber casu
            del lines

            indexS1 = obstypes.index('S1') + 2
            indexS2 = obstypes.index('S2') + 2
            indexS5 = obstypes.index('S5') + 2

            for sat in satellite:
                print ("pocitam pre: " + sat)
                t = []
                signals1 = []
                signals2 = []
                signals5 = []
                sinele = []
                a = []

                for k, d in enumerate(data):
                    # print d
                    # print k, sats[k]

                    if sat not in sats[k]:
                        continue
                    # print "pocitam pre: " + sat, k

                    for i in d:  # jedna sekunda pre x druzic
                        # print i
                        sv = i[0]
                        obssek = i[1]
                        #print obssek
                        if sat not in sv:
                            continue

                        for j in range(len(navtime)):
                            navsv = int(navvalues[j][0])

                            if int(sat[1:]) != navsv:
                                continue

                            navdatetime = time.strptime(str(navtime[j]), "%Y-%m-%d %H:%M:%S")  # delenie casu
                            navtimesecond = datetime.datetime(navdatetime[0], navdatetime[1], navdatetime[2],
                                                              navdatetime[3], navdatetime[4], navdatetime[5])
                            dat_nav = datetime.date(int(navdatetime[0]), int( navdatetime[1]), int( navdatetime[2]))  # datum prvej observacie
                            navsek = (navdatetime[3]*60*60) + (navdatetime[4]*60) + navdatetime[5]
                            if dat != dat_nav:
                                continue
                            Dt = obssek-navsek
                            #print navsv, navsek, obssek
                            if Dt < 0:
                                break
                            if Dt >= platnos_spravy:    # platonst spravy 7200 sekund = 2 hodiny uveden v input_data
                                continue

                            sattelite = poloha_druzice.fvypocet_poloha(navvalues[j], Dt)  # vypocet polohy druzice

                            azi, ele = fComputeAziEle(stanice[c][0], [sattelite[0], sattelite[1], sattelite[2]])   # vypocet elevacneho uhla, azimutu

                            # zobrazenie podkladovych map je v WGS 84 Web Mercator (EPSG:3857) tj. suradnice su prepocitane tiez a
                            # tym musi byt aj azimut prepocitany
                            latst, lonst, hel = fXYZ_to_LatLonH(stanice[c][0][0], stanice[c][0][1], stanice[c][0][2])
                            surstx, sursty = fwgs4326towgs3857(latst, lonst)
                            latsv, lonsv, hel2 = fXYZ_to_LatLonH(sattelite[0], sattelite[1], sattelite[2])
                            sursvx, sursvy = fwgs4326towgs3857(latsv, lonsv)

                            azi2 = fCalculateAzimuth([surstx, sursty], [sursvx, sursvy])

                            # print azi, ele
                            ele_min = stanice[c][1][0]
                            ele_max = stanice[c][1][1]
                            azi_min = stanice[c][2][0]
                            azi_max = stanice[c][2][1]

                            # filtrovanie vhodnych signalov podla azimutu a elevacneho uhla
                            if ele_min <= ele <= ele_max and azi_min <= azi2 <= azi_max:

                                t.append(i[1])
                                sinele.append(m.sin(m.radians(ele)))
                                signals1.append(m.pow(10, float(i[indexS1]) / 20))
                                signals2.append(m.pow(10, float(i[indexS2]) / 20))
                                signals5.append(m.pow(10, float(i[indexS5]) / 20))

                if len(sinele) < pocet_observacii:            # filter poctu observacii
                    print "pre ", sat, " malo observacii - vypocet neprebehne", len(sinele)
                    continue

                #print max(sinele), len(signals1), len(signals2)

                sinelexx = [[], []]
                s1 = [[], []]
                s2 = [[], []]
                s5 = [[], []]
                t2 = [[], []]
                sinele3 = [[], []]

                for l in range(len(sinele)-1):
                    if sinele[l] - sinele[l+1] > 0:
                        sinelexx[0].append(sinele[l])

                    elif sinele[l] - sinele[l+1] < 0:
                        sinelexx[1].append(sinele[l])
                ind = -1
                for q in range(len(sinelexx)):
                    if len(sinelexx[q]) != 1:
                        continue
                    ind = sinele.index(sinelexx[q][0])

                if ind == -1:
                    s1[0] = signals1
                    s2[0] = signals2
                    s5[0] = signals5
                    sinele3[0] = sinele
                    t2[0] = t
                else:
                    s1[0] = signals1[:ind + 1]
                    s2[0] = signals2[:ind + 1]
                    s5[0] = signals5[:ind + 1]
                    sinele3[0] = sinele[:ind + 1]
                    t2[0] = t[:ind + 1]

                    s1[1] = signals1[ind + 1:]
                    s2[1] = signals2[ind + 1:]
                    s5[1] = signals5[ind + 1:]
                    sinele3[1] = sinele[ind + 1:]
                    t2[1] = t[ind + 1:]

                for w in range(len(sinelexx)):

                    if len(sinele3[w]) < pocet_observacii:
                        print "pre ", sat, " malo observacii (" + str(len(sinele3[w])) + ") - vypocet neprebehne"
                        continue

                    print "pre " + sat + " pocet observacii: " + str(len(sinele3[w]))

                    SNR1 = np.array(s1[w])  # zoznam na np.array
                    SNR2 = np.array(s2[w])
                    SNR5 = np.array(s5[w])
                    x = np.array(sinele3[w])
                    timenum = np.array(t2[w])

                    alt_prob = []
                    i = vyska_min
                    while i < vyska_max:
                        alt_prob.append(i)
                        i = i + krok
                    # print alt_prob

                    f = []
                    for i in range(len(alt_prob)):      # prepocet vysky na frekvencu pre dlzky vln L1, L2, L5
                        f.append((2 * alt_prob[i]) / 0.190)
                    freqS1 = np.array(f)

                    f = []
                    for i in range(len(alt_prob)):
                        f.append((2 * alt_prob[i]) / 0.244)
                    freqS2 = np.array(f)

                    f = []
                    for i in range(len(alt_prob)):
                        f.append((2 * alt_prob[i]) / 0.255)
                    freqS5 = np.array(f)

                    modelSNR1, residuals1, SNR1_vector = model.modelSNR(timenum, SNR1, deg)
                    modelSNR1_1, residuals1_1, SNR1_vector_1 = model.modelSNR(x, SNR1_vector, deg)

                    modelSNR2, residuals2, SNR2_vector = model.modelSNR(timenum, SNR2, deg)
                    modelSNR2_1, residuals2_1, SNR2_vector_1 = model.modelSNR(x, SNR2_vector, deg)

                    if sat in satellites_S5:
                        modelSNR5, residuals5, SNR5_vector = model.modelSNR(timenum, SNR5, deg)
                        modelSNR5_1, residuals5_1, SNR5_vector_1 = model.modelSNR(x, SNR5_vector, deg)

                    try:  # skuska imporu kniznice Astropy, ak  nenaimportuje skusi kniznicu scipy ak ani ta sa nenaimportuje vypocet skonci.
                        from astropy.stats import LombScargle

                        powerS1 = LombScargle(SNR1_vector_1, x).power(freqS1)
                        powerS2 = LombScargle(SNR2_vector_1, x).power(freqS2)
                        ind2S1 = np.ndarray.argmax(powerS1)
                        ind2S2 = np.ndarray.argmax(powerS2)

                        if sat in satellites_S5:
                            powerS5 = LombScargle(SNR5_vector_1, x).power(freqS5)
                            ind2S5 = np.ndarray.argmax(powerS5)
                            vyska3.append(alt_prob[ind2S5])
                            doy_list2.append(dat)
                    except:
                        print sys.exc_info()[1]
                        try:
                            import scipy.signal as signal

                            powerS1 = signal.lombscargle(timenum, SNR1_vector_1, freqS1)  # vypocet periodogramu
                            powerS2 = signal.lombscargle(timenum, SNR2_vector_1, freqS2)  # vypocet periodogramu
                            ind2S1 = np.ndarray.argmax(powerS1)
                            ind2S2 = np.ndarray.argmax(powerS2)

                            if sat in satellites_S5:
                                powerS5 = signal.lombscargle(timenum, SNR5_vector_1, freqS5)  # vypocet periodogramu
                                ind2S5 = np.ndarray.argmax(powerS5)
                                vyska3.append(alt_prob[ind2S5])
                                doy_list2.append(dat)

                                """
                                plt.subplot(3, 1, 3)
                                plt.plot(alt_prob, powerS5, color='green', label='SNR5')
                                plt.xlabel('vyska anteny')
                                plt.ylabel('volt')"""

                        except:
                            print sys.exc_info()[1]
                            print "nepodarilo sa importovat ani jednu kniznicu potrebnu pre vypocet. Konec!"

                    if sat in satellites_S5:
                        print alt_prob[ind2S1], alt_prob[ind2S2], alt_prob[ind2S5]

                    print alt_prob[ind2S1], alt_prob[ind2S2]

                    vyska1.append(alt_prob[ind2S1])
                    doy_list1.append(dat)
                    vyska2.append(alt_prob[ind2S2])
                    if dat in SHMU[0] and dat in doy_list1:
                        sneh_ind = SHMU[0].index(dat)
                        SHMUdata[0].append(SHMU[0][sneh_ind])
                        sneh = (SHMU[1][sneh_ind] / 100.0)
                        SHMUdata[1].append(SHMU[1][sneh_ind])

                    """
                    plt.subplot(3, 1, 1)
                    plt.plot(alt_prob, powerS1, color='red', label='SNR1')
                    plt.ylabel('volt')
    
                    plt.subplot(3, 1, 2)
                    plt.plot(alt_prob, powerS2, color='blue', label='SNR2')
                    plt.ylabel('volt')
    
                    plt.show()
                    #plt.savefig("peiodogram_c"+str(doy_eph)+".png")
                    #plt.close()"""

        if len(doy_list1) == 0:
            print 'Ziadny vypocet vysky. Skontroluj data'
            os.system("pause")
            sys.exit()
        # priemer vysok
        mean1 = np.mean(vyska1)
        mean2 = np.mean(vyska2)
        mean3 = np.mean(vyska3)
        #print mean1, mean2, mean3
        # print vyska1
        # print vyska2
        # print vyska3

        # graf 1 vypocitane vysky a referencna vyska
        m1 = [stanice[c][3], stanice[c][3]]
        dm = [min(doy_list1), max(doy_list1)]
        plt.plot(doy_list1, vyska1, 'o', markersize=10, color="blue")
        plt.plot(doy_list1, vyska2, '^', markersize=10, color="red")

        if sat in satellites_S5:
            plt.plot(doy_list2,vyska3, 'p', markersize=10, color="green")

        plt.plot(dm, m1, color="black")
        plt.savefig(cesta_vystup + stanice[c][4] + '_ref_vyska.png')
        #plt.show()

        vyska_r1 = []
        for v in vyska1:
            vyska_r1.append((kame[3] - v) * 100)

        vyska_r2 = []
        for v in vyska2:
            vyska_r2.append((kame[3] - v) * 100)

        vyska_r3 = []
        for v in vyska3:
            vyska_r3.append((kame[3] - v) * 100)

        # vypocet RMSE pre jednotlive signaly
        rmse_val = reader.rmse(np.array(vyska_r1), np.array(SHMUdata[1]))
        print("rms error L1 is: " + str(rmse_val))
        rmse_val = reader.rmse(np.array(vyska_r2), np.array(SHMUdata[1]))
        print("rms error L2 is: " + str(rmse_val))
        rmse_val = reader.rmse(np.array(vyska_r3), np.array(SHMUdata[1]))
        print("rms error L5 is: " + str(rmse_val))

        # graf pre odhadovanu vysku snehu
        plt.plot(doy_list1, vyska_r1, 'o', markersize=10, color="blue")
        plt.plot(doy_list1, vyska_r2, '^', markersize=10, color="red")
        if sat in satellites_S5:
            plt.plot(doy_list2, vyska_r3, 'p', markersize=10, color="green")
        plt.plot(SHMUdata[0], SHMUdata[1])
        plt.savefig(cesta_vystup + stanice[c][4] + '_rozdiel.png')
        #plt.show()

        c = c+1
except:
    print
    print sys.exc_info()[0]
    import traceback
    print traceback.format_exc()
finally:
    print
    #raw_input()



