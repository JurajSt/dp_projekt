import os, sys
from datetime import date
import numpy as np
from numpy import *
from input_data import *
import math as m
import model
from azi_ele_coor import *
import matplotlib.pyplot as plt
from dataSHMU import *
import reader
import ogr

try:
    import matplotlib.pyplot as plt
except:
    print sys.exc_info()[1]
    print "nie je nainstalovanapotrebna kniznica. Pre fungovanie aplikacie treba doinstalovat"

try:
    ephfiles = np.sort(next(os.walk(eph_path))[2])
    #obsfiles = np.sort(next(os.walk(obs_path))[2])

    doy_list1 = []
    doy_list2 = []
    SHMUdata = [[], []]
    vyska1 = []
    vyska2 = []
    vyska3 = []

    deg = 6   # stupen polynomu pri odstranovani trendu

    # datum zaciatku a konca vypoctu
    dat_start = date(rok_start, mes_start, den_start)
    dat_end = date(rok_konec, mes_konec, den_konec)
    delta = dat_end-dat_start
    pocet_dni = int((delta.total_seconds())/86400)

    row = 1
    col = 0
    c = 0
    while c < len(stanice):

        SHMU = fSHMU(c) # spracovanie snehovych dat

        for f in ephfiles:

            if 'sp3' != f.split('.')[1]:
                continue

            e = eph_path + f # cesta k SP3 suboru
            print
            print e

            sp3, date_e = reader.fReadSP3(e)   # nacitanie presnych efemerid
            doy_eph = reader.fdayofyear(int(date_e[0]), int(date_e[1]), int(date_e[2]))  # prepocet datumu na den v roku
            dat = date(int(date_e[0]), int(date_e[1]), int(date_e[2])) # datum prveho zaznamu sp3 suboru
            # print int(date_e[0]), int(date_e[1]), int(date_e[2])

            if dat > dat_end: # kontrola ci sa datum zaznamu sa nachadza v zadanom intervale datumu
                break
            if dat_start > dat:
                continue

            if 0 < doy_eph < 10:            # zistenie nazvu pre observacne data pomocou doy
                namefile = '00' + str(doy_eph) + '0'
            elif 10 <= doy_eph < 100:
                namefile = '0' + str(doy_eph) + '0'
            else:
                namefile = str(doy_eph) + '0'

            rok = date_e[0][2:]
            obsfile = obs_path + stanice[c][4] + namefile + '.' + rok + 'o'
            if not os.path.exists(obsfile):
                print 'ziadny observacny subor'
                continue
            # else: continue
            r_obs = open(obsfile, "r")
            lines = r_obs.readlines()
            r_obs.close()

            print obsfile
            print dat

            #parsrovanie observacnych dat
            header, version, headlines, headlength, obstimes, sats, numallsvs, numsvs, date_o = reader.scan(lines)
            data, obstypes = reader.processBlocks(lines, header, obstimes, headlines, headlength, sats, numallsvs, satellite)

            del lines

            indexS1 = obstypes.index('S1') + 2
            indexS2 = obstypes.index('S2') + 2
            indexS5 = obstypes.index('S5') + 2

            # svs = []
            den_hodiny = []
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
                        if sat not in sv:
                            continue

                        lg = reader.fInterpolateSatelliteXYZ(sp3, sv, i[1])  # vypocet lagrangeovej interpolacie
                        azi, ele = fComputeAziEle(stanice[c][0], [lg[0], lg[1], lg[2]])   # vypocet elevacneho uhla, azimutu

                        # zobrazenie podkladovych map je v WGS 84 Web Mercator (EPSG:3857) tj. suradnice su prepocitane tiez a
                        # tym musi byt aj azimut prepocitany
                        latst, lonst, hel = fXYZ_to_LatLonH(stanice[c][0][0], stanice[c][0][1], stanice[c][0][2])
                        surstx, sursty = fwgs4326towgs3857(latst, lonst)
                        latsv, lonsv, hel2 = fXYZ_to_LatLonH(lg[0], lg[1], lg[2])
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
                            #print ele
                            break

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
                ind = 0
                for q in range(len(sinelexx)):
                    if len(sinelexx[q]) != 1:
                        continue
                    ind = sinele.index(sinelexx[q][0])

                if ind == 0:
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

                    try:        # skuska imporu kniznice Astropy, ak  nenaimportuje skusi kniznicu scipy ak ani ta sa nenaimportuje vypocet skonci.
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

                    print alt_prob[ind2S1], alt_prob[ind2S2]#, alt_prob[ind2S5]

                    vyska1.append(alt_prob[ind2S1])
                    doy_list1.append(dat)
                    vyska2.append(alt_prob[ind2S2])

                    if dat in SHMU[0] and dat in doy_list1:
                        sneh_ind = SHMU[0].index(dat)
                        SHMUdata[0].append(SHMU[0][sneh_ind])
                        sneh = (SHMU[1][sneh_ind] / 100.0)
                        SHMUdata[1].append(SHMU[1][sneh_ind])
                        # print SHMUdata[1], SHMU[1][sneh_ind]

                    """
                    plt.subplot(3, 1, 1)
                    plt.plot(alt_prob, powerS1, color='red', label='SNR1')
                    plt.ylabel('volt')
    
                    plt.subplot(3, 1, 2)
                    plt.plot(alt_prob, powerS2, color='blue', label='SNR2')
                    plt.ylabel('volt')
    
                    plt.subplot(3, 1, 3)
                    plt.plot(alt_prob, powerS5, color='green', label='SNR5')
                    plt.xlabel('vyska anteny')
                    plt.ylabel('volt')
                    plt.show()"""
                    #plt.savefig("peiodogram_c"+str(doy_eph)+".png")
                    #plt.close()


        if  len(doy_list1) == 0:
            print 'Ziadny vypocet vysky. Skontroluj data'
            os.system("pause")
            sys.exit()
        # priemer vysok
        mean1 = np.mean(vyska1)
        mean2 = np.mean(vyska2)
        if sat in satellites_S5:
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
        plt.savefig(cesta_vystup+ stanice[c][4]+ '_ref_vyska.png')
        #plt.show()

        vyska_r1 = []
        for v in vyska1:
            vyska_r1.append((kame[3] - v) * 100)

        vyska_r2 = []
        for v in vyska2:
            vyska_r2.append((kame[3] - v) * 100)

        if sat in satellites_S5:
            vyska_r3 = []
            for v in vyska3:
                vyska_r3.append((kame[3] - v) * 100)


        # vypocet RMSE pre jednotlive signaly
        print
        try:
            rmse_val = reader.rmse(np.array(vyska_r1), np.array(SHMUdata[1]))
            print("rms chyba L1 is: " + str(rmse_val))
            rmse_val = reader.rmse(np.array(vyska_r2), np.array(SHMUdata[1]))
            print("rms chyba L2 is: " + str(rmse_val))
            if sat in satellites_S5:
                rmse_val = reader.rmse(np.array(vyska_r3), np.array(SHMUdata[1]))
                print("rms chyba L5 is: " + str(rmse_val))
        except ValueError as e:
            print "ValueError: " + e
            print "datum observacie nie je v datach v meranej vyske snehu"

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

