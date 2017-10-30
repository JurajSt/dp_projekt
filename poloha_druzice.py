import math

def fvypocet_poloha(value, Dt):
    # konstanty
    omega = 0.000072921151467  # uhlova rychlost rotacie zeme
    mi = 398600500000000.0  # geocentricka gravitacna konstanta
    pi = math.pi

    n0 = math.sqrt((mi / math.pow(math.pow(value[11],2),3)))
    n = n0 + value[6]
    M = value[7] + n * Dt
##    print float(sprava[1].split()[3])
##    E = M
    a = M
    b = 1

    while round(a,5) != round(b,5):
        E0= a - ((a-value[9] * math.sin(a) - M)/(1 - value[9] * math.cos(a)))
        E1 = E0 - ((E0-value[9] * math.sin(E0) - M)/(1 - value[9] * math.cos(E0)))
        del(a,b)                # mazem z dovodu prepisovania premennej
        a = E1
        b = E0
##        print a, b

    cosv = (math.cos(E1)- value[9]) / (1-value[9] * math.cos(E1))
    if E1<0:
        E1 = (2*pi)+E1
    if E1>pi:
        v =(2*pi)-math.acos(cosv)
    else:
        v = math.acos(cosv)
    w = value[18]

    cos2x = math.pow(math.cos(v+w),2) - math.pow(math.sin(v+w),2)
    sin2x = 2 * math.sin(v+w) * math.cos(v+w)

    du = value[8] * cos2x + value[10] * sin2x
    u = v + w + du

    dr = value[17] * cos2x + value[5] * sin2x
    r = math.pow(value[11],2) * (1 - value[9] * math.cos(E1)) + dr

    x_p = r * math.cos(u)
    y_p = r * math.sin(u)
    r_kontrola = math.sqrt(math.pow(x_p,2) + math.pow(y_p,2))

##    print "r", r,"r_konnt", r_kontrola
##    if str(r) == str(r_kontrola):                   # taktiez musim dat str() inak padmienka nefunguje
##            print "vypocet je zatial v poriadku"    # podminka kontroluje spravnost vypoctu
##            print r, r_kontrola
##    elif str(r) != str(r_kontrola):                 # taktiez musim dat str inak padmienka nefunguje
##        print "vo vypocte je chyba"                 # bez str() urci dne rovnake hodnoty ako nerovne
##        print r, r_kontrola                         # a skript skonci
##        os.system("pause")
##        sys.exit()

    di = value[13] * cos2x + value[15] * sin2x
    i = value[16] + value[20] * Dt + di
    O = value[14] + (value[19] - omega) * Dt - omega * value[12]

    X = x_p * math.cos(O) - y_p * math.sin(O) * math.cos(i)
    Y = x_p * math.sin(O) + y_p * math.cos(O) * math.cos(i)
    Z = y_p * math.sin(i)

##    pole_XYZ.append(sprava[0][0:2].replace(" ",""))
##    pole_XYZ.append(cas_spravy)
##    pole_XYZ.append(X)
##    pole_XYZ.append(Y)
##    pole_XYZ.append(Z)

    return X, Y, Z