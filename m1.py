import numpy as np
import analyza_eph
import analyza_nav
import sys, os
from input_data import *
import fresnelZone

try:
    try:
        ephfiles = np.sort(next(os.walk(eph_path))[2])
    except StopIteration:
        ephfiles = []
    try:
        obsfiles = np.sort(next(os.walk(obs_path))[2])
    except StopIteration:
        print 'Zadne observacne data!!! Konec'
        sys.exit()
    try :
        navfiles = np.sort(next(os.walk(nav_path))[2])
    except StopIteration:
        navfiles = []

    if len(obsfiles) == 0 and len(ephfiles) == 0:
        print 'Zadne navigacne data ani presne efemeridy!!! Konec'
        sys.exit()

    if len(ephfiles) == 0:

        if analyza1 == 'a':
            analyza_nav.fAnalyza()
            fresnelZone.fFresnelZ(fresnel_azi,fresnel_ele)
            import kruh         # vyokona sa kod v skripte kruh

        else:
            import m2

    if len(navfiles) == 0:

        if analyza1 == 'a':
            #analyza_eph.fAnalyza()
            #fresnelZone.fFresnelZ(fresnel_azi, fresnel_ele)
            import kruh
        else:
            import m3



    os.system("pause")
except:
    print
    print sys.exc_info()[0]
    import traceback
    print traceback.format_exc()
finally:
    print


