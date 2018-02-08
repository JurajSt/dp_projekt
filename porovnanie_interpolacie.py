import rinex_reader as reader
import poloha_druzice
import CubicSplineInterpolation
import LagrangeInterpolation
import numpy as np
import matplotlib.pyplot as plt

#PG24  14763.334976 -21041.066684   6317.119555 0
#PG24  14669.539566 -20119.871238   8978.917668 900
#PG24  14522.923576 -18919.148151  11484.588392 1800
#PG24  14361.477832 -17450.283830  13790.621311 2700
#PG24  14221.379379 -15733.092760  15857.042668 3600
x = [14763334.976, 14669539.566, 14522923.576, 14361477.832, 14221379.379]
y = [-21041066.684, -20119871.238, -18919148.151, -17450283.830, -15733092.760 ]
z = [6317119.555, 8978917.668, 11484588.392, 13790621.311, 15857042.668]
t = [0, 900, 1800, 2700, 3600]

xs = np.arange(t[0], t[-1], 1)
Dt = 0
point = []
navig = 'data/bordel/ganp2740_1.15n'
navdata = reader.rinexnav(navig)
navvalues = navdata.values
while Dt <=3600:
    sattelite = poloha_druzice.fvypocet_poloha(navvalues[0],Dt)
    #print str(sattelite[0]-x[i])+"\t"+str(sattelite[1]-y[i])+"\t"+str(+sattelite[2]-z[i])
    point.append(sattelite)
    Dt = Dt+1
int = CubicSplineInterpolation.fCubicSplineInterpolation(x,y,z,t)
lg = LagrangeInterpolation.fLagrangeInterpolation(x, y, z, t)

xn = []
yn = []
zn = []

xe = []
ye = []
ze = []

xl = []
yl = []
zl = []
for i in range(0,3600):
    #print i, str(point[i][0] - int[0](i)) + "\t" + str(point[i][1] - int[1](i)) + "\t" + str(point[i][2] - int[2](i)) + \
    #         "\t\t" + str(point[i][0] - lg[0](i)) + "\t" + str(point[i][1] - lg[1](i)) + "\t" + str(point[i][2] - lg[2](i))
    xn.append(point[i][0])
    yn.append(point[i][1])
    zn.append(point[i][2])
    xe.append(int[0](i))
    ye.append(int[1](i))
    ze.append(int[2](i))
    xl.append(int[0](i))
    yl.append(int[1](i))
    zl.append(int[2](i))


#plt.plot(x, y,'o', label='data from navig message')
plt.plot(xn, yn, 'b', label='data from navig message')
#plt.plot(int[0](xs), int[1](xs),'b', label="data from interpolation")
plt.plot(xe, ye, 'r', label='data from eph CSL')
plt.plot(xl, yl, 'g', label='data from eph LagrI')
plt.legend(loc='lower left', ncol=2)
plt.show()