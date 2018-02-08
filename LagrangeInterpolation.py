from scipy.interpolate import lagrange
import numpy as np
import matplotlib.pyplot as plt

#x = [-13764.905816, -13765.792108, -13686.601303, -13566.283173, -13443.272722]
#y = [22307.596354, 21682.713765, 20767.079339, 19564.712410, 18088.589738]
#z = [3468.56267, 6266.923297, 8955.906822, 11488.553232, 13820.663301]
#t = [0, 900, 1800, 2700, 3600]

def fLagrangeInterpolation(x,y,z, t):
    xs = np.arange(t[0], t[-1], 1)
    lx = lagrange(t,x)
    ly = lagrange(t,y)
    lz = lagrange(t,z)


    return [lx(xs), ly(xs), lz(xs), xs]

#plt.plot(t, x, 'o', label='data')
#plt.plot(xs, lx(xs), label="S")
#plt.legend(loc='lower right')
#plt.show()