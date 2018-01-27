from scipy.interpolate import CubicSpline
import numpy as np
import matplotlib.pyplot as plt

#x = [-13764.905816, -13765.792108, -13686.601303, -13566.283173, -13443.272722]
#y = [22307.596354, 21682.713765, 20767.079339, 19564.712410, 18088.589738]
#z = [3468.56267, 6266.923297, 8955.906822, 11488.553232, 13820.663301]
#t = [0, 900, 1800, 2700, 3600]

#x = [-13764.905816, -13686.601303]
#y = [22307.596354, 20767.079339]
#z = [3468.56267, 8955.906822]
#t = [0, 1800]

def fCubicSplineInterpolation(x, y, z, t):
    xs = np.arange(t[0], t[-1], 1)
    cx = CubicSpline(t, x)
    cy = CubicSpline(t, y)
    cz = CubicSpline(t, z)


    return [cx, cy, cz]

#plt.figure(figsize=(6.5, 4))
#plt.plot(x, y, 'o', label='data')
#plt.plot(cx(xs), cy(xs), label="S")
#plt.xlim(0, 3600)
#plt.legend(loc='lower left', ncol=2)
#plt.show()
#print cx(900), cy(900), cz(900)