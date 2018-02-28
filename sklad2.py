import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()


ax1.plot(np.arange(100), np.arange(100))  #dummy data

ax1_ticks = ax1.get_xticks()
ax2_scale = ax1_ticks / 14.0 # scale the original axis ticks by 14

ax1.set_xlabel("Pounds")

ax2.set_xticks(ax2_scale)
ax2.set_xlabel("Stone")

plt.show()

new_tick_locations = np.array([.2, .5, .9])

def tick_function(X):
    V = 1/(1+X)
    return ["%.3f" % z for z in V]

s = tick_function(new_tick_locations)
