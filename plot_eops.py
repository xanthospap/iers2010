#! /usr/bin/python

import fileinput
import matplotlib.pyplot as plt

t = []
xp = []; yp = []
dut = []; lod = []
dX = []; dY = []
Rt = []
Rxp = []; Ryp = []
Rdut = []; Rlod = []
RdX = []; RdY = []

for line in fileinput.input():
    try:
        fc = [ float(x) for x in line.split() ]
        t.append  (fc[0])
        xp.append (fc[1])
        yp.append (fc[2])
        dut.append(fc[3])
        lod.append(fc[4])
        dX.append (fc[5])
        dY.append (fc[6])
    except:
        fc = [ float(x) for x in line.split()[1:] ]
        Rt.append  (fc[0])
        Rxp.append (fc[1])
        Ryp.append (fc[2])
        Rdut.append(fc[3])
        Rlod.append(fc[4])
        RdX.append (fc[5])
        RdY.append (fc[6])
        

fig, axs = plt.subplots(2, 2, sharex=True)
axs[0, 0].plot(t, xp, Rt, Rxp, '+')
axs[0, 0].set_title(r'$x$pole [sec]')
axs[0, 1].plot(t, yp, 'tab:orange', Rt, Ryp, '+')
axs[0, 1].yaxis.set_label_position("right")
axs[0, 1].yaxis.tick_right()
axs[0, 1].set_title(r'$y$pole [sec]')
axs[1, 0].plot(t, dut, 'tab:green', Rt, Rdut, '+')
axs[1, 0].set_title(r'$\Delta UT$ [sec]')
axs[1, 1].plot(t, lod, 'tab:red', Rt, Rlod, '+')
axs[1, 1].set_title(r'$\Delta LOD$ [sec]')
axs[1, 1].yaxis.set_label_position("right")
axs[1, 1].yaxis.tick_right()

#for ax in axs.flat:
#    ax.set(xlabel='x-label', ylabel='y-label')

plt.show()
