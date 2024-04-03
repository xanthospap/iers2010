#! /usr/bin/python

import matplotlib.pyplot as plt
import sys

def resolve(fn):
    t = []
    xp = []; yp = []
    dut = []; lod = []
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            fc = [ float(x) for x in line.split() ]
            t.append  (fc[0])
            xp.append (fc[1])
            yp.append (fc[2])
            dut.append(fc[3])
            try:
                lod.append(fc[4])
            except:
                pass
    return t, xp, yp, dut, lod


t1, xp1, yp1, dut1, dlod1 = resolve(sys.argv[1])
t2, xp2, yp2, dut2, dlod2 = resolve(sys.argv[2])

fig, axs = plt.subplots(2, 2, sharex=True)
axs[0, 0].plot(t1, xp1)
axs[0, 0].set_title(r'$\Delta x_p$ [$\mu as$]')

axs[0, 1].plot(t1, yp1)
axs[0, 1].yaxis.tick_right()
axs[0, 1].set_title(r'$\Delta y_p$ [$\mu as$]')

axs[1, 0].plot(t1, dut1)
axs[1, 0].set_title(r'$\Delta UT1$ [$\mu s$]')

axs[1, 1].plot(t1, dlod1)
axs[1, 1].set_title(r'$\Delta LOD$ [$\mu s$]')
axs[1, 1].yaxis.set_label_position("right")
axs[1, 1].yaxis.tick_right()

plt.savefig('argv1.png', bbox_inches='tight')

fig, axs = plt.subplots(2, 2, sharex=True)
axs[0, 0].plot(t2, xp2)
axs[0, 0].set_title(r'$\Delta x_p$ [$\mu as$]')

axs[0, 1].plot(t2, yp2)
axs[0, 1].yaxis.tick_right()
axs[0, 1].set_title(r'$\Delta y_p$ [$\mu as$]')

axs[1, 0].plot(t2, dut2)
axs[1, 0].set_title(r'$\Delta UT1$ [$\mu s$]')

if dlod2 != []:
    axs[1, 1].plot(t2, dlod2)
    axs[1, 1].set_title(r'$\Delta LOD$ [$\mu s$]')
    axs[1, 1].yaxis.set_label_position("right")
    axs[1, 1].yaxis.tick_right()

plt.savefig('argv2.png', bbox_inches='tight')

for i in range(len(t1)):
    assert(t1[i] == t2[i])

fig, axs = plt.subplots(2, 2, sharex=True)
axs[0, 0].plot(t1, [ x[0]-x[1] for x in zip(xp1,xp2) ])
axs[0, 0].set_title(r'$\Delta x_p$ [$\mu as$]')

axs[0, 1].plot(t2, [ x[0]-x[1] for x in zip(yp1,yp2) ])
axs[0, 1].yaxis.tick_right()
axs[0, 1].set_title(r'$\Delta y_p$ [$\mu as$]')

axs[1, 0].plot(t2, [ x[0]-x[1] for x in zip(dut1,dut2) ])
axs[1, 0].set_title(r'$\Delta UT1$ [$\mu s$]')

if dlod2 != [] and dlod1 != []:
    axs[1, 1].plot(t2, dlod2)
    axs[1, 1].set_title(r'$\Delta LOD$ [$\mu s$]')
    axs[1, 1].yaxis.set_label_position("right")
    axs[1, 1].yaxis.tick_right()

plt.savefig('argv3.png', bbox_inches='tight')
