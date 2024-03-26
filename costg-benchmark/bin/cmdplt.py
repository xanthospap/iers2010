#! /usr/bin/python

import sys
import math
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

t = []
a1x = []; a1y=[]; a1z=[];
a2x = []; a2y=[]; a2z=[];

if __name__ == "__main__":
    for line in sys.stdin:
        if line[0] != '#':
            print(line)
            imjd, sec, x1, y1, z1, x2, y2, z2 = [ float(x) for x in line.split() ]
            t.append(sec)
            a1x.append(x1); a1y.append(y1); a1z.append(z1);
            a2x.append(x2); a2y.append(y2); a2z.append(z2);

    ## do nothing on empty input
    if len(t) <= 1: sys.exit(1)

## Get statistics of differences
    ar = [z[0]-z[1] for z in zip(a1x,a2x)]
    sx = stats.describe(ar)
    xstr = r'X: max={:.2e} mean={:+.2e}$\pm${:.3e}'.format(max(abs(x) for x in sx.minmax), sx.mean, math.sqrt(sx.variance))
    ar = [z[0]-z[1] for z in zip(a1y,a2y)]
    sy = stats.describe(ar)
    ystr = r'Y: max={:.2e} mean={:+.2e}$\pm${:.3e}'.format(max(abs(x) for x in sy.minmax), sy.mean, math.sqrt(sy.variance))
    ar = [z[0]-z[1] for z in zip(a1z,a2z)]
    sz = stats.describe(ar)
    zstr = r'Z: max={:.2e} mean={:+.2e}$\pm${:.3e}'.format(max(abs(x) for x in sz.minmax), sz.mean, math.sqrt(sz.variance))

## Plot differences (this-COST-G) per component
    plt.rcParams["font.family"] = "monospace"

    ax1 = plt.subplot(311)
    plt.plot(t, [z[0]-z[1] for z in zip(a1x,a2x)])
    plt.text(t[0], sx.minmax[0], xstr)
    plt.tick_params('x', labelsize=6)
    plt.tick_params('x', labelbottom=False)
    
    ax2 = plt.subplot(312, sharex=ax1)
    plt.plot(t, [z[0]-z[1] for z in zip(a1y,a2y)])
    plt.text(t[0], sy.minmax[0], ystr)
    plt.tick_params('x', labelbottom=False)
    plt.ylabel('Gravitational acceleration discrepancies $[m/s^2]$')
    
    ax3 = plt.subplot(313, sharex=ax1)
    plt.plot(t, [z[0]-z[1] for z in zip(a1z,a2z)])
    plt.text(t[0], sz.minmax[0], zstr)
    plt.tick_params('x', labelsize=6)
    plt.xlabel('Seconds of Day')
    
    plt.show()

## Plot norm differences
    def norm(x,y,z): return math.sqrt(x*x+y*y+z*z)
    nd = [ norm(a1x[i], a1y[i], a1z[i]) - norm(a2x[i], a2y[i], a2z[i]) for i in range(len(a1x)) ]
    sz = stats.describe(nd)
    nstr = r'$\delta a$: max={:.2e} mean={:+.2e}$\pm${:.3e}'.format(max(abs(x) for x in sz.minmax), sz.mean, math.sqrt(sz.variance))
    fig, ax = plt.subplots()
    ax.plot(t, nd, linewidth=2.0)
    ax.set_xlabel('Seconds of Day')
    ax.set_ylabel('Gravitational acceleration discrepancies $[m/s^2]$')
    plt.show()
