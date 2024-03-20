#! /usr/bin/python

import sys
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

## Get statistics of differences
    sx = stats.describe([z[0]-z[1] for z in zip(a1x,a2x)])
    sy = stats.describe([z[0]-z[1] for z in zip(a1y,a2y)])
    sz = stats.describe([z[0]-z[1] for z in zip(a1z,a2z)])

## Plot differences (this-COST-G)
    ax1 = plt.subplot(311)
    plt.plot(t, [z[0]-z[1] for z in zip(a1x,a2x)])
    plt.tick_params('x', labelsize=6)
# share x only
    ax2 = plt.subplot(312, sharex=ax1)
    plt.plot(t, [z[0]-z[1] for z in zip(a1y,a2y)])
# make these tick labels invisible
    plt.tick_params('x', labelbottom=False)
# share x and y
    ax3 = plt.subplot(313, sharex=ax1)
    plt.plot(t, [z[0]-z[1] for z in zip(a1z,a2z)])
    plt.tick_params('x', labelbottom=False)
    plt.show()
