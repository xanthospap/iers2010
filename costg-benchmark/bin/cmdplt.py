#! /usr/bin/python

import os
import sys
import math
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description=
    'Plot acceleration diffs w.r.t costg results.',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Apr, 2024'''))

parser.add_argument(
    '-s',
    '--saveas',
    metavar='SAVEAS',
    dest='save_as',
    default=None,
    required=False,
    help='Save the file using this file(name)')

parser.add_argument(
    '--quite',
    action='store_true',
    dest='quiet',
    help='Do not show plot(s)')

if __name__ == "__main__":
    ## parse cmd
    args = parser.parse_args()

    t = []
    a1x = []; a1y=[]; a1z=[];
    a2x = []; a2y=[]; a2z=[];

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

    major_ticks = np.arange(0, t[-1], 3600*3)
    minor_ticks = np.arange(0, t[-1], 3600*1)

    ax1 = plt.subplot(311)
    ax1.set_xticks(major_ticks)
    ax1.set_xticks(minor_ticks, minor=True)
    ax1.grid(axis='x', which='minor', alpha=0.2)
    ax1.grid(axis='x', which='major', alpha=0.5)
    plt.plot(t, [z[0]-z[1] for z in zip(a1x,a2x)])
    plt.text(t[0], sx.minmax[0], xstr)
    plt.tick_params('x', labelsize=6)
    plt.tick_params('x', labelbottom=False)
    
    ax2 = plt.subplot(312, sharex=ax1)
    ax2.set_xticks(major_ticks)
    ax2.set_xticks(minor_ticks, minor=True)
    ax2.grid(axis='x', which='minor', alpha=0.2)
    ax2.grid(axis='x', which='major', alpha=0.5)
    plt.plot(t, [z[0]-z[1] for z in zip(a1y,a2y)])
    plt.text(t[0], sy.minmax[0], ystr)
    plt.tick_params('x', labelbottom=False)
    plt.ylabel('Gravitational acceleration discrepancies $[m/s^2]$')
    
    ax3 = plt.subplot(313, sharex=ax1)
    ax3.set_xticks(major_ticks)
    ax3.set_xticks(minor_ticks, minor=True)
    ax3.grid(axis='x', which='minor', alpha=0.2)
    ax3.grid(axis='x', which='major', alpha=0.5)
    plt.plot(t, [z[0]-z[1] for z in zip(a1z,a2z)])
    plt.text(t[0], sz.minmax[0], zstr)
    plt.tick_params('x', labelsize=6)
    plt.xlabel('Seconds of Day')
    
    plt.show()

    if args.save_as:
        fn = os.path.join(args.save_as + 'raw_diffs.svg')
        plt.savefig(fn, format='svg', dpi=300)

## Plot norm differences
    def norm(x,y,z): return math.sqrt(x*x+y*y+z*z)
    nd = [ norm(a1x[i], a1y[i], a1z[i]) - norm(a2x[i], a2y[i], a2z[i]) for i in range(len(a1x)) ]
    sz = stats.describe(nd)
    nstr = r'$\delta a$: max={:.2e} mean={:+.2e}$\pm${:.3e}'.format(max(abs(x) for x in sz.minmax), sz.mean, math.sqrt(sz.variance))

    ax1 = plt.subplot(211)
    ax1.plot(t, nd, linewidth=2.0)
    #ax1.set_xlabel('Seconds of Day')
    #ax1.set_ylabel('Gravitational acceleration discrepancies $[m/s^2]$')
    ax1.set_ylabel(r'$|\ddot{r}_cost| - |\ddot{r}|$ $[m/s^2]$')
    ax1.set_xticks(major_ticks)
    ax1.set_xticks(minor_ticks, minor=True)
    ax1.grid(axis='x', which='minor', alpha=0.2)
    ax1.grid(axis='x', which='major', alpha=0.5)
    plt.tick_params('x', labelbottom=False)
    
    nd = [ norm(a1x[i]-a2x[i], a1y[i]-a2y[i], a1z[i]-a2z[i]) for i in range(len(a1x)) ]
    sz = stats.describe(nd)
    ax2 = plt.subplot(212)
    ax2.plot(t, nd, linewidth=2.0)
    ax2.set_xlabel('Seconds of Day')
    #ax2.set_ylabel('Gravitational acceleration discrepancies $[m/s^2]$')
    ax2.set_ylabel(r'$|\ddot{\delta r}|$ $[m/s^2]$')
    ax2.set_xticks(major_ticks)
    ax2.set_xticks(minor_ticks, minor=True)
    ax2.grid(axis='x', which='minor', alpha=0.2)
    ax2.grid(axis='x', which='major', alpha=0.5)
    plt.tick_params('x', labelsize=6)

    plt.show()
    
    if args.save_as:
        fn = os.path.join(args.save_as + 'norm_diffs.svg')
        # plt.savefig(fn, format='svg', dpi=300)
