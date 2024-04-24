#! /usr/bin/python


## high quality plots e.g. for LaTeX
EXPORT_PGF = False
if EXPORT_PGF:
    import matplotlib
    matplotlib.use("pgf")

import os
import sys
import math
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
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

if EXPORT_PGF:
    plt.style.use('bmh')
    matplotlib.rcParams.update({
        "pgf.texsystem": "pdflatex",
        'font.family': 'serif',
        'font.size': 6,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'axes.titlesize': 10,
        'axes.labelsize': 8,
        'lines.linewidth': 1,
        'lines.markersize' : 5,
        'xtick.labelsize' : 5,
        'ytick.labelsize' : 5
    })

if __name__ == "__main__":
    ## parse cmd
    args = parser.parse_args()

    t = []
    a1x = []; a1y=[]; a1z=[];
    a2x = []; a2y=[]; a2z=[];
    title=None

    for line in sys.stdin:
        if line[0] != '#':
            print(line)
            imjd, sec, x1, y1, z1, x2, y2, z2 = [ float(x) for x in line.split() ]
            t.append(sec)
            a1x.append(x1); a1y.append(y1); a1z.append(z1);
            a2x.append(x2); a2y.append(y2); a2z.append(z2);
        else:
            if line.startswith('#title'):
                title = line.replace('#title','').strip()

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

    fig, axs = plt.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    
    axs[0].plot(t, [z[0]-z[1] for z in zip(a1x,a2x)])
    axs[0].text(t[0], sx.minmax[0], xstr)
    
    axs[1].plot(t, [z[0]-z[1] for z in zip(a1y,a2y)])
    axs[1].text(t[0], sy.minmax[0], ystr)
    
    axs[2].xaxis.set_major_locator(MultipleLocator(3600*3))
    axs[2].xaxis.set_major_formatter(lambda x, pos: str(int(x/3600e0))+'h')
    axs[2].xaxis.set_minor_locator(MultipleLocator(3600))
    axs[2].plot(t, [z[0]-z[1] for z in zip(a1z,a2z)])
    axs[2].text(t[0], sz.minmax[0], zstr)
    
    plt.xlabel('Hours of Day')
    fig.text(0.0, 0.5, r'$\delta \ddot{r}_x$, $\delta \ddot{r}_y$ and $\delta \ddot{r}_z$ in $[m/sec^2]$', va='center', rotation='vertical')

##  x-grids on
    axs[0].xaxis.grid(True, which='major')
    axs[1].xaxis.grid(True, which='major')
    axs[2].xaxis.grid(True, which='major')
    axs[0].xaxis.grid(True, which='minor')
    axs[1].xaxis.grid(True, which='minor')
    axs[2].xaxis.grid(True, which='minor')

    if title: fig.suptitle(title)
    plt.tight_layout()
    
    if not EXPORT_PGF: plt.show()

    if args.save_as:
        fn = os.path.join(args.save_as + 'raw_diffs')
        if EXPORT_PGF:
            fn += '.pgf'
            plt.savefig(fn, format='pgf')
        else:
            plt.savefig(fn+'.jpg')

## Plot norm differences
    def norm(x,y,z): return math.sqrt(x*x+y*y+z*z)
    nd = [ norm(a1x[i], a1y[i], a1z[i]) - norm(a2x[i], a2y[i], a2z[i]) for i in range(len(a1x)) ]
    sz = stats.describe(nd)
    nstr = r'$\delta a$: max={:.2e} mean={:+.2e}$\pm${:.3e}'.format(max(abs(x) for x in sz.minmax), sz.mean, math.sqrt(sz.variance))

    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    
    axs[0].plot(t, nd, linewidth=2.0)
    axs[0].set_ylabel(r'$|\ddot{r}_cost| - |\ddot{r}|$ $[m/s^2]$')
    
    nd = [ norm(a1x[i]-a2x[i], a1y[i]-a2y[i], a1z[i]-a2z[i]) for i in range(len(a1x)) ]
    sz = stats.describe(nd)
    axs[1].plot(t, nd, linewidth=2.0)
    axs[1].set_ylabel(r'$|\ddot{\delta r}|$ $[m/s^2]$')
    
    axs[0].xaxis.grid(True, which='major')
    axs[1].xaxis.grid(True, which='major')
    axs[0].xaxis.grid(True, which='minor')
    axs[1].xaxis.grid(True, which='minor')
    plt.xlabel('Hours of Day')

    if title: fig.suptitle(title)
    plt.tight_layout()

    if not EXPORT_PGF: plt.show()
    
    if args.save_as:
        fn = os.path.join(args.save_as + 'norm_diffs')
        if EXPORT_PGF:
            fn += '.pgf'
            plt.savefig(fn, format='pgf')
        else:
            plt.savefig(fn+'.jpg')
