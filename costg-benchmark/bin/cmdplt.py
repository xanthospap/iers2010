#! /usr/bin/python

import os
import sys
import math
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
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
    '--scale',
    metavar='SCALE',
    dest='scale',
    default=1e0,
    required=False,
    type=float,
    help='Scale plot')

parser.add_argument(
    '--no-diff',
    action='store_true',
    dest='nodif',
    required=False,
    help='Plot acceleration results instead of differences.')

parser.add_argument(
    '--with-diff',
    action='store_true',
    dest='withdif',
    required=False,
    help='Plot acceleration results and differences (same plot).')

parser.add_argument(
    '--plot-norm',
    action='store_true',
    dest='pnorm',
    required=False,
    help='Plot acceleration differences norm on another plot.')

parser.add_argument(
    '--quite',
    action='store_true',
    dest='quiet',
    help='Do not show plot(s)')

if __name__ == "__main__":
    args = parser.parse_args()
    scale = args.scale

    t = []
    a1x = []; a1y=[]; a1z=[];
    a2x = []; a2y=[]; a2z=[];
    title=None

    for line in sys.stdin:
        if line[0] != '#' and not line.startswith('mjd'):
            try:
                imjd, sec, x1, y1, z1, x2, y2, z2 = [ float(x) for x in line.split(' ') ]
                t.append(sec)
                a1x.append(x1); a1y.append(y1); a1z.append(z1);
                a2x.append(x2); a2y.append(y2); a2z.append(z2);
            except:
                pass
        else:
            if line.startswith('#title'):
                title = line.replace('#title','').strip()

    ## do nothing on empty input
    if len(t) <= 1: sys.exit(1)

## Get statistics of differences
    ar = [z[0]-z[1] for z in zip(a1x,a2x)]
    sx = stats.describe(ar)
    xstr = r"$\delta \ddot{{r}}_x$: max={:.1e} mean={:+.1e}$\pm${:.2e}".format(max(abs(x) for x in sx.minmax), sx.mean, math.sqrt(sx.variance))
    ar = [z[0]-z[1] for z in zip(a1y,a2y)]
    sy = stats.describe(ar)
    ystr = r'$\delta \ddot{{r}}_y$: max={:.1e} mean={:+.1e}$\pm${:.2e}'.format(max(abs(x) for x in sy.minmax), sy.mean, math.sqrt(sy.variance))
    ar = [z[0]-z[1] for z in zip(a1z,a2z)]
    sz = stats.describe(ar)
    zstr = r'$\delta \ddot{{r}}_z$: max={:.1e} mean={:+.1e}$\pm${:.2e}'.format(max(abs(x) for x in sz.minmax), sz.mean, math.sqrt(sz.variance))

    fig, axs = plt.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    group_labels = []

## Plot acceleration results
    if (args.nodif or args.withdif):
        axs[0].plot(t, [scale*z for z in a1x])
        axs[0].plot(t, [scale*z for z in a2x])
        axs[0].text(t[0], sx.minmax[0]*scale, xstr)
        
        axs[1].plot(t, [scale*z for z in a1y])
        axs[1].plot(t, [scale*z for z in a2y])
        axs[1].text(t[0], sy.minmax[0]*scale, ystr)
        
        axs[2].xaxis.set_major_locator(MultipleLocator(3600*3))
        axs[2].xaxis.set_major_formatter(lambda x, pos: str(int(x/3600e0))+'h')
        axs[2].xaxis.set_minor_locator(MultipleLocator(3600))
        axs[2].plot(t, [scale*z for z in a1z])
        axs[2].plot(t, [scale*z for z in a2z])
        axs[2].text(t[0], sz.minmax[0]*scale, zstr)
        
        plt.xlabel('Hours of Day')
        fig.text(0.0, 0.5, r'$\delta \ddot{r}_x$, $\delta \ddot{r}_y$ and $\delta \ddot{r}_z$ in $[m/sec^2]\times$'+'{:.1e}'.format(1/scale), va='center', rotation='vertical')
        
        group_labels.append("acc. group 1")
        group_labels.append("acc. group 2")

## Plot differences (this-COST-G) per component
    if (args.withdif or (not args.nodif)):
        axs[0].plot(t, [scale*(z[0]-z[1]) for z in zip(a1x,a2x)])
        axs[0].text(t[0], sx.minmax[0]*scale, xstr)
        
        axs[1].plot(t, [scale*(z[0]-z[1]) for z in zip(a1y,a2y)])
        axs[1].text(t[0], sy.minmax[0]*scale, ystr)
        
        axs[2].xaxis.set_major_locator(MultipleLocator(3600*3))
        axs[2].xaxis.set_major_formatter(lambda x, pos: str(int(x/3600e0))+'h')
        axs[2].xaxis.set_minor_locator(MultipleLocator(3600))
        axs[2].plot(t, [scale*(z[0]-z[1]) for z in zip(a1z,a2z)])
        axs[2].text(t[0], sz.minmax[0]*scale, zstr)
        
        plt.xlabel('Hours of Day')
        fig.text(0.0, 0.5, r'$\delta \ddot{r}_x$, $\delta \ddot{r}_y$ and $\delta \ddot{r}_z$ in $[m/sec^2]\times$'+'{:.1e}'.format(1/scale), va='center', rotation='vertical')
        
        group_labels.append("acc. differences")

##  x-grids on
    axs[0].xaxis.grid(True, which='major')
    axs[1].xaxis.grid(True, which='major')
    axs[2].xaxis.grid(True, which='major')
    axs[0].xaxis.grid(True, which='minor')
    axs[1].xaxis.grid(True, which='minor')
    axs[2].xaxis.grid(True, which='minor')

    if title: fig.suptitle(title)
    # Legend
    if len(group_labels) >= 2 :
        # Get colours for current style
        colours = plt.rcParams['axes.prop_cycle'].by_key()['color']
        # Set up handles (the bits that are drawn in the legend)
        handles = []
        for idx in range(len(group_labels)):
            colour = colours[idx]
            handles.append(Patch(edgecolor=colour, facecolor=colour, fill=True))
            # Actually create our figure legend, using the handles and labels
            fig.legend(handles=handles, labels=group_labels, loc="lower center", ncol=len(group_labels))
        fig.patch.set_facecolor('white')
        plt.tight_layout()
        # Remove space from the right handside to make space for the legend
        plt.subplots_adjust(right=0.9)
    else:
        plt.tight_layout()
    
    plt.show()

    if args.save_as:
        fn = os.path.join(args.save_as + 'raw_diffs')
        plt.savefig(fn+'.jpg')

## Plot norm differences
    if (args.pnorm):
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
        
        axs[1].xaxis.set_major_locator(MultipleLocator(3600*3))
        axs[1].xaxis.set_major_formatter(lambda x, pos: str(int(x/3600e0))+'h')
        axs[1].xaxis.set_minor_locator(MultipleLocator(3600))
        
        axs[0].xaxis.grid(True, which='major')
        axs[1].xaxis.grid(True, which='major')
        axs[0].xaxis.grid(True, which='minor')
        axs[1].xaxis.grid(True, which='minor')
        plt.xlabel('Hours of Day')

        if title: fig.suptitle(title)
        plt.tight_layout()
        
        if args.save_as:
            fn = os.path.join(args.save_as + 'norm_diffs')
            plt.savefig(fn+'.jpg')
