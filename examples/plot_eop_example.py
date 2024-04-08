#! /usr/bin/python

EXPORT_PGF = True

if EXPORT_PGF:
    import matplotlib
    matplotlib.use("pgf")

import sys
import fileinput
import matplotlib.pyplot as plt
import numpy as np

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

def set_size(width_pt, fraction=1, subplots=(1, 1)):
    """Set figure dimensions to sit nicely in our document.

    Parameters
    ----------
    width_pt: float
            Document width in points
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    golden_ratio = (5**.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[0] / subplots[1])

    return (fig_width_in, fig_height_in)


traw=[]; xpraw=[]; ypraw=[]; dut1raw=[]; lodraw=[];
tcor=[]; xpcor=[]; ypcor=[]; dut1cor=[]; lodcor=[];

for line in fileinput.input():
    t,xp,yp,dut1,lod=[float(x) for x in line.split()[1:]]
    if line.startswith('[RAW]'):
        traw.append(t); xpraw.append(xp); ypraw.append(yp); dut1raw.append(dut1); lodraw.append(lod)
    else:
        tcor.append(t); xpcor.append(xp); ypcor.append(yp); dut1cor.append(dut1); lodcor.append(lod)

if EXPORT_PGF:
    fig, axs = plt.subplots(2, 2, figsize=set_size(418,1,(2,2)))
else:
    fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(traw, xpraw, tcor, xpcor)
axs[0, 0].set_title(r'$x_p$')
axs[0, 0].set(ylabel=r'$x_p$ in $[arcsec]$')
axs[0, 0].legend([r'${x_p}_{iers}$', r'${x_p}_{corrected}$'])
#
axs[0, 1].plot(traw, ypraw, tcor, ypcor)
axs[0, 1].set_title(r'$y_p$')
axs[0, 1].set(ylabel=r'$y_p$ in $[arcsec]$')
axs[0, 1].legend([r'${y_p}_{iers}$', r'${y_p}_{corrected}$'])
#
axs[1, 0].plot(traw, dut1raw, tcor, dut1cor)
axs[1, 0].set_title(r'$\Delta UT1$')
axs[1, 0].set(ylabel=r'$\Delta UT1$ in $[sec]$')
axs[1, 0].legend([r'${\Delta UT1}_{iers}$', r'${\Delta UT1}_{corrected}$'])
#
axs[1, 1].plot(traw, lodraw, tcor, lodcor)
axs[1, 1].set_title(r'$LOD$')
axs[1, 1].set(ylabel=r'$LOD$ in $[sec]$')
axs[1, 1].legend([r'${LOD}_{iers}$', r'${LOD}_{corrected}$'])

fig.tight_layout()

if EXPORT_PGF:
    plt.savefig('eop_interpolation.pgf', format='pgf')
else:
    plt.show()
