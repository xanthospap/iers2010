#! /usr/bin/python

import sys
import matplotlib
matplotlib.use("pgf")

import matplotlib.pyplot as plt
import numpy as np

# Use the seborn style
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

def eopts(fn, type):
    tp = type.upper()
    vmjd=[]; vx=[]; vy=[]; vut1=[]; vlod=[];
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            if line.lstrip().startswith(tp):
                mjd, xp, yp, ut1, lod = [ float(x) for x in line.split()[1:] ]
                vmjd.append(mjd); 
                vx.append(xp);vy.append(yp);vut1.append(ut1);vlod.append(lod);
    return vmjd, vx, vy, vut1, vlod

if __name__ == "__main__":

    ## resolve EOPS
    t1, xp1, yp1, ut11, lod1 = eopts(sys.argv[1], "OCTIDE")
    t2, xp2, yp2, ut12, lod2 = eopts(sys.argv[1], "LIBRTN")

    fig, axs = plt.subplots(2, 2, figsize=set_size(418,1,(2,2)))
    axs[0, 0].plot(t1, xp1, t2, xp2)
    axs[0, 0].set_title(r'$\Delta x_p$')
    axs[0, 0].set(ylabel=r'$\Delta x_p$ in $[\mu arcsec]$')
    axs[0, 0].legend([r'$\Delta x_p$ ocean tide', r'$\Delta x_p$ libration'])
    #
    axs[0, 1].plot(t1, yp1, t2, yp2)
    axs[0, 1].set_title(r'$\Delta y_p$')
    axs[0, 1].set(ylabel=r'$\Delta y_p$ in $[\mu arcsec]$')
    axs[0, 1].legend([r'$\Delta y_p$ ocean tide', r'$\Delta y_p$ libration'])
    #
    axs[1, 0].plot(t1, ut11, t2, ut12)
    axs[1, 0].set_title(r'$\Delta UT1$')
    axs[1, 0].set(ylabel=r'$\Delta UT1$ in $[\mu sec]$', xlabel="MJD")
    axs[1, 0].legend([r'$\Delta UT1$ ocean tide', r'$\Delta UT1$ libration'])
    #
    axs[1, 1].plot(t1, lod1, t2, lod2)
    axs[1, 1].set_title(r'$\Delta LOD$')
    axs[1, 1].set(ylabel=r'$\Delta LOD$ in $[\mu sec]$', xlabel="MJD")
    axs[1, 1].legend([r'$\Delta LOD$ ocean tide', r'$\Delta LOD$ libration'])

    # plt.show()
    fig.tight_layout()
    plt.savefig('eop_variations.pgf', format='pgf')
