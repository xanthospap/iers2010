#! /usr/bin/python

import os
import sys
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import argparse

parser = argparse.ArgumentParser(
    description="Plot state integration diffs w.r.t reference results.",
    epilog=(
        """National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Apr, 2024"""
    ),
)

parser.add_argument(
    "-i",
    "--costg01",
    metavar="COSTG_01FILE",
    dest="costg_fn",
    default=None,
    required=False,
    help="Reference (interpolated) EOP values from costg benchmark (01earthRotation_interpolatedEOP.txt).",
)

parser.add_argument(
    "--diffs",
    action="store_true",
    dest="plot_diffs",
    help="Plot diffs instead of results.",
)


def parseRefEops(fn):
    retlst = []
    with open(fn, "r") as fin:
        for line in fin.readlines():
            if len(line) > 5:
                try:
                    # gps time [mjd], xp, yp, s' [rad], dUT1, LOD [seconds], X, Y, s [rad]
                    # t, xp, yp, sp, dut1, lod, Xcip, Ycip, s = [float(x) for x in line.split()]
                    retlst.append([float(x) for x in line.split()])
                except:
                    pass
    return retlst


def arcseconds(radians):
    return radians * (180.0 * 3600.0) / np.pi


def match(l1, l2):
    matches = []

    def grep(data, t, tol=1e-12):
        return [sublist for sublist in data if abs(sublist[0] - t) < tol]

    for ref in l1:
        cmp = grep(l2, ref[0])
        if cmp != []:
            cmp = cmp[0]
            matches.append([ref[0]] + [x[0] - x[1] for x in zip(ref[1:], cmp[1:])])
            print(f"Ref:{ref}")
            print(f"Cmp:{cmp}")
            print(f"Diffs:{matches[-1]}")
    return matches


if __name__ == "__main__":
    args = parser.parse_args()

    dps = 0.3

    liblst = []
    for line in sys.stdin:
        try:
            # gps time [mjd], xp, yp, s' [rad], dUT1, LOD [seconds], X, Y, s [rad]
            # t, xp, yp, sp, dut1, lod, Xcip, Ycip, s = [float(x) for x in line.split()]
            l = line.split()
            liblst.append(
                [int(l[0]) + float(l[1]) / 86400e0] + [float(x) for x in l[2:]]
            )
        except:
            pass

    csglst = parseRefEops(args.costg_fn)

    rdata = csglst
    ldata = liblst
    if args.plot_diffs:
        rdata = match(rdata, ldata)
        ldata = rdata

    ylabel_ang = r"[arc$^{\prime\prime}$]"
    ylabel_clk = r"[$^{\prime\prime}$]"
    with PdfPages("foo.pdf") as pdf:
        fig, ax = plt.subplots(2, 1, sharex=True)
        ax[0].scatter(
            [x[0] for x in rdata],
            [arcseconds(x[1]) for x in rdata],
            s=dps,
            label="_nolegend_",
        )
        ax[0].scatter(
            [x[0] for x in ldata],
            [arcseconds(x[1]) for x in ldata],
            s=dps,
            label="_nolegend_",
        )
        ax[1].scatter(
            [x[0] for x in rdata],
            [arcseconds(x[2]) for x in rdata],
            s=dps,
            label="reference",
        )
        ax[1].scatter(
            [x[0] for x in ldata],
            [arcseconds(x[2]) for x in ldata],
            s=dps,
            label="library",
        )
        ax[0].grid(True)
        ax[1].grid(True)
        ax[0].set_ylabel(ylabel_ang)
        ax[1].set_ylabel(ylabel_ang)
        ax[1].set_xlabel("Epoch")
        ax[0].set_title(r"$x_p$, $y_p$")
        if not args.plot_diffs:
            ax[1].legend()
        fig.subplots_adjust(hspace=0)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(2, 1, sharex=True)
        ax[0].scatter([x[0] for x in rdata], [arcseconds(x[3]) for x in rdata], s=dps)
        ax[0].scatter([x[0] for x in ldata], [arcseconds(x[3]) for x in ldata], s=dps)
        ax[1].scatter([x[0] for x in rdata], [arcseconds(x[8]) for x in rdata], s=dps)
        ax[1].scatter([x[0] for x in ldata], [arcseconds(x[8]) for x in ldata], s=dps)
        ax[0].grid(True)
        ax[1].grid(True)
        ax[0].set_ylabel(ylabel_ang)
        ax[1].set_ylabel(ylabel_ang)
        ax[1].set_xlabel("Epoch")
        ax[0].set_title(r"$s$, $sp$")
        fig.subplots_adjust(hspace=0)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(2, 1, sharex=True)
        ax[0].scatter([x[0] for x in rdata], [arcseconds(x[6]) for x in rdata], s=dps)
        ax[0].scatter([x[0] for x in ldata], [arcseconds(x[6]) for x in ldata], s=dps)
        ax[1].scatter([x[0] for x in rdata], [arcseconds(x[7]) for x in rdata], s=dps)
        ax[1].scatter([x[0] for x in ldata], [arcseconds(x[7]) for x in ldata], s=dps)
        ax[0].grid(True)
        ax[1].grid(True)
        ax[0].set_ylabel(ylabel_ang)
        ax[1].set_ylabel(ylabel_ang)
        ax[1].set_xlabel("Epoch")
        ax[0].set_title(r"$X_{cip}$, $Y_{cip}$")
        fig.subplots_adjust(hspace=0)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(2, 1, sharex=True)
        ax[0].scatter([x[0] for x in rdata], [x[4] for x in rdata], s=dps)
        ax[0].scatter([x[0] for x in ldata], [x[4] for x in ldata], s=dps)
        ax[1].scatter([x[0] for x in rdata], [x[5] for x in rdata], s=dps)
        ax[1].scatter([x[0] for x in ldata], [x[5] for x in ldata], s=dps)
        ax[0].grid(True)
        ax[1].grid(True)
        ax[0].set_ylabel(ylabel_clk)
        ax[1].set_ylabel(ylabel_clk)
        ax[1].set_xlabel("Epoch")
        ax[0].set_title("dUT1, LOD")
        fig.subplots_adjust(hspace=0)
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
