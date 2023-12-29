#! /usr/bin/python

import argparse
import ftplib
import datetime

## GFZ AOD1B ftp site
gfz_ftp = "isdcftp.gfz-potsdam.de"
generic_dir = "

## t: date or datetime
## tidal: 0->nontidal, 1->tidal
## rl: RL (int)
## 
def product_name(**kwargs):
    if kwargs['tidal'] == 0:
        return "AOD1B_DATE_X_RL.asc.gz".
            replace("DATE", "{:}".format(kwarg['t'].strftime("%Y-%m-%d"))).
            replace("RL", "{:02d}".format(kwargs['rl'))
    

#isdcftp.gfz-potsdam.de
#/grace/Level-1B/GFZ/AOD/RL06/2023
#AOD1B_2023-12-27_X_06.asc.gz
#AOD1B YYYY-MM-DD S RL.EXT.gz (non-tidal)
#
#ftp://isdcftp.gfz-potsdam.de/grace/Level-1B/GFZ/AOD/RL06/TIDES
#

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Download AOD1B Product files.',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Dec, 2024'''))

parser.add_argument(
    '-f',
    '--from',
    metavar='FROM',
    dest='tfrom',
    default=None,
    required=True,
    help='Starting date')
