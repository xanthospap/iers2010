#! /usr/bin/python

import os, sys
import argparse
import ftplib
import datetime
import gzip
import shutil
from argparse import RawTextHelpFormatter

## GFZ AOD1B ftp site
gfz_ftp = "isdcftp.gfz-potsdam.de"

## t: date or datetime
## tidal: 0->nontidal, 1->tidal
## rl: RL (int)
## 
def product_name(**kwargs):
    if kwargs['tidal'] == 0:
        fn = 'AOD1B_DATE_X_RL.asc.gz'.replace('DATE', '{:}'.format(kwargs['t'].strftime('%Y-%m-%d'))).replace('RL', '{:02d}'.format(kwargs['rl']))
        uri = '/grace/Level-1B/GFZ/AOD/RLNN/YYYY'.replace('YYYY', '{:}'.format(kwargs['t'].strftime('%Y'))).replace('NN', '{:02d}'.format(kwargs['rl']))
        url = gfz_ftp ##'isdcftp.gfz-potsdam.de'
    return { 'fn': fn, 'dir': uri, 'url': url }

def ftpget(server, fdir, fn, outdir, saveas=None):
    saveas = fn if saveas is None else saveas
    saveas = os.path.join(outdir, saveas)
    ftp = ftplib.FTP(server)
    ftp.login()
    ftp.cwd(fdir)
    #ftp.retrlines('LIST')
    ftp.retrbinary("RETR {}".format(fn), open(saveas, 'wb').write)
    ftp.quit()
    return saveas

#isdcftp.gfz-potsdam.de
#/grace/Level-1B/GFZ/AOD/RL06/2023
#AOD1B_2023-12-27_X_06.asc.gz
#AOD1B YYYY-MM-DD S RL.EXT.gz (non-tidal)
#
#ftp://isdcftp.gfz-potsdam.de/grace/Level-1B/GFZ/AOD/RL06/TIDES
#

parser = argparse.ArgumentParser(
    formatter_class=RawTextHelpFormatter,
    description=
    'Download AOD1B Product file(s).',
    epilog=(r'''National Technical University of Athens,
    Dionysos Satellite Observatory
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Dec, 2024'''))

parser.add_argument(
    '-f',
    '--from',
    metavar='START_DATE',
    dest='tfrom',
    default=None,
    required=True,
    type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'),
    help='Starting date as YYYY-MM-DD.')

parser.add_argument(
    '-t',
    '--to',
    metavar='END_DATE',
    dest='tend',
    default=None,
    required=False,
    type=lambda s: datetime.datetime.strptime(s, '%Y-%m-%d'),
    help='Ending date as YYYY-MM-DD (exlusive).')

parser.add_argument(
    '-n',
    '--rl',
    metavar='RLnn',
    dest='rlnn',
    default=6,
    type=int,
    required=False,
    help='RL number (e.g. 6 for RL06).')

parser.add_argument(
    '-o',
    '--output-dir',
    metavar='OUT_DIR',
    dest='out_dir',
    default=os.getcwd(),
    required=False,
    help='Directory where the downloaded file(s) will be saved.')

parser.add_argument(
    '--tidal',
    dest='is_tidal',
    action='store_true',
    help='Download tidal product (default is non-tidal).')

parser.add_argument(
    '--decompress',
    dest='dcmp',
    action='store_true',
    help='Decompress the downloaded file(s).')

if __name__ == "__main__":
    args = parser.parse_args()
# if ending date not give, set it to starting date plus one day
    args.tend = args.tend if args.tend is not None else args.tfrom + datetime.timedelta(days=1)
# download products, one file/one day at a time
    t = args.tfrom
    while t < args.tend:
        dct = product_name(**{'t':t, 'rl': args.rlnn, 'tidal': args.is_tidal})
# if file already there, skip download
        if os.path.isfile(os.path.join(args.out_dir, dct['fn'])) or os.path.isfile(os.path.join(args.out_dir, dct['fn'].replace('.gz',''))):
            print('Note: skipping download; target file already present')
            saveas = os.path.join(args.out_dir, dct['fn']) if os.path.isfile(os.path.join(args.out_dir, dct['fn'])) else os.path.join(args.out_dir, dct['fn'].replace('.gz',''))
        else:
            saveas = ftpget(dct['url'], dct['dir'], dct['fn'], args.out_dir)
# double-check we downloaded the target file
        if not os.path.isfile(saveas) and not os.path.isfile(saveas.replace('.gz', '')):
            print('ERROR. Failed to fetch product file {:} from {:}@{:}'.format(dct['fn'], dct['dir'], dct['url']))
            sys.exit(9)
# do we need to decompress ?
        if saveas.endswith('.gz') and args.dcmp:
            ns = saveas[0:-3]
            with gzip.open(saveas, 'rb') as fin:
                with open(ns, 'wb') as fout:
                    shutil.copyfileobj(fin, fout)
            saveas = ns
# next day
        t += datetime.timedelta(days=1)
