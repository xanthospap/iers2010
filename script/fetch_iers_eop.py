#! /usr/bin/python3

##
## Script to download EOP C04 (14/20) file(s) from IERS
##
import os
import sys
import shutil
import argparse
import urllib.request

URL = "https://hpiers.obspm.fr/iers/eop/"

available_files = {
    "eopc04.12h.1984-now": {"format": "c04", "series": 20, "interval": 12, "dir": "eopc04"},
    "eopc04.1962-now": {"format": "c04", "series": 20, "interval": 24, "dir": "eopc04"},
    "eopc04.62-now": {"format": "c04", "series": 14, "interval": 24, "dir": "eopc04_14"}
}

def choose_eop(**kwargs):
    erno = 0
    for istr in [ "format", "series", "interval" ]:
        if istr not in kwargs:
            print("ERROR. Argument \"{:}\" not provided!".format(istr), 
                file=sys.stderr)
            erno += 1
    if erno:
        raise RuntimeError("ERROR No file downloaded\n")

    def match_dicts(d1,d2):
        for e in d1.items():
            if e[0] in d2 and d2[e[0]] != e[1]:
                return False
        return True

    for f in available_files.items():
        fname = f[0]
        dct   = f[1]
        if match_dicts(kwargs, dct):
            return os.path.join(os.path.join(URL, dct["dir"]), fname)

    return None

def download(**kwargs):
    target = choose_eop(**kwargs)
    if not target:
        raise RuntimeError("ERROR. Failed to match EOP file to input options\n")
    if "save_as" in kwargs and kwargs["save_as"] is None:
        save_as = os.path.join(os.getcwd(), os.path.basename(target))
    else:
        save_as = kwargs["save_as"] if "save_as" in kwargs else os.path.join(os.getcwd(), os.path.basename(target))
    with urllib.request.urlopen(target) as remote, open(save_as, 'w') as local:
        html = remote.read().decode('utf-8')
        print(html, file=local)

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Download IERS C04 file(s) from IERS',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Jul, 2023'''))

parser.add_argument(
    '-o',
    '--output',
    metavar='OUTPUT',
    dest='save_as',
    default=None,
    required=False,
    help='Save the downloaded file using this file(name)')

parser.add_argument(
    '--format',
    metavar='FORMAT',
    dest='format',
    required=False,
    default='c04',
    help='The format of the EOP file to download. Currently only \'c04\' is an available option.')

parser.add_argument(
    '--series',
    metavar='SERIES',
    dest='series',
    required=False,
    default=20,
    type=int,
    help='The series of the EOP file to download. Available options are: \'14\' and \'20\'.')

parser.add_argument(
    '--interval',
    metavar='INTERVAL',
    dest='interval',
    required=False,
    default=24,
    type=int,
    help='The data interval of the EOP file to download in hours. Available options are: \'12\' and \'24\'.')

parser.add_argument(
    '--verbose',
    action='store_true',
    dest='verbose',
    help='Verbose mode on')

if __name__ == "__main__":
    ## parse cmd
    args = parser.parse_args()
    download(**vars(args))
