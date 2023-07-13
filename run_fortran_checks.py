#! /usr/bin/python3

import os
import sys
import argparse
import re
import subprocess

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Run validatation programs against IERS2010 Fortran library',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    July, 2023'''))

parser.add_argument(
    '--cpp-dir',
    metavar='CPP_DIR',
    dest='cpp_dir',
    default=os.path.abspath(os.getcwd()),
    required=False,
    help='Directory with C++ test executables')

parser.add_argument(
    '--for-dir',
    metavar='FOR_DIR',
    dest='for_dir',
    default=os.path.abspath(os.getcwd()),
    required=False,
    help='Directory with FORTRAN subroutines')

parser.add_argument(
    '--data-dir',
    metavar='DATA_DIR',
    dest='data_dir',
    default=os.path.abspath(os.getcwd()),
    required=False,
    help='Directory with test data')

parser.add_argument(
    '--markdown',
    action='store_true',
    dest='markdown',
    help='Print results in Markdown format (matrix)')

parser.add_argument(
    '--verbose',
    action='store_true',
    dest='verbose',
    help='Verbose mode on')

# Fortran subroutines
for_sr = ['FCNNUT.F']
def for_sr2fc(for_sr): return 'py'+os.path.splitext(for_sr.lower())[0]
def for_sr2py(for_sr): return os.path.splitext(for_sr.lower())[0] + '.py'

# C++ test progs
cpp_sr = ['test-fcnnut.out']

if __name__ == '__main__':

    ## parse cmd
    args = parser.parse_args()

    verboseprint = print if args.verbose else lambda *a, **k: None

    for s in zip(for_sr, cpp_sr):
# Compile Fortran subroutine
        sr = s[0]
        fsrc = os.path.join(args.for_dir, sr)
        if not os.path.isfile(fsrc):
            printf("ERROR. Failed to locate file {:}".format(fsrc))
            sys.exit(1)
        cfsrc = for_sr2fc(sr)
        exe = 'f2py3 -c {:} -m {:}'.format(fsrc, cfsrc).split()
        verboseprint('Running command: [{:}]'.format(exe))
        result = subprocess.run(exe, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=False)
## run the C++ test program; STDOUT to '.tmp'
        sr = s[1]
        cprog = os.path.join(args.cpp_dir, sr)
        exe = '{:}'.format(cprog).split()
        verboseprint('Running command: [{:}]'.format(exe))
        with open('.tmp', 'w') as fout:
            result = subprocess.run(exe, stdout=fout, stderr=subprocess.STDOUT, check=False)
## run the python script
        sr = s[0]
        pprog = os.path.join(args.for_dir, for_sr2py(sr))
        exe = '{:} {:}'.format(pprog, '.tmp').split()
        verboseprint('Running command: [{:}]'.format(exe))
        subprocess.run(exe, check=False)
