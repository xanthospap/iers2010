#! /usr/bin/python

import os
import sys
import argparse
import re
import subprocess

def find_progs(prog_dir):
    return [f for f in os.listdir(prog_dir) if re.match(r'sofa-.*\.out', f)]

def filter_results(fn):
    argl=[];numtestsl=[];numfailsl=[];maxerrorl=[];statusl=[];
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            if line.startswith("Function") or line.startswith("Program") or line.startswith("---"):
                pass
            else:
# pn06/   dpsi  10000   3347 +2.078386880e-12 FAILED
                l = line.split()
                function = l[0]
                if len(l) == 6:
                    argl.append(l[1])
                    offset = 2
                else:
                    argl.append("")
                    offset = 1
                numtestsl.append(int(l[offset]))
                numfailsl.append(int(l[offset+1]))
                maxerrorl.append(float(l[offset+2]))
                statusl.append(l[offset+3])
    return function, [{'arg':e[0], 'num_tests': e[1], 'num_fails':e[2], 'max_error':e[3], 'status':e[4]} for e in zip(argl,numtestsl,numfailsl,maxerrorl,statusl)]

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass

parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=
    'Run validatation programs against SOFA library',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    May, 2023'''))

parser.add_argument(
    '--progs-dir',
    metavar='PROGS_DIR',
    dest='progs_dir',
    default=os.path.abspath(os.getcwd()),
    required=False,
    help='Directory with SOFA test executables')

parser.add_argument(
    '--plots-dir',
    metavar='PLOTS_DIR',
    dest='plots_dir',
    required=False,
    default='figures',
    help='Directory to save plots at.')

parser.add_argument(
    '--factor',
    metavar='FACTOR',
    dest='factor',
    type=float,
    required=False,
    default=1e12,
    help='Factor to be used when plotting, e.g. 1e9')

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

if __name__ == '__main__':

    ## parse cmd
    args = parser.parse_args()

    verboseprint = print if args.verbose else lambda *a, **k: None

    progs_no_path = find_progs(args.progs_dir)
    if len(progs_no_path) == 0:
        print('No validation programs found! Nothing to do ...')
        sys.exit(0)

    temp_fn = ".tmp"
    results = {}
    for prog in progs_no_path:
        exe = os.path.join(args.progs_dir, prog)
        #for vp in ["xy06", "s06", "c2ixys", "era00", "sp00", "pom00"]:
        #    if vp in prog:
        verboseprint('Running command: [{:}]'.format(exe))
        ftmp = open(temp_fn, "w")
        result = subprocess.run(exe, stdout=ftmp, check=False)
        ftmp.close()
        fun,fargs = filter_results(temp_fn)
        results[fun] = fargs

    if args.markdown:
        print('{:10s}|{:10s}|{:10s}|{:10s}|{:10s}|{:10s}'.format("function", "argument", "num tests", "num fails", "max error", "status"))
        print('{:10s}|{:10s}|{:10s}|{:10s}|{:10s}|{:10s}'.format("-"*10,"-"*10,"-"*10,"-"*10,"-"*10,"-"*10))
    for f,fargs in results.items():
        for farg in fargs:
            if args.markdown:
                print('{:10s}|{:10s}|{:10d}|{:10d}|{:.3e}|{:10s}'.format(f, farg['arg'], farg['num_tests'], farg['num_fails'], farg['max_error'], farg['status']))
            else:
                print('{:6s} {:6s} {:8d} {:8d} {:.3e} {:10s}'.format(f, farg['arg'], farg['num_tests'], farg['num_fails'], farg['max_error'], farg['status']))
