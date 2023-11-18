#!/usr/bin/python3

import os
import sys
import subprocess
import argparse
import ftplib
import importlib.util
import json

class myFormatter(argparse.ArgumentDefaultsHelpFormatter,
                  argparse.RawTextHelpFormatter):
    pass


parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description='Run validatation programs for the iers2010 library',
    epilog=('''National Technical University of Athens,
    Dionysos Satellite Observatory\n
    Send bug reports to:
    Xanthos, xanthos@mail.ntua.gr
    Nov, 2023'''))

parser.add_argument(
    '--progs-dir',
    metavar='PROGS_DIR',
    dest='progs_dir',
    default=os.path.abspath(os.path.join(os.getcwd(), 'test')),
    required=False,
    help='Directory with test executables (top-level test directory)')

parser.add_argument(
    '--data-dir',
    metavar='DATA_DIR',
    dest='data_dir',
    default=os.path.abspath(os.getcwd()),
    required=False,
    help='Directory with test data')

parser.add_argument(
    '--reference-results',
    metavar='REFERENCE_RESULTS_FILE',
    dest='ref_respy',
    default=os.path.join(
        os.path.abspath(
            os.getcwd()),
        'test/reference_results.py'),
    required=False,
    help='File with reference test results, for comparing against')

parser.add_argument(
    '--verbose',
    action='store_true',
    dest='verbose',
    help='Verbose mode on')

test_subdirs = ['unit_tests']

if __name__ == '__main__':
    # parse cmd
    args = parser.parse_args()

# verbose print
    verboseprint = print if args.verbose else lambda *a, **k: None

# run tests recursively in each-subdir
    for sdir in test_subdirs:
        path = os.path.join(args.progs_dir, sdir)
        if not os.path.isdir(path):
            print('ERROR Expected to find dir {:} but couldn\'t!'.format(path), file=sys.stderr)
            sys.exit(1)
        for file in os.listdir(path):
            if file.endswith('.out'):
                exe = os.path.join(path, file)
                verboseprint('Running command: {:}'.format([exe]))
# assert a 0 exit code
                subprocess.run([exe], stdout=sys.stderr, check=True)
                print('Test {:} OK!'.format(exe))
            elif file.endswith('.json'):
                with open(os.path.join(path, file), 'r') as jin: data = json.load(jin)
                for prog, dct in data.items():
                    exe = '{:}'.format(dct['cxx']) 
                    cmdargs = ['{:}'.format(dct['flags']), '{:}'.format(dct['name']), '-I{:}'.format(dct['incp']), ]
                    verboseprint('Running command {:}'.format(' '.join([exe]+cmdargs)))
                    result = subprocess.run([exe] + cmdargs, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=False)
                    if result.returncode != int(dct['exit']):
                        print('ERROR Expected a return code {:} and got {:}; exe {:}'.format(
                            dct['exit'], result.returncode, prog), file=sys.stderr)
                        sys.exit(2)
                    else:
                        print('Test {:} OK! (failed to compile)'.format(prog))
