#!/usr/bin/python3

import os
import sys
import subprocess
import argparse
import ftplib
import importlib.util
import urllib.request
import json

def get_data_files(data_files_dct):
    for d in data_files_dct:
        if not os.path.isfile(d["local"]):
            try:
                print("Downloading data file {:} to {:}".format(d["url"], d["local"]))
                urllib.request.urlretrieve(d["url"], d["local"])
            except:
                print("ERROR. Failed to download data file {:}".format(d["url"]), file=sys.stderr)
                sys.exit(1)

def prog_needs_args(prog_path, special_progs_dct):
    prog_name = os.path.basename(prog_path)
    for d in special_progs_dct:
        if prog_name == d['prog']:
            return True
    return False

def check_file_vs_str(file, str):
    with open(file, 'r') as fin:
        fstr = fin.read()
        return (fstr == str)

def run_progs_with_args(special_progs_dct):
    for d in special_progs_dct:
        if not os.path.isfile(os.path.join(d['path'],d['prog'])):
            print('ERROR Failed to find executable {:}'.format(d['prog']), file=sys.stderr)
            sys.exit(1)
        exe = '{:}'.format(os.path.join(d['path'],d['prog'])) 
        cmdargs = ['{:}'.format(x) for x in d['args']]
        print('Running command {:}'.format(' '.join([exe]+cmdargs)))
        with open('.tmp.result', 'w') as fout:
            result = subprocess.run([exe] + cmdargs, stdout=fout, stderr=subprocess.STDOUT, check=False)
        if result.returncode != int(d['exit']):
            print('ERROR Expected exit code {:} and got {:}; program: {:}'.format(d['exit'], result.returncode, exe), file=sys.stderr)
            sys.exit(1)
        if 'results' in d:
            if not check_file_vs_str('.tmp.result', d['results']):
                print('ERROR. Failed test {:}; expected string did not match output'.format(exe), file=sys.stderr)
                sys.exit(1)

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
    dest='reffn',
    default=os.path.join(
        os.path.abspath(
            os.getcwd()),
        'data/reference_results.py'),
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

# run tests with command line arguments
    with open(args.reffn, "rb") as source_file:
        code = compile(source_file.read(), args.reffn, "exec")
    exec(code)
    get_data_files(data_files)
    run_progs_with_args(special_progs)

# run tests recursively in each-subdir
    for sdir in test_subdirs:
        path = os.path.join(args.progs_dir, sdir)
        if not os.path.isdir(path):
            print('ERROR Expected to find dir {:} but couldn\'t!'.format(path), file=sys.stderr)
            sys.exit(1)
        for file in os.listdir(path):
            if file.endswith('.out') and not prog_needs_args(file, special_progs):
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
