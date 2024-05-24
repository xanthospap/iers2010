#!/usr/bin/python3

import os
import sys
import subprocess
import argparse
import ftplib
import importlib.util
import urllib.request
import json
import gzip
import shutil

def get_data_files(data_files_dct):
    for d in data_files_dct:
        if not os.path.isfile(d["local"]):
            try:
                print("Downloading data file {:} to {:}".format(d["url"], d["local"]))
                urllib.request.urlretrieve(d["url"], d["local"])
            except:
                print("ERROR. Failed to download data file {:}".format(d["url"]), file=sys.stderr)
                sys.exit(1)
            if 'actions' in d:
                actions = d['actions']
                for action in actions:
                    if action == 'decompress':
                        with gzip.open(d['local'], 'rb') as f_in:
                            with open(d['local'][0:-3], 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)

def prog_needs_args(prog_path, special_progs_dct):
    prog_name = os.path.basename(prog_path)
    for d in special_progs_dct:
        if prog_name == d['prog']:
            return True
    return False

def report(prog_name, exit_status, comment):
    print('{:50s} {:7s} {:20s}'.format(prog_name, exit_status, comment))

def check_file_vs_str(file, vstr):
    with open(file, 'r') as fin:
        fstr = fin.read().replace(' ', '').replace('\n','')
    return (fstr == vstr.replace(' ','').replace('\n',''))

def check_file_vs_str_verbose(file, vstr):
    with open(file, 'r') as fin:
        fstr = fin.read().replace(' ', '').replace('\n','')
    rstr = vstr.replace(' ','').replace('\n','')
    for s in zip(fstr,rstr):
        print('[{:}/{:}]'.format(s[0],s[1]))
        if s[0] != s[1]: return

def check_w_precision(fn1, fn2str, cols, precision):
    ## just for fun, collect and report actual precision
    accuracy = {}
    for c in cols: accuracy[c] = 0e0
    with open(fn1, 'r') as f1in: f1lines = f1in.readlines()
    f2lines = fn2str.splitlines()
    line_nr = 1
    for l12 in zip(f1lines, f2lines):
        l1 = l12[0].split()
        l2 = l12[1].split()
        assert(len(l1) == len(l2))
        for c in range(0, len(l1)):
            if c not in cols:
                if l1[c] != l2[c]:
                    print('Result Missmatch! line/col is {:}/{:}'.format(line_nr,c+1))
                assert(l1[c] == l2[c])
        for k in zip(cols,precision):
            if abs(float(l1[k[0]]) - float(l2[k[0]])) > k[1]:
                print('Result Missmatch! line/col is {:}/{:} diff {:.15e} > {:.3e}'.format(line_nr,k[0],abs(float(l1[k[0]]) - float(l2[k[0]])),k[1]))
            if abs(float(l1[k[0]]) - float(l2[k[0]])) > accuracy[k[0]]:
                accuracy[k[0]] = abs(float(l1[k[0]]) - float(l2[k[0]]))
            assert(abs(float(l1[k[0]]) - float(l2[k[0]])) < k[1])
        line_nr += 1
    rep = '[NOTE] Actual/Measured accuracy by column(s): '
    for k,v in accuracy.items():
        rep += ' {:}/{:.15e}'.format(k,v)
    return rep

def run_progs_with_args(special_progs_dct):
    for d in special_progs_dct:
        comment_str = ''
        if not os.path.isfile(os.path.join(d['path'],d['prog'])):
            print('ERROR Failed to find executable {:}@{:}'.format(d['prog'], d['path']), file=sys.stderr)
            sys.exit(1)
        exe = '{:}'.format(os.path.join(d['path'],d['prog'])) 
        cmdargs = ['{:}'.format(x) for x in d['args']]
## run prrogram and store results in '.tmp.result'
        with open('.tmp.result', 'w') as fout:
            result = subprocess.run([exe] + cmdargs, stdout=fout, stderr=subprocess.DEVNULL, check=False)
## if exit is the special keyword 'nzero'
        if d['exit'] == 'nzero':
            if result.returncode == 0:
                print('ERROR Expected exit code other than zero and got {:}; program: {:}'.format(result.returncode, exe), file=sys.stderr)
                sys.exit(1)
            else:
                comment_str += 'Expecting runtime error; successefuly failed'
## check exit code
        elif result.returncode != int(d['exit']):
            print('ERROR Expected exit code {:} and got {:}; program: {:}'.format(d['exit'], result.returncode, exe), file=sys.stderr)
            sys.exit(1)
## check result string vs reference result
        if 'results' in d:
            if 'testwp' in d:
                treport = check_w_precision('.tmp.result', d['results'], d['testwp']['cols'], d['testwp']['precision'])
                print(treport)
                report(d['prog'], 'OK', 'precision ok')
            else:
                if not check_file_vs_str('.tmp.result', d['results']):
                    print('ERROR. Failed test {:}; expected string did not match output'.format(exe), file=sys.stderr)
                    sys.exit(1)
                report(d['prog'], 'OK', '{:}'.join(cmdargs))
        else:
            report(d['prog'], 'OK', '{:}'.join(cmdargs)+comment_str)


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

# Report header    
    print('{:50s} {:7s} {:20s}'.format('Program Name', 'Status', 'Comment'))
    print('{:50s} {:7s} {:20s}'.format('-'*50, '-'*7, '-'*20))

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
# choose all files that end with .out but are not listed in the dictionary
# reference_results::special_progs
        for file in os.listdir(path):
            if file.endswith('.out') and not prog_needs_args(file, special_progs):
                exe = os.path.join(path, file)
                verboseprint('Running command: {:}'.format([exe]))
# assert a 0 exit code
                subprocess.run([exe], stdout=sys.stderr, check=True)
                report(file, 'OK', '')
# encountered the json file - these are mocks, should not compile
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
                        #print('Test {:} OK! (failed to compile)'.format(prog))
                        report(prog, 'OK', 'Mock test; successefuly failed to compile')
