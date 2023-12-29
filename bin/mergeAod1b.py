#! /usr/bin/python

##
## TODO Verify the merged AOD1B file, i.e. make a harsh check on its integrity 
##      and validity.
##

import os, sys
import argparse
import datetime
import copy
import re
import shutil
from argparse import RawTextHelpFormatter

def readAod1bHeader(fn):
    dct = {}
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            if line.strip() == "END OF HEADER":
                break
            else:
                l = line.split(':', 1)
                key = l[0].rstrip()
                val = l[1].strip()
                dct[key] = val
    return dct

def catAod1bData(fn, fout):
    with open(fn, 'r') as fin:
        while True:
            line = fin.readline()
            if line.strip() == "END OF HEADER":
                break
        s = fin.read()
    fout.write(s)

class Aod1b:
    def __init__(self, fn):
        self.fn = fn
        self.dct = readAod1bHeader(fn)
    def timeEpoch(self):
        return datetime.datetime.strptime(self.dct['TIME EPOCH (GPS TIME)'], "%Y-%m-%d %H:%M:%S")
    def firstObs(self):
        val = self.dct['TIME FIRST OBS(SEC PAST EPOCH)']
        d = re.match(r'[0-9]*.[0-9]* *\(([0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2})', val)[1]
        return datetime.datetime.strptime(d, "%Y-%m-%d %H:%M:%S")
    def lastObs(self):
        val = self.dct['TIME LAST OBS(SEC PAST EPOCH)']
        d = re.match(r'[0-9]*.[0-9]* *\(([0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2})', val)[1]
        return datetime.datetime.strptime(d, "%Y-%m-%d %H:%M:%S")
    def numDataSets(self):
        key = self.dct['NUMBER OF DATA SETS']
        return int(key)
    def numDataRecords(self):
        key = self.dct['NUMBER OF DATA RECORDS']
        return int(key)
    def headersMatch(self, other):
        for k,v in self.dct.items():
            if k not in ['TIME FIRST OBS(SEC PAST EPOCH)', 'TIME LAST OBS(SEC PAST EPOCH)', 'PRODUCT START CREATE TIME(UTC)', 'PRODUCT END CREATE TIME(UTC)', 'FILESIZE (BYTES)', 'FILENAME']:
                if v != other.dct[k]:
                    print('Warning! Mismatch in AOD1B headers! Field is {:}'.format(k))
                    return False
        return True
    def merge(self, other):
        if not self.headersMatch(other):
            print('Cannot merger file {:} and {:}; Header mismatch!'.format(self.fn, other.fn))
            return None
        first, second = self, other
        if first.firstObs() > second.firstObs():
            first, second = other, self
        assert(first.lastObs() < second.firstObs())
        cpy = copy.deepcopy(self)
        microsec = (second.lastObs() - first.timeEpoch()).total_seconds() * 1e6
        cpy.dct['TIME LAST OBS(SEC PAST EPOCH)'] = '{:.6f} ({:})'.format(microsec*1e-6, second.lastObs().strftime('%Y-%m-%d %H:%M:%S'))
        cpy.dct['NUMBER OF DATA RECORDS'] = ' {:}'.format(first.numDataRecords() + second.numDataRecords())
        cpy.dct['NUMBER OF DATA SETS'] = ' {:}'.format(first.numDataSets() + second.numDataSets())
        return cpy
    def isPriorTo(self, other):
        return self.lastObs() < other.firstObs();
    def __lt__(self, other): return self.isPriorTo(other)
    def dump(self, ofile=sys.stdout):
        for k,v in self.dct.items():
            print('{:30s}: {:}'.format(k,v), file=ofile)

def sortAod1bFileList(file_list):
    Aod1bFileList = []
    for f in file_list: Aod1bFileList.append(Aod1b(f))
    sortedAod1b = sorted(Aod1bFileList)
    return [ f.fn for f in sortedAod1b ]

parser = argparse.ArgumentParser(
    formatter_class=RawTextHelpFormatter,
    description=
    'Merge AOD1B Product file(s).',
    epilog=(r'''National Technical University of Athens,
    Dionysos Satellite Observatory
    Send bug reports to:
    Xanthos Papanikolaou, xanthos@mail.ntua.gr
    Dec, 2024'''))

parser.add_argument(
    '-f',
    '--files',
    metavar='AOD1B_FILES',
    dest='files',
    default=None,
    required=True,
    nargs = '+',
    help='White-space seperated list of AOD1B files to merge.')

parser.add_argument(
    '-o',
    '--output-dir',
    metavar='OUT_DIR',
    dest='out_dir',
    default=os.getcwd(),
    required=False,
    help='Directory where the downloaded file(s) will be saved.')

if __name__ == "__main__":
    args = parser.parse_args()
    tmpout = '.aod1bdata'
    
    sfiles = sortAod1bFileList(args.files)
    with open(tmpout, 'w') as fout:
        merged = Aod1b(sfiles[0])
        catAod1bData(sfiles[0], fout)
        for f in sfiles[1:]:
            merged = merged.merge(Aod1b(f))
            catAod1bData(f, fout)

    outf = os.path.join(args.out_dir, sfiles[0].replace('.asc', '.merged.asc'))
    with open(outf, 'w') as fout: 
        merged.dump(fout)
        with open(tmpout, 'r') as fin:
            fout.write(fin.read())
