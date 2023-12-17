#! /usr/bin/python
# -*- coding: utf-8 -*-

##
## Read the file tables6_iers2010.csv and extract the (three) tables
## as ready-to-use C++ style arrays of type Step2TidesCoeffs.
## The output of this script can be directly used in the
## iers2010_step2.cpp source code file.
##

import sys
delimeter = ';'

def iers_line_to_cpp(line):
    l = line.split(delimeter)
# Doodson Number;τ;s;h;p;N’;ps;l;l’;F;D;Ω;Re(δkf)*1e-5;Im(δkf)*1e-5 ;In-Phase Amp * 1e-12;Out-Of-Phase Amp * 1e-12
# 0              1 2 3 4 5  6  7 8  9 101112           13            14                   15
    assert( len(l) == 16 )
    return ('{{/*{:}*/ {:}, {:+.8e}, {:+.8e} }},'.format(l[0], ','.join(['{:2d}'.format(int(x)) for x in l[1:7]]), float(l[14]), float(l[15]))) 

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('ERROR Must provide a CSV file!', file=sys.stderr)
        sys.exit(1)

    num_freqs = 0
    with open(sys.argv[1], 'r') as fin:
        line  = fin.readline()
# assert(header)
        if line.replace(' ','').strip() != 'DoodsonNumber;τ;s;h;p;N’;ps;l;l’;F;D;Ω;Re(δkf)*1e-5;Im(δkf)*1e-5;In-PhaseAmp*1e-12;Out-Of-PhaseAmp*1e-12':
            print('ERROR Failed to validate header line in CSV file!', file=sys.stderr)
            print('Expected:DoodsonNumber;τ;s;h;p;N’;ps;l;l’;F;D;Ω;Re(δkf)*1e-5;Im(δkf)*1e-5;In-PhaseAmp*1e-12;Out-Of-PhaseAmp*1e-12')
            print('Found   :{:}'.format(line.replace(' ',''))) 
            sys.exit(1)
        for line in fin.readlines():
            print(iers_line_to_cpp(line))
            num_freqs += 1
    print('Number of records: {:}'.format(num_freqs))
