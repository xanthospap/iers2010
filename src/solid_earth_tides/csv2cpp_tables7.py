#! /usr/bin/python
# -*- coding: utf-8 -*-

##
## Read the file tables7_iers2010.csv and extract the (three) tables
## as ready-to-use C++ style arrays of type 
##

import sys
delimeter = ';'

def iers_line_to_cpp(line):
    l = line.split(delimeter)
# Frequency;Doodson;τ;s;h;p;N’;ps;l;l’;F;D;Ω;Δrf_ip;Δrf_op;Δtf_ip;Δtf_op
# 0         1       2 3 4 5 6  7  8 9  10111213     14     15     16
    assert( len(l) == 17 )
    return ('{{/*{:}*/ {{{:}}}, {:+.4e}, {:+.4e}, {:+.4e}, {:+.4e} }},'.format(l[1], ','.join(['{:2d}'.format(int(x)) for x in l[2:8]]), float(l[13]), float(l[14]), float(l[15]), float(l[16]))) 

if __name__ == '__main__':

    if len(sys.argv) != 2:
        print('ERROR Must provide a CSV file!', file=sys.stderr)
        sys.exit(1)

    num_freqs = 0
    with open(sys.argv[1], 'r') as fin:
        line  = fin.readline()
# assert(header)
        if line.replace(' ','').strip() != 'Frequency;Doodson;τ;s;h;p;N’;ps;l;l’;F;D;Ω;Δrf_ip;Δrf_op;Δtf_ip;Δtf_op':
            print('ERROR Failed to validate header line in CSV file!', file=sys.stderr)
            print('Expected:Frequency;Doodson;τ;s;h;p;N’;ps;l;l’;F;D;Ω;Δrf_ip;Δrf_op;Δtf_ip;Δtf_op')
            print('Found   :{:}'.format(line.replace(' ',''))) 
            sys.exit(1)
        for line in fin.readlines():
            print(iers_line_to_cpp(line))
            num_freqs += 1
    print('Number of records: {:}'.format(num_freqs))
