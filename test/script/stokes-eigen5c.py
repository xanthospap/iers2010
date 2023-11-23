#! /usr/bin/python

import math
import datetime
import sys

## ---------------------------------------------------------------------------
## Dummy script to compute C(4,0) and S(4,0) for given epochs using the 
## EIGEN-5C model.
## This script is meant to assist testing of parsing ICGEM files.
## ---------------------------------------------------------------------------

def cs40(t):
#gfct   4    0 0.539988144712e-06 0.000000000000e+00 0.1095e-10 0.0000e+00 20041001
#dot    4    0 0.470000000000e-11 0.000000000000e+00 0.0000e+00 0.0000e+00
    t0 = datetime.datetime.strptime('20041001', '%Y%m%d')
    diff = t - t0
    dt = (diff.days + diff.seconds/86400e0) / 365.25e0
    C = 0.539988144712e-06 + 0.470000000000e-11 * dt
    S = 0e0
    return C,S

for year in range(1995, 2025):
    c,s = cs40(datetime.datetime(year,1,1))
    print("C = {:+.15e} S = {:+.15e}".format(c,s))
