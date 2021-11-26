#! /usr/bin/python

import sys
import datetime
import julian

for line in sys.stdin:
    if line.startswith('*  '):
        l = line.split()
        y, m, d, h, mn = [ int(x) for x in l[1:6] ]
        d = datetime.datetime(y, m, d, h, mn, int(float(l[6])))
        jd = julian.to_jd(d)
    elif line.startswith('PL'):
        l = line.split()
        x,y,z = [float(x)*1e3 for x in l[1:4] ]
        print("{:.3f} {:.9f} {:.9f} {:.9f} {:.9f}".format(2451545e0, jd-2451545e0, x, y, z))
