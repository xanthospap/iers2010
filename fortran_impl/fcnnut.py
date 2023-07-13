#! /usr/bin/python3

import sys
import numpy
import pyfcnnut

## f2py -c FCNNUT.F -m pyfcnnut
## results in microas

x = numpy.array(0, dtype=numpy.float64)
y = numpy.array(0, dtype=numpy.float64)
dx = numpy.array(0,dtype=numpy.float64)
dy = numpy.array(0,dtype=numpy.float64)
max_error = [0,0,0,0]

num_tests = 0
with open(sys.argv[1], 'r') as fin:
    for line in fin.readlines():
        if line[0] != '#':
            mjd, xref, yref, dxref, dyref = [ float(_) for _ in line.split() ]
            pyfcnnut.fcnnut(mjd, x, y, dx, dy)
            # print("{:+.1e} {:+.1e} {:+.1e} {:+.1e}".format(x-xref, y-yref, dx-dxref, dy-dyref))
            error = [x-xref, y-yref, dx-dxref, dy-dyref]
            for i,e in enumerate(error):
                if abs(e) > abs(max_error[i]): max_error[i] = e
            num_tests += 1

print("{:20s}| dX={:+.1e} dY={:+.1e} dsX={:+.1e} dsY={:+.1e} [\u03BCasec]".format("fcnnut", *max_error))
