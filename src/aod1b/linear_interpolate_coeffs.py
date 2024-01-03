#! /usr/bin/python

##
##  Dummy python script to perform linear interpoation of AOD1B coeffs.
##  This script is used to test the libiers2010
##

import sys, os
import datetime
import re

def extract(n,m,fn,ctype):
    lst = []
    with open(fn, 'r') as fin:
        line = fin.readline()
        while line:
            if line.strip() == "END OF HEADER":
                break
            line = fin.readline()
        line = fin.readline()
        while line:
# DATA SET 01:  16471 COEFFICIENTS FOR 2023-12-25 00:00:00 OF TYPE atm
            f = re.match("DATA SET ([0-9]*): *([0-9]*) *COEFFICIENTS FOR ([0-9-]*) ([0-9:]*) OF TYPE ([a-z]*)", line.strip())
            set_nr = int(f[1])
            num_lines = int(f[2])
            t = datetime.datetime.strptime(' '.join([f[3],f[4]]), "%Y-%m-%d %H:%M:%S")
            tp = f[5]
            Clm = Slm = None
            if tp.lower() == ctype:
                for i in range(num_lines):
                    line = fin.readline()
                    l = line.split()
                    if [int(x) for x in l[0:2]] == [n,m]:
                        Clm = float(l[2])
                        Slm = float(l[3])
                lst.append((t,Clm,Slm))
            else:
                for i in range(num_lines):
                    line = fin.readline()
            line = fin.readline()
    return lst

def fdays(t1,t2):
   difference = t2 - t1
   return difference.total_seconds() / datetime.timedelta(days=1).total_seconds()

def interpolate(lst, t):
    for i,t0 in enumerate(lst[:-1]):
        if t >= t0[0] and t < lst[i+1][0]:
            t1 = lst[i][0]
            t2 = lst[i+1][0]
            c1 = lst[i][1]
            c2 = lst[i+1][1]
            s1 = lst[i][2]
            s2 = lst[i+1][2]
            t2i = fdays(t,t2)
            t21 = fdays(t1,t2)
            ti1 = fdays(t1,t)
            return c1 * t2i/t21 + c2 * ti1/t21, s1 * t2i/t21 + s2 * ti1/t21
    raise RuntimeError

if __name__ == "__main__":
    data = extract(10,9,sys.argv[1],'atm')
    t = data[0][0]
    while t < data[-1][0] + datetime.timedelta(hours=3):
        try:
            c,s = interpolate(data, t)
            print('{:}.000000000 {:+.15e} {:+.15e}'.format(t.strftime("%Y/%m/%d %H:%M:%S"), c,s))
        except: pass
        t += datetime.timedelta(seconds=180)
