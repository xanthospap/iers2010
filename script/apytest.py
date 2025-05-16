#!/usr/bin/env python
# -*- coding: utf-8 -*-

##
##  This script is used to create test/validation data for the 
##  test/unit_test/i2c_astropy family of test programs.
##
##  It uses astropy to transform between celestial and terrestrial RFs. It 
##  parses an Sp3 file (which except from position records should also have 
##  velocities) and transforms the states from ITRF to GCRF, assuming time 
##  in TAI.
##

from collections.abc import Callable

import numpy as np
import pandas as pd

from astropy.time import Time, TimeDelta
from astropy.coordinates import ITRS, GCRS, CartesianRepresentation, CartesianDifferential
from astropy import units as u

import sys
import datetime

DELTA = TimeDelta(19., format="sec")  # GPS time is 19s behind TAI

def pt2mjd(t): return Time(t).mjd

def xtractSatelliteState(fn, sat):
    entries = []
    with open(fn, 'r') as fin:
        for line in fin.readlines():
            #try:
            if line.startswith('* '):
                tstr = ' '.join( [str(int(x)) for x in line.split()[1:6]]+["{:.6f}".format(float(line.split()[6]))] )
                t = pt2mjd(datetime.datetime.strptime(tstr, '%Y %m %d %H %M %S.%f'))
            elif line.startswith('P'+sat.upper()):
                x, y, z = [float(x)*1e+3 for x in line.split()[1:4]]
            elif line.startswith('V'+sat.upper()):
                vx, vy, vz = [float(x)*1e-1 for x in line.split()[1:4]]
                entries.append([t, np.array([x, y, z]), np.array([vx,vy,vz])])
            else:
                pass
            #except:
            #    pass
    return entries

def eci_to_ecef(
    epoch: Time, pos: CartesianRepresentation, vel: CartesianDifferential
) -> tuple[np.ndarray, np.ndarray]:
    """Convert ECI (GCRS) to ECEF (ITRS) using astropy and IERS data."""
    gcrs = GCRS(pos.with_differentials(vel), obstime=epoch)
    itrs = gcrs.transform_to(ITRS(obstime=epoch))
    pos = itrs.cartesian.xyz.to_value()
    vel = itrs.cartesian.differentials['s'].d_xyz.to_value()
    return pos, vel

def ecef_to_eci(
    epoch: Time, pos: CartesianRepresentation, vel: CartesianDifferential
) -> tuple[np.ndarray, np.ndarray]:
    """Convert ECEF (ITRS) to ECI (GCRS) using astropy and IERS data."""
    itrs = ITRS(pos.with_differentials(vel), obstime=epoch)
    gcrs = itrs.transform_to(GCRS(obstime=epoch))
    pos = gcrs.cartesian.xyz.to_value()
    vel = gcrs.cartesian.differentials['s'].d_xyz.to_value()
    return pos, vel

if __name__ == "__main__":
    data = xtractSatelliteState(sys.argv[1], sys.argv[2])
    epoch = Time([d[0] for d in data], format="mjd", scale="tai")
    vel = CartesianDifferential([d[2][0] for d in data], [d[2][1] for d in data], [d[2][2] for d in data], unit=u.m / u.s)
    pos = CartesianRepresentation([d[1][0] for d in data], [d[1][1] for d in data], [d[1][2] for d in data], unit=u.m)

    itrs = ITRS(pos.with_differentials(vel), obstime=epoch)
    gcrs = itrs.transform_to(GCRS(obstime=epoch))
    pos = gcrs.cartesian.xyz.to_value()
    vel = gcrs.cartesian.differentials['s'].d_xyz.to_value()

    ip = itrs.cartesian.xyz.to_value()
    iv = itrs.cartesian.differentials['s'].d_xyz.to_value()
    for i, t in enumerate(epoch):
        print(f"{t.mjd:.15f} {ip[0][i]:.9f} {ip[1][i]:.9f} {ip[2][i]:.9f} {iv[0][i]*1e3:.12f} {iv[1][i]*1e3:.12f} {iv[2][i]*1e3:.12f} {pos[0][i]:.6f} {pos[1][i]:.6f} {pos[2][i]:.6f} {vel[0][i]*1e3:.9f} {vel[1][i]*1e3:.9f} {vel[2][i]*1e3:.9f}")
