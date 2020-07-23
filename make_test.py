#! /usr/bin/python

from shutil import copyfile
import os
import random
import fileinput
import datetime
import math
import jdcal

def replace_inplace(filename, text2cut, text2paste):
  with fileinput.FileInput(filename, inplace=True, backup='.bak') as file:
    for line in file:
     print(line.replace(text2cut, text2paste), end='')


cpp_test = "test/test_iers2010.cpp"
cpp_test_am = "test/test_iers2010.cpp.am"
if not os.path.isfile(cpp_test_am):
  print("[ERROR] Failed to find cpp test file \"{:}\"".format(cpp_test_am))
else:
  copyfile(cpp_test_am, cpp_test)

# for_test_am = "fortran_impl/MAIN.F.am"
# for_test = "fortran_impl/MAIN.F"
# if not os.path.isfile(for_test_am):
#   print("[ERROR] Failed to find fortran test file \"{:}\"".format(for_test_am))
# else:
#   copyfile(for_test_am, for_test)

## generate random parameters
julian_day = random.uniform(2442413.5e0, 2460676.5e0)
julian_century = (julian_day - 2451545e0)/36525e0
mjd = julian_day - 2400000.5e0
gdt = jdcal.jd2gcal(2400000.5, mjd)
tf = "{:04d}-{:02d}-{:02d}".format(gdt[0],gdt[1],gdt[2])
pdt = datetime.datetime.strptime(tf, "%Y-%m-%d")
year = pdt.year
fractional_doy = float(pdt.strftime("%j")) + gdt[3]
dlat = random.uniform(-90e0, 90e0)
dlon = random.uniform(-180e0, 180e0)
dhgt = random.uniform(-1e3, 9e3)
tmpr = random.uniform(-30e0, 45e0) + 273.15e0
elev = random.uniform(0e0, 90e0)
pres = random.uniform(0e0, 2e3)
wvpr = random.uniform(0e0, 30e0)
lwvl = random.uniform(4.5e-1, 6e-1)
hfac = random.uniform(0e0, 1e0)
wfac = random.uniform(0e0, 1e0)
grbl = random.randint(0, 1)

## replace fundarg input parameter (julian century)
replace_inplace(cpp_test, "fundarg_inp", "{:20.15f}".format(julian_century))

## replace pmsdnut2 input parameter (modified julian day)
replace_inplace(cpp_test, "pmsdnut2_inp", "{:20.15f}".format(mjd))

## replace utlibr input parameter (modified julian day)
replace_inplace(cpp_test, "utlibr_inp", "{:20.15f}".format(mjd))

## replace utlibr input parameter (modified julian day)
replace_inplace(cpp_test, "fcnnut_inp", "{:20.15f}".format(mjd))

## replace arg2 input parameter (year, fractional doy)
in_string = "{:4d}, {:20.15f}".format(year, fractional_doy)
replace_inplace(cpp_test, "arg2_inp", in_string)

## replace cnmtx input parameter (modified julian day)
replace_inplace(cpp_test, "cnmtx_inp", "{:20.15f}".format(mjd))

## replace ortho_eop input parameter (modified julian day)
replace_inplace(cpp_test, "ortho_eop_inp", "{:20.15f}".format(mjd))

## replace rg_zont2 input parameter (julian century)
replace_inplace(cpp_test, "rg_zont2_inp", "{:20.15f}".format(julian_century))

## replace fcul_a input parameter (latitude, height, temperature, elevation)
in_string = "{:20.15f}, {:10.4f}e0, {:10.3f}e0, {:5.2f}e0".format(dlat, dhgt, tmpr, elev)
replace_inplace(cpp_test, "fcul_a_inp", in_string)

## replace fcul_b input parameter (latitude, height, doy, elevation)
in_string = "{:20.15f}, {:10.4f}e0, {:20.15f}e0, {:5.2f}e0".format(dlat, dhgt, fractional_doy, elev)
replace_inplace(cpp_test, "fcul_b_inp", in_string)

## replace fcul_zd_hpa input parameter (latitude, height, pressure, water vapor, 
##                                     laser wavelength)
in_string = "{:20.15f}, {:10.4f}e0, {:10.5f}e0, {:10.5f}e0, {:10.5f}e0".format(dlat, dhgt, pres, wvpr, lwvl)
replace_inplace(cpp_test, "fcul_zd_hpa_inp", in_string)

## replace gmf input parameter (mjd, latitude, longtitude, height, zenith distance)
in_string = "{:20.15f}, {:20.15f}, {:20.15f}, {:10.3f}e0, {:15.10f}e0".format(mjd, math.radians(dlat),
  math.radians(dlon), dhgt, math.radians(90e0-elev))
replace_inplace(cpp_test, "gmf_inp", in_string)

## replace vmf1 input parameter (hydro coef., wet coef., mjd, latitude, zenith distance)
in_string = "{:20.15f}, {:20.15f}, {:20.15f}, {:20.15f}, {:15.10f}e0".format(
  hfac, wfac, mjd, math.radians(dlat), math.radians(90e0-elev))
replace_inplace(cpp_test, "vmf1_inp", in_string)

## replace vmf1_ht input parameter (hydro coef., wet coef., mjd, latitude, height, zenith distance)
in_string = "{:20.15f}, {:20.15f}, {:20.15f}, {:20.15f}, {:15.5f}e0, {:15.10f}e0".format(
  hfac, wfac, mjd, math.radians(dlat), dhgt, math.radians(90e0-elev))
replace_inplace(cpp_test, "vmf1_ht_inp", in_string)

## replace gpt input parameter (mjd, latitude, longtitude, height)
in_string = "{:20.15f}, {:20.15f}, {:20.15f}, {:10.3f}e0".format(mjd, math.radians(dlat),
  math.radians(dlon), dhgt)
replace_inplace(cpp_test, "gpt_inp", in_string)

## replace gpt2 input parameter (mjd, latitude, longtitude, height, nstat, True/False)
in_string = "{:20.15f}, {:20.15f}, {:20.15f}, {:10.3f}e0, {:1d}".format(mjd, math.radians(dlat),
  math.radians(dlon), dhgt, grbl)
replace_inplace(cpp_test, "gpt2_inp", in_string)
