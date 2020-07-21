#! /usr/bin/python

from shutil import copyfile
import os
import random
import fileinput
import datetime
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

for_test_am = "fortran_impl/MAIN.F.am"
for_test = "fortran_impl/MAIN.F"
if not os.path.isfile(for_test_am):
  print("[ERROR] Failed to find fortran test file \"{:}\"".format(for_test_am))
else:
  copyfile(for_test_am, for_test)

julian_day = random.uniform(2442413.5e0, 2460676.5e0)
julian_century = (julian_day - 2451545e0)/36525e0
mjd = julian_day - 2400000.5e0
gdt = jdcal.jd2gcal(2400000.5, mjd)
tf = "{:04d}-{:02d}-{:02d}".format(gdt[0],gdt[1],gdt[2])
pdt = datetime.datetime.strptime(tf, "%Y-%m-%d")
year = pdt.year
fractional_doy = float(pdt.strftime("%j")) + gdt[3]

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
