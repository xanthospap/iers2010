#! /usr/bin/python

import os
import sys
import re
import subprocess

##
## Script should be executed inside the 'sofa_unit_tests' folder
##

for fn in os.listdir():
    if re.match(r'sofa-[a-z0-9A-Z]*\.out', fn):
        print("> Running Unit Test {:}".format(fn))
        print("{}".format(''.join(['-']*80)))
        subprocess.run(["./{:}".format(fn)])
        print("")