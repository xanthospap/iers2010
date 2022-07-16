#! /usr/bin/python

import os
import sys
import re
import subprocess

##
## Script should be executed inside the 'sofa_unit_tests' folder
##
if os.path.dirname(os.path.realpath(__file__)) != 'sofa_unit_tests':
    dn = os.path.dirname(os.path.realpath(__file__))
    print('ERROR. Script should be run inside the \'sofa_unit_tests\' folder; current dir: {:}'.format(dn), file=sys.stderr)

for fn in os.listdir():
    if re.match(r'sofa-[a-z0-9A-Z]*\.out', fn):
        print("> Running Unit Test {:}".format(fn))
        print("{}".format(''.join(['-']*80)))
        subprocess.run(["./{:}".format(fn)])
        print("")