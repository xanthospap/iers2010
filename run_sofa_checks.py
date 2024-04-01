#! /usr/bin/python3

import os
import sys
import argparse
import subprocess
import glob

test_dir = 'test/sofa'

if __name__ == '__main__':
    if not os.path.isdir(test_dir):
        print("Failed to find testing directory \'test/sofa\'", file=sys.sdterr)
        sys.exit(1)

    print("Running test suite against SOFA library")

    listing = os.listdir(test_dir)
    for prog in [ os.path.join(test_dir, fn) for fn in listing if fn.endswith('.out') ]:
        subprocess.run([prog], stdout=sys.stderr, check=True)

    sys.exit(0)
