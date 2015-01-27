#! /usr/bin/env python

"""
  This script takes as cmd two files, containing results from running
  test-cases for the libiers10++ library, and compares output for
  each function.
  The input files are format-specific. Run the test.e programs in
  test/test.e and test/fortran/test.e to produce them.
"""

import os, sys
import re

##  print usage and exit
def usage () :
	print 'Invalid cmds; Runs as cmpcvf.py <result_file_1> <result_file_2>'
	sys.exit (0)

##  given an input file and a string, search the input file to find the string.
##+ The actual string searched for, is 'Function <string>' at the start of line
##+ using regular expressions. If the string is found, the input buffre is set 
##+ to the line the string was found and 0 is returned; else 1 is returned.
def gotostr (ifile,strng) :
	ifile.seek (0)
	string = '^ *Function '+strng+' *$'
	line = ifile.readline ()
	while True :
		result = re.match (string,line)
		if result:
			return 0
		line = ifile.readline ()
		if not line:
			return 1


## check number of command line arguments
if len (sys.argv) != 3 :
	print 'ERROR!'
	usage ()

## try opening input files
try :
	fin1 = open (sys.argv[1], 'r')
	try :
		fin2 = open (sys.argv[2], 'r')
	except:
		print 'ERROR. Could not open file',sys.argv[2]
		sys.exit (1)
except:
	print 'ERROR. Could not open file',sys.argv[1]
	sys.exit (1)

## get a list of all functions in the first file
function_list = []
for line in fin1 :
	gr = re.search (r"^ *Function ",line)
	if gr :
		function_list.append ( line.split()[1] )
fin1.seek (0)

for f in function_list:
	i = gotostr (fin1,f)


## close the input files
fin1.close ()
fin2.close ()
