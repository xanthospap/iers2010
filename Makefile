# --------------------------------------------------------------
# MAKEFILE FOR PROJECT IERS2010
# --------------------------------------------------------------

export package := libiers10++
export version := 1.0

include Makefile.conf

all:
	$(MAKE) -C src $@

test:
	(cd test; g++ -Wall -std=c++11 -L../lib/ -I../inc/ test.cpp -liers2010)
