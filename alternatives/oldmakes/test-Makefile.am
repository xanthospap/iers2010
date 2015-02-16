check_PROGRAMS    = test
dist_test_SOURCES = test.cpp
test_CXXFLAGS     = -I$(top_srcdir)/src -std=c++11
test_LDADD        = ../libiers10++.a

check_SCRIPTS     = greptest.sh
TESTS             = $(check_SCRIPTS)

greptest.sh:
	echo './test | grep "Library Status : OK"' > greptest.sh
	chmod +x greptest.sh

CLEANFILES        = greptest.sh
