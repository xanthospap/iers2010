check_PROGRAMS    = test
dist_test_SOURCES = test.cpp
EXTRA_DIST        = hardisp_results_onsa \
					hardisp_results_reyk \
					reyk.blq \
					onsa.blq
test_CXXFLAGS     = -I$(top_srcdir)/src -L$(top_srcdir)/src -std=c++11
test_LDADD        = $(top_srcdir)/src/libiers10++.a

check_SCRIPTS     = greptest.sh
TESTS             = $(check_SCRIPTS)

EXTRA_PROGRAMS            = test-hardisp
dist_test_hardisp_SOURCES = test-hardisp.cpp
test_hardisp_CXXFLAGS     = -I$(top_srcdir)/src -L$(top_srcdir)/src -std=c++11
test_hardisp_LDADD        = $(top_srcdir)/src/libiers10++.a

greptest.sh:
	echo './test | grep "Library Status : OK"' > greptest.sh
	chmod +x greptest.sh

CLEANFILES        = greptest.sh
