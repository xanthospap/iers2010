noinst_PROGRAMS = testIers2010 \
testBlqIn \
testHardisp

check_PROGRAMS = testFundarg \
testPmsdnut2 \
testUtlibr \
testFcnnut \
testDehanttideinel \
testArg2 \
testCnmtx \
testOrthoeop \
testRgZont2 \
testFcula \
testFculb \
testFculZdhPa

MCXXFLAGS = \
-std=c++17 \
-g \
-pg \
-Wall \
-Wextra \
-Werror \
-pedantic \
-W \
-Wshadow \
-Wdisabled-optimization \
-DDEBUG

testIers2010_SOURCES  = test_iers2010.cpp
testIers2010_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testIers2010_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testBlqIn_SOURCES  = test_blq.cpp
testBlqIn_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testBlqIn_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testHardisp_SOURCES  = test_hardisp.cpp
testHardisp_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testHardisp_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testFundarg_SOURCES  = test_fundarg.cpp
testFundarg_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testFundarg_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testPmsdnut2_SOURCES  = test_pmsdnut2.cpp
testPmsdnut2_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testPmsdnut2_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testUtlibr_SOURCES  = test_utlibr.cpp
testUtlibr_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testUtlibr_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testFcnnut_SOURCES  = test_fcnnut.cpp
testFcnnut_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testFcnnut_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testDehanttideinel_SOURCES  = test_dehanttideinel.cpp
testDehanttideinel_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testDehanttideinel_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testArg2_SOURCES  = test_arg2.cpp
testArg2_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testArg2_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testCnmtx_SOURCES  = test_cnmtx.cpp
testCnmtx_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testCnmtx_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testOrthoeop_SOURCES  = test_orthoeop.cpp
testOrthoeop_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testOrthoeop_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testRgZont2_SOURCES  = test_rgzont2.cpp
testRgZont2_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testRgZont2_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testFcula_SOURCES  = test_fcula.cpp
testFcula_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testFcula_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testFculb_SOURCES  = test_fculb.cpp
testFculb_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testFculb_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a

testFculZdhPa_SOURCES  = test_fculzdhpa.cpp
testFculZdhPa_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src
testFculZdhPa_LDADD    = $(top_srcdir)/src/libiers2010.la /usr/local/lib/libggdatetime.a
#check_PROGRAMS    = test
#dist_test_SOURCES = test.cpp
#EXTRA_DIST        = hardisp_results_onsa \
#					hardisp_results_reyk \
#					reyk.blq \
#					onsa.blq
#test_CXXFLAGS     = -I$(top_srcdir)/src -L$(top_srcdir)/src -std=c++11
#test_LDADD        = $(top_srcdir)/src/libiers10++.a
#
#check_SCRIPTS     = greptest.sh
#TESTS             = $(check_SCRIPTS)
#
#EXTRA_PROGRAMS            = test-hardisp
#dist_test_hardisp_SOURCES = test-hardisp.cpp
#test_hardisp_CXXFLAGS     = -I$(top_srcdir)/src -L$(top_srcdir)/src -std=c++11
#test_hardisp_LDADD        = $(top_srcdir)/src/libiers10++.a
#
#greptest.sh:
#	echo './test | grep "Library Status : OK"' > greptest.sh
#	chmod +x greptest.sh
#
#CLEANFILES        = greptest.sh
