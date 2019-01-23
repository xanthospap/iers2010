lib_LTLIBRARIES = libiers2010.la
bin_PROGRAMS = hardisp

## Only use the following flags for debuging purposes
libiers2010_la_CXXFLAGS = \
	-std=c++17 \
	-g \
	-pg \
	-Wall \
	-Wextra \
	-Werror \
	-pedantic \
	-W \
	-Wshadow \
	-Winline \
	-Wdisabled-optimization \
	-DDEBUG

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

dist_include_HEADERS = \
	iers2010.hpp

dist_libiers2010_la_SOURCES = \
	fundarg.cpp \
	pmsdnut2.cpp \
	utlibr.cpp \
	fcnnut.cpp \
	arg2.cpp \
	cal2jd.cpp \
	dat.cpp \
	st1idiu.cpp \
	st1isem.cpp \
	st1l1.cpp \
	step2diu.cpp \
	step2lon.cpp \
	dehanttideinel.cpp \
	etutc.cpp \
	tdfrph.cpp \
	spline.cpp \
	eval.cpp \
	recurs.cpp \
	shells.cpp \
	admint.cpp \
	cnmtx.cpp \
	ortho_eop.cpp \
	rg_zont2.cpp \
	fcul_a.cpp \
	fcul_b.cpp \
	fcul_zd_hpa.cpp \
	gmf.cpp \
	vmf1.cpp \
	vmf1_ht.cpp \
	gpt.cpp \
	gpt2.cpp

hardisp_SOURCES  = hardisp.cpp
hardisp_CXXFLAGS = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
hardisp_LDADD    = $(top_srcdir)/src/libiers2010.la
