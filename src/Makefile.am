bin_PROGRAMS               = hardisp
hardisp_SOURCES            = hardisp.cpp
hardisp_CXXFLAGS           = -std=c++11 -g -Wall -Wextra -pedantic -W -Wshadow -Winline -Wdisabled-optimization -I$(srcdir) -L$(srcdir)
hardisp_LDADD              = libiers10++.a
pkglib_LIBRARIES           = libiers10++.a
dist_pkginclude_HEADERS    = iers2010.hpp hardisp.hpp dehanttideinel.hpp
libiers10___a_CXXFLAGS     = -std=c++11 -g -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -Wdisabled-optimization
dist_pkgdata_DATA          = ../gpt2_5.grd
dist_libiers10___a_SOURCES = fundarg.cpp \
							 pmsdnut2.cpp \
							 utlibr.cpp \
							 fcnnut.cpp \
							 arg2.cpp \
							 vmf1_ht.cpp \
							 vmf1.cpp \
							 gmf.cpp \
							 dehanttideinel.cpp \
							 rg_zont2.cpp \
							 cnmtx.cpp \
							 ortho_eop.cpp \
							 apg.cpp \
							 gpt.cpp \
							 gpt2.cpp \
							 fcul_a.cpp \
							 fcul_b.cpp \
							 fcul_zd_hpa.cpp \
							 st1idiu.cpp \
							 st1isem.cpp \
							 st1l1.cpp \
							 step2diu.cpp \
							 step2lon.cpp \
							 cal2jd.cpp \
							 dat.cpp \
							 etutc.cpp \
							 eval.cpp \
							 recurs.cpp \
							 shells.cpp \
							 utils.cpp \
							 spline.cpp \
							 tdfrph.cpp \
							 admint.cpp
dist_libiers10___a_EXTRA   = gpt2-unset.cpp
EXTRA_DIST                 = gpt2-unset.cpp
nodist_libiers10___a_SOURCES = gpt2.cpp
gpt2.cpp: gpt2-unset.cpp
	cp gpt2-unset.cpp gpt2.cpp
	$(SED) -i "s|/foo/bar/gpt2_5.grd|$(abs_top_builddir)/gpt2_5.grd|g" $@

## TODO Remove the gpt2-unset.cpp file
