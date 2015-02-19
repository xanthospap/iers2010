# --------------------------------------------------------------
# MAKEFILE FOR LIBIERS10++ LIBRARY
# PART 2 / MAKEFILE FOR SRC
# --------------------------------------------------------------

##
##  WARNING
##  Variable builddir should be overwritten by the top-level
##+ Makefile;
##

##  Include the Makefile configuration file
## ===================================================
include ../Makefile.conf

## Source related folders (for convinience)
## ===================================================
dhtdir      = $(srcdir)/dehanttideinel
hspdir      = $(srcdir)/hardisp

## Compiler options
## ===================================================
CXX         = g++
DEFS        = -DHAVE_CONFIG_H $(STD_C11) $(OPT_LEVEL) $(USE_EXTERNAL_CONSTS) $(QUICK_EXIT)
CXXFLAGS    = -c -Wall $(DEFS) -g -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -Wdisabled-optimization
LDFLAGS     = -Wall $(DEFS)
LIBS        = 

LIB       = $(srcdir)/../$(package)v$(version).a

SRC       = fundarg.cpp  pmsdnut2.cpp  utlibr.cpp fcnnut.cpp arg2.cpp \
            dehanttideinel.cpp rg_zont2.cpp cnmtx.cpp ortho_eop.cpp apg.cpp \
            gpt.cpp gpt2.cpp st1idiu.cpp st1isem.cpp st1l1.cpp step2diu.cpp \
            step2lon.cpp cal2jd.cpp dat.cpp etutc.cpp eval.cpp recurs.cpp \
            shells.cpp utils.cpp spline.cpp tdfrph.cpp admint.cpp vmf1_ht.cpp \
			vmf1.cpp gmf.cpp

OBJ_     := $(SRC:%.cpp=%.o)
OBJ       = $(addprefix $(builddir)/, $(OBJ_))

HEAD      = $(srcdir)/iers2010.hpp

all: mkdir $(LIB)

mkdir:
	mkdir -p $(builddir)

install:
	$(INSTALL) -d $(libdir)
	$(INSTALL) -d $(includedir)
	$(INSTALL_PROGRAM) -m 0744 $(LIB) $(libdir)
	$(INSTALL_PROGRAM) -m 0744 $(HEAD) $(includedir)

clean:
	-rm $(OBJ) 
	-rm $(LIB)
	-rm -rf $(builddir)

$(LIB): $(OBJ)
	ar cru $(LIB) $(OBJ)
	ranlib $(LIB)
	cp $(LIB) ../$(package).a

$(builddir)/fundarg.o: $(srcdir)/fundarg.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/pmsdnut2.o: $(srcdir)/pmsdnut2.cpp ../obj/fundarg.o $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/utlibr.o: $(srcdir)/utlibr.cpp ../obj/fundarg.o $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/fcnnut.o: $(srcdir)/fcnnut.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/arg2.o: $(srcdir)/arg2.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/rg_zont2.o: $(srcdir)/rg_zont2.cpp $(HEAD) ../obj/fundarg.o
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/dehanttideinel.o: $(srcdir)/dehanttideinel.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/cnmtx.o: $(srcdir)/cnmtx.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/ortho_eop.o: $(srcdir)/ortho_eop.cpp $(srcdir)/cnmtx.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/gpt.o: $(srcdir)/gpt.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/gpt2.o: $(srcdir)/gpt2.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/apg.o: $(srcdir)/apg.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/vmf1_ht.o: $(srcdir)/vmf1_ht.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/vmf1.o: $(srcdir)/vmf1.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/gmf.o: $(srcdir)/gmf.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/st1idiu.o: $(dhtdir)/st1idiu.cpp $(dhtdir)/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/st1isem.o: $(dhtdir)/st1isem.cpp $(dhtdir)/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/st1l1.o: $(dhtdir)/st1l1.cpp $(dhtdir)/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/step2diu.o: $(dhtdir)/step2diu.cpp $(dhtdir)/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/step2lon.o: $(dhtdir)/step2lon.cpp $(dhtdir)/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/cal2jd.o: $(dhtdir)/cal2jd.cpp $(dhtdir)/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/dat.o: $(dhtdir)/dat.cpp $(dhtdir)/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/admint.o: $(hspdir)/admint.cpp $(hspdir)/hardisp.hpp $(hspdir)/tdfrph.cpp $(hspdir)/spline.cpp $(hspdir)/eval.cpp $(hspdir)/shells.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/etutc.o: $(hspdir)/etutc.cpp $(hspdir)/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/eval.o: $(hspdir)/eval.cpp $(hspdir)/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/recurs.o: $(hspdir)/recurs.cpp $(hspdir)/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/shells.o: $(hspdir)/shells.cpp $(hspdir)/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/utils.o: $(hspdir)/utils.cpp $(hspdir)/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/spline.o: $(hspdir)/spline.cpp $(hspdir)/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
$(builddir)/tdfrph.o: $(hspdir)/tdfrph.cpp $(hspdir)/hardisp.hpp $(hspdir)/etutc.cpp $(hspdir)/utils.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
