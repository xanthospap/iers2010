# --------------------------------------------------------------
# MAKEFILE FOR LIBIERS10++ LIBRARY
# PART 2 / MAKEFILE FOR SRC
# @configure_input@
# --------------------------------------------------------------

C11       = -std=c++11

CXXFLAGS  = -c $(C11) -Wall $(OPT) $(USE_EXTERNAL_CONSTS) $(QUICK_EXIT)
LDFLAGS   = $(C11) -Wall $(OPT) $(USE_EXTERNAL_CONSTS) $(QUICK_EXIT)

LIB       = ../$(package)v$(version).a

SRC       = fundarg.cpp  pmsdnut2.cpp  utlibr.cpp fcnnut.cpp arg2.cpp \
            dehanttideinel.cpp rg_zont2.cpp cnmtx.cpp ortho_eop.cpp apg.cpp \
            gpt.cpp gpt2.cpp st1idiu.cpp st1isem.cpp st1l1.cpp step2diu.cpp \
            step2lon.cpp cal2jd.cpp dat.cpp etutc.cpp eval.cpp recurs.cpp \
            shells.cpp utils.cpp spline.cpp tdfrph.cpp admint.cpp

OBJ_     := $(SRC:%.cpp=%.o)
OBJ       = $(addprefix ../obj/, $(OBJ_))

HEAD      = iers2010.hpp

all: mkdir $(LIB)

mkdir:
	mkdir -p ../obj

install:
	install -d $(libdir)
	install -d $(includedir)
	install -m 0744 $(LIB) $(libdir)
	install -m 0744 $(HEAD) $(includedir)

clean:
	-rm $(OBJ) 
	-rm $(LIB)
	-rm -rf ../obj

$(LIB): $(OBJ)
	ar rc $(LIB) $(OBJ)

../obj/fundarg.o: fundarg.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/pmsdnut2.o: pmsdnut2.cpp ../obj/fundarg.o $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/utlibr.o: utlibr.cpp ../obj/fundarg.o $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/fcnnut.o: fcnnut.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/arg2.o: arg2.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/rg_zont2.o: rg_zont2.cpp $(HEAD) ../obj/fundarg.o
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/dehanttideinel.o: dehanttideinel.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/cnmtx.o: cnmtx.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/ortho_eop.o: ortho_eop.cpp cnmtx.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/gpt.o: gpt.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/gpt2.o: gpt2.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/apg.o: apg.cpp $(HEAD)
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/st1idiu.o: dehanttideinel/st1idiu.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/st1isem.o: dehanttideinel/st1isem.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/st1l1.o: dehanttideinel/st1l1.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/step2diu.o: dehanttideinel/step2diu.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/step2lon.o: dehanttideinel/step2lon.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/cal2jd.o: dehanttideinel/cal2jd.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/dat.o: dehanttideinel/dat.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/admint.o: hardisp/admint.cpp hardisp/hardisp.hpp hardisp/tdfrph.cpp hardisp/spline.cpp hardisp/eval.cpp hardisp/shells.cpp 
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/etutc.o: hardisp/etutc.cpp hardisp/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/eval.o: hardisp/eval.cpp hardisp/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/recurs.o: hardisp/recurs.cpp hardisp/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/shells.o: hardisp/shells.cpp hardisp/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/utils.o: hardisp/utils.cpp hardisp/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/spline.o: hardisp/spline.cpp hardisp/hardisp.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<
../obj/tdfrph.o: hardisp/tdfrph.cpp hardisp/hardisp.hpp hardisp/etutc.cpp hardisp/utils.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
