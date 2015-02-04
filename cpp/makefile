# --------------------------------------------------------------
# MAKEFILE FOR IERS2010 LIBRARY
# --------------------------------------------------------------

# --------------------------------------------------------------
#
CFLAGS   = -c $(C11) -Wall $(OPT) $(USE_EXTERNAL_CONSTS) $(QUICK_EXIT)
LFLAGS   = $(C11) -Wall $(OPT) $(USE_EXTERNAL_CONSTS) $(QUICK_EXIT)
#
LIB      = ../lib/libiers2010.a
#
SRC      = fundarg.cpp  pmsdnut2.cpp  utlibr.cpp fcnnut.cpp arg2.cpp dehanttideinel.cpp \
           rg_zont2.cpp cnmtx.cpp ortho_eop.cpp apg.cpp gpt.cpp \
           st1idiu.cpp st1isem.cpp st1l1.cpp step2diu.cpp step2lon.cpp cal2jd.cpp dat.cpp \
           etutc.cpp eval.cpp recurs.cpp shells.cpp utils.cpp spline.cpp tdfrph.cpp admint.cpp

#
OBJ_    := $(SRC:%.cpp=%.o)
OBJ      = $(addprefix ../obj/, $(OBJ_))
#
HEAD     = iers2010.hpp
#
# --------------------------------------------------------------

all: print mkdirs $(LIB) link donem

mkdirs:
	-@mkdir ../lib 2>/dev/null
	-@mkdir ../obj 2>/dev/null
	-@mkdir ../inc 2>/dev/null

clean:
	-rm $(OBJ) 2>/dev/null

delete:
	-rm $(OBJ) 2>/dev/null
	-rm $(LIB) 2>/dev/null
	-rm $(CURDIR)/../inc/iers2010.hpp 2>/dev/null
	-rm $(CURDIR)/../inc/dehanttideinel.hpp 2>/dev/null

link:
	@echo "... Linking header file(s) ..."
	@/bin/ln -nsf $(CURDIR)/iers2010.hpp $(CURDIR)/../inc/iers2010.hpp
	@/bin/ln -nsf $(CURDIR)/dehanttideinel/dehanttideinel.hpp $(CURDIR)/../inc/dehanttideinel.hpp

print:
	@echo ""
	@echo "********************************************************************************"
	@echo " Compiling IERS-2010 library."
	@echo "--------------------------------------------------------------------------------"

donem:
	@echo "--------------------------------------------------------------------------------"
	@echo " done!"
	@echo " Created library: $(LIB)"
	@echo "********************************************************************************"
	@echo ""

$(LIB): $(OBJ)
	ar rc $(LIB) $(OBJ)

../obj/fundarg.o: fundarg.cpp $(HEAD)
	$(CXX) $(CFLAGS) -o ../obj/fundarg.o fundarg.cpp
../obj/pmsdnut2.o: pmsdnut2.cpp ../obj/fundarg.o $(HEAD)
	$(CXX) $(CFLAGS) -o ../obj/pmsdnut2.o pmsdnut2.cpp
../obj/utlibr.o: utlibr.cpp ../obj/fundarg.o $(HEAD)
	$(CXX) $(CFLAGS) -o ../obj/utlibr.o utlibr.cpp
../obj/fcnnut.o: fcnnut.cpp $(HEAD)
	$(CXX) $(CFLAGS) -o ../obj/fcnnut.o fcnnut.cpp
../obj/arg2.o: arg2.cpp $(HEAD)
	$(CXX) $(CFLAGS) -o ../obj/arg2.o arg2.cpp
../obj/rg_zont2.o: rg_zont2.cpp $(HEAD) ../obj/fundarg.o
	$(CXX) $(CFLAGS) -o ../obj/rg_zont2.o rg_zont2.cpp
../obj/dehanttideinel.o: dehanttideinel.cpp $(HEAD)
	$(CXX) $(CFLAGS) -o ../obj/dehanttideinel.o dehanttideinel.cpp
../obj/cnmtx.o: cnmtx.cpp $(HEAD)
	$(CXX) $(CFLAGS) -o ../obj/cnmtx.o cnmtx.cpp
../obj/ortho_eop.o: ortho_eop.cpp cnmtx.cpp $(HEAD)
	$(CXX) $(CFLAGS) -o ../obj/ortho_eop.o ortho_eop.cpp
../obj/gpt.o: gpt.cpp $(HEAD)
	$(CXX) $(CFLAGS) -o ../obj/gpt.o gpt.cpp
../obj/apg.o: apg.cpp $(HEAD)
	$(CXX) $(CFLAGS) -o ../obj/apg.o apg.cpp
../obj/st1idiu.o: dehanttideinel/st1idiu.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CFLAGS) -o ../obj/st1idiu.o dehanttideinel/st1idiu.cpp
../obj/st1isem.o: dehanttideinel/st1isem.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CFLAGS) -o ../obj/st1isem.o dehanttideinel/st1isem.cpp
../obj/st1l1.o: dehanttideinel/st1l1.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CFLAGS) -o ../obj/st1l1.o dehanttideinel/st1l1.cpp
../obj/step2diu.o: dehanttideinel/step2diu.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CFLAGS) -o ../obj/step2diu.o dehanttideinel/step2diu.cpp
../obj/step2lon.o: dehanttideinel/step2lon.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CFLAGS) -o ../obj/step2lon.o dehanttideinel/step2lon.cpp
../obj/cal2jd.o: dehanttideinel/cal2jd.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CFLAGS) -o ../obj/cal2jd.o dehanttideinel/cal2jd.cpp
../obj/dat.o: dehanttideinel/dat.cpp dehanttideinel/dehanttideinel.hpp
	$(CXX) $(CFLAGS) -o ../obj/dat.o dehanttideinel/dat.cpp
../obj/admint.o: hardisp/admint.cpp hardisp/hardisp.hpp hardisp/tdfrph.cpp hardisp/spline.cpp hardisp/eval.cpp hardisp/shells.cpp 
	$(CXX) $(CFLAGS) -o ../obj/admint.o hardisp/admint.cpp
../obj/etutc.o: hardisp/etutc.cpp hardisp/hardisp.hpp
	$(CXX) $(CFLAGS) -o ../obj/etutc.o hardisp/etutc.cpp
../obj/eval.o: hardisp/eval.cpp hardisp/hardisp.hpp
	$(CXX) $(CFLAGS) -o ../obj/eval.o hardisp/eval.cpp
../obj/recurs.o: hardisp/recurs.cpp hardisp/hardisp.hpp
	$(CXX) $(CFLAGS) -o ../obj/recurs.o hardisp/recurs.cpp
../obj/shells.o: hardisp/shells.cpp hardisp/hardisp.hpp
	$(CXX) $(CFLAGS) -o ../obj/shells.o hardisp/shells.cpp
../obj/utils.o: hardisp/utils.cpp hardisp/hardisp.hpp
	$(CXX) $(CFLAGS) -o ../obj/utils.o hardisp/utils.cpp
../obj/spline.o: hardisp/spline.cpp hardisp/hardisp.hpp
	$(CXX) $(CFLAGS) -o ../obj/spline.o hardisp/spline.cpp
../obj/tdfrph.o: hardisp/tdfrph.cpp hardisp/hardisp.hpp hardisp/etutc.cpp hardisp/utils.cpp
	$(CXX) $(CFLAGS) -o ../obj/tdfrph.o hardisp/tdfrph.cpp
