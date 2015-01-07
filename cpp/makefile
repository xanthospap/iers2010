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
SRC      = fundarg.cpp  pmsdnut2.cpp  utlibr.cpp fcnnut.cpp arg2.cpp
#
OBJ_    := $(SRC:%.cpp=%.o)
OBJ      = $(addprefix ../obj/, $(OBJ_))
#
HEAD     = iers2010.hpp
#
# --------------------------------------------------------------

all: print mkdirs $(LIB) link clean donem

mkdirs:
	-@mkdir ../lib
	-@mkdir ../obj
	-@mkdir ../inc

clean:
	-rm $(OBJ)

delete:
	-rm $(OBJ)
	-rm $(LIB)
	-rm $(CURDIR)/../inc/iers2010.hpp

link:
	@echo "... Linking header file(s) ..."
	@/bin/ln -nsf $(CURDIR)/iers2010.hpp $(CURDIR)/../inc/iers2010.hpp

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
