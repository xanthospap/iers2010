FC = gfortran 
FF = -c -Wall -O 

IERS_OBS = PMSDNUT2.o \
	UTLIBR.o \
	FUNDARG.o \
	FCNNUT.o \
	CAL2JD.o \
	DAT.o \
	NORM8.o \
	SPROD.o \
	ST1IDIU.o \
	ST1ISEM.o \
	ST1L1.o \
	STEP2DIU.o \
	STEP2LON.o \
	ZERO_VEC8.o \
	DEHANTTIDEINEL.o \
	CNMTX.o \
	ORTHO_EOP.o \
	RG_ZONT2.o \
	FCUL_A.o \
	FCUL_B.o \
	FCUL_ZD_HPA.o \
	GMF.o \
	VMF1.o \
	GPT.o \
	GPT2.o \
	ARG2.o

IERS_OBS_D0 = PMSDNUT2_D0.o \
	UTLIBR_D0.o \
	FUNDARG.o \
	FCNNUT.o \
	CAL2JD.o \
	DAT.o \
	NORM8.o \
	SPROD.o \
	ST1IDIU.o \
	ST1ISEM.o \
	ST1L1.o \
	STEP2DIU.o \
	STEP2LON.o \
	ZERO_VEC8.o \
	DEHANTTIDEINEL.o \
	CNMTX.o \
	ORTHO_EOP.o \
	RG_ZONT2.o \
	FCUL_A.o \
	FCUL_B.o \
	FCUL_ZD_HPA.o \
	GMF.o \
	VMF1.o \
	GPT.o \
	GPT2.o \
	ARG2.o

all: test_iers test_iers_d0

test_iers: MAIN.F $(IERS_OBS)
	$(FC) -o $@ MAIN.F $(IERS_OBS)
test_iers_d0: MAIN.F $(IERS_OBS_D0)
	$(FC) -o $@ MAIN.F $(IERS_OBS_D0)
testFundarg: test_fundarg.f $(IERS_OBS)
	$(FC) -o $@ test_fundarg.f $(IERS_OBS)
testPmsdnut2: test_pmsdnut2.f $(IERS_OBS)
	$(FC) -o $@ test_pmsdnut2.f $(IERS_OBS)
testPmsdnut2_d0: test_pmsdnut2.f $(IERS_OBS_D0)
	$(FC) -o $@ test_pmsdnut2.f $(IERS_OBS_D0)
testUtlibr: test_utlibr.f $(IERS_OBS)
	$(FC) -o $@ test_utlibr.f $(IERS_OBS)
testUtlibr_d0: test_utlibr.f $(IERS_OBS_D0)
	$(FC) -o $@ test_utlibr.f $(IERS_OBS_D0)
testFcnnut: test_fcnnut.f $(IERS_OBS)
	$(FC) -o $@ test_fcnnut.f $(IERS_OBS)
testDehanttideinel: test_dehanttideinel.f $(IERS_OBS)
	$(FC) -o $@ test_dehanttideinel.f $(IERS_OBS)
testArg2: test_arg2.f $(IERS_OBS)
	$(FC) -o $@ test_arg2.f $(IERS_OBS)
testCnmtx: test_cnmtx.f $(IERS_OBS)
	$(FC) -o $@ test_cnmtx.f $(IERS_OBS)
testOrthoeop: test_orthoeop.f $(IERS_OBS)
	$(FC) -o $@ test_orthoeop.f $(IERS_OBS)
testRgZont2: test_rgzont2.f $(IERS_OBS)
	$(FC) -o $@ test_rgzont2.f $(IERS_OBS)

FUNDARG.o : FUNDARG.F
	$(FC) $(FF) -o $@ FUNDARG.F
PMSDNUT2.o : PMSDNUT2.F
	$(FC) $(FF) -o $@ PMSDNUT2.F
PMSDNUT2_D0.o : PMSDNUT2_D0.F
	$(FC) $(FF) -o $@ PMSDNUT2_D0.F
UTLIBR.o : UTLIBR.F
	$(FC) $(FF) -o $@ UTLIBR.F
UTLIBR_D0.o : UTLIBR_D0.F
	$(FC) $(FF) -o $@ UTLIBR_D0.F
FCNNUT.o : FCNNUT.F FUNDARG.o
	$(FC) $(FF) -o $@ FCNNUT.F
CAL2JD.o : CAL2JD.F
	$(FC) $(FF) -o $@ CAL2JD.F
DAT.o : DAT.F 
	$(FC) $(FF) -o $@ DAT.F
NORM8.o : NORM8.F 
	$(FC) $(FF) -o $@ NORM8.F
SPROD.o : SPROD.F 
	$(FC) $(FF) -o $@ SPROD.F
ST1IDIU.o : ST1IDIU.F 
	$(FC) $(FF) -o $@ ST1IDIU.F
ST1ISEM.o : ST1ISEM.F 
	$(FC) $(FF) -o $@ ST1ISEM.F
ST1L1.o : ST1L1.F 
	$(FC) $(FF) -o $@ ST1L1.F
STEP2DIU.o : STEP2DIU.F 
	$(FC) $(FF) -o $@ STEP2DIU.F
STEP2LON.o : STEP2LON.F 
	$(FC) $(FF) -o $@ STEP2LON.F
ZERO_VEC8.o : ZERO_VEC8.F 
	$(FC) $(FF) -o $@ ZERO_VEC8.F
DEHANTTIDEINEL.o : DEHANTTIDEINEL.F CAL2JD.o DAT.o NORM8.o SPROD.o ST1IDIU.o ST1ISEM.o ST1L1.o STEP2DIU.o STEP2LON.o ZERO_VEC8.o
	$(FC) $(FF) -o $@ DEHANTTIDEINEL.F
CNMTX.o : CNMTX.F
	$(FC) $(FF) -o $@ CNMTX.F
ORTHO_EOP.o : ORTHO_EOP.F CNMTX.o
	$(FC) $(FF) -o $@ ORTHO_EOP.F
RG_ZONT2.o : RG_ZONT2.F FUNDARG.o
	$(FC) $(FF) -o $@ RG_ZONT2.F
FCUL_A.o : FCUL_A.F
	$(FC) $(FF) -o $@ FCUL_A.F
FCUL_B.o : FCUL_B.F
	$(FC) $(FF) -o $@ FCUL_B.F
FCUL_ZD_HPA.o : FCUL_ZD_HPA.F
	$(FC) $(FF) -o $@ FCUL_ZD_HPA.F
GMF.o : GMF.F
	$(FC) $(FF) -o $@ GMF.F
VMF1.o : VMF1.F
	$(FC) $(FF) -o $@ VMF1.F
GPT.o : GPT.F
	$(FC) $(FF) -o $@ GPT.F
GPT2.o : GPT2.F
	$(FC) $(FF) -o $@ GPT2.F