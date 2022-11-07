FC = gfortran 
FF = -c -Wall -O 

IERS_OBS = PMSDNUT2.o \
	UTLIBR.o \
	FUNDARG.o \
	FCNNUT.o \
	CNMTX.o \
	ORTHO_EOP.o \
	RG_ZONT2.o \
	FCUL_A.o \
	FCUL_B.o \
	FCUL_ZD_HPA.o \
	GMF.o \
	VMF1.o \
	VMF1_HT.o \
	GPT.o \
	GPT2.o \
	ARG2.o

IERS_OBS_D0 = PMSDNUT2_D0.o \
	UTLIBR_D0.o \
	FUNDARG.o \
	FCNNUT.o \
	CNMTX.o \
	ORTHO_EOP.o \
	RG_ZONT2.o \
	FCUL_A.o \
	FCUL_B.o \
	FCUL_ZD_HPA.o \
	GMF.o \
	VMF1.o \
	VMF1_HT.o \
	GPT.o \
	GPT2.o \
	ARG2.o

all: test_iers \
	test_iers_d0 \
	testFundarg \
	testPmsdnut2 \
	testPmsdnut2_d0 \
	testUtlibr \
	testUtlibr_d0 \
	testFcnnut \
	testArg2 \
	testCnmtx \
	testOrthoeop \
	testRgZont2 \
	testFcula \
	testFculb \
	testFculZdhPa \
	testGmf \
	testVmf1 \
	testVmf1Ht \
	testGpt \
	testGpt2 \
	validateFcnnut \
	makeTestRgZont2 \
	makeTestCnmtx \
	makeTestOrthoEop

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
validateFcnnut: validate_fcnnut.f $(IERS_OBS)
	$(FC) -o $@ validate_fcnnut.f $(IERS_OBS)
testArg2: test_arg2.f $(IERS_OBS)
	$(FC) -o $@ test_arg2.f $(IERS_OBS)
testCnmtx: test_cnmtx.f $(IERS_OBS)
	$(FC) -o $@ test_cnmtx.f $(IERS_OBS)
testOrthoeop: test_orthoeop.f $(IERS_OBS)
	$(FC) -o $@ test_orthoeop.f $(IERS_OBS)
testRgZont2: test_rgzont2.f $(IERS_OBS)
	$(FC) -o $@ test_rgzont2.f $(IERS_OBS)
testFcula: test_fcula.f $(IERS_OBS)
	$(FC) -o $@ test_fcula.f $(IERS_OBS)
testFculb: test_fculb.f $(IERS_OBS)
	$(FC) -o $@ test_fculb.f $(IERS_OBS)
testFculZdhPa: test_fcul_zdh_pa.f $(IERS_OBS)
	$(FC) -o $@ test_fcul_zdh_pa.f $(IERS_OBS)
testGmf: test_gmf.f $(IERS_OBS)
	$(FC) -o $@ test_gmf.f $(IERS_OBS)
testVmf1: test_vmf1.f $(IERS_OBS)
	$(FC) -o $@ test_vmf1.f $(IERS_OBS)
testVmf1Ht: test_vmf1_ht.f $(IERS_OBS)
	$(FC) -o $@ test_vmf1_ht.f $(IERS_OBS)
testGpt: test_gpt.f $(IERS_OBS)
	$(FC) -o $@ test_gpt.f $(IERS_OBS)
testGpt2: test_gpt2.f $(IERS_OBS)
	$(FC) -o $@ test_gpt2.f $(IERS_OBS)
makeTestRgZont2: make_test_rgzont2.f $(IERS_OBS)
	$(FC) -o $@ make_test_rgzont2.f $(IERS_OBS)
makeTestCnmtx: make_test_cnmtx.f $(IERS_OBS)
	$(FC) -o $@ make_test_cnmtx.f $(IERS_OBS)
makeTestOrthoEop: make_test_orthoeop.f $(IERS_OBS)
	$(FC) -o $@ make_test_orthoeop.f $(IERS_OBS)

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
VMF1_HT.o : VMF1_HT.F
	$(FC) $(FF) -o $@ VMF1_HT.F
GPT.o : GPT.F
	$(FC) $(FF) -o $@ GPT.F
GPT2.o : GPT2.F
	$(FC) $(FF) -o $@ GPT2.F
