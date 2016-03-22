      PROGRAM MAIN

C     RESULT REPORTING
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT2 = "(F20.15)"

C     ----------------------------------------------
C     Variables and data for FUNDARG
C     ----------------------------------------------
      DOUBLE PRECISION T, LR, LPR, FR, DR, OMR
      DOUBLE PRECISION    L,  LP,  F,  D,  OM
      DATA T, LR, LPR, FR, DR, OMR/0.07995893223819302D0,
     .   2.291187512612069099D0,
     .   6.212931111003726414D0,
     .   3.658025792050572989D0,
     .   4.554139562402433228D0,
     .   -0.5167379217231804489D0/

C     ----------------------------------------------
C     Variables and data for PMSDNUT2
C     ----------------------------------------------
      DOUBLE PRECISION PM(2), PMR(2)
      DATA (PMR(I),I=1,2) /24.83144238273364834D0,
     . -14.09240692041837661D0/

C     ----------------------------------------------
C     Variables and data for UTLIBR
C     ----------------------------------------------
      DOUBLE PRECISION  dUT1,  dLOD
      DOUBLE PRECISION dUT1AR, dUT1BR, dLODAR, dLODBR
      DATA dUT1AR, dUT1BR, dLODAR, dLODBR/
     . 2.441143834386761746D0,
     . -2.655705844335680244D0,
     . -14.78971247349449492D0,
     . 27.39445826599846967D0/

C     ----------------------------------------------
C     Variables and data for FCNNUT
C     ----------------------------------------------
      DOUBLE PRECISION  FCX,  FCY,  FCDX,  FCDY
      DOUBLE PRECISION FCXR, FCYR, FCDXR, FCDYR
      DATA FCXR, FCYR, FCDXR, FCDYR/
     . -176.8012290066270680D0,
     . -93.51855308903756736D0,
     . 3.745573770491803067D0,
     . 3.745573770491803067D0/

C     ----------------------------------------------
C     Variables and data for RG_ZONT2
C     ----------------------------------------------
      DOUBLE PRECISION ZDUTR, ZDLODR, ZDOMEGAR
      DOUBLE PRECISION ZDUT,  ZDLOD,  ZDOMEGA
      DATA ZDUTR, ZDLODR, ZDOMEGAR/
     . 7.983287678576557467E-002,
     . 5.035331113978199288E-005,
     . -4.249711616463017E-014/

C     ----------------------------------------------
C     Variables and data for VMF1
C     ----------------------------------------------
      DOUBLE PRECISION VMF1HR, VMF1WR, VMF1H, VMF1W
      DATA VMF1HR, VMF1WR/ 3.424342122738070593D0,
     . 3.448299714692572238D0/

C     ----------------------------------------------
C     Variables and data for GPT2 
C     ----------------------------------------------
      DOUBLE PRECISION GPT_P, GPT_T, GPT_DT, GPT_E,
     . GPT_AH, GPT_AW, GPT_UNDU

C     CALL FUNDARG
      CALL FUNDARG(T, L, LP, F, D, OM)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> FUNDARG Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(L-LR)
      WRITE(*, FMT1) DABS(LP-LPR)
      WRITE(*, FMT1) DABS(F-FR)
      WRITE(*, FMT1) DABS(D-DR)
      WRITE(*, FMT1) DABS(OM-OMR)

C     CALL PMSDNUT2
      T = 54335D0
      CALL PMSDNUT2(T, PM)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> PMSDNUT2 Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(PM(1)-PMR(1))
      WRITE(*, FMT1) DABS(PM(2)-PMR(2))

C     CALL UTLIBR
      T = 44239.1D0
      CALL UTLIBR(T, dUT1, dLOD)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> UTLIBR Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(dUT1-dUT1AR)
      WRITE(*, FMT1) DABS(dLOD-dLODAR)
      T = 55227.4D0
      CALL UTLIBR(T, dUT1, dLOD)
C     REPORT RESULTS
      WRITE(*, FMT1) DABS(dUT1-dUT1BR)
      WRITE(*, FMT1) DABS(dLOD-dLODBR)

C     CALL FCNNUT
      T = 54790D0
      CALL FCNNUT(T, FCX, FCY, FCDX, FCDY)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> FCNNUT Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(FCX-FCXR)
      WRITE(*, FMT1) DABS(FCY-FCYR)
      WRITE(*, FMT1) DABS(FCDX-FCDXR)
      WRITE(*, FMT1) DABS(FCDY-FCDYR)

C     CALL RG_ZONT2
      T = .07995893223819302D0
      CALL RG_ZONT2(T, ZDUT, ZDLOD, ZDOMEGA)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> RG_ZONT2 Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(ZDUT-ZDUTR)
      WRITE(*, FMT1) DABS(ZDLOD-ZDLODR)
      WRITE(*, FMT1) DABS(ZDOMEGA-ZDOMEGAR)

C     CALL VMF1
      CALL VMF1(0.00127683D0, 0.00060955D0, 55055D0,
     .  0.6708665767D0, 1.278564131D0, VMF1H, VMF1W)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> VMF1 Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(VMF1H-VMF1HR)
      WRITE(*, FMT1) DABS(VMF1W-VMF1WR)

C     CALL GPT2
      CALL GPT2(56141.d0,0.8412486994612668D0,0.28571039855147173D0,
     .  156.D0, 1, 0, GPT_P, GPT_T, GPT_DT, GPT_E, GPT_AH,
     .  GPT_AW, GPT_UNDU)
C     REPORT RESULTS
C      PRINT *, '----------------------------------------'
C      PRINT *, '> GPT2 Results: (Test A)'
C      PRINT *, '----------------------------------------'
C      WRITE(*, FMT2) GPT_P
C      WRITE(*, FMT1) GPT_T
C      WRITE(*, FMT1) GPT_DT
C      WRITE(*, FMT1) GPT_E
C      WRITE(*, FMT1) GPT_AH
C      WRITE(*, FMT1) GPT_AW
C      WRITE(*, FMT1) GPT_UNDU
C     CALL GPT2
C      CALL GPT2(56141.d0,0.8412486994612668D0,0.28571039855147173D0,
C     .  156.D0, 1, 1, GPT_P, GPT_T, GPT_DT, GPT_E, GPT_AH,
C     .  GPT_AW, GPT_UNDU)
C     REPORT RESULTS
C      PRINT *, '----------------------------------------'
C      PRINT *, '> GPT2 Results: (Test B)'
C      PRINT *, '----------------------------------------'
C      WRITE(*, FMT2) GPT_P
C      WRITE(*, FMT1) GPT_T
C      WRITE(*, FMT1) GPT_DT
C      WRITE(*, FMT1) GPT_E
C      WRITE(*, FMT1) GPT_AH
C      WRITE(*, FMT1) GPT_AW
C      WRITE(*, FMT1) GPT_UNDU

C     ALL DONE
      END
