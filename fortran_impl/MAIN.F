      PROGRAM MAIN

C     RESULT REPORTING
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"

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

C     ALL DONE
      END
