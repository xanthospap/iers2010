      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

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

C     CALL FCNNUT
      CALL FCNNUT(54790D0, FCX, FCY, FCDX, FCDY)

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
