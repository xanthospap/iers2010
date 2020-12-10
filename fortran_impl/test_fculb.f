      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for FCUL_B
C     ----------------------------------------------
      DOUBLE PRECISION FCULB

C     CALL FUNDARG
      FCULB = FCUL_B(30.67166667D0, 2075D0, 224D0, 15D0)

C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> FCUL_B Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(FCULB-3.800758725284345996D0)

C     ALL DONE
      END
