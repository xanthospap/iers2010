      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for FCUL_A
C     ----------------------------------------------
      DOUBLE PRECISION FCULA

C     CALL FUNDARG
      FCULA = FCUL_A(30.67166667D0, 2075D0, 300.15D0, 15D0)

C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> FCUL_A Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(FCULA-3.800243667312344087D0)

C     ALL DONE
      END
