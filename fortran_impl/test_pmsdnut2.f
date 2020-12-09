      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for PMSDNUT2
C     ----------------------------------------------
      DOUBLE PRECISION PM(2), PMR(2)
      DATA (PMR(I),I=1,2) /24.83144238273364834D0,
     . -14.09240692041837661D0/

      CALL PMSDNUT2(54335D0, PM)

C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> PMSDNUT2 Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(PM(1)-PMR(1))
      WRITE(*, FMT1) DABS(PM(2)-PMR(2))

C     ALL DONE
      END
