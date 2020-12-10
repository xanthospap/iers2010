      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for GMF
C     ----------------------------------------------
      DOUBLE PRECISION GMFH, GMFW

      CALL GMF(55055D0, 0.6708665767D0, -1.393397187D0, 844.715D0,
     . 1.278564131D0, GMFH, GMFW)
     
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> GMF Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(GMFH-3.425245519339138678D0)
      WRITE(*, FMT1) DABS(GMFW-3.449589116182419257D0)

C     ALL DONE
      END
