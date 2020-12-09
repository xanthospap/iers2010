      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

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

C     CALL UTLIBR
      CALL UTLIBR(44239.1D0, dUT1, dLOD)

C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> UTLIBR Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(dUT1-dUT1AR)
      WRITE(*, FMT1) DABS(dLOD-dLODAR)

      CALL UTLIBR(55227.4D0, dUT1, dLOD)
C     REPORT RESULTS
      WRITE(*, FMT1) DABS(dUT1-dUT1BR)
      WRITE(*, FMT1) DABS(dLOD-dLODBR)

C     ALL DONE
      END
