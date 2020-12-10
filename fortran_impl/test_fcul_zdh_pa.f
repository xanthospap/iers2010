      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for FCUL_ZDH_PA
C     ----------------------------------------------
      DOUBLE PRECISION FCUL_ZTD, FCUL_ZHD, FCUL_ZWD

C     CALL FCUL_ZD_HPA
      CALL FCULZD_HPA(30.67166667D0, 2010.344D0, 798.4188D0,
     . 14.322D0, 0.532D0, FCUL_ZTD, FCUL_ZHD, FCUL_ZWD)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> FCULZD_HPA  Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(FCUL_ZTD-1.935225924846803114D0)
      WRITE(*, FMT1) DABS(FCUL_ZHD-1.932992176591644462D0)
      WRITE(*, FMT1) DABS(FCUL_ZWD-0.2233748255158703871D-02)

C     ALL DONE
      END
