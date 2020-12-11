      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for GPT
C     ----------------------------------------------
      DOUBLE PRECISION PRES, TEMP, UNDU

      CALL GPT(55055D0, 0.6708665767D0, -1.393397187D0, 812.546D0,
     .  PRES, TEMP, UNDU)
     
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> GPT Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(PRES-918.0710638757363995D0)
      WRITE(*, FMT1) DABS(TEMP-19.31914181012882992D0)
      WRITE(*, FMT1) DABS(UNDU+42.19185643717770517D0)

C     ALL DONE
      END
