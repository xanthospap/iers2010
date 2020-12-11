      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for GPT2 
C     ----------------------------------------------
      DOUBLE PRECISION GPT_P(1), GPT_T(1), GPT_DT(1), GPT_E(1),
     . GPT_AH(1), GPT_AW(1), GPT_UNDU(1)
      DOUBLE PRECISION GPT21LON(1),GPT21LAT(1),GPT21HEL(1)

      CALL GPT2(56141.d0,0.8412486994612668D0,0.28571039855147173D0,
     .  156.D0, 1, 0, GPT_P, GPT_T, GPT_DT, GPT_E, GPT_AH,
     .  GPT_AW, GPT_UNDU)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> GPT2 Results: (Test A)'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(GPT_P-1002.56D0)
      WRITE(*, FMT1) DABS(GPT_T-22.12D0)
      WRITE(*, FMT1) DABS(GPT_DT+6.53D0)
      WRITE(*, FMT1) DABS(GPT_E-15.63D0)
      WRITE(*, FMT1) DABS(GPT_AH-0.0012647D0)
      WRITE(*, FMT1) DABS(GPT_AW-0.0005726D0)
      WRITE(*, FMT1) DABS(GPT_UNDU-44.06D0)
      
      CALL GPT2(56141.d0,0.8412486994612668D0,0.28571039855147173D0,
     .  156.D0, 1, 1, GPT_P, GPT_T, GPT_DT, GPT_E, GPT_AH,
     .  GPT_AW, GPT_UNDU)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> GPT2 Results: (Test B)'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(GPT_P-1003.49D0)
      WRITE(*, FMT1) DABS(GPT_T-11.95D0)
      WRITE(*, FMT1) DABS(GPT_DT+5.47D0)
      WRITE(*, FMT1) DABS(GPT_E-9.58D0)
      WRITE(*, FMT1) DABS(GPT_AH-0.0012395D0)
      WRITE(*, FMT1) DABS(GPT_AW-0.0005560D0)
      WRITE(*, FMT1) DABS(GPT_UNDU-44.06D0)

C     ALL DONE
      END
