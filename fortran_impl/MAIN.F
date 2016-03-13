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

C     ALL DONE
      END
