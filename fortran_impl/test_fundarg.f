      PROGRAM MAIN

C     RESULT REPORTING
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A30,E22.15,A15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT2 = "(A22,I3,A12)"

C     ----------------------------------------------
C     Variables and data for FUNDARG
C     ----------------------------------------------
      DOUBLE PRECISION T, LR, LPR, FR, DR, OMR
      DOUBLE PRECISION    L,  LP,  F,  D,  OM
      DOUBLE PRECISION IJC(24)
      INTEGER I

C     Test MJD for validation
      IJC = (/ 0D0, 0.1D0, 0.01D0, 0.001D0, 0.0001D0, 0.00001D0,
     . 0.000001D0, 0.0000001D0, 0.00000001D0, 0.000000001D0, 
     . 0.0000000001D0, 0.00000000001D0, 0.000000000001D0,
     . 0.0000000000001D0, 0.00000000000001D0, 0.000000000000001D0,
     . 0.1234567890001D0, 1.1234567890001D0, 2.1234567890001D0,
     . 4.1234567890001D0, 5.1234567890001D0, 1.9999999990001D0,
     . 1.9999999990000D0, 1.9999999999999D0/)

C    Print the input JC array (C-style)
      WRITE (*,'(A25)') 'const double ijc[] = {'
      DO 11 I=1,23
        WRITE(*, '(E20.15,A1)') IJC(I),','
   11 CONTINUE
      WRITE (*,'(E20.15,A5)') IJC(24),'};'

C     call FUNDARG and print C-pseudocode for validation
      DO 10 I=1,24
        CALL FUNDARG(IJC(I), L, LP, F, D, OM)
        WRITE(*, FMT2) 'iers2010::fundarg(ijc[', I-1, 
     . '], fundarg);'
        WRITE(*, FMT1) 'assert(std::abs(fundarg[0]-(',L,
     .  '))<PRECISION);'
        WRITE(*, FMT1) 'assert(std::abs(fundarg[1]-(',LP,
     .  '))<PRECISION);'
        WRITE(*, FMT1) 'assert(std::abs(fundarg[2]-(',F,
     .  '))<PRECISION);'
        WRITE(*, FMT1) 'assert(std::abs(fundarg[3]-(',D,
     .  '))<PRECISION);'
        WRITE(*, FMT1) 'assert(std::abs(fundarg[4]-(',OM,
     .  '))<PRECISION);'
   10 CONTINUE

C     ALL DONE
      END
