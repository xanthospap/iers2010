      PROGRAM MAIN

C     RESULT REPORTING
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A30,E30.20,A15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT2 = "(A18,E30.20,A11)"

C     ----------------------------------------------
C     Variables and data for FUNDARG
C     ----------------------------------------------
      DOUBLE PRECISION T, LR, LPR, FR, DR, OMR
      DOUBLE PRECISION    L,  LP,  F,  D,  OM
      DOUBLE PRECISION IJC(24)
      INTEGER I, NUM_TESTS, IMJD, MJD_MAX, MJD_MIN
      DOUBLE PRECISION RND, MJD, JD
      
      MJD_MIN = 38761D0
      MJD_MAX = 62502D0
      NUM_TESTS = 100

C     Test MJD for validation
      IJC = (/ 0D0, 0.1D0, 0.01D0, 0.001D0, 0.0001D0, 0.00001D0,
     . 0.000001D0, 0.0000001D0, 0.00000001D0, 0.000000001D0, 
     . 0.0000000001D0, 0.00000000001D0, 0.000000000001D0,
     . 0.0000000000001D0, 0.00000000000001D0, 0.000000000000001D0,
     . 0.1234567890001D0, 1.1234567890001D0, 2.1234567890001D0,
     . 4.1234567890001D0, 5.1234567890001D0, 1.9999999990001D0,
     . 1.9999999990000D0, 1.9999999999999D0/)

C     CALL FUNDARG( 0.34953976531531522065D0, L,LP,F,D,OM)
C     CALL FUNDARG( -0.34953976531531522065D0, L,LP,F,D,OM)

C
C     call FUNDARG and print C-pseudocode for validation
C     (predefined tests)
C
      DO 10 I=1,24
        CALL FUNDARG(IJC(I), L, LP, F, D, OM)
        WRITE(*, FMT2) 'iers2010::fundarg(', IJC(I),
     . ', fundarg);'
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

C
C     Random Tests
C     
      DO I=1,NUM_TESTS
        IMJD = IRAND(0)
        IF (IMJD.gt.62502) IMJD = IMJD / 62502 + 38761
        IF (IMJD.lt.38761) IMJD = 38761 + IMJD / 1000
        CALL RANDOM_NUMBER(RND)
        MJD = IMJD + RND
        JD = MJD + 2400000.5D0
        T = (JD - 2451545D0)/36525D0
      
        CALL FUNDARG(T, L, LP, F, D, OM)
        WRITE(*, FMT2) 'iers2010::fundarg(', T,
     . ', fundarg);'
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

      END DO

C     ALL DONE
      END
