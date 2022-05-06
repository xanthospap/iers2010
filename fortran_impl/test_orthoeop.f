      PROGRAM MAIN

C     RESULT REPORTING
      CHARACTER(LEN=*), PARAMETER  :: FMT2 = "(A20,E30.20,A16)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A22,E30.20,A15)"

C     ----------------------------------------------
C     Variables and data for ORTHOEOP 
C     ----------------------------------------------
      DOUBLE PRECISION EOP(3)
      INTEGER I, NUM_TESTS, IMJD, MJD_MAX, MJD_MIN
      DOUBLE PRECISION RND, MJD
      
      MJD_MIN = 38761D0
      MJD_MAX = 62502D0
      NUM_TESTS = 1000
      
      DO I=1,NUM_TESTS
        IMJD = IRAND(0)
        IF (IMJD.gt.62502) IMJD = IMJD / 62502 + 38761
        IF (IMJD.lt.38761) IMJD = 38761 + IMJD / 1000
        CALL RANDOM_NUMBER(RND)
        MJD = IMJD + RND
      
        CALL ORTHO_EOP(MJD, EOP)
        WRITE(*, FMT2) 'iers2010::ortho_eop(', MJD,
     . ', dx, dy, dut1);'
        WRITE(*, FMT1) 'assert(std::abs(dx-(',EOP(1),
     .  '))<PRECISION);'
        WRITE(*, FMT1) 'assert(std::abs(dy-(',EOP(2),
     .  '))<PRECISION);'
        WRITE(*, FMT1) 'assert(std::abs(dut1-(',EOP(3),
     .  '))<PRECISION);'

      END DO

C     ALL DONE
      END
