      PROGRAM MAIN

C     ----------------------------------------------
C     Variables and data for UTLIBR
C     ----------------------------------------------
      DOUBLE PRECISION  dUT1,  dLOD
      INTEGER I, NUM_TESTS, IMJD, MJD_MAX, MJD_MIN
      DOUBLE PRECISION RND, MJD, JD
      
      MJD_MIN = 38761D0
      MJD_MAX = 62502D0
      NUM_TESTS = 1000

C
C     Random Tests
C     
      DO I=1,NUM_TESTS
        IMJD = IRAND(0)
        IF (IMJD.gt.62502) IMJD = IMJD / 62502 + 38761
        IF (IMJD.lt.38761) IMJD = 38761 + IMJD / 1000
        CALL RANDOM_NUMBER(RND)
        MJD = IMJD + RND
C       JD = MJD + 2400000.5D0
C       T = (JD - 2451545D0)/36525D0
      
        CALL UTLIBR(MJD, dUT1, dLOD)
        WRITE(*, '(A17,E30.20,X,A14)') 'iers2010::utlibr(', MJD,
     . ', dut1, dlod);'
        WRITE(*, '(A22,E30.20,A15)') 'assert(std::abs(dut1-(',dUT1,
     .  '))<PRECISION1);'
        WRITE(*, '(A22,E30.20,A15)') 'assert(std::abs(dlod-(',dLOD,
     .  '))<PRECISION2);'

      END DO

C     ALL DONE
      END
