      PROGRAM MAIN

      INTEGER NUM_TESTS, IMJD
      DOUBLE PRECISION PM(2), PMR(2), MJD_MIN, MJD_MAX, MJD
      REAL RND
      
      MJD_MIN = 38761D0
      MJD_MAX = 62502D0
      NUM_TESTS = 1000
      
      DO I=1,NUM_TESTS
        IMJD = IRAND(0)
        IF (IMJD.gt.62502) IMJD = IMJD / 62502 + 38761
        IF (IMJD.lt.38761) IMJD = 38761 + IMJD / 1000
        CALL RANDOM_NUMBER(RND)
        MJD = IMJD + RND
      
        CALL PMSDNUT2(MJD, PM)
        WRITE(*,'(A19,E25.20,A6)') 'iers2010::pmsdnut2(',MJD,',x,y);'
        WRITE(*,'(A19,E30.20,A14)') 'assert(std::abs(x-(',
     .    PM(1),'))<PRECISION);'
        WRITE(*,'(A19,E30.20,A14)') 'assert(std::abs(y-(',
     .    PM(2),'))<PRECISION);'

      END DO

C     ALL DONE
      END
