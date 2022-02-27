      PROGRAM MAIN

      INTEGER NUM_TESTS, IMJD
      DOUBLE PRECISION PM(2), PMR(2), MJD_MIN, MJD_MAX, MJD
      REAL RND
      
      MJD_MIN = 38761D0
      MJD_MAX = 62502D0
      NUM_TESTS = 100
      
      DO I=1,NUM_TESTS
      IMJD = IRAND(0)
      IF (IMJD.gt.62502) IMJD = IMJD / 62502 + 38761
      IF (IMJD.lt.38761) IMJD = 38761 + IMJD / 1000
      CALL RANDOM_NUMBER(RND)
      MJD = IMJD + RND
      CALL PMSDNUT2(MJD, PM)
      WRITE(*,'(A17I3A15A30)') 'SUBCASE("Example ',I,': Checking for ',
     $ 'descripancies  > PRECISION") {'
      WRITE(*, '(A4F25.15A8F20.15A8F20.15A3)') 'mjd=', MJD, 
     $ 'e0; dxr=', PM(1), 'e0; dyr=', PM(2), 'e0;'
      WRITE(*,'(A22)') 'pmsdnut2(mjd, dx, dy);'
      WRITE(*,'(A38)') 'REQUIRE(std::abs(dx-dxr) < PRECISION);'
      WRITE(*,'(A38)') 'REQUIRE(std::abs(dy-dyr) < PRECISION);'
      WRITE(*,'(A4)') '}'
      END DO

C     ALL DONE
      END
