      PROGRAM MAIN

C     RESULT REPORTING
      CHARACTER(LEN=*), PARAMETER  :: FMT2 = "(A15,I4,A1,E30.20,A10)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A23,I2,A3,E30.20,A15)"

C     ----------------------------------------------
C     Variables and data for ARG2
C     ----------------------------------------------
      DOUBLE PRECISION ANGLE(11)
      INTEGER I, NUM_TESTS, J
      INTEGER IYEAR, IDOY
      DOUBLE PRECISION FDOY, RND

      NUM_TESTS = 100
      
      DO I=1,NUM_TESTS
        CALL RANDOM_NUMBER(RND)
        IYEAR = 1973 + FLOOR(80 * RND)
        IF (IYEAR.gt.2050) IYEAR = 2050
        IF (IYEAR.lt.1974) IYEAR = 1974

        CALL RANDOM_NUMBER(RND)
        IDOY = FLOOR( RND * 366 )
        IF (IDOY.gt.365) IDOY = 365
        IF (IDOY.gt.400) IDOY = 0

        CALL RANDOM_NUMBER(RND)
        FDOY = IDOY + RND
        
        CALL ARG2(IYEAR, FDOY, ANGLE)
        WRITE(*, FMT2) 'iers2010::arg2(',IYEAR,',',FDOY,
     . ', angles);'
        DO J=1,11
          WRITE(*, FMT1) 'assert(std::abs(angles[',J-1,']-(',ANGLE(J),
     .    '))<PRECISION);'
        END DO

      END DO

C     ALL DONE
      END
