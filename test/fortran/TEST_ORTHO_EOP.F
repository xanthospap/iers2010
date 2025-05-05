      PROGRAM MAIN

      DOUBLE PRECISION T, MJD, MIN_MJD, MAX_MJD, EOP(3)
      INTEGER NUM_TESTS, I,J

      NUM_TESTS = 1000
      MIN_MJD = 47892D0
      MAX_MJD = 66154D0

      DO I=1,NUM_TESTS-1
        CALL RANDOM_NUMBER(T)
        MJD = (MAX_MJD-MIN_MJD)*T + MIN_MJD
        CALL ORTHO_EOP(MJD, EOP)
        WRITE(*,'(E30.22,X,E25.17,X,E25.17,X,E25.17)') 
     .  MJD, EOP(1), EOP(2), EOP(3)
      END DO

      END
