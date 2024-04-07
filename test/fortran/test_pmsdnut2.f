      PROGRAM MAIN

      DOUBLE PRECISION T, MJD, MIN_MJD, MAX_MJD, PM(2)
      INTEGER NUM_TESTS, I,J

C     NUM_TESTS = 100
C     MIN_MJD = 45700D0
C     MAX_MJD = 62502D0

C     DO I=1,NUM_TESTS-1
C       CALL RANDOM_NUMBER(T)
C       MJD = (MAX_MJD-MIN_MJD)*T + MIN_MJD
C       CALL PMSDNUT2(MJD, PM)
C       WRITE(*,'(E30.22,X,E25.17,X,E25.17)') 
C    .  MJD, PM(1), PM(2)
C     END DO

*     12 January 2013
      MIN_MJD = 56304D0
      MAX_MJD = 56310D0

      T = MIN_MJD
      DO WHILE (T .LE. MAX_MJD)
        CALL PMSDNUT2(T,PM)
        WRITE(*,'(E30.22,X,E25.17,X,E25.17)') 
     .  T, PM(1), PM(2)
        T = T + .01
      END DO

      END
