      PROGRAM MAIN

      DOUBLE PRECISION T, MJD, MIN_MJD, MAX_MJD, EOP(3)
      INTEGER I,J

*     12 January 2013
      MIN_MJD = 56304D0
      MAX_MJD = 56310D0

      T = MIN_MJD
      DO WHILE (T .LE. MAX_MJD)
        CALL ORTHO_EOP(T, EOP)
        WRITE(*,'(E30.22,X,E25.17,X,E25.17,X,E25.17)') 
     .  T, EOP(1), EOP(2), EOP(3)
        T = T + .01
      END DO

      END
