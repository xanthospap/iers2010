      PROGRAM MAIN

      DOUBLE PRECISION T, MJD, MIN_MJD, MAX_MJD, EOP(3), PM(2), DUT1, DLOD
      INTEGER I,J

      MIN_MJD = 59945D0
      MAX_MJD = 59950D0

      T = MIN_MJD
      DO WHILE (T .LE. MAX_MJD)
        CALL ORTHO_EOP(T, EOP)
        CALL PMSDNUT2(T, PM)
        CALL UTLIBR(T, DUT1, DLOD)
        WRITE(*,'(E30.22,X,E25.17,X,E25.17,X,E25.17,X,E25.17,X,E25.17,
     .  X,E25.17,X,E25.17)') 
     .  T, EOP(1), EOP(2), EOP(3), PM(1), PM(2), DUT1, DLOD
        T = T + .05
      END DO

      END
