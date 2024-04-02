      PROGRAM MAIN

      DOUBLE PRECISION T, MJD, MIN_MJD, MAX_MJD, EOP(3), PM(2), DUT1, DLOD
      DOUBLE PRECISION IEOP(4)
      INTEGER I,J

*     12 January 2013
      MIN_MJD = 56304D0
      MAX_MJD = 56310D0

      T = MIN_MJD
      DO WHILE (T .LE. MAX_MJD)
        CALL ORTHO_EOP(T, EOP)
        CALL PMSDNUT2(T, PM)
        CALL UTLIBR(T, DUT1, DLOD)
        CALL PMUT1_OCEANS(T,EOP(1),EOP(2),EOP(3),DLOD)
        CALL PM_GRAVI(T,EOP(1),EOP(2))
        WRITE(*,'(E30.22,X,E25.17,X,E25.17,X,E25.17,X,E25.17,X,E25.17,
     .  X,E25.17,X,E25.17)') 
     .  T, EOP(1), EOP(2), EOP(3), PM(1), PM(2), DUT1, DLOD
        T = T + .01
      END DO

      T = 56304D0
      CALL ORTHO_EOP(T, EOP)
      WRITE(*,'(A20,F10.3,X,F10.3,X,F10.3)') 
     .'12 January 2013', EOP(1), EOP(2), EOP(3) 
      CALL PMUT1_OCEANS(T,EOP(1),EOP(2),EOP(3),DLOD)
      WRITE(*,'(A20,F10.3,X,F10.3,X,F10.3)') 
     .'12 January 2013', EOP(1), EOP(2), EOP(3) 

      T = 47100.01D0
      CALL ORTHO_EOP(T, EOP)
      WRITE(*,'(A20,F10.3,X,F10.3,X,F10.3)') 
     .'Test Case', EOP(1), EOP(2), EOP(3) 


      END
