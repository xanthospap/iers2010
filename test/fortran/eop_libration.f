      PROGRAM MAIN

      DOUBLE PRECISION T, MJD, MIN_MJD, MAX_MJD
      DOUBLE PRECISION COR_X,COR_Y,COR_UT1,COR_LOD
      INTEGER I,J

*     12 January 2013
      MIN_MJD = 56304D0
      MAX_MJD = 56310D0

      T = MIN_MJD
      DO WHILE (T .LE. MAX_MJD)
*       Libration variations on polar motion
        CALL PM_GRAVI(T,COR_X,COR_Y)
        CALL UTLIBR(T,COR_UT1,COR_LOD)
        WRITE(*,'(E30.22,X,E25.17,X,E25.17,X,E25.17,X,E25.17)') 
     .  T,COR_X,COR_Y,COR_UT1,COR_LOD
        T = T + .01
      END DO

      END
