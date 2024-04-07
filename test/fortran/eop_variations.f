      PROGRAM MAIN

*
*  Compute EOP variations due to:
*  a) Ocen tides and
*  b) Libration
*  based on interp.f and UTLIBR.
*  Units are microarcseconds and microseconds.

      DOUBLE PRECISION T, MJD, MIN_MJD, MAX_MJD, DT
      DOUBLE PRECISION COR_X,COR_Y,COR_UT1,COR_LOD
      INTEGER I,J

*     12 January 2013
      MIN_MJD = 56304D0
      MAX_MJD = 56310D0

      T = MIN_MJD
      DO WHILE (T .LE. MAX_MJD)
        DT = T - MIN_MJD
*       Ocean tide variations
        CALL PMUT1_OCEANS(T,COR_X,COR_Y,COR_UT1,COR_LOD)
        WRITE(*,'(A10,E30.22,X,E25.17,X,E25.17,X,E25.17,X,E25.17)') 
     .  'OCTIDE', T,COR_X,COR_Y,COR_UT1,COR_LOD
*       Libration variations on polar motion
        CALL PM_GRAVI(T,COR_X,COR_Y)
        CALL UTLIBR(T,COR_UT1,COR_LOD)
        WRITE(*,'(A10,E30.22,X,E25.17,X,E25.17,X,E25.17,X,E25.17)') 
     .  'LIBRTN', T,COR_X,COR_Y,COR_UT1,COR_LOD
        T = T + .01
      END DO

      END
