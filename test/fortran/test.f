      PROGRAM MAIN
      INTEGER YEAR,I
      DOUBLE PRECISION T,L,LP,F,D,OM,X,DX,Y,DY,PM(2)
      DOUBLE PRECISION DUT1,DLOD,IANGLE(11),ANGLE(11),DAY

*     ------------------------------------------------------------------------
*     FUNDARG
*     ------------------------------------------------------------------------
      WRITE (*,*) 'Function FUNADRG'
      WRITE (*,*) 'Abs. differences in radians'
      T = 0.07995893223819302D0
      CALL FUNDARG ( T, L, LP, F, D, OM )
      WRITE (*,10) '|dl|  = ',DABS (2.291187512612069099D0-L)
      WRITE (*,10) '|dlp| = ',DABS (6.212931111003726414D0-LP)
      WRITE (*,10) '|df|  = ',DABS (3.658025792050572989D0-F)
      WRITE (*,10) '|dd|  = ',DABS (4.554139562402433228D0-D)
      WRITE (*,10) '|dom| = ',DABS (-0.5167379217231804489D0-OM)
 10   FORMAT (A,F20.18)

*     ------------------------------------------------------------------------
*     FCNNUT
*     ------------------------------------------------------------------------
      WRITE (*,*) 'Function FCNNUT'
      WRITE (*,*) 'Abs. differences in microarcseconds:'
      T = 54790D0
      CALL FCNNUT ( T,X,Y,DX,DY )
      WRITE (*,11) '|dx|  = ',DABS (-176.8012290066270680D0-X)
      WRITE (*,11) '|dy|  = ',DABS (-93.51855308903756736D0-Y)
      WRITE (*,11) '|ddx| = ',DABS (3.745573770491803067D0-DX)
      WRITE (*,11) '|ddy| = ',DABS (3.745573770491803067D0-DY)
 11   FORMAT (A,F20.18)

*     ------------------------------------------------------------------------
*     PMSDNUT2
*     ------------------------------------------------------------------------
      WRITE (*,*) 'Function PMSDNUT2'
      WRITE (*,*) 'Abs. differences in microarcseconds:'
      T = 54335D0
      CALL PMSDNUT2 (T, PM)
      WRITE (*,12) '|pm(1)|  = ',DABS (24.83144238273364834D0-PM(1))
      WRITE (*,12) '|pm(2)|  = ',DABS (-14.09240692041837661D0-PM(2))
 12   FORMAT (A,F20.18)

*     ------------------------------------------------------------------------
*     UTLIBR
*     ------------------------------------------------------------------------
      WRITE (*,*) 'Function UTLIBR (2 tests)'
      WRITE (*,*) 'Abs. differences in mus and mus/day:'
      T = 44239.1D0
      CALL UTLIBR (T,DUT1,DLOD)
      WRITE (*,*) '->test A.'
      WRITE (*,13) '|ddut1|  = ',DABS (2.441143834386761746D0-DUT1)
      WRITE (*,13) '|ddut1|  = ',DABS (-14.78971247349449492D0-DLOD)
      T = 55227.4D0
      CALL UTLIBR (T,DUT1,DLOD)
      WRITE (*,*) '->test B.'
      WRITE (*,13) '|ddut1|  = ',DABS (-2.655705844335680244D0-DUT1)
      WRITE (*,13) '|ddut1|  = ',DABS (27.39445826599846967D0-DLOD)
 13   FORMAT (A,F20.18)

*     ------------------------------------------------------------------------
*     ARG2
*     ------------------------------------------------------------------------
      WRITE (*,*) 'Function ARG2'
      WRITE (*,*) 'Abs. differences in radians:'
      ANGLE(1)  = 2.849663065753787805D0
      ANGLE(2)  = 6.28318080000000023D0
      ANGLE(3)  = 4.926040134021299366D0
      ANGLE(4)  = 1.608450491115348768D0
      ANGLE(5)  = 2.375021572352622456D0
      ANGLE(6)  = 0.4746414933980958040D0
      ANGLE(7)  = 3.908159227647345801D0
      ANGLE(8)  = 2.551018561669245344D0
      ANGLE(9)  = 5.041990012540757959D0
      ANGLE(10) = 4.206816878908014701D0
      ANGLE(11) = 1.608463638294885811D0
      YEAR = 2008
      DAY = 311.5D0
      CALL ARG2 (YEAR,DAY,IANGLE)
      DO I=1,11
        WRITE (*,14) '|dangle(',I,')|  = ', DABS (ANGLE(I)-IANGLE(I))
      ENDDO
 14   FORMAT (A,I2,A,F20.18)

*     ------------------------------------------------------------------------
*     RGZONT2
*     ------------------------------------------------------------------------
      WRITE (*,*) 'Function RG_ZONT2'
      WRITE (*,*) 'Abs. differences in (see results):'
      T = .07995893223819302D0;
      CALL RG_ZONT2 (T, X, DX, Y);
      WRITE (*,15) '|ddut|   = ',DABS(7.983287678576557467D-002
     .  -X),' seconds'
      WRITE (*,15) '|ddlod|  = ',DABS(5.035331113978199288D-005
     .  -DX),' seconds/day'
      WRITE (*,15) '|ddomega|= ',DABS(-4.249711616463017D-014
     .  -Y),' radians/second'
 15   FORMAT (A,F20.18,A)

*     ------------------------------------------------------------------------
*     ORTHO_EOP
*     ------------------------------------------------------------------------
      WRITE (*,*) 'Function ORTHO_EOP'
      WRITE (*,*) 'Abs. differences in microarcseconds:'
      T = 47100D0
      CALL ORTHO_EOP (T,ANGLE);
      WRITE (*,16) '|deop(1)|   = ',DABS(-162.8386373279636530D0
     . -ANGLE(1))
      WRITE (*,16) '|deop(2)|   = ',DABS(117.7907525842668974D0
     . -ANGLE(2))
      WRITE (*,16) '|deop(3)|   = ',DABS(-23.39092370609808214D0
     . -ANGLE(3))
 16   FORMAT (A,F20.18)

      END PROGRAM MAIN