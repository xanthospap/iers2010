      SUBROUTINE PM_GRAVI (rjd,cor_x,cor_y)
C
C    This subroutine provides, in time domain, the diurnal
C    lunisolar effet on polar motion (")
C    
C    N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
C
C    These corrections should be added to "average"
C    EOP values to get estimates of the instantaneous values.
C
C     PARAMETERS ARE :
C     rjd      - epoch of interest given in mjd
C     cor_x    - tidal correction in x (sec. of arc)
C     cor_y    - tidal correction in y (sec. of arc)
C
C     coded by Ch. Bizouard (2002)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER nlines
      PARAMETER(nlines=10)
      DOUBLE PRECISION ARG(6)    ! Array of the tidal arguments   
      REAL*4 XCOS(nlines),XSIN(nlines),YCOS(nlines),YSIN(nlines)
      INTEGER NARG(nlines,6)
      
      halfpi = 1.5707963267948966d0
      secrad=2.d0*halfpi/(180.d0*3600.d0)	

c  Diurnal lunisolar tidal terms present in x (microas),y(microas)      
c  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
       data( 
     & NARG(j,1),NARG(j,2),NARG(j,3),NARG(j,4),NARG(j,5),NARG(j,6),
     & XSIN(j),XCOS(j),YSIN(j),YCOS(j),j=1,nlines)/    
     & 1,-1, 0,-2, 0,-1,    -.44,   .25,   -.25,  -.44,
     & 1,-1, 0,-2, 0,-2,   -2.31,  1.32,  -1.32, -2.31,
     & 1, 1, 0,-2,-2,-2,    -.44,   .25,   -.25,  -.44,
     & 1, 0, 0,-2, 0,-1,   -2.14,  1.23,  -1.23, -2.14,
     & 1, 0, 0,-2, 0,-2,  -11.36,  6.52,  -6.52,-11.36,
     & 1,-1, 0, 0, 0, 0,     .84,  -.48,    .48,   .84,
     & 1, 0, 0,-2, 2,-2,   -4.76,  2.73,  -2.73, -4.76,
     & 1, 0, 0, 0, 0, 0,   14.27, -8.19,   8.19, 14.27,
     & 1, 0, 0, 0, 0,-1,    1.93, -1.11,   1.11,  1.93,
     & 1, 1, 0, 0, 0, 0,     .76,  -.43,    .43,   .76/
 
      T = (rjd - 51544.5D0)/36525.0D0  ! julian century

C Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
C et leur derivee temporelle 

      ARG(1) = (67310.54841d0 +
     .        (876600d0*3600d0 + 8640184.812866d0)*T +
     .         0.093104d0*T**2 -
     .         6.2d-6*T**3)*15.0d0 + 648000.0d0
      ARG(1)=dmod(ARG(1),1296000d0)*secrad 
   

      ARG(2) = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2
     .  + 1717915923.2178d0*T + 485868.249036d0
      ARG(2) = DMOD(ARG(2),1296000d0)*secrad
      

      ARG(3) = -0.00001149d0*T**4 - 0.000136d0*T**3
     .  -  0.5532d0*T**2
     .  + 129596581.0481d0*T + 1287104.79305d0
      ARG(3) = DMOD(ARG(3),1296000d0)*secrad

          
      ARG(4) = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2
     .  + 1739527262.8478d0*T + 335779.526232d0
      ARG(4) = DMOD(ARG(4),1296000d0)*secrad

    
      ARG(5) = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2
     .  + 1602961601.2090d0*T + 1072260.70369d0
      ARG(5) = DMOD(ARG(5),1296000d0)*secrad

  
      ARG(6) = -0.00005939d0*T**4 + 0.007702d0*T**3
     .  + 7.4722d0*T**2
     .  - 6962890.2665d0*T + 450160.398036d0
      ARG(6) = DMOD(ARG(6),1296000d0)*secrad


C CORRECTIONS

	cor_x  = 0.d0
	cor_y  = 0.d0

 	do j=1,nlines
 	
	ag  = 0.d0
		do i=1,6
 		ag  = ag  + dble(narg(j,i))*ARG(i)
		enddo
	ag=dmod(ag,4.d0*halfpi)

        cor_x =cor_x+dble(XCOS(j))*dcos(ag)+dble(XSIN(j))*dsin(ag)
        cor_y =cor_y+dble(YCOS(j))*dcos(ag)+dble(YSIN(j))*dsin(ag) 

        enddo
  
      cor_x = cor_x !* 1.0d-6   ! arcsecond (")
      cor_y = cor_y !* 1.0d-6   ! arcsecond (")
 
      RETURN

      END
