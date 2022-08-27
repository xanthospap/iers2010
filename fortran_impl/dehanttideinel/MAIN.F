      IMPLICIT NONE
      
      INTEGER YR,MONTH,DAY,EX,I,SECONDS
      DOUBLE PRECISION XSTA(3),XSUN(3),XMON(3),DXTIDE(3),FHR

C     TEST CASE 1; units [m]
      XSTA(1) = 4075578.385D0
      XSTA(2) =  931852.890D0
      XSTA(3) = 4801570.154D0
      XSUN(1) = 137859926952.015D0
      XSUN(2) = 54228127881.4350D0
      XSUN(3) = 23509422341.6960D0
      XMON(1) = -179996231.920342D0
      XMON(2) = -312468450.131567D0
      XMON(3) = -169288918.592160D0
      YR      = 2009
      MONTH   = 4
      DAY     = 13
      FHR     = 0.00D0
      EX      = 1
      WRITE(*,'(A17I3A15A30)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > PRECISION") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A8F15.4A12F15.4A12F15.4A3)') 'xsta(0)=', XSTA(1), 
     $ 'e0; xsta(1)=', XSTA(2), 'e0; xsta(2)=', XSTA(3), 'e0;'
      WRITE(*,'(A8F20.4A12F20.4A12F20.4A3)') 'xsun(0)=', XSUN(1), 
     $ 'e0; xsun(1)=', XSUN(2), 'e0; xsun(2)=', XSUN(3) , 'e0;'
      WRITE(*,'(A8F20.7A12F20.7A12F20.7A3)') 'xmon(0)=', XMON(1), 
     $ 'e0; xmon(1)=', XMON(2), 'e0; xmon(2)=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A9F20.10A13F20.10A13F20.10A3)') 'xtide(0)=',
     $ DXTIDE(1), 'e0; xtide(1)=', DXTIDE(2), 'e0; xtide(2)=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A46)') 'res=dehanttideinel_from_utc(xsta,xsun,xmon,t);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(0)-xtide(0))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(1)-xtide(1))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(2)-xtide(2))<PRECISION);'
      WRITE(*,'(A4)') '}'

C     gradually increase fractional hours, 1hour every step
      DO I=1,23,3
      FHR = I
      EX = EX + 1
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A17I3A15A30)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > PRECISION") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A8F15.4A12F15.4A12F15.4A3)') 'xsta(0)=', XSTA(1), 
     $ 'e0; xsta(1)=', XSTA(2), 'e0; xsta(2)=', XSTA(3), 'e0;'
      WRITE(*,'(A8F20.4A12F20.4A12F20.4A3)') 'xsun(0)=', XSUN(1), 
     $ 'e0; xsun(1)=', XSUN(2), 'e0; xsun(2)=', XSUN(3) , 'e0;'
      WRITE(*,'(A8F20.7A12F20.7A12F20.7A3)') 'xmon(0)=', XMON(1), 
     $ 'e0; xmon(1)=', XMON(2), 'e0; xmon(2)=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A9F20.10A13F20.10A13F20.10A3)') 'xtide(0)=',
     $ DXTIDE(1), 'e0; xtide(1)=', DXTIDE(2), 'e0; xtide(2)=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A46)') 'res=dehanttideinel_from_utc(xsta,xsun,xmon,t);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(0)-xtide(0))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(1)-xtide(1))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(2)-xtide(2))<PRECISION);'
      WRITE(*,'(A4)') '}'
      END DO

C     TEST CASE 2; units [m]
      XSTA(1) =  1112189.660D0
      XSTA(2) = -4842955.026D0
      XSTA(3) =  3985352.284D0
      XSUN(1) = -54537460436.2357D0
      XSUN(2) =  130244288385.279D0
      XSUN(3) =  56463429031.5996D0
      XMON(1) =  300396716.912D0
      XMON(2) =  243238281.451D0
      XMON(3) =  120548075.939D0
      YR      = 2012
      MONTH   = 7
      DAY     = 13
      FHR     = 0.00D0
      EX = EX + 1
      WRITE(*,'(A17I3A15A30)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > PRECISION") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A8F15.4A12F15.4A12F15.4A3)') 'xsta(0)=', XSTA(1), 
     $ 'e0; xsta(1)=', XSTA(2), 'e0; xsta(2)=', XSTA(3), 'e0;'
      WRITE(*,'(A8F20.4A12F20.4A12F20.4A3)') 'xsun(0)=', XSUN(1), 
     $ 'e0; xsun(1)=', XSUN(2), 'e0; xsun(2)=', XSUN(3) , 'e0;'
      WRITE(*,'(A8F20.7A12F20.7A12F20.7A3)') 'xmon(0)=', XMON(1), 
     $ 'e0; xmon(1)=', XMON(2), 'e0; xmon(2)=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A9F20.10A13F20.10A13F20.10A3)') 'xtide(0)=',
     $ DXTIDE(1), 'e0; xtide(1)=', DXTIDE(2), 'e0; xtide(2)=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A46)') 'res=dehanttideinel_from_utc(xsta,xsun,xmon,t);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(0)-xtide(0))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(1)-xtide(1))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(2)-xtide(2))<PRECISION);'
      WRITE(*,'(A4)') '}'

C     gradually increase fractional hours, 1hour every step
      DO I=2,23,3
      FHR = I
      EX = EX + 1
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A17I3A15A30)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > PRECISION") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A8F15.4A12F15.4A12F15.4A3)') 'xsta(0)=', XSTA(1), 
     $ 'e0; xsta(1)=', XSTA(2), 'e0; xsta(2)=', XSTA(3), 'e0;'
      WRITE(*,'(A8F20.4A12F20.4A12F20.4A3)') 'xsun(0)=', XSUN(1), 
     $ 'e0; xsun(1)=', XSUN(2), 'e0; xsun(2)=', XSUN(3) , 'e0;'
      WRITE(*,'(A8F20.7A12F20.7A12F20.7A3)') 'xmon(0)=', XMON(1), 
     $ 'e0; xmon(1)=', XMON(2), 'e0; xmon(2)=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A9F20.10A13F20.10A13F20.10A3)') 'xtide(0)=',
     $ DXTIDE(1), 'e0; xtide(1)=', DXTIDE(2), 'e0; xtide(2)=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A46)') 'res=dehanttideinel_from_utc(xsta,xsun,xmon,t);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(0)-xtide(0))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(1)-xtide(1))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(2)-xtide(2))<PRECISION);'
      WRITE(*,'(A4)') '}'
      END DO
      
C     TEST CASE 3; units [m]
      XSTA(1) =  1112200.5696D0 
      XSTA(2) = -4842957.8511D0 
      XSTA(3) =  3985345.9122D0 
      XSUN(1) =  100210282451.6279D0 
      XSUN(2) =  103055630398.3160D0 
      XSUN(3) =  56855096480.4475D0 
      XMON(1) =  369817604.4348D0 
      XMON(2) =  1897917.5258D0 
      XMON(3) =  120804980.8284D0 
      YR      = 2015
      MONTH   = 7
      DAY     = 15
      FHR     = 0.00D0
      EX = EX + 1
      WRITE(*,'(A17I3A15A30)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > PRECISION") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A8F15.4A12F15.4A12F15.4A3)') 'xsta(0)=', XSTA(1), 
     $ 'e0; xsta(1)=', XSTA(2), 'e0; xsta(2)=', XSTA(3), 'e0;'
      WRITE(*,'(A8F20.4A12F20.4A12F20.4A3)') 'xsun(0)=', XSUN(1), 
     $ 'e0; xsun(1)=', XSUN(2), 'e0; xsun(2)=', XSUN(3) , 'e0;'
      WRITE(*,'(A8F20.7A12F20.7A12F20.7A3)') 'xmon(0)=', XMON(1), 
     $ 'e0; xmon(1)=', XMON(2), 'e0; xmon(2)=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A9F20.10A13F20.10A13F20.10A3)') 'xtide(0)=',
     $ DXTIDE(1), 'e0; xtide(1)=', DXTIDE(2), 'e0; xtide(2)=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A46)') 'res=dehanttideinel_from_utc(xsta,xsun,xmon,t);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(0)-xtide(0))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(1)-xtide(1))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(2)-xtide(2))<PRECISION);'
      WRITE(*,'(A4)') '}'

C     DO I=3,23,3
C     FHR = I + 0.05 * I
C     IF (FHR.GT.24) FHR = 23.99D0
      DO I=0,86399,3600
      FHR = I / 86400D0
      EX = EX + 1
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A17I3A15A30)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > PRECISION") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A8F15.4A12F15.4A12F15.4A3)') 'xsta(0)=', XSTA(1), 
     $ 'e0; xsta(1)=', XSTA(2), 'e0; xsta(2)=', XSTA(3), 'e0;'
      WRITE(*,'(A8F20.4A12F20.4A12F20.4A3)') 'xsun(0)=', XSUN(1), 
     $ 'e0; xsun(1)=', XSUN(2), 'e0; xsun(2)=', XSUN(3) , 'e0;'
      WRITE(*,'(A8F20.7A12F20.7A12F20.7A3)') 'xmon(0)=', XMON(1), 
     $ 'e0; xmon(1)=', XMON(2), 'e0; xmon(2)=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A9F20.10A13F20.10A13F20.10A3)') 'xtide(0)=',
     $ DXTIDE(1), 'e0; xtide(1)=', DXTIDE(2), 'e0; xtide(2)=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A46)') 'res=dehanttideinel_from_utc(xsta,xsun,xmon,t);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(0)-xtide(0))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(1)-xtide(1))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(2)-xtide(2))<PRECISION);'
      WRITE(*,'(A4)') '}'
      END DO

C     TEST CASE 4; units [m]
      XSTA(1) = 1112152.8166D0 
      XSTA(2) = -4842857.5435D0 
      XSTA(3) = 3985496.1783D0 
      XSUN(1) = 8382471154.1312895D0 
      XSUN(2) = 10512408445.356153D0 
      XSUN(3) = -5360583240.3763866D0 
      XMON(1) = 380934092.93550891D0 
      XMON(2) = 2871428.1904491195D0 
      XMON(3) = 79015680.553570181D0 
      YR      = 2017
      MONTH   = 1
      DAY     = 15
      FHR     = 0.00D0
      EX = EX + 1
      WRITE(*,'(A17I3A15A30)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > PRECISION") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A8F15.4A12F15.4A12F15.4A3)') 'xsta(0)=', XSTA(1), 
     $ 'e0; xsta(1)=', XSTA(2), 'e0; xsta(2)=', XSTA(3), 'e0;'
      WRITE(*,'(A8F20.4A12F20.4A12F20.4A3)') 'xsun(0)=', XSUN(1), 
     $ 'e0; xsun(1)=', XSUN(2), 'e0; xsun(2)=', XSUN(3) , 'e0;'
      WRITE(*,'(A8F20.7A12F20.7A12F20.7A3)') 'xmon(0)=', XMON(1), 
     $ 'e0; xmon(1)=', XMON(2), 'e0; xmon(2)=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A9F20.10A13F20.10A13F20.10A3)') 'xtide(0)=',
     $ DXTIDE(1), 'e0; xtide(1)=', DXTIDE(2), 'e0; xtide(2)=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A46)') 'res=dehanttideinel_from_utc(xsta,xsun,xmon,t);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(0)-xtide(0))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(1)-xtide(1))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(2)-xtide(2))<PRECISION);'
      WRITE(*,'(A4)') '}'

C     DO I=4,23,3
C     FHR = I + 0.1*I
C     IF (FHR.GT.24) FHR = 23.99D0
      DO I=0,86399,3600
      FHR = I / 86400D0
      EX = EX + 1
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A17I3A15A30)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > PRECISION") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A8F15.4A12F15.4A12F15.4A3)') 'xsta(0)=', XSTA(1), 
     $ 'e0; xsta(1)=', XSTA(2), 'e0; xsta(2)=', XSTA(3), 'e0;'
      WRITE(*,'(A8F20.4A12F20.4A12F20.4A3)') 'xsun(0)=', XSUN(1), 
     $ 'e0; xsun(1)=', XSUN(2), 'e0; xsun(2)=', XSUN(3) , 'e0;'
      WRITE(*,'(A8F20.7A12F20.7A12F20.7A3)') 'xmon(0)=', XMON(1), 
     $ 'e0; xmon(1)=', XMON(2), 'e0; xmon(2)=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A9F20.10A13F20.10A13F20.10A3)') 'xtide(0)=',
     $ DXTIDE(1), 'e0; xtide(1)=', DXTIDE(2), 'e0; xtide(2)=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A46)') 'res=dehanttideinel_from_utc(xsta,xsun,xmon,t);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(0)-xtide(0))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(1)-xtide(1))<PRECISION);'
      WRITE(*,'(A50)') 'REQUIRE(std::abs(res(2)-xtide(2))<PRECISION);'
      WRITE(*,'(A4)') '}'
      END DO
      END
