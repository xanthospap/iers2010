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
      WRITE(*,'(A17I3A15A25)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > 1e-9") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A9F15.4A13F15.4A13F15.4A3)') 'xsta.x()=', XSTA(1), 
     $ 'e0; xsta.y()=', XSTA(2), 'e0; xsta.z()=', XSTA(3), 'e0;'
      WRITE(*,'(A9F20.4A13F20.4A13F20.4A3)') 'xsun.x()=', XSUN(1), 
     $ 'e0; xsun.y()=', XSUN(2), 'e0; xsun.z()=', XSUN(3) , 'e0;'
      WRITE(*,'(A9F20.7A13F20.7A13F20.7A3)') 'xmon.x()=', XMON(1), 
     $ 'e0; xmon.y()=', XMON(2), 'e0; xmon.z()=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A10F20.10A14F20.10A14F20.10A3)') 'xtide.x()=',
     $ DXTIDE(1), 'e0; xtide.y()=', DXTIDE(2), 'e0; xtide.z()=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A37)') 'res=dehanttideinel(xsta,xsun,xmon,t);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.x()-xtide.x()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.y()-xtide.y()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.z()-xtide.z()) < 1e-9);'
      WRITE(*,'(A4)') '}'

      DO I=1,23,3
      FHR = I
      EX = EX + 1
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A17I3A15A25)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > 1e-9") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A9F15.4A13F15.4A13F15.4A3)') 'xsta.x()=', XSTA(1), 
     $ 'e0; xsta.y()=', XSTA(2), 'e0; xsta.z()=', XSTA(3), 'e0;'
      WRITE(*,'(A9F20.4A13F20.4A13F20.4A3)') 'xsun.x()=', XSUN(1), 
     $ 'e0; xsun.y()=', XSUN(2), 'e0; xsun.z()=', XSUN(3) , 'e0;'
      WRITE(*,'(A9F20.7A13F20.7A13F20.7A3)') 'xmon.x()=', XMON(1), 
     $ 'e0; xmon.y()=', XMON(2), 'e0; xmon.z()=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A10F20.10A14F20.10A14F20.10A3)') 'xtide.x()=',
     $ DXTIDE(1), 'e0; xtide.y()=', DXTIDE(2), 'e0; xtide.z()=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A37)') 'res=dehanttideinel(xsta,xsun,xmon,t);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.x()-xtide.x()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.y()-xtide.y()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.z()-xtide.z()) < 1e-9);'
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
      WRITE(*,'(A17I3A15A25)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > 1e-9") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A9F15.4A13F15.4A13F15.4A3)') 'xsta.x()=', XSTA(1), 
     $ 'e0; xsta.y()=', XSTA(2), 'e0; xsta.z()=', XSTA(3), 'e0;'
      WRITE(*,'(A9F20.4A13F20.4A13F20.4A3)') 'xsun.x()=', XSUN(1), 
     $ 'e0; xsun.y()=', XSUN(2), 'e0; xsun.z()=', XSUN(3) , 'e0;'
      WRITE(*,'(A9F20.7A13F20.7A13F20.7A3)') 'xmon.x()=', XMON(1), 
     $ 'e0; xmon.y()=', XMON(2), 'e0; xmon.z()=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A10F20.10A14F20.10A14F20.10A3)') 'xtide.x()=',
     $ DXTIDE(1), 'e0; xtide.y()=', DXTIDE(2), 'e0; xtide.z()=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A37)') 'res=dehanttideinel(xsta,xsun,xmon,t);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.x()-xtide.x()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.y()-xtide.y()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.z()-xtide.z()) < 1e-9);'
      WRITE(*,'(A4)') '}'

      DO I=2,23,3
      FHR = I
      EX = EX + 1
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A17I3A15A25)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > 1e-9") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A9F15.4A13F15.4A13F15.4A3)') 'xsta.x()=', XSTA(1), 
     $ 'e0; xsta.y()=', XSTA(2), 'e0; xsta.z()=', XSTA(3), 'e0;'
      WRITE(*,'(A9F20.4A13F20.4A13F20.4A3)') 'xsun.x()=', XSUN(1), 
     $ 'e0; xsun.y()=', XSUN(2), 'e0; xsun.z()=', XSUN(3) , 'e0;'
      WRITE(*,'(A9F20.7A13F20.7A13F20.7A3)') 'xmon.x()=', XMON(1), 
     $ 'e0; xmon.y()=', XMON(2), 'e0; xmon.z()=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A10F20.10A14F20.10A14F20.10A3)') 'xtide.x()=',
     $ DXTIDE(1), 'e0; xtide.y()=', DXTIDE(2), 'e0; xtide.z()=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A37)') 'res=dehanttideinel(xsta,xsun,xmon,t);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.x()-xtide.x()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.y()-xtide.y()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.z()-xtide.z()) < 1e-9);'
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
      WRITE(*,'(A17I3A15A25)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > 1e-9") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A9F15.4A13F15.4A13F15.4A3)') 'xsta.x()=', XSTA(1), 
     $ 'e0; xsta.y()=', XSTA(2), 'e0; xsta.z()=', XSTA(3), 'e0;'
      WRITE(*,'(A9F20.4A13F20.4A13F20.4A3)') 'xsun.x()=', XSUN(1), 
     $ 'e0; xsun.y()=', XSUN(2), 'e0; xsun.z()=', XSUN(3) , 'e0;'
      WRITE(*,'(A9F20.7A13F20.7A13F20.7A3)') 'xmon.x()=', XMON(1), 
     $ 'e0; xmon.y()=', XMON(2), 'e0; xmon.z()=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A10F20.10A14F20.10A14F20.10A3)') 'xtide.x()=',
     $ DXTIDE(1), 'e0; xtide.y()=', DXTIDE(2), 'e0; xtide.z()=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A37)') 'res=dehanttideinel(xsta,xsun,xmon,t);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.x()-xtide.x()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.y()-xtide.y()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.z()-xtide.z()) < 1e-9);'
      WRITE(*,'(A4)') '}'

C     DO I=3,23,3
C     FHR = I + 0.05 * I
C     IF (FHR.GT.24) FHR = 23.99D0
      DO I=0,86399,3600
      FHR = I / 86400D0
      EX = EX + 1
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A17I3A15A25)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > 1e-9") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A9F15.4A13F15.4A13F15.4A3)') 'xsta.x()=', XSTA(1), 
     $ 'e0; xsta.y()=', XSTA(2), 'e0; xsta.z()=', XSTA(3), 'e0;'
      WRITE(*,'(A9F20.4A13F20.4A13F20.4A3)') 'xsun.x()=', XSUN(1), 
     $ 'e0; xsun.y()=', XSUN(2), 'e0; xsun.z()=', XSUN(3) , 'e0;'
      WRITE(*,'(A9F20.7A13F20.7A13F20.7A3)') 'xmon.x()=', XMON(1), 
     $ 'e0; xmon.y()=', XMON(2), 'e0; xmon.z()=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A10F20.10A14F20.10A14F20.10A3)') 'xtide.x()=',
     $ DXTIDE(1), 'e0; xtide.y()=', DXTIDE(2), 'e0; xtide.z()=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A37)') 'res=dehanttideinel(xsta,xsun,xmon,t);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.x()-xtide.x()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.y()-xtide.y()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.z()-xtide.z()) < 1e-9);'
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
      WRITE(*,'(A17I3A15A25)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > 1e-9") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A9F15.4A13F15.4A13F15.4A3)') 'xsta.x()=', XSTA(1), 
     $ 'e0; xsta.y()=', XSTA(2), 'e0; xsta.z()=', XSTA(3), 'e0;'
      WRITE(*,'(A9F20.4A13F20.4A13F20.4A3)') 'xsun.x()=', XSUN(1), 
     $ 'e0; xsun.y()=', XSUN(2), 'e0; xsun.z()=', XSUN(3) , 'e0;'
      WRITE(*,'(A9F20.7A13F20.7A13F20.7A3)') 'xmon.x()=', XMON(1), 
     $ 'e0; xmon.y()=', XMON(2), 'e0; xmon.z()=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A10F20.10A14F20.10A14F20.10A3)') 'xtide.x()=',
     $ DXTIDE(1), 'e0; xtide.y()=', DXTIDE(2), 'e0; xtide.z()=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A37)') 'res=dehanttideinel(xsta,xsun,xmon,t);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.x()-xtide.x()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.y()-xtide.y()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.z()-xtide.z()) < 1e-9);'
      WRITE(*,'(A4)') '}'

C     DO I=4,23,3
C     FHR = I + 0.1*I
C     IF (FHR.GT.24) FHR = 23.99D0
      DO I=0,86399,3600
      FHR = I / 86400D0
      EX = EX + 1
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A17I3A15A25)') 'SUBCASE("Example ',EX,': Checking for ',
     $ 'descripancies  > 1e-9") {'
      CALL DEHANTTIDEINEL(XSTA,YR,MONTH,DAY,FHR,XSUN,XMON,DXTIDE)
      WRITE(*,'(A9F15.4A13F15.4A13F15.4A3)') 'xsta.x()=', XSTA(1), 
     $ 'e0; xsta.y()=', XSTA(2), 'e0; xsta.z()=', XSTA(3), 'e0;'
      WRITE(*,'(A9F20.4A13F20.4A13F20.4A3)') 'xsun.x()=', XSUN(1), 
     $ 'e0; xsun.y()=', XSUN(2), 'e0; xsun.z()=', XSUN(3) , 'e0;'
      WRITE(*,'(A9F20.7A13F20.7A13F20.7A3)') 'xmon.x()=', XMON(1), 
     $ 'e0; xmon.y()=', XMON(2), 'e0; xmon.z()=', XMON(3), 'e0;'
      WRITE(*,'(A25I4A9I2A17I2A11F10.0A3)') 't=datetime<seconds>(year(',
     $ YR, '), month(', MONTH,'), day_of_month(', DAY,'), seconds(', 
     $ FHR*3600D0, '));'
      WRITE(*,'(A10F20.10A14F20.10A14F20.10A3)') 'xtide.x()=',
     $ DXTIDE(1), 'e0; xtide.y()=', DXTIDE(2), 'e0; xtide.z()=', 
     $ DXTIDE(3), 'e0;'
      WRITE(*,'(A37)') 'res=dehanttideinel(xsta,xsun,xmon,t);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.x()-xtide.x()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.y()-xtide.y()) < 1e-9);'
      WRITE(*,'(A45)') 'REQUIRE(std::abs(res.z()-xtide.z()) < 1e-9);'
      WRITE(*,'(A4)') '}'
      END DO
      END
