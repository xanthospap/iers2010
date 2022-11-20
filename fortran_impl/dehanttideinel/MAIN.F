      PROGRAM MAIN

      DOUBLE PRECISION PRN, FHR
      INTEGER NUM_TESTS, I,J, YEAR, MONTH, DAY, SECDAY
      DOUBLE PRECISION XSTA1(3), XSTA2(3), XSTA3(3), XSTA4(3)
      DOUBLE PRECISION XSTA5(3), XSTA6(3), XSTA7(3), XSTA8(3), XSTA9(3)
      DOUBLE PRECISION XSTA0(3), XSUN(3), XMON(3)
      DOUBLE PRECISION DXTIDE(3)

      NUM_TESTS = 10

C Coordinates of sites, use IDS log files for 
C DORIS beacons
C NY-ALESUND II
      XSTA0(1) =  1201300.046D0
      XSTA0(2) =  251874.435D0
      XSTA0(3) =  6238000.324D0
C HOFN                                              
      XSTA1(1) =  2679718.473D0
      XSTA1(2)= -727978.273D0
      XSTA1(3)= 5722772.163D0
C KAUAI                                             
      XSTA2(1) = -5543981.353D0
      XSTA2(2)=-2054584.146D0
      XSTA2(3)= 2387480.974D0
C FUTUNA                                            
      XSTA3(1) = -6178323.900D0
      XSTA3(2)= -202689.617D0
      XSTA3(3)=-1566023.046D0
C BELGRANO                                          
      XSTA4(1) =  1106046.627D0
      XSTA4(2)= -763739.010D0
      XSTA4(3)=-6214243.195D0
C TERRE ADELIE                                      
      XSTA5(1) = -1940878.515D0
      XSTA5(2)= 1628473.041D0
      XSTA5(3)=-5833723.413D0
C NOUMEA                                            
      XSTA6(1) = -5743533.479D0
      XSTA6(2)= 1380501.922D0
      XSTA6(3)=-2397907.623D0
C ASCENSION                                         
      XSTA7(1) =  6121154.081D0
      XSTA7(2)=-1563976.723D0
      XSTA7(3)= -872606.019D0
C MAHE                                              
      XSTA8(1) =  3597847.792D0
      XSTA8(2)= 5240876.363D0
      XSTA8(3)= -516790.797D0
C DIONYSOS                                          
      XSTA9(1) =  4595212.468D0
      XSTA9(2)= 2039473.691D0
      XSTA9(3)= 3912617.891D0

C Site coordinates are constant; write them down as a C++ vector
      WRITE(*,'(A50)') 'std::vector<Eigen::Matrix<double,3,1>> xstavec;'
      WRITE(*,'(A45,F15.4,A1,F15.4,A1,F15.4,A3)') 
     . 'xstavec.push_back(Eigen::Matrix<double,3,1>(',
     . XSTA0(1), ',', XSTA0(2), ',', XSTA0(3), '));'
      WRITE(*,'(A45,F15.4,A1,F15.4,A1,F15.4,A3)') 
     . 'xstavec.push_back(Eigen::Matrix<double,3,1>(',
     . XSTA1(1), ',', XSTA1(2), ',', XSTA1(3), '));'
      WRITE(*,'(A45,F15.4,A1,F15.4,A1,F15.4,A3)') 
     . 'xstavec.push_back(Eigen::Matrix<double,3,1>(',
     . XSTA2(1), ',', XSTA2(2), ',', XSTA2(3), '));'
      WRITE(*,'(A45,F15.4,A1,F15.4,A1,F15.4,A3)') 
     . 'xstavec.push_back(Eigen::Matrix<double,3,1>(',
     . XSTA3(1), ',', XSTA3(2), ',', XSTA3(3), '));'
      WRITE(*,'(A45,F15.4,A1,F15.4,A1,F15.4,A3)') 
     . 'xstavec.push_back(Eigen::Matrix<double,3,1>(',
     . XSTA4(1), ',', XSTA4(2), ',', XSTA4(3), '));'
      WRITE(*,'(A45,F15.4,A1,F15.4,A1,F15.4,A3)') 
     . 'xstavec.push_back(Eigen::Matrix<double,3,1>(',
     . XSTA5(1), ',', XSTA5(2), ',', XSTA5(3), '));'
      WRITE(*,'(A45,F15.4,A1,F15.4,A1,F15.4,A3)') 
     . 'xstavec.push_back(Eigen::Matrix<double,3,1>(',
     . XSTA6(1), ',', XSTA6(2), ',', XSTA6(3), '));'
      WRITE(*,'(A45,F15.4,A1,F15.4,A1,F15.4,A3)') 
     . 'xstavec.push_back(Eigen::Matrix<double,3,1>(',
     . XSTA7(1), ',', XSTA7(2), ',', XSTA7(3), '));'
      WRITE(*,'(A45,F15.4,A1,F15.4,A1,F15.4,A3)') 
     . 'xstavec.push_back(Eigen::Matrix<double,3,1>(',
     . XSTA8(1), ',', XSTA8(2), ',', XSTA8(3), '));'
      WRITE(*,'(A45,F15.4,A1,F15.4,A1,F15.4,A3)') 
     . 'xstavec.push_back(Eigen::Matrix<double,3,1>(',
     . XSTA9(1), ',', XSTA9(2), ',', XSTA9(3), '));'

C We need (ECEF) coordinates for the Sun and Moon. Have some reference ones,
C and we'll add pseudo-random noise to them, in the order of +- 100km
      XSUN(1) = 137859926952.015D0
      XSUN(2) = 54228127881.4350D0 
      XSUN(3) = 23509422341.6960D0
      XMON(1) = -179996231.920342D0
      XMON(2) = -312468450.131567D0
      XMON(3) = -169288918.592160D0

C Iterate for the number of test cases 
      DO I=1,NUM_TESTS
      WRITE(*, '(A30,I3)') '// Test case #',I

C Create a random date
      call random_number(RN)
      YEAR  = 1995 + FLOOR(40*RN)
      call random_number(RN)
      MONTH = 1 + FLOOR(12*RN)
      call random_number(RN)
      DAY = 1 + FLOOR(30*RN)
      call random_number(RN)
      SECDAY = FLOOR(86401*RN)
C Write as a C++ datetime
      WRITE(*,'(A35)') 't=dso::datetime<dso::nanoseconds>('
      WRITE(*,'(A10,I4,A15,I2,A25,I2,A20,I6,A20)') 'dso::year(',
     . YEAR,'), dso::month(', MONTH, '), dso::day_of_month(', DAY, 
     . '),dso::nanoseconds(',SECDAY,'*1000000000L));'

C add some random km [0-1e3] to the Sun and Moon positions
      call random_number(RN)
      XSUN(1) = XSUN(1) + (-0.5+RN) * 1D5
      call random_number(RN)
      XSUN(2) = XSUN(2) + (-0.5+RN) * 1D5
      call random_number(RN)
      XSUN(3) = XSUN(3) + (-0.5+RN) * 1D5
      call random_number(RN)
      XMON(1) = XMON(1) + (-0.5+RN) * 1D5
      call random_number(RN)
      XMON(2) = XMON(2) + (-0.5+RN) * 1D5
      call random_number(RN)
      XMON(3) = XMON(3) + (-0.5+RN) * 1D5
      WRITE(*, '(A35,E30.22,A1,E30.22,A1,E30.22,A3)') 
     . 'xsun = Eigen::Matrix<double,3,1>(',XSUN(1), ',', XSUN(2), 
     . ',', XSUN(3), ');'
      WRITE(*, '(A35,E30.22,A1,E30.22,A1,E30.22,A3)') 
     . 'xmon = Eigen::Matrix<double,3,1>(',XMON(1), ',', XMON(2), 
     . ',', XMON(3), ');'
      
C Call the C++ function
      WRITE(*, '(A55)') 
     . 'iers2010::dehanttideinel(t,xsun,xmon,xstavec,xcorvec);'

C Call the Fortran function, for each station, and compare results
      CALL DEHANTTIDEINEL(XSTA0,YEAR,MONTH,DAY,SECDAY/3600D0,
     . XSUN,XMON,DXTIDE)
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[0](0)-(',
     . DXTIDE(1),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[0](1)-(',
     . DXTIDE(2),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[0](2)-(',
     . DXTIDE(3),'))<PRECISION);'
      CALL DEHANTTIDEINEL(XSTA1,YEAR,MONTH,DAY,SECDAY/3600D0,
     . XSUN,XMON,DXTIDE)
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[1](0)-(',
     . DXTIDE(1),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[1](1)-(',
     . DXTIDE(2),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[1](2)-(',
     . DXTIDE(3),'))<PRECISION);'
      CALL DEHANTTIDEINEL(XSTA2,YEAR,MONTH,DAY,SECDAY/3600D0,
     . XSUN,XMON,DXTIDE)
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[2](0)-(',
     . DXTIDE(1),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[2](1)-(',
     . DXTIDE(2),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[2](2)-(',
     . DXTIDE(3),'))<PRECISION);'
      CALL DEHANTTIDEINEL(XSTA3,YEAR,MONTH,DAY,SECDAY/3600D0,
     . XSUN,XMON,DXTIDE)
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[3](0)-(',
     . DXTIDE(1),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[3](1)-(',
     . DXTIDE(2),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[3](2)-(',
     . DXTIDE(3),'))<PRECISION);'
      CALL DEHANTTIDEINEL(XSTA4,YEAR,MONTH,DAY,SECDAY/3600D0,
     . XSUN,XMON,DXTIDE)
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[4](0)-(',
     . DXTIDE(1),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[4](1)-(',
     . DXTIDE(2),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[4](2)-(',
     . DXTIDE(3),'))<PRECISION);'
      CALL DEHANTTIDEINEL(XSTA5,YEAR,MONTH,DAY,SECDAY/3600D0,
     . XSUN,XMON,DXTIDE)
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[5](0)-(',
     . DXTIDE(1),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[5](1)-(',
     . DXTIDE(2),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[5](2)-(',
     . DXTIDE(3),'))<PRECISION);'
      CALL DEHANTTIDEINEL(XSTA6,YEAR,MONTH,DAY,SECDAY/3600D0,
     . XSUN,XMON,DXTIDE)
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[6](0)-(',
     . DXTIDE(1),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[6](1)-(',
     . DXTIDE(2),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[6](2)-(',
     . DXTIDE(3),'))<PRECISION);'
      CALL DEHANTTIDEINEL(XSTA7,YEAR,MONTH,DAY,SECDAY/3600D0,
     . XSUN,XMON,DXTIDE)
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[7](0)-(',
     . DXTIDE(1),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[7](1)-(',
     . DXTIDE(2),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[7](2)-(',
     . DXTIDE(3),'))<PRECISION);'
      CALL DEHANTTIDEINEL(XSTA8,YEAR,MONTH,DAY,SECDAY/3600D0,
     . XSUN,XMON,DXTIDE)
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[8](0)-(',
     . DXTIDE(1),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[8](1)-(',
     . DXTIDE(2),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[8](2)-(',
     . DXTIDE(3),'))<PRECISION);'
      CALL DEHANTTIDEINEL(XSTA9,YEAR,MONTH,DAY,SECDAY/3600D0,
     . XSUN,XMON,DXTIDE)
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[9](0)-(',
     . DXTIDE(1),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[9](1)-(',
     . DXTIDE(2),'))<PRECISION);'
      WRITE(*,'(A32,F15.12,A20)') 'assert(std::abs(xcorvec[9](2)-(',
     . DXTIDE(3),'))<PRECISION);'
      
      WRITE(*, '(A20,I3,A15)') 'printf("Test case #',I-1, 
     . ' passed\n");'
      END DO

      END
