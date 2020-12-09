      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for DEHANT (test 1)
C     ----------------------------------------------
      DOUBLE PRECISION XSTA1(3), XSUN1(3), XMON1(3)
      DATA (XSTA1(I),I=1,3) /4075578.385D0,
     . 931852.890D0, 4801570.154D0/
      DATA (XSUN1(I),I=1,3) /137859926952.015D0,
     . 54228127881.4350D0, 23509422341.6960D0/
      DATA (XMON1(I),I=1,3) /-179996231.920342D0,
     . -312468450.131567D0, -169288918.592160D0/

C     ----------------------------------------------
C     Variables and data for DEHANT (test 2)
C     ----------------------------------------------
      DOUBLE PRECISION XSTA2(3), XSUN2(3), XMON2(3)
      DATA (XSTA2(I),I=1,3) /1112189.660D0,
     . -4842955.026D0, 3985352.284D0/
      DATA (XSUN2(I),I=1,3) /-54537460436.2357D0,
     . 130244288385.279D0, 56463429031.5996D0/
      DATA (XMON2(I),I=1,3) /300396716.912D0,
     . 243238281.451D0, 120548075.939D0/

C     ----------------------------------------------
C     Variables and data for DEHANT (test 3)
C     ----------------------------------------------
      DOUBLE PRECISION XSTA3(3), XSUN3(3), XMON3(3)
      DATA (XSTA3(I),I=1,3) /1112200.5696D0,
     . -4842957.8511D0, 3985345.9122D0/
      DATA (XSUN3(I),I=1,3) /100210282451.6279D0,
     . 103055630398.3160D0, 56855096480.4475D0/
      DATA (XMON3(I),I=1,3) /369817604.4348D0,
     . 1897917.5258D0, 120804980.8284D0/

C     ----------------------------------------------
C     Variables and data for DEHANT (test 4)
C     ----------------------------------------------
      DOUBLE PRECISION XSTA4(3), XSUN4(3), XMON4(3),
     . DXTIDE(3)
      DATA (XSTA4(I),I=1,3) /1112152.8166D0,
     . -4842857.5435D0, 3985496.1783D0/
      DATA (XSUN4(I),I=1,3) /8382471154.1312895D0,
     . 10512408445.356153D0, -5360583240.3763866D0/
      DATA (XMON4(I),I=1,3) /380934092.93550891D0,
     . 2871428.1904491195D0, 79015680.553570181D0/

C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> DEHANTTIDEINEL Results'
      PRINT *, '----------------------------------------'

      PRINT *, '// Test Case A'
      CALL DEHANTTIDEINEL(XSTA1,2009,4,13,0D0,XSUN1,XMON1,DXTIDE)
      WRITE(*, FMT1) DABS(DXTIDE(1)-0.7700420357108125891D-01)
      WRITE(*, FMT1) DABS(DXTIDE(2)-0.6304056321824967613D-01)
      WRITE(*, FMT1) DABS(DXTIDE(3)-0.5516568152597246810D-01)

      PRINT *, '// Test Case B'
      CALL DEHANTTIDEINEL(XSTA2,2012,7,13,0D0,XSUN2,XMON2,DXTIDE)
      WRITE(*, FMT1) DABS(DXTIDE(1)+0.2036831479592075833D-01)
      WRITE(*, FMT1) DABS(DXTIDE(2)-0.5658254776225972449D-01)
      WRITE(*, FMT1) DABS(DXTIDE(3)+0.7597679676871742227D-01)
      
      PRINT *, '// Test Case C'
      CALL DEHANTTIDEINEL(XSTA3,2015,7,15,0D0,XSUN3,XMON3,DXTIDE)
      WRITE(*, FMT1) DABS(DXTIDE(1)-0.00509570869172363845D0)
      WRITE(*, FMT1) DABS(DXTIDE(2)-0.0828663025983528700D0)
      WRITE(*, FMT1) DABS(DXTIDE(3)+0.0636634925404189617D0)
     
      PRINT *, '// Test Case D'
      CALL DEHANTTIDEINEL(XSTA4,2017,1,15,0D0,XSUN4,XMON4,DXTIDE)
      WRITE(*, FMT1) DABS(DXTIDE(1)-.0050957086917236384D0)
      WRITE(*, FMT1) DABS(DXTIDE(2)-.082866302598352870D0)
      WRITE(*, FMT1) DABS(DXTIDE(3)+.063663492540418962D0)


C     ALL DONE
      END
