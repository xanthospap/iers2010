      PROGRAM MAIN

C     RESULT REPORTING
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT2 = "(F20.15)"

C     ----------------------------------------------
C     Variables and data for FUNDARG
C     ----------------------------------------------
      DOUBLE PRECISION T, LR, LPR, FR, DR, OMR
      DOUBLE PRECISION    L,  LP,  F,  D,  OM
      DATA T, LR, LPR, FR, DR, OMR/0.07995893223819302D0,
     .   2.291187512612069099D0,
     .   6.212931111003726414D0,
     .   3.658025792050572989D0,
     .   4.554139562402433228D0,
     .   -0.5167379217231804489D0/

C     ----------------------------------------------
C     Variables and data for PMSDNUT2
C     ----------------------------------------------
      DOUBLE PRECISION PM(2), PMR(2)
      DATA (PMR(I),I=1,2) /24.83144238273364834D0,
     . -14.09240692041837661D0/

C     ----------------------------------------------
C     Variables and data for UTLIBR
C     ----------------------------------------------
      DOUBLE PRECISION  dUT1,  dLOD
      DOUBLE PRECISION dUT1AR, dUT1BR, dLODAR, dLODBR
      DATA dUT1AR, dUT1BR, dLODAR, dLODBR/
     . 2.441143834386761746D0,
     . -2.655705844335680244D0,
     . -14.78971247349449492D0,
     . 27.39445826599846967D0/

C     ----------------------------------------------
C     Variables and data for FCNNUT
C     ----------------------------------------------
      DOUBLE PRECISION  FCX,  FCY,  FCDX,  FCDY
      DOUBLE PRECISION FCXR, FCYR, FCDXR, FCDYR
      DATA FCXR, FCYR, FCDXR, FCDYR/
     . -176.8012290066270680D0,
     . -93.51855308903756736D0,
     . 3.745573770491803067D0,
     . 3.745573770491803067D0/

C     ----------------------------------------------
C     Variables and data for RG_ZONT2
C     ----------------------------------------------
      DOUBLE PRECISION ZDUTR, ZDLODR, ZDOMEGAR
      DOUBLE PRECISION ZDUT,  ZDLOD,  ZDOMEGA
      DATA ZDUTR, ZDLODR, ZDOMEGAR/
     . 7.983287678576557467E-002,
     . 5.035331113978199288E-005,
     . -4.249711616463017E-014/

C     ----------------------------------------------
C     Variables and data for VMF1
C     ----------------------------------------------
      DOUBLE PRECISION VMF1HR, VMF1WR, VMF1H, VMF1W
      DATA VMF1HR, VMF1WR/ 3.424342122738070593D0,
     . 3.448299714692572238D0/

C     ----------------------------------------------
C     Variables and data for GPT2 
C     ----------------------------------------------
      DOUBLE PRECISION GPT_P(1), GPT_T(1), GPT_DT(1), GPT_E(1),
     . GPT_AH(1), GPT_AW(1), GPT_UNDU(1)
      DOUBLE PRECISION GPT21LON(1),GPT21LAT(1),GPT21HEL(1)

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
C     Variables and data for FCUL_ZD_HPA
C     ----------------------------------------------
      DOUBLE PRECISION FCUL_ZTD, FCUL_ZHD, FCUL_ZWD

C     CALL FUNDARG
      CALL FUNDARG(T, L, LP, F, D, OM)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> FUNDARG Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(L-LR)
      WRITE(*, FMT1) DABS(LP-LPR)
      WRITE(*, FMT1) DABS(F-FR)
      WRITE(*, FMT1) DABS(D-DR)
      WRITE(*, FMT1) DABS(OM-OMR)

C     CALL PMSDNUT2
      T = 54335D0
      CALL PMSDNUT2(T, PM)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> PMSDNUT2 Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(PM(1)-PMR(1))
      WRITE(*, FMT1) DABS(PM(2)-PMR(2))

C     CALL UTLIBR
      T = 44239.1D0
      CALL UTLIBR(T, dUT1, dLOD)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> UTLIBR Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(dUT1-dUT1AR)
      WRITE(*, FMT1) DABS(dLOD-dLODAR)
      T = 55227.4D0
      CALL UTLIBR(T, dUT1, dLOD)
C     REPORT RESULTS
      WRITE(*, FMT1) DABS(dUT1-dUT1BR)
      WRITE(*, FMT1) DABS(dLOD-dLODBR)

C     CALL FCNNUT
      T = 54790D0
      CALL FCNNUT(T, FCX, FCY, FCDX, FCDY)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> FCNNUT Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(FCX-FCXR)
      WRITE(*, FMT1) DABS(FCY-FCYR)
      WRITE(*, FMT1) DABS(FCDX-FCDXR)
      WRITE(*, FMT1) DABS(FCDY-FCDYR)

C     CALL RG_ZONT2
      T = .07995893223819302D0
      CALL RG_ZONT2(T, ZDUT, ZDLOD, ZDOMEGA)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> RG_ZONT2 Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(ZDUT-ZDUTR)
      WRITE(*, FMT1) DABS(ZDLOD-ZDLODR)
      WRITE(*, FMT1) DABS(ZDOMEGA-ZDOMEGAR)

C     CALL VMF1
      CALL VMF1(0.00127683D0, 0.00060955D0, 55055D0,
     .  0.6708665767D0, 1.278564131D0, VMF1H, VMF1W)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> VMF1 Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(VMF1H-VMF1HR)
      WRITE(*, FMT1) DABS(VMF1W-VMF1WR)

C     CALL GPT2
      CALL GPT2(56141.d0,0.8412486994612668D0,0.28571039855147173D0,
     .  156.D0, 1, 0, GPT_P, GPT_T, GPT_DT, GPT_E, GPT_AH,
     .  GPT_AW, GPT_UNDU)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> GPT2 Results: (Test A)'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) GPT_P
      WRITE(*, FMT1) GPT_T
      WRITE(*, FMT1) GPT_DT
      WRITE(*, FMT1) GPT_E
      WRITE(*, FMT1) GPT_AH
      WRITE(*, FMT1) GPT_AW
      WRITE(*, FMT1) GPT_UNDU
      CALL GPT2(56141.d0,0.8412486994612668D0,0.28571039855147173D0,
     .  156.D0, 1, 1, GPT_P, GPT_T, GPT_DT, GPT_E, GPT_AH,
     .  GPT_AW, GPT_UNDU)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> GPT2 Results: (Test B)'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) GPT_P
      WRITE(*, FMT1) GPT_T
      WRITE(*, FMT1) GPT_DT
      WRITE(*, FMT1) GPT_E
      WRITE(*, FMT1) GPT_AH
      WRITE(*, FMT1) GPT_AW
      WRITE(*, FMT1) GPT_UNDU

C     CALL DEHANTTIDEINEL
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> DEHANTTIDEINEL Results: (Test A)'
      PRINT *, '----------------------------------------'
      CALL DEHANTTIDEINEL(XSTA1,2009,4,13,0D0,XSUN1,XMON1,DXTIDE)
      WRITE(*, FMT2) DABS(DXTIDE(1)-0.7700420357108125891D-01)
      WRITE(*, FMT2) DABS(DXTIDE(2)-0.6304056321824967613D-01)
      WRITE(*, FMT2) DABS(DXTIDE(3)-0.5516568152597246810D-01)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> DEHANTTIDEINEL Results: (Test B)'
      PRINT *, '----------------------------------------'
      CALL DEHANTTIDEINEL(XSTA2,2012,7,13,0D0,XSUN2,XMON2,DXTIDE)
      WRITE(*, FMT2) DABS(DXTIDE(1)+0.2036831479592075833D-01)
      WRITE(*, FMT2) DABS(DXTIDE(2)-0.5658254776225972449D-01)
      WRITE(*, FMT2) DABS(DXTIDE(3)+0.7597679676871742227D-01)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> DEHANTTIDEINEL Results: (Test C)'
      PRINT *, '----------------------------------------'
      CALL DEHANTTIDEINEL(XSTA3,2015,7,15,0D0,XSUN3,XMON3,DXTIDE)
      WRITE(*, FMT2) DABS(DXTIDE(1)-0.00509570869172363845D0)
      WRITE(*, FMT2) DABS(DXTIDE(2)-0.0828663025983528700D0)
      WRITE(*, FMT2) DABS(DXTIDE(3)+0.0636634925404189617D0)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> DEHANTTIDEINEL Results: (Test D)'
      PRINT *, '----------------------------------------'
      CALL DEHANTTIDEINEL(XSTA4,2017,1,15,0D0,XSUN4,XMON4,DXTIDE)
      WRITE(*, FMT2) DABS(DXTIDE(1)-.0050957086917236384D0)
      WRITE(*, FMT2) DABS(DXTIDE(2)-.082866302598352870D0)
      WRITE(*, FMT2) DABS(DXTIDE(3)+.063663492540418962D0)

C     CALL FCUL_ZD_HPA
      CALL FCULZD_HPA(30.67166667D0, 2010.344D0, 798.4188D0,
     . 14.322D0, 0.532D0, FCUL_ZTD, FCUL_ZHD, FCUL_ZWD)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> FCULZD_HPA  Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT2) DABS(FCUL_ZTD-1.935225924846803114D0)
      WRITE(*, FMT2) DABS(FCUL_ZHD-1.932992176591644462D0)
      WRITE(*, FMT2) DABS(FCUL_ZWD-0.2233748255158703871D-02)

C     CALL GPT2
C     GPT21LAT(1)=0.8412486994612668D0
C     GPT21LON(1)=0.28571039855147173e0
C     GPT21HEL(1)=156.D0
C     CALL GPT2(56141.D0,GPT2LAT,GPT2LON,GPT2HEL,1,0, 
C    . GPT_P,GPT_T,GPT_DT,GPT_E,GPT_AH,GPT_AW,GPT_UNDU)
C     REPORT RESULTS
C     PRINT *, '----------------------------------------'
C     PRINT *, '> GPT2 Results: (Test a)'
C     PRINT *, '----------------------------------------'
C     WRITE(*, FMT2) DABS(GPT_P(1)-1002.56D0)
C     WRITE(*, FMT2) DABS(GPT_T(1)-22.12D0)
C     WRITE(*, FMT2) DABS(GPT_DT(1)+6.53D0)
C     WRITE(*, FMT2) DABS(GPT_E(1)-15.63D0)
C     WRITE(*, FMT2) DABS(GPT_AH(1)-0.0012647D0)
C     WRITE(*, FMT2) DABS(GPT_AW(1)-0.0005726D0)
C     WRITE(*, FMT2) DABS(GPT_UNDU(1)-44.06D0)
C     GPT21LAT(1)=0.8412486994612668D0
C     GPT21LON(1)=0.28571039855147173e0
C     GPT21HEL(1)=156.D0
C     CALL GPT2(56141.D0,GPT2LAT,GPT2LON,GPT2HEL,1,1, 
C    . GPT_P,GPT_T,GPT_DT,GPT_E,GPT_AH,GPT_AW,GPT_UNDU)
C     REPORT RESULTS
C     PRINT *, '----------------------------------------'
C     PRINT *, '> GPT2 Results: (Test b)'
C     PRINT *, '----------------------------------------'
C     WRITE(*, FMT2) DABS(GPT_P(1)-1003.49D0)
C     WRITE(*, FMT2) DABS(GPT_T(1)-11.95D0)
C     WRITE(*, FMT2) DABS(GPT_DT(1)+5.47D0)
C     WRITE(*, FMT2) DABS(GPT_E(1)-9.58D0)
C     WRITE(*, FMT2) DABS(GPT_AH(1)-0.0012395D0)
C     WRITE(*, FMT2) DABS(GPT_AW(1)-0.0005560D0)
C     WRITE(*, FMT2) DABS(GPT_UNDU(1)-44.06D0)

C     ALL DONE
      END
