      PROGRAM MAIN

C     RESULT REPORTING
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(A22,E22.15,A15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT2 = "(A22,I3,A25)"

C     ----------------------------------------------
C     Variables and data for FCNNUT
C     ----------------------------------------------
      DOUBLE PRECISION  FCX,  FCY,  FCDX,  FCDY
      DOUBLE PRECISION FCXR, FCYR, FCDXR, FCDYR
      DOUBLE PRECISION IMJD(24)
      INTEGER I

C     Test MJD for validation
      IMJD = (/ 54790D0, 45700D0, 45700.0000000001D0, 
     . 45700.1234500000D0, 45700.9999900000D0, 50082.123456789D0,
     . 50082.9999999999D0, 50083.0000000009D0, 50083.123456789D0,
     . 53735.8888888888D0, 53735.9888888888D0, 53735.9D0, 537360D0,
     . 53736.0000000001D0, 537360.0000000011D0, 537360.0000000111D0,
     . 537360.0000001111D0, 537360.0000111111D0, 54101D0, 54102D0,
     . 56293.0D0, 56293.9990000000D0, 56293.9999999999D0, 
     . 56793.0D0 /)

C    Print the input MJD array (C-style)
      WRITE (*,'(A25)') 'const double imjd[] = {'
      DO 11 I=1,23
        WRITE(*, '(E20.15,A1)') IMJD(I),','
   11 CONTINUE
      WRITE (*,'(E20.15,A5)') IMJD(24),'};'

C     call FCNNUT and print C-pseudocode for validation
      DO 10 I=1,24
        CALL FCNNUT(IMJD(I), FCX, FCY, FCDX, FCDY)
        WRITE(*, FMT2) 'iers2010::fcnnut(imjd[', I-1, 
     . '], fcx, fcy, fcdx, fcdy);'
        WRITE(*, FMT1) 'assert(std::abs(fcx-(',FCX,'))<PRECISION);'
        WRITE(*, FMT1) 'assert(std::abs(fcy-(',FCY,'))<PRECISION);'
        WRITE(*, FMT1) 'assert(std::abs(fcdx-(',FCDX,'))<PRECISION);'
        WRITE(*, FMT1) 'assert(std::abs(fcdy-(',FCDY,'))<PRECISION);'
   10 CONTINUE

C     ALL DONE
      END
