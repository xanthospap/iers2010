      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for CNMTX 
C     ----------------------------------------------
      INTEGER I
      DOUBLE PRECISION H(12), H_REF(12)
      DATA (H_REF(I),I=1,12) /
     .    15.35873641938967360D0,
     .    9.784941251812741214D0,
     .    -5.520740128266865554D0,
     .    3.575314211234633888D0,
     .    -13.93717453496387648D0,
     .    -9.167400321705855504D0,
     .    5.532815475865292321D0,
     .    9.558741883500834646D0,
     .    -10.22541212627272600D0,
     .    0.8367570529461261231D0,
     .    1.946355176475630611D0,
     .    -13.55702062247304696D0/

C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> CNMTX Results'
      PRINT *, '----------------------------------------'

      CALL CNMTX(54964.0D0, H)
      DO 500 I = 1, 12
        WRITE(*, FMT1) DABS(H_REF(I)-H(I))
  500 CONTINUE

C     ALL DONE
      END
