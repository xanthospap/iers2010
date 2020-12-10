      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F25.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for RG_ZONT2
C     ----------------------------------------------
      DOUBLE PRECISION DUTR, DLODR, DOMEGAR
      DOUBLE PRECISION DUT,  DLOD,  DOMEGA
      DATA DUTR, DLODR, DOMEGAR/
     . 7.983287678576557467D-2,
     . 5.035331113978199288D-5,
     . -4.249711616463017D-14/

C     CALL RG_ZONT2
      CALL RG_ZONT2(.07995893223819302D0, DUT, DLOD, DOMEGA)
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> RG_ZONT2 Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(DUT-DUTR)
      WRITE(*, FMT1) DABS(DLOD-DLODR)
      WRITE(*, FMT1) DABS(DOMEGA-DOMEGAR)

C     ALL DONE
      END
