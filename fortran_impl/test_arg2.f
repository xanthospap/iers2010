      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for ARG2
C     ----------------------------------------------
      INTEGER I
      DOUBLE PRECISION ANGLE_REF(11), ANGLE(11)
      DATA (ANGLE_REF(I),I=1,11) /
     .   2.849663065753787805D0,
     .   6.28318080000000023D0,
     .   4.926040134021299366D0,
     .   1.608450491115348768D0,
     .   2.375021572352622456D0,
     .   0.4746414933980958040D0,
     .   3.908159227647345801D0,
     .   2.551018561669245344D0,
     .   5.041990012540757959D0,
     .   4.206816878908014701D0,
     .   1.608463638294885811D0/

C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> ARG2 Results'
      PRINT *, '----------------------------------------'

      CALL ARG2(2008, 311.5D0, ANGLE)
      DO 500 I = 1, 11
        WRITE(*, FMT1) DABS(ANGLE_REF(I)-ANGLE(I))
  500 CONTINUE

C     ALL DONE
      END
