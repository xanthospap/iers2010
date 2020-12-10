      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for ORTHOEOP 
C     ----------------------------------------------
      DOUBLE PRECISION EOP(3), EOP_REF(3)
      DATA (EOP_REF(I),I=1,3) /
     .    -162.8386373279636530D0,
     .    117.7907525842668974D0,
     .    -23.39092370609808214D0/

      CALL ORTHO_EOP(47100D0, EOP)

C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> ORTHO_EOP RESULTS:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(EOP_REF(1)-EOP(1))
      WRITE(*, FMT1) DABS(EOP_REF(2)-EOP(2))
      WRITE(*, FMT1) DABS(EOP_REF(3)-EOP(3))

C     ALL DONE
      END
