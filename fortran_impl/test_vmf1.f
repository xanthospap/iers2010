      PROGRAM MAIN

C     RESULT REPORTING
C     CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(F18.15)"
      CHARACTER(LEN=*), PARAMETER  :: FMT1 = "(E12.6)"

C     ----------------------------------------------
C     Variables and data for VMF1
C     ----------------------------------------------
      DOUBLE PRECISION VMF1H, VMF1W

      CALL VMF1(0.00127683D0, 0.00060955D0, 55055D0, 0.6708665767D0, 
     .  1.278564131D0, VMF1H, VMF1W)
     
C     REPORT RESULTS
      PRINT *, '----------------------------------------'
      PRINT *, '> VMF1 Results:'
      PRINT *, '----------------------------------------'
      WRITE(*, FMT1) DABS(VMF1H-3.424342122738070593D0)
      WRITE(*, FMT1) DABS(VMF1W-3.448299714692572238D0)

C     ALL DONE
      END
