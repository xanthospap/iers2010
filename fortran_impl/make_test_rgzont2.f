      PROGRAM MAIN

C     RESULT REPORTING
      DOUBLE PRECISION T, DUT,  DLOD,  DOMEGA
      INTEGER NUM_TESTS, I

      NUM_TESTS = 100

      WRITE (*,'(A46)') 'struct Result {double t,dut,dlod,domega;}R[]={'

      DO I=1,NUM_TESTS-1
        CALL RANDOM_NUMBER(T)
        CALL RG_ZONT2(T, DUT, DLOD, DOMEGA)
C     REPORT RESULTS
        WRITE(*,'(A1,E27.22,A1,E25.17,A1,E25.17,A1,E25.17,A2)') 
     .  '{',T,',',DUT,',',DLOD,',',DOMEGA,'},'
      END DO

      I=NUM_TESTS
      CALL RANDOM_NUMBER(T)
      CALL RG_ZONT2(T, DUT, DLOD, DOMEGA)
        WRITE(*,'(A1,E27.22,A1,E25.17,A1,E25.17,A1,E25.17,A3)') 
     .  '{',T,',',DUT,',',DLOD,',',DOMEGA,'}};'

      WRITE (*,'(A16,I5,A8)') 'for (int i=0; i<',I,'; i++) {'
      WRITE (*,'(A45)') '  iers2010::rg_zont2(R[i].t,dut,dlod,domega);'
      WRITE (*,'(A47)') 
     . '  assert(std::abs(dut-R[i].dut)<PRECISION_DUT);'
      WRITE (*,'(A50)') 
     . '  assert(std::abs(dlod-R[i].dlod)<PRECISION_DLOD);'
      WRITE (*,'(A56)') 
     . '  assert(std::abs(domega-R[i].domega)<PRECISION_DOMEGA);'
      WRITE(*,'(A1)') '}'


C End of program
      END
