      PROGRAM MAIN

C     RESULT REPORTING
      DOUBLE PRECISION T, MJD, MIN_MJD, MAX_MJD, EOP(3)
      INTEGER NUM_TESTS, I,J

      NUM_TESTS = 100
      MIN_MJD = 45700D0
      MAX_MJD = 62502D0

      WRITE (*,'(A46)') 'struct Result {double t,dx,dy,dut;} R[]={'

      DO I=1,NUM_TESTS-1
        CALL RANDOM_NUMBER(T)
        MJD = (MAX_MJD-MIN_MJD)*T + MIN_MJD
        CALL ORTHO_EOP(MJD, EOP)
C     REPORT RESULTS
        WRITE(*,'(A1,E30.22,A1,E25.18,A1,E25.17,A1,E25.17,A2)') 
     .  '{',MJD,',', EOP(1),',',EOP(2),',',EOP(3),'},'
      END DO
    
      I=NUM_TESTS
      CALL RANDOM_NUMBER(T)
      MJD = (MAX_MJD-MIN_MJD)*T + MIN_MJD
      CALL ORTHO_EOP(MJD, EOP)
      WRITE(*,'(A1,E30.22,A1,E25.18,A1,E25.17,A1,E25.17,A4)') 
     .  '{',MJD,',', EOP(1),',',EOP(2),',',EOP(3),'}};'
      
      WRITE (*,'(A16,I5,A8)') 'for (int i=0; i<',I,'; i++) {'
      WRITE (*,'(A45)') '  iers2010::ortho_eop(R[i].t,dx,dy,dut);'
      WRITE (*,'(A48)') 
     . '    assert(std::abs(dx-R[i].dx)<PRECISION_MAS);'
      WRITE (*,'(A48)') 
     . '    assert(std::abs(dy-R[i].dy)<PRECISION_MAS);'
      WRITE (*,'(A48)') 
     . '    assert(std::abs(dut-R[i].dut)<PRECISION_MS);'
      WRITE (*,'(A1)')  '}'

C End of program
      END