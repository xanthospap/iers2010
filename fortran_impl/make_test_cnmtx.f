      PROGRAM MAIN

C     RESULT REPORTING
      DOUBLE PRECISION T, MJD, MIN_MJD, MAX_MJD, H(12)
      INTEGER NUM_TESTS, I,J

      NUM_TESTS = 100
      MIN_MJD = 45700D0
      MAX_MJD = 62502D0

      WRITE (*,'(A46)') 'struct Result {double t, h[12];} R[]={'

      DO I=1,NUM_TESTS-1
        CALL RANDOM_NUMBER(T)
        MJD = (MAX_MJD-MIN_MJD)*T + MIN_MJD
        CALL CNMTX(MJD,H)
C     REPORT RESULTS
        WRITE(*,'(A1,E30.22,A3)') '{',MJD,',{'
        DO J=1,11
            WRITE(*,'(E25.17,A1)') H(J), ','
        END DO
        J = 12
        WRITE(*,'(E25.17,A3)') H(J), '}},'
      END DO

      I=NUM_TESTS
      CALL RANDOM_NUMBER(T)
      MJD = (MAX_MJD-MIN_MJD)*T + MIN_MJD
      CALL CNMTX(MJD,H)
C     REPORT RESULTS
      WRITE(*,'(A1,E30.22,A3)') '{',MJD,',{'
      DO J=1,11
          WRITE(*,'(E25.17,A1)') H(J), ','
      END DO
      J = 12
      WRITE(*,'(E25.17,A5)') H(J), '}}};'


      WRITE (*,'(A16,I5,A8)') 'for (int i=0; i<',I,'; i++) {'
      WRITE (*,'(A35)') '  iers2010::oeop::cnmtx(R[i].t,h);'
      WRITE (*,'(A35)') '  for (int j=0; j<12; j++) {'
      WRITE (*,'(A48)') 
     . '    assert(std::abs(h[j]-R[i].h[j])<PRECISION);'
      WRITE (*,'(A4)')  '  }'
      WRITE (*,'(A1)')  '}'

C End of program
      END