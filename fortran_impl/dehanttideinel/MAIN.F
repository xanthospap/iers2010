      PROGRAM Main
      IMPLICIT NONE

C     Declare variables
      INTEGER, PARAMETER :: MAX_RECORDS = 2
      INTEGER :: YR, MONTH, DAY, IOSTAT, I
      DOUBLE PRECISION :: FHR, DUT1
      DOUBLE PRECISION :: XSTA(3), XSUN(3), XMON(3), DXTIDE(3)
      
      CHARACTER(LEN=256) :: FILENAME
      CHARACTER(LEN=512) :: LINE  ! Read the whole line
      INTEGER :: UNIT_IN

C     Ask for input file
      WRITE(*,*) 'Enter input filename:'
      READ(*,*) FILENAME

C     Open input file
      OPEN(UNIT=UNIT_IN, FILE=FILENAME, STATUS='OLD',
     &     ACTION='READ', IOSTAT=IOSTAT, FORM='FORMATTED')

      IF (IOSTAT .NE. 0) THEN
         WRITE(*,*) 'ERROR: Could not open file: ', TRIM(FILENAME)
         STOP
C     ELSE
C        WRITE(*,*) 'SUCCESS: Opened file: ', TRIM(FILENAME)
      END IF

      WRITE(*,*) 'Starting to read data...'

      I = 0
      DO WHILE (I .LT. MAX_RECORDS)
         READ(UNIT_IN, '(A)', IOSTAT=IOSTAT) LINE
         IF (IOSTAT .NE. 0) THEN
            WRITE(*,*) 'ERROR: EOF or read error detected on line ', I
            EXIT
         ENDIF

C        WRITE(*,*) 'Raw input line:', TRIM(LINE)

         READ(LINE, *, IOSTAT=IOSTAT) YR, MONTH, DAY, FHR, DUT1, 
     &    XSTA(1), XSTA(2), XSTA(3), XSUN(1), XSUN(2), XSUN(3), 
     &    XMON(1), XMON(2), XMON(3)

         IF (IOSTAT .NE. 0) THEN
            WRITE(*,*) 'ERROR: Failed to parse line:', TRIM(LINE)
            EXIT
         ENDIF

C        WRITE(*,*) 'Parsed values:', YR, MONTH, DAY, FHR,
C    &    XSTA(1), XSTA(2), XSTA(3), XSUN(1), XSUN(2), XSUN(3),
C    &    XMON(1), XMON(2), XMON(3)

C        WRITE(*,*) 'Calling DEHANTTIDEINEL...'
         CALL DEHANTTIDEINEL(XSTA, YR, MONTH, DAY, FHR, XSUN, XMON,
     &        DXTIDE)

         WRITE(*,*) YR, MONTH, DAY, FHR, DUT1, XSTA(1), XSTA(2), 
     &      XSTA(3), XSUN(1), XSUN(2), XSUN(3), XMON(1), XMON(2), 
     &      XMON(3), DXTIDE(1), DXTIDE(2), DXTIDE(3)

C        WRITE(*,*) 'Iteration:', I
         I = I + 1
      END DO

C     Close input file
      CLOSE(UNIT_IN)

      WRITE(*,*) 'Processing complete.'

      END PROGRAM Main
