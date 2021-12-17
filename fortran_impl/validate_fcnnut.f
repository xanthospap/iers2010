      PROGRAM MAIN

      implicit none
      double precision :: mjd, x, y, dx, dy
      CHARACTER(100) :: num1char
      CHARACTER(LEN=100), PARAMETER :: FMT1 = "(A15,F20.15,X,F20.15,X,
     .   F20.15,X,F20.15)"
      

C     READ(*,*) mjd
      CALL GET_COMMAND_ARGUMENT(1,num1char)
      READ(num1char,*) mjd

      CALL FCNNUT(mjd, x, y, dx, dy)

C     REPORT RESULTS
      WRITE(*, FMT1) 'fcnnut', x, y, dx, dy
      
      END
