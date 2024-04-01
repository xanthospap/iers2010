      implicit none
      
      integer i,N,Nmax
      parameter(Nmax=10000)
      double precision tmp,tmin,tmax,dt,t,X,Y,dX,dY
      double precision DATE(NMAX),XC(NMAX),XS(NMAX),SX(NMAX)
      
!     Example of using fcnnut routine: producing X coordinate
!     between two mjd dates (tmin and tmax) with step dt and
!     using the C04 FCN model.

!     First: read the table-asc-c04.txt file to get the
!     FCN amplitudes and their epochs.
      
      open(unit=1,file='table-asc-c04.txt')
      i=1
      read(1,*)
10    read(1,*,end=11) tmp,DATE(i),XC(i),XS(i),SX(i)
      i=i+1
      goto 10
11    close(1)
      N=i-1

!     Second: call for fcnnut routine for the epochs
!     you want and print the output.

      tmin=51544.D0  ! departure date in mjd
      tmax=53000.D0  ! end date in mjd
!     dt=10.D0       ! step in days
      dt=10.03D0 
      
      t=tmin
      do while(t.le.tmax)
      call fcnnut(DATE,XC,XS,SX,N,t,X,Y,dX,dY)  ! output in microas
!     print *,t,X
      write(*,'(F20.5,1X,F18.12,1X,F18.12,1X,F18.12,1X,F18.12)')
     . t,X,Y,dX,dY
      t=t+dt
      end do
      
      end
      
      

      subroutine fcnnut(DATE,XC,XS,SX,N,mjd,X,Y,dX,dY)
      
!------------------------------------------------------------------------
!
!     Empirical model of FCN
!
!     Input:
!            DATE, table of mjd epochs of length N
!                  or col. 2 of table-asc.txt
!            XC, table of real parts of length N (microas)
!                  or col. 3 of table-asc.txt
!            XS, table of imag parts of length N (microas)
!                  or col. 4 of table-asc.txt
!            SX, table of uncertainty of length N (microas)
!                  or col. 5 of table-asc.txt
!            N, length of the tables
!            mjd, mjd epoch for which you want the FCN
!
!     Output: X, Y, CIP offsets (microas)
!             dX, dY, their uncertainties (microas)
!
!     No subroutine needed
!
!     Author: sebastien.lambert@obspm.fr
!
!------------------------------------------------------------------------

      implicit none
      
!     Input variables

      integer N
      double precision DATE(N),XC(N),XS(N),YC(N),YS(N),SX(N),SY(N)
      double precision mjd

!     Output variables
      
      double precision X,Y    ! FCN contributions
      double precision dX,dY  ! uncertainties on X, Y
      
!     Internal variables

      integer i,j
      double precision pi,per,phi,mpe
      double precision axc,axs,ayc,ays
      double precision daxc,daxs,dayc,days,dt,t

      pi=3.14159265358979323846D0
      
!     Mean prediction error

      mpe=0.3D0                         ! microas per day

!     FCN parameters

      per=-430.21D0                      ! period in days
      phi=(2.D0*pi/per)*(mjd-51544.5D0)  ! phase in rad
              
!     Amplitudes extracted from the table

      do i=1,N
         YC(i)=XS(i)
         YS(i)=-XC(i)
         SY(i)=SX(i)
      end do

!     Prediction of the amplitude at the input date

      if (mjd.le.DATE(1)) then
         axc=XC(1)
         axs=XS(1)
         ayc=YC(1)
         ays=YS(1)
      else if (mjd.ge.DATE(N)) then
         axc=XC(N)
         axs=XS(N)
         ayc=YC(N)
         ays=YS(N)
      else
         do i=1,N-1
            if (mjd.ge.DATE(i).and.mjd.lt.DATE(i+1)) then
               t=mjd-DATE(i)
               dt=DATE(i+1)-DATE(i)
               daxc=XC(i+1)-XC(i)
               daxs=XS(i+1)-XS(i)
               dayc=YC(i+1)-YC(i)
               days=YS(i+1)-YS(i)
               axc=XC(i)+(daxc/dt)*t
               axs=XS(i)+(daxs/dt)*t
               ayc=YC(i)+(dayc/dt)*t
               ays=YS(i)+(days/dt)*t
            end if
         end do
      end if
      
!     Computation of X and Y

      X=axc*dcos(phi)-axs*dsin(phi)  ! microas
      Y=ayc*dcos(phi)-ays*dsin(phi)  ! microas
      
!     Prediction of the uncertainty at the input date

      if (mjd.le.DATE(1)) then
         axc=SX(1)+mpe*(DATE(1)-mjd)
         axs=SX(1)+mpe*(DATE(1)-mjd)
         ayc=SY(1)+mpe*(DATE(1)-mjd)
         ays=SY(1)+mpe*(DATE(1)-mjd)
      else if (mjd.ge.DATE(N)) then
         axc=SX(N)+mpe*(mjd-DATE(N))
         axs=SX(N)+mpe*(mjd-DATE(N))
         ayc=SY(N)+mpe*(mjd-DATE(N))
         ays=SY(N)+mpe*(mjd-DATE(N))
      else
         do i=1,N-1
            if (mjd.ge.DATE(i).and.mjd.lt.DATE(i+1)) then
               t=mjd-DATE(i)
               dt=DATE(i+1)-DATE(i)
               daxc=SX(i+1)-SX(i)
               daxs=SX(i+1)-SX(i)
               dayc=SY(i+1)-SY(i)
               days=SY(i+1)-SY(i)
               axc=dabs(SX(i)+(daxc/dt)*t)
               axs=dabs(SX(i)+(daxs/dt)*t)
               ayc=dabs(SY(i)+(dayc/dt)*t)
               ays=dabs(SY(i)+(days/dt)*t)
            end if
         end do
      end if
      
!     Computation of the uncertainties

      dX=axc+axs  ! microas
      dY=ayc+ays  ! microas
      
      end
