
C     **************************************************
C     * nucl_sub                                       *
C     **************************************************

C     WRITTEN BY JaeGun Jung, November 2007

C     This subroutine prepare to implement differential equations. 
C     for solving nucleation. The original code is aero_dms.f in
C     previous version of codes.

      SUBROUTINE nucl_sub(ygas,pn_p,tstart,tend,pn,inucl)

      IMPLICIT NONE

C-----INPUTS------------------------------------------------------------
 
C     inucl -  a nucleation flag
C     ygas - gas phase species mixing ratio as ppt
C     tstart - start time of nucleation calculation
C     tend - end time of nucleation calculation

C-----OUTPUTS-----------------------------------------------------------

C     pn_p - new particle concentrations
C     ygas - gas phase species mixing ratio changed

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer inucl ! a nucleation flag
      real ygas(ngas)
      real pn_p(nsect)
      real tstart
      real tend
      real pn(nsect)

C-----VARIABLES DECLARATIONS--------------------------------------------

      integer i ! a counter
      integer nok, nbad ! number of ok and bad in odeint subroutine
      real yscal(keqna)
      real ylast(keqna)
      real maxrelerr
      real dtlast
      real dtnext
      real timelast
      real y(keqna)
      real yprime(keqna)
      real time
      real dttry
      real yout(keqna)

C-----EXTERNAL FUNCTIONS------------------------------------------------

      external derivnucl, rkqs

C-----CODE--------------------------------------------------------------

cdbg      print*,'boxmass=',boxmass
cdbg      pause

      ! Converting Gc's to ygas's, kg to ppt
      ygas(mgsvi)=1.0e+12*28.9/(boxmass*gmw(srtso4))*Gc(srtso4)
      ygas(mgnh3)=1.0e+12*28.9/(boxmass*gmw(srtnh3))*Gc(srtnh3)

      ! Initialization
      y(1)=pn_p(1) ! the smallest size bin
      y(2)=ygas(mgsvi) 
      y(3)=ygas(mgnh3) 

C-----Runge Kutta with time step checking (rkqc from Numerical Recipies)
  
      time = tstart
      dttry=1.e-3  ! hr
cdbg      maxrelerr=1.0e-6
      maxrelerr=5.e-2
cdbg      maxrelerr=1.0e+3
      time=tstart
      timelast=tstart
      call derivnucl(keqna,time,y,yprime,pn,ygas,inucl,0.0) 
          ! keqna is defined in aervaria.inc.
          ! It is a number of equation will be calculated. 
      do while (time .lt. tend)
        call rkqs(y,yprime,keqna,time,dttry,maxrelerr,yscal,dtlast,
     &              dtnext,derivnucl,tend,pn,ygas,inucl)
crkqs        call odeint(y,keqna,time,tend,maxrelerr,dttry,dttry*1e-4,nok,
crkqs     &              nbad,derivnucl,tend,rkqs,pn,ygas,inucl)
crkqs        call stiff(y,yprime,keqna,time,dttry,maxrelerr,yscal,dtlast,
crkqs     &              dtnext,derivnucl,tend,pn,ygas,inucl)
        dttry=dtnext   
        if (time .lt.tend) then
          timelast=time
          do i=1, keqna
            ylast(i)=y(i)
          enddo
        else
          do i=1, keqna
            y(i)=(y(i)-ylast(i))/(time-timelast)*(tend-timelast)
     &        +ylast(i)
          enddo
          time=tend
        endif
      enddo

C-----Calculated y values by adaptive runge kutta method go to output values.
 
      !Save calculated values
      if (y(1).ge.0.0) then
        pn_p(1)=y(1)
      else
        pn_p(1)=0.0
      endif

      if (y(2).ge.0.0) then
        ygas(mgsvi)=y(2)
      else
        ygas(mgsvi)=0.0
      endif

      if (y(3).ge.0.0) then
        ygas(mgnh3)=y(3)
      else
        ygas(mgnh3)=0.0
      endif

      if (icond_test .eq. 1) then
        Gc(srtso4)=boxmass*ygas(mgsvi)*1.0e-12*100.0/28.9
      else
        Gc(srtso4)=boxmass*ygas(mgsvi)*1.0e-12*gmw(srtso4)/28.9
      endif
      Gc(srtnh3)=boxmass*ygas(mgnh3)*1.0e-12*gmw(srtnh3)/28.9

      RETURN
      END
