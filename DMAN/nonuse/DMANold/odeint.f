
C     *************************************************
C     * odeint                                        *
C     *************************************************

C     WRITTEN BY JaeGun Jung, December 2007

C     This subroutine is a Runge-Kutta driver with adaptive stepsize
C     control. Integrate the starting values ystart(nvar) from x1 to
C     x2 with accuracy eps, storing intermediate results in the com-
C     mon block /path/. h1 should be set as a guessed first stepsize,
C     hmin as the minimum allowed stepsize (can be zero). On output
C     nok and nbad are the number of good and bad (but retried and
C     fixed) steps taken, and ystart is replaced by values at the 
C     end of the integration interval. derivnucl is the subroutine
C     for calculating the right-hand side derivative, while rkqs is
C     the name of the stepper routine to be used. /path/ contains its
C     own information about how often an intermediate value is to be
C     stored.     

C-----INPUTS------------------------------------------------------------     

C     (variable) - (description)
C     (variable) - (description)

C-----OUTPUTS-----------------------------------------------------------

C     (variable) - (description)
C     (variable) - (description)

      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,
     & derivnucl,end,rkqs,pn,ygas,inucl)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer inucl
      real end
      real pn(nsect)
      real ygas(ngas)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      real eps,h1,hmin,x1,x2,ystart(nvar),TINY
      integer i,kmax,kount,nstp

      parameter (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.0e-3)
cbdg      parameter (MAXSTP=100,NMAX=50,KMAXX=200,TINY=1.0e-30)
                           !Increase TINY when it fails 12/04/07 jgj

      real dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     &     yp(NMAX,KMAXX),yscal(NMAX)

C-----EXTERNAL FUNCTIONS------------------------------------------------

      external derivnucl, rkqs

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----COMMON VARIABLES--------------------------------------------------

      common /path/ kmax, kount, dxsav, xp, yp
           ! User storage for intermediate results. Preset dxsav and kmax.
C-----CODE--------------------------------------------------------------

      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      kmax=KMAXX
      dxsav=h1
      do i=1, nvar
        y(i)=ystart(i)
      enddo
      if (kmax.gt.0) xsav=x-2.*dxsav         ! Assures storage of first step.
      do nstp=1, MAXSTP                      ! Take at most MAXSTP steps.
        call derivnucl(nvar,x,y,dydx,pn,ygas,inucl,0.0)
        do i=1, nvar  ! Scaling used to monitor accuracy. This general-purpose
          yscal(i)=eps*(abs(y(i))+abs(h*dydx(i)))+TINY  ! can be modified if need be.
cdbg          yscal(i)=eps*abs(y(i))+TINY  ! can be modified if need be.
        enddo
        if (kmax.gt.0) then
          if (abs(x-xsav).gt.abs(dxsav)) then! Store intermediate results.
            if (kount.lt.kmax-1) then
              kount=kount+1
              xp(kount)=x
              do i=1, nvar
                yp(i,kount)=y(i)
              enddo
              xsav=x
            endif
          endif 
        endif
        if ((x+h-x2)*(x+h-x1).gt.0.) h=x2-x!If Stepsize can overshoot, decrease 
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivnucl,
     & end,pn,ygas,inucl)
        if (hdid.eq.h) then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if ((x-x2)*(x2-x1).ge.0) then      ! Are we done?
          do i=1, nvar
            ystart(i)=y(i)
          enddo
          if (kmax.ne.0) then
            kount=kount+1                  ! Save final step.
            xp(kount)=x
            do i=1, nvar
              yp(i,kount)=y(i)
            enddo
          endif
          return                           ! Normal exit.
        endif
        if (abs(hnext).lt.hmin) then
          pause 'stepsize smaller than minimum in odeint'
          print*,'abs(hnext)=',abs(hnext),' hmin=',hmin
          print*,'nok=',nok,' nbad=',nbad
          print*,'xp and yps'
          do kount=1, kmax
            print*,xp(kount),(yp(i,kount),i=1,nvar)
          enddo
        endif
        h=hnext
      enddo
      pause 'too many steps in odeint'
      print*,'nok=',nok,' nbad=',nbad
      print*,'xp and yps'
      do kount=1, kmax
        print*,xp(kount),(yp(i,kount),i=1,nvar)
      enddo
     
      


      RETURN
      END
