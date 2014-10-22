
C     **************************************************
C     * rkqs                                           *
C     **************************************************

C     RE-WRITTEN BY JaeGun Jung, November 2007

C     This subroutine solve differential equation using Runge
C     Kutta adaptive time step, which allows controlling time 
C     step based on truncation, error between y(h) and 2*y(h/2). 

C     Referrence: NUMERICAL RECIPES IN FORTRAN 77
C     Chapter 16.2 Adaptive Stepsize Control for Runge-Kutta

C     Fifth-order Runge-Kutta step with monitoring of local trucation
C     error to ensure accuracy and adjust stepsize. Input are the
C     dependent variable vector y's and its derivative dydx's at the
C     starting value of the independent variable x. Also input are the 
C     stepsize to be attempted htry, the required accuracy eps, and
C     the vector yscal against which the error is scaled. On output,
C     y and x are replaced by their new values, hdid is the stepsize
C     that was actually accomplished, and hnext is the estimated next
C     stepsize. derivs (here derivnucl) is the user-supplied subroutine
C     that computes the right-hand side derivatives.

C                            Last change:  JGJ    19 Jul 04  

C-----INPUTS------------------------------------------------------------

C     Description of VARIABLES 
C     ========================

C     YSCAL = Expecting truncation error with HTRY
C     HTRY = Time step to try
C     HNEXT = Time step that will be changed after
C     YERR = Truncation error between y(h) and 2*y(h/2)

C     YSCAL = RELTOL * YTEMP + 1.0E-3 (NOW).
C     Now it is impossible to control yscal based on DYDX
C     because YSCAL is just truncation error. Truncation
C     error is not related with DYDX. And smallest particles
C     from nucleation are more than 1.0e+8 as maximum rate.
C     So the control with DYDT is less useful for nucleation
 
C                                             COMMENT by JGJ

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE RKQS(Y,DYDX,N,x,HTRY,reltol,YSCAL,HDID,HNEXT,derivnucl,
     & end,pn,ygas,inucl)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------
   
      integer n       !number of equations
      real y(n)       !y values in y = f(x)
      real dydx(n)    !derivatives
      real x          !time of main routine
      real htry       !dx trying
      real reltol     !maximum relative error
      real yscal(n)
      real hdid       !dx tried
      real hnext      !dx to try
      real end        !end time
      real pn(nsect)
      real ygas(ngas)
      integer inucl   ! a nucleation flag

C-----VARIABLES DECLARATIONS--------------------------------------------

      integer nmax
      integer i
      real eps
      real errmax,h,htemp,xnew
      real safety
      real pgrow
      real pshrnk
      real errcon
      real deltat ! time step in this subroutine
      real ltf    ! left time fraction
      real ltfset ! left time fraction setup

      parameter (nmax=50)

      real yerr(nmax)
      real YTEMP(NMAX)

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      PARAMETER (SAFETY=0.9,ERRCON=1.89E-4)

      PARAMETER (PGROW=-0.20,PSHRNK=-0.25)

      parameter (ltfset=1.0e-3) !left time fraction setup to prevent
                                !an infinite loop
      parameter (eps=1.0e-3)    !a minimum allowable relative error
                                !eqn (16.6.20) in the referrence

C-----EXTERNAL FUNCTIONS------------------------------------------------

      external derivnucl

C-----CODE--------------------------------------------------------------

      !Initialization
      h=htry

 1    call rk4nucl(y,dydx,n,x,h,ytemp,yerr,derivnucl,pn,ygas,inucl)
      errmax=0.
      DO 11 I=1,N
cdbg        YSCAL(I)=reltol*ytemp(i)+1.0e-3
        YSCAL(I)=reltol*ytemp(i)+1.0e-40
cdbg        YSCAL(I)=max(eps,abs(reltol*ytemp(i)))
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
 11   continue
codeint      errmax=errmax/eps               ! Scale relative to required tolerance
                                      ! It is added by including odeint
      do i=1,n                        !dbg
        if (yscal(i).eq.0.) then      !dbg
          write(*,*)'yerr yscal, i', yerr(i), yscal(i), i !dbg
          write(*,*)'ytemp(i)',ytemp(i) !dbg
        endif                         !dbg
      enddo                           !dbg

      if (errmax .gt. 1.) then
        htemp=safety*h*(errmax**pshrnk)
        h=sign(max(abs(htemp),0.1*abs(h)),h)
        if (x+h.gt.end) h=end-x ! tmg (08/15/03)
        xnew=x+h
        if (xnew.eq.x) then
cdbg          deltat=(end-x) !time step here
cdbg          ltf=h/deltat !left time fraction
cdbg          if (ltf.le.ltfset) then ! Current h is less than 0.1 % of deltat
cdbg            x=end
cdbg            do i=1,n
cdbg              if(y(i).lt.0.0)y(i)=0.0 !set all y positive, physical properties
cdbg            enddo
cdbg            RETURN
cdbg          else
cdbg            if (ltf.le.ltfset*50) then ! Current h is less than 5 % of deltat
cdbg              write(*,*)'WARNING:'
cbdg              write(*,*)'Current h is less than or equal to 5 % of 
cbdg&deltat'                        !ltf is less than or equal to ltfset*50
cdbg              write(*,*)'h=',h
cdbg              x=end
cdbg              do i=1,n
cdbg                if(y(i).lt.0.0)y(i)=0.0 !set all y positive
cdbg              enddo
cdbg              RETURN
cdbg            else
              write(*,*)'stepsize underflow in rkqc' 
              write(*,*)'dydx=',dydx
              write(*,*)'y=',y
              write(*,*)'ltf=',ltf
              write(*,*)'h=',h
cbdg            endif
cdbg          endif
        endif 
        goto 1
      else
        IF(ERRMAX.GT.ERRCON)THEN
          HNEXT=SAFETY*H*(ERRMAX**PGROW)
        ELSE
          HNEXT=5.*H
        ENDIF
        hdid=h
        x=x+h
        DO 12 I=1,N
          Y(I)=YTEMP(I)
 12     CONTINUE
      endif 


      RETURN
      END

C=======================================================================

C     *************************************************
C     * rk4nucl                                       *
C     *************************************************

C     RE-WRITTEN BY JaeGun Jung, November 2007

C     This subroutine implements Runge Kutta method for solving ODEs

      SUBROUTINE RK4NUCL(Y,DYDX,N,X,H,YOUT,yerr,derivnucl,pn,ygas,inucl)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer inucl
      integer n
      real Y(N),DYDX(N),YOUT(N),yerr(n)
      real x
      real h
      real pn(nsect)
      real ygas(ngas)
cdbg      DIMENSION Y(N),DYDX(N),YOUT(N),yerr(n)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i
      integer nmax

      PARAMETER (NMAX=50)

      real ak2(nmax),ak3(nmax),ak4(nmax),ak5(nmax),ak6(nmax),
     &     ytemp(nmax),a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,
     &     b52,b53,b54,b61,b62,b63,b64,b65,c1,c3,c4,dc1,dc3,
     &     dc4,dc5,dc6
      real c6

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(a2=.2,a3=.3,a4=.6,a5=1.,a6=.875,b21=.2,b31=3./40.,
     &     b32=9./40.,b41=.3,b42=-.9,b43=1.2,b51=-11./54.,b52=2.5,
     &     b53=-70./27.,b54=35./27.,b61=1631./55296.,b62=175./512.,
     &     b63=575./13824.,b64=44275./110592.,b65=253./4096.,
     &     c1=37./378.,c3=250./621.,c4=125./594.,c6=512./1771.,
     &     dc1=c1-2825./27648.,dc3=c3-18575./48384.,
     &     dc4=c4-13525./55296.,dc5=-277./14336.,dc6=c6-.25)

C-----EXTERNAL FUNCTIONS------------------------------------------------

      external derivnucl

C-----CODE--------------------------------------------------------------

      DO 11 I=1,N
        YTEMP(I)=Y(I)+b21*H*DYDX(I) 
11    CONTINUE
      CALL DERIVNUCL(n,X+a2*h,YTEMP,ak2,pn,ygas,inucl,a2*h)
      DO 12 I=1,N
        YTEMP(I)=Y(I)+H*(b31*DYDX(I)+b32*ak2(i))
12    CONTINUE
      CALL DERIVNUCL(n,X+a3*h,YTEMP,ak3,pn,ygas,inucl,a3*h) 
      DO 13 I=1,N
        YTEMP(I)=Y(I)+H*(b41*DYDX(I)+b42*ak2(i)+b43*ak3(i))
13    CONTINUE
      CALL DERIVNUCL(n,X+a4*H,YTEMP,ak4,pn,ygas,inucl,a4*H)
      DO 14 I=1,N
        YTEMP(I)=Y(I)+H*(b51*DYDX(I)+b52*ak2(i)+b53*ak3(i)+
     &         b54*ak4(i))
14    CONTINUE
      call derivnucl(n,x+a5*h,ytemp,ak5,pn,ygas,inucl,a5*h)
      do 15 i=1,n
        YTEMP(i)=Y(i)+H*(b61*DYDX(i)+b62*ak2(i)+b63*ak3(i)+
     &           b64*ak4(i)+b65*ak5(i))
15    continue
      call derivnucl(n,x+a6*h,ytemp,ak6,pn,ygas,inucl,a6*h)
      do 16 i=1,n
         yout(i)=y(i)+h*(c1*dydx(i)+c3*ak3(i)+c4*ak4(i)+
     &          c6*ak6(i))
16    continue
      do 17 i=1,n
         yerr(i)=h*(dc1*dydx(i)+dc3*ak3(i)+dc4*ak4(i)+dc5*ak5(i)
     &          +dc6*ak6(i))  ! truncation err= delta of Eq(16.2.6)
17    continue


      RETURN
      END
