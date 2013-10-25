
C     *************************************************
C     * stiff                                         *
C     *************************************************

C     WRITTEN BY JaeGun Jung, Novemeber 2007

C     This subroutine uses fourth-order Rosenbrock step for integrating
C     stiff o.d.e.'s with monitoring of local truncation error to adjust
C     stepsize. Input are the dependent variable vector y's and its deriv-
C     ative dydx's at the starting value of the independent variable x. 
C     Also input are the stepsize to be attempted htry, the required accu-
C     racy eps, and the vector yscal's against which the error is scaled.
C     On output, y and x are replaced by their new values, hdid is the st-
C     epsize that was actually accomplished, and hnext is the estimated 
C     next stepsize. derivs (here derivnucl) is a user-supplied subroutine
C     that computes the derivatives of the right-hand side with respect to
C     x, while jacobn (a fixed name) is a user-supplied subroutine that 
C     computes the Jacobi matrix of derivatives of the right-hand side with
C     respect to the components of y. Parameters: NMAX is the maximum value
C     of n; GROW and SHRNK are the largest and smallest factors by which
C     stepsize can change in one step; ERRCON = (GROW/SAFETY)**(1/PGROW)
C     and handles the case when errmax is close to zero.   

C     Cited code in Chapter 16.6 of Numerical Recipes in Fortran 77

C     Routines called: derivnucl, jacobn, lubksb, ludcmp

C-----INPUTS------------------------------------------------------------     

C     (variable) - (description)
C     (variable) - (description)

C-----OUTPUTS-----------------------------------------------------------

C     (variable) - (description)
C     (variable) - (description)

      SUBROUTINE stiff(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivnucl,
     & end,pn,ygas,inucl)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer n,MAXTRY
      real y(n), dydx(n), x, htry, eps, yscal(n), hdid, hnext
      real end
      real pn(nsect)
      real ygas(ngas)
      integer inucl

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer NMAX
      
      parameter (NMAX=50)

      integer i, j, jtry, indx(NMAX)
      real SAFETY, GROW, PGROW, SHRNK, PSHRNK, ERRCON, GAM, A21, A31, 
     & A32, A2X, A3X, C21, C31, C32, C41, C42, C43, B1, B2, B3, B4, B5,
     & E1, E2, E3, E4, C1X, C2X, C3X, C4X
      real d, errmax, h, xsav, a(NMAX,NMAX), dfdx(NMAX), dfdy(NMAX,NMAX)
     & , dysav(NMAX), err(NMAX), g1(NMAX), g2(NMAX), g3(NMAX), g4(NMAX)
     & , ysav(NMAX)

C-----EXTERNAL FUNCTIONS------------------------------------------------

      external derivnucl

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (SAFETY=0.9, GROW=1.5, PGROW=-0.25, SHRNK=0.5,
     & PSHRNK=-1./3., ERRCON=0.1296, MAXTRY=40)

      !The Kaps-Rentrop parameters
      parameter (GAM=1./2., A21=2., A31=48./25., A32=6./25., C21=-8.,
     & C31=372./25., C32=12./5., C41=-112./125., C42=-54./125.,
     & C43=-2./5., B1=19./9., B2=1./2., B3=25./108., B4=125./108.,
     & E1=17./54., E2=7./36., E3=0., E4=125./108., C1X=1./2.,
     & C2X=-3./2.,C3X=121./50., C4X=29./250., A2X=1., A3X=3./5.)

      !The Shampine parameters
cKR      parameter (GAM=0.231, A21=2.0, A31=4.52470820736, A32=
cKR     & 4.16352878860, C21=-5.07167533877, C31=6.02015272865, C32=
cKR     & 0.159750684673, C41=-1.856343618677, C42=-8.50538085819, 
cKR     & C43= -2.08407513602, B1=1.282612945268, E1=-2.30215540292, 
cKR     & E2=-3.07363448539, E3=0.873280801802, E4=1.282612945268, 
cKR     & C1X=GAM, C2X=-0.396296677520e-01, C3X=0.550778939579, 
cKR     & C4X=-0.553509845700e-01, A2X=0.462, A3X=0.880208333333)

C-----CODE--------------------------------------------------------------

      !Save initial values
      xsav= x
      do j=1, n
        ysav(i)=y(i)
        dysav(i)=dydx(i)
      enddo
   
      call jacobn(xsav,ysav,dfdx,dfdy,n,nmax,pn,ysav,inucl)
        ! The user must supply this subroutine to return the n-by-n matrix
        ! dfdy and the vector dfdx.
    
      !Set stepsize to the initial trial value.
      h=htry

      !Set up the matrix l - GAM*hf'
      do jtry=1, MAXTRY
        do i=1, n
          do j=1, n
            a(i,j)=-dfdy(i,j)
          enddo
          a(i,i)=1./(GAM*h)+a(i,i)
        enddo
        call ludcmp(a,n,NMAX,indx,d) !LU decomposition of the matrix.
        do i=1, n                    !Set up right hand side for g1.
          g1(i)=dysav(i)+h*C1X*dfdx(i) 
        enddo
        call lubksb(a,n,NMAX,indx,d) !Solve for g1.
        do i=1, n                    !Compute intermediate values of y and x.
          y(i)=ysav(i)+A21*g1(i)
        enddo
        x=xsav+A2X*h
        call derivnucl(n,x,y,dydx,pn,ygas,inucl)  
                                     !Compute dydx at the intermediate values.
        do i=1, n                    !Set up right-hand side for g2.
          g2(i)=dydx(i)+h*C2X*dfdx(i)+C21*g1(i)/h
        enddo   
        call lubksb(a,n,NMAX,indx,g2)!Solve for g2.
        do i=1, n                    !Compute intermediate values of y and x.
          y(i)=ysav(i)+A31*g1(i)+A32*g2(i)
        enddo
        x= xsav+A3X*h
        call derivnucl(n,x,y,dydx,pn,ygas,inucl)        
                                     !Compute dydx at the intermediate values.
        do i=1, n                    !Set up right-hand side for g3.
          g3(i)=dydx(i)+h*C3X*dfdx(i)+(C31*g1(i)+C32*g2(i)+C43*g3(i))/h
        enddo
        call lubksb(a,n,NMAX,indx,g3)!Solve for g3.
        do i=1, n                    !Set up right-hand side for g4.
          g4(i)=dydx(i)+h*C4X*dfdx(i)+(C41*g1(i)+C42*g2(i)+C43*g3(i))/h
        enddo
        call lubksb(a,n,NMAX,indx,g4)!Solve for g4.
        do i=1, n         !Get fourth-order estimate of y and error estimates.
          y(i)=ysav(i)+B1*g1(i)+B2*g2(i)+B3*g3(i)+B4*g4(i)
          err(i)=E1*g1(i)+E2*g2(i)+E3*g3(i)+E4*g4(i)
        enddo
        x=xsav+h
        if(x.eq.xsav)pause'stepsize not significant in stiff'
        errmax=0.                    !Evaluate accuracy.
        do i=1, n
          errmax=max(errmax,abs(err(i)/yscal(i)))
        enddo
        errmax=errmax/eps            !Scale relative to required tolerance.
        if(errmax.le.1.)then
          hdid=h
          if(errmax.gt.ERRCON)then
            hnext=SAFETY*h*errmax*PGROW
          else
            hnext=GROW*h
          endif
          RETURN
        else
          hnext=SAFETY*h*errmax**PSHRNK
          h=sign(max(abs(hnext),SHRNK*abs(h)),h)
        endif
      enddo                          !Go back and re-try step.
      pause 'exceeded MAXTRY in stiff'



      END
