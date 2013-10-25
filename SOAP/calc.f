      subroutine calc(nsol,sx,smw,scsat,sctot,mwpre,cpre,
     +                znum,idiag,iout,igrdchm,ichm,jchm,kchm,tempk)
      implicit none
c
c     CALC calls the HYBRD solver to calculate secondary organic
c     aerosol phase compistion within the SOAP algorithm
c
c     Called by SOAP
c
      include 'soap.com'
c
c Variable declarations
c
      integer LR, NPAR
      parameter (LR=(NSOAP*(NSOAP+1))/2)
      parameter (NPAR=NSOAP*3+2)
c
      integer nsol
      real sx(nsol), smw(nsol), scsat(nsol), sctot(nsol)
      real mwpre, cpre, tempk
      integer znum,idiag,iout,igrdchm,ichm,jchm,kchm

      real par(NPAR), fjac(NSOAP,NSOAP), fvec(NSOAP), diag(NSOAP)
      real qtf(NSOAP), wa1(NSOAP), wa2(NSOAP), wa3(NSOAP), wa4(NSOAP)
      real r(LR)
c
      real xtol, epsfcn, factor, sum, randn, sum2
      integer i, seed, jcont, nfev, info
      integer maxfev, ml, mu, mode, nprint, ldfjac, ierr
c
c SPFCN evaluates the function defining the solution for HYBRD
c
      external spfcn
c
c Entry point
c
c
c INITIALIZE VARIABLES
c
      seed=1
      jcont=1
c
c input parameters for the HYBRD solver
c xtol is the relative error tolerance
c recommended setting is sqrt(spmpar(1)), or tighter
c
      xtol = 5.0e-5
      maxfev = 200 
      epsfcn = 0.
      mode = 1
      nprint = 0
      factor = 100.
      ldfjac = NSOAP
      ml = nsol - 1
      mu = nsol - 1
c
c CREATE PARAMETER ARRAY TO PASS TO FUNCTION
c
      do i=1,nsol
         par(i)=scsat(i)
         par(nsol+i)=sctot(i)
         par(2*nsol+i)=smw(i)
      enddo
      par(3*nsol+1)=cpre
      par(3*nsol+2)=mwpre
      go to 3600
c
c LOOP THAT ASSIGNS RANDOM NUMBERS TO SX FOR SECOND GUESS AND BEYOND
c
3500  jcont=jcont+1
      if (mod(jcont,25).eq.0) then
         write (idiag,'(a,2i4,a,4i4)') 
     &        ' SOAP info = ', info, jcont,
     &           ' grd,i,j,k = ', igrdchm, ichm, jchm, kchm 
      elseif (jcont.gt.200) then
         write (iout,'(//,a)') ' ERROR in CALC:'
         write (iout,'(/,a)') ' ERROR: SOAP unable to converge'
         write (iout,'(a,i3,a,3i4)') ' Grid = ', igrdchm,
     &                    ' cell(i,j,k) = ', ichm, jchm, kchm
         write (iout,'(a,f8.1)') ' Temperature = ', tempk
         write (iout,'(a5,2a15)') 
     &                      ' spec','total ug/m3','mole fraction'
         write (iout,'(i5,1p2e15.3)') (i,sctot(i),sx(i),i=1,nsol)
         call camxerr()
      end if 
      sum = 0.0
      do i=1,nsol
         call rndnum(seed,randn)
         sx(i) = randn
         sum = sum+randn
      enddo
      if (pflag.eq.1 .and. jcont.lt.90) then
         sum = sum / (real(jcont/10+1)/10.)
      endif
      do i=1,nsol
         sx(i) = sx(i)/sum
      enddo
c
c CALL SUBROUTINE HYBRD TO CALCULATE SX
c
3600   call hybrd(spfcn,nsol,sx,fvec,xtol,par,maxfev,ml,mu,epsfcn,diag,
     *            mode,factor,nprint,info,nfev,fjac,ldfjac,r,LR,
     *            qtf,wa1,wa2,wa3,wa4,NPAR)
c
c check for problems in HYBRD
c
      if (info.eq.0) go to 9000
      if (info.ne.1) go to 3500
c
c check for significant errors in sx
c
      sum = 0.0
      sum2 = 0.0
      ierr = 0
      do i=1,nsol
         if (sx(i).lt.-xtol) ierr = 1
         sx(i) = amax1(sx(i),0.0)
         if (sx(i).gt.1.0+xtol) ierr = 1
         sum = sum+sx(i)
         sum2 = sum2 + (sctot(i)-sx(i)*scsat(i))/smw(i)
      enddo
      if (pflag.eq.1) then
        sum2 = sum2 + cpre/mwpre
        sum = sum + cpre/mwpre/sum2
      endif
      if (ierr.ne.0 .or. abs(sum-1.0).gt.nsol*xtol) go to 3500
c
c clean up any trivial errors in sx and exit
c
      do i=1,nsol
         sx(i) = sx(i)/sum
      enddo
c
      znum=jcont
      return
c
 9000 continue
      write (iout,'(//,a)') ' ERROR in CALC:'
      write (iout,'(/,2a)') ' ERROR: illegal input to HYBRD',
     &                       ' solver in SOAP'
      call camxerr()
      end
