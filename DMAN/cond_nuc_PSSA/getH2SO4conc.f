
C     **************************************************
C     *  getH2SO4conc                                  *
C     **************************************************

C     WRITTEN BY Jeff Pierce, May 2007

C     This subroutine uses newtons method to solve for the steady state 
C     H2SO4 concentration when nucleation is occuring.

C     It solves for H2SO4 in 0 = P - CS*H2SO4 - M(H2SO4)

C     where P is the production rate of H2SO4, CS is the condensation sink
C     and M(H2SO4) is the loss of mass towards making new particles.

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     H2SO4rate - H2SO4 generation rate [kg box-1 s-1]
C     CS - condensation sink [s-1]
C     NH3conc - ammonium in box [kg box-1]
Cxxxxx     prev - logical flag saying if a previous guess should be used or not
Cxxxxx     gasConc_prev - the previous guess [kg/box] (not used if prev is false)

C-----OUTPUTS-----------------------------------------------------------

C     gasConc - gas H2SO4 [kg/box]

      SUBROUTINE getH2SO4conc(H2SO4rate,CS,NH3conc,gasConc, dmappt)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision H2SO4rate
      double precision CS
      double precision NH3conc
      double precision gasConc
      logical prev
      double precision gasConc_prev

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      double precision fn, rnuc ! nucleation rate [# cm-3 s-1] and critical radius [nm]
	double precision fn1, rnuc1 ! nucleation rate [# cm-3 s-1] and critical radius [nm]
      double precision res      ! d[H2SO4]/dt, need to find the solution where res = 0
      double precision massnuc     ! mass being removed by nucleation [kg s-1 box-1]
      double precision gasConc1 ! perturbed gasConc
      double precision gasConc_hi, gasConc_lo
      double precision res1     ! perturbed res
      double precision res_new  ! new guess for res
      double precision dresdgasConc ! derivative for newtons method
      double precision Gci(icomp-1)      !array to carry gas concentrations
      logical nflg              !says if nucleation occured
      double precision H2SO4min !minimum H2SO4 concentration in parameterizations (molec/cm3)
      double precision pi
      integer iter,iter1
      integer itermax !Set maximum iteration by jgj 
	double precision CSeps ! low limit for CS
	double precision max_H2SO4conc !maximum H2SO4 concentration in parameterizations (kg/box)
	double precision nh3ppt   !ammonia concentration in ppt
	double precision dmappt   !dimethyl amine concentration in ppt
      double precision maxres     !max res 50 ppt/hr = 16.7 kg box-1 s-1
                                  !by jgj 11/25/07
 
C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------


C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654)
      parameter(H2SO4min=1.D4) !molecules cm-3
	parameter(CSeps=1.0d-20)
      parameter(maxres=16.7) 
      parameter(itermax=1000) !the max iteration

C-----CODE--------------------------------------------------------------

      do i=1,icomp-1
         Gci(i)=0.d0
      enddo
      Gci(srtnh4)=NH3conc
cdbg      print*,'At the beginning of getH2SO4conc'
cdbg      print*,'H2SO4rate=',H2SO4rate,'gasConc=',gasConc
cdbg      print*,'CS=',CS,'CSeps=',CSeps
cdbg      print*,'H2SO4rate/CS=',H2SO4rate/CS

	! make sure CS doesn't equal zero
c	CS = max(CS,CSeps)

	! some specific stuff for napari vs. vehk
	if ((bin_nuc.eq.1).or.(tern_nuc.eq.1))then
                nh3ppt = (1.0e+21*8.314)*Gci(srtnh4)*temp/(pres*boxvol
     &                   *gmw(srtnh4))
cjgj		!nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12
cjgj		nh3ppt = nh3ppt_o
		if ((nh3ppt.gt.100.0).and.(tern_nuc.eq.1))then
cjgj			max_H2SO4conc=1.0D9*boxvol/1000.d0*98.d0/6.022d23
			max_H2SO4conc=1.0D11*boxvol/1000.d0*98.d0/6.022d23
                !The bounded max_H2SO4conc cannot make res below zero within
                !given condition sometimes. Set it large to be converged.
                !by jgj 01/10/2008
		elseif (bin_nuc.eq.1)then
			max_H2SO4conc=1.0D11*boxvol/1000.d0*98.d0/6.022d23
		else
			max_H2SO4conc = 1.0D100
		endif
	else
		max_H2SO4conc = 1.0D100
	endif	

C Checks for when condensation sink is very small
	if (CS.gt.CSeps) then
		gasConc = H2SO4rate/CS
	else
		if ((bin_nuc.eq.1).or.(tern_nuc.eq.1)) then
			gasConc = max_H2SO4conc
		else
			print*,'condesation sink too small in getH2SO4conc.f'
			STOP
		endif
	endif

	gasConc = min(gasConc,max_H2SO4conc)
      Gci(srtso4) = gasConc
      call getNucRate(Gci,fn,rnuc,nflg,CS,dmappt)

      if (nflg) then ! nucleation occured
         gasConc_lo = H2SO4min*boxvol/(1000.d0/98.d0*6.022d23) !convert to kg/box

C     Test to see if gasConc_lo gives a res < 0 (this means ANY nucleation is too high)
         Gci(srtso4) = gasConc_lo*1.000001d0
         call getNucRate(Gci,fn1,rnuc1,nflg,CS,dmappt)
         if (nflg) then
            massnuc = 4.d0/3.d0*pi*(rnuc1*1.d-9)**3*1350.*fn1*boxvol*
c            massnuc = 4.d0/3.d0*pi*(rnuc1*1.d-9)**3*1800.*fn1*boxvol*
     &           98.d0/96.d0
cdbg            print*,'res',res
cdbg            print*,'H2SO4rate',H2SO4rate
cdbg            print*,'CS',CS
cdbg            print*,'CS*gasConc_lo',CS*gasConc_lo
Cjrp            print*,'mnuc',mnuc
            res = H2SO4rate - CS*gasConc_lo - massnuc
            if (res.lt.0.d0) then ! any nucleation too high
               print*,'nucleation cuttoff'
               gasConc = gasConc_lo*1.000001 ! make gasConc too low for nucleation to occur
               return
            endif
         endif
         
         gasConc_hi = gasConc ! we know this must be the upper limit (since no nucleation)
         !take density of nucleated particle to be 1350 kg/m3
         massnuc = 4.d0/3.d0*pi*(rnuc*1.d-9)**3*1350.*fn*boxvol*
c         massnuc = 4.d0/3.d0*pi*(rnuc*1.d-9)**3*1800.*fn*boxvol*
     &        98.d0/96.d0
Cjrp         print*,'H2SO4rate',H2SO4rate,'CS*gasConc',CS*gasConc,
Cjrp     &        'mnuc',mnuc
         res = H2SO4rate - CS*gasConc - massnuc

	! check to make sure that we can get solution
            if (res.gt.maxres) then
              print*,'res = H2SO4rate - CS*gasConc - massnuc'
              print*,'res',res,' H2SO4rate',H2SO4rate,' CS',CS
              print*,'gasConc',gasConc,' massnuc',massnuc
	      print*,'gas production rate too high in getH2SO4conc.f'
	      return
c	      STOP
            endif

         iter = 0
Cjrp         print*, 'iter',iter
Cjrp         print*,'gasConc_lo',gasConc_lo,'gasConc_hi',gasConc_hi
Cjrp         print*,'res',res
         do while (abs(res) .gt. 1.D-8)
            iter = iter+1
            if (res .lt. 0.d0) then ! H2SO4 concentration too high, must reduce
               gasConc_hi = gasConc ! old guess is new upper bound
            elseif (res .gt. 0.d0) then ! H2SO4 concentration too low, must increase
               gasConc_lo = gasConc ! old guess is new lower bound
            endif
Cjrp            print*, 'iter',iter
Cjrp            print*,'gasConc_lo',gasConc_lo,'gasConc_hi',gasConc_hi
            gasConc = sqrt(gasConc_hi*gasConc_lo) ! take new guess as logmean
            Gci(srtso4) = gasConc
            call getNucRate(Gci,fn,rnuc,nflg,CS,dmappt)
            massnuc = 4.d0/3.d0*pi*(rnuc*1.d-9)**3*1350.*fn*boxvol*
c 	      massnuc = 4.d0/3.d0*pi*(rnuc*1.d-9)**3*1800.*fn*boxvol*
     &           98.d0/96.d0
            res = H2SO4rate - CS*gasConc - massnuc
         if (iter.gt.itermax) then
           print*,'iter',iter
           print*,'res = H2SO4rate - CS*gasConc - massnuc'
           print*,'res',res,' H2SO4rate',H2SO4rate,' CS',CS
           print*,'gasConc',gasConc,' massnuc',massnuc
	   print*,'gasConc too high in getH2SO4conc.f'
           print*,'nh3ppt=',nh3ppt
           print*,'fn',fn,'rnuc',rnuc
           STOP
         endif
         enddo

cdbg         print*,'IN getH2SO4conc'
cdbg         print*,'fn',fn
cdbg         print*,'H2SO4rate',H2SO4rate
cdbg         print*,'massnuc',massnuc,'CS*gasConc',CS*gasConc

      else ! nucleation didn't occur
      endif
cjgj      print*,'gasConc=',gasConc
cdbg      print*,'At the end of getH2SO4conc'
cdbg      print*,'H2SO4rate=',H2SO4rate,'gasConc=',gasConc

      return
      end
