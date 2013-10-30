
C     **************************************************
C     *  so4cond_oxd                                   *
C     **************************************************

C     WRITTEN BY Peter Adams, June 2000
C     Based on thermocond.f

C     This is a version for oxidation of aqueous chemistry!

C     This routine determines the condensational driving force for
C     mass transfer of sulfate between gas and aerosol phases.  It then calls
C     a mass- and number-conserving algorithm for condensation (or
C     evaporation) of aerosol.

C     An adaptive time step is used to prevent numerical difficulties.
C     To account for the changing gas phase concentration of sulfuric
C     acid, its decrease during a condensational time step is well-
C     approximated by an exponential decay with a constant, sK (Hz).
C     sK is calculated from the mass and number distribution of the
C     aerosol.  Not only does this approach accurately take into account
C     the changing sulfuric acid concentration, it is also used to
C     predict (and limit) the final sulfuric acid concentration.  This
C     approach is more accurate and faster (allows longer condensational
C     time steps) than assuming a constant supersaturation of sulfate.

C     Ammonia condensation is added, the species that are not taken account
C     are skipped in two j(species) loops. The gas specise that are tiny
C     set the flag rather than going to end of subroutine. It is because
C     multicompent species are dealt in this subroutine compared to this
C     subroutine's original code. 12/03/2007 by JaeGun Jung

C     A correction factor, corfactor is added after tested with analytic 
C     solution. 12/11/2007 by JaeGun Jung

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell
C     dt - total model time step to be taken (s)

C-----OUTPUTS-----------------------------------------------------------

ccondtest      SUBROUTINE so4cond(dt,ygas)
cdbg      SUBROUTINE so4cond(Nki,Mki,Gci,Nkf,Mkf,Gcf,dt,xkDMAN)
cnogas      SUBROUTINE so4cond(Nki,Mki,Gci,Nkf,Mkf,Gcf,dt,moxid
cnogas     &                  ,ichm,jchm,kchm)
      SUBROUTINE so4cond_oxd(Nki,Mki,Nkf,Mkf,dt,moxid,iact,ichm,jchm
     &                      ,kchm,xkDMAN)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'
      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

cdbg      real ygas(ngas) ! gas species mixing ratio
cnogas      double precision Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
cnogas      double precision Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      double precision Nki(ibins), Mki(ibins, icomp)
      double precision Nkf(ibins), Mkf(ibins, icomp)
      real dt ! timestep
      double precision xkDMAN(ibins+1) ! xk as inputs
      real*4 moxid(ibins, icomp-1) 
                   ! sulfate produced by aqueous chemistry [=] ug/m3
      integer iact ! a starting activation bin
      integer ichm ! i coordinate in PMCAMx
      integer jchm ! j coordinate in PMCAMx
      integer kchm ! k coordinate in PMCAMx

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer j,k,jj,kk        !counters
      integer Gcflag(icomp)    !flags for checking insignificant Gc
cdbg      integer itr, itrts       !iterations for master loop and timestep
cdbg      integer nsteps, floor, ceil
cdbg      integer insteps          !a counter for nsteps

      real time                !amount of time (s) that has been simulated
      real cdt                 !internal, adaptive time step
      real mu                  !viscosity of air (kg/m s)
      real mfp                 !mean free path of air molecule (m)
      real Kn                  !Knudsen number of particle
      real Dpk(ibins)          !diameter (m) of particles in bin k
      real density             !density (kg/m3) of particles
      real Di                  !diffusivity of gas in air (m2/s)
      real beta(icomp-1)       !correction for non-continuum
      real ratio               !ammonia mass ratio with respect to sulfate
      real alphanh3       ! ammonia accomodation coef. in case of calculation
      real corfactor!Top-hat method in tmconds produce a retardation of growth
                    !corfactor is decided by comparing an analytical solution.
      real gasfrac  !gas reduction rate is regulated depending on timestep ,jgj
      real cvt

      double precision dp(ibins, icomp-1)  !Driving force for condensation (Pa)
      double precision tau(ibins)          !condensation parameter (see cond.f)
      double precision atau(ibins, icomp)  !same as tau, but all species
      double precision atauc(ibins, icomp) !same as atau, but for const dp
      double precision tj(icomp-1), tk(ibins,icomp-1)  
                                !factors used for calculating tau
      double precision sK(icomp)!exponential decay const for gas, jgj
      double precision R        !constants
      double precision zeta13   !from Eqn B30 of Tzivion et al.(1989)
      double precision mp       !particle mass (kg)
      double precision Nko(ibins), Mko(ibins, icomp), Gco(icomp-1) 
      double precision tdt      !the value 2/3
                                !output of cond routine
      double precision mi(icomp)!initial and final aerosol masses (updates Gc)
      double precision mf(icomp)!initial and final aerosol masses (updates Gc)
      double precision tr       ! used to calculate time step limits
      double precision mc, ttr
      double precision Neps     !value below which Nk is insignificant
      double precision cthresh  !determines minimum gas conc. for cond.
cneg_numb_cons      double precision Ntotf, Ntoto, dNerr  !used to track number cons.

      double precision Mknh3max !Maximum allowable NH3
      double precision taumax   !Maxmum tau for ammonia condensation
      double precision moxd(ibins) !moxid/Nact, NOT USING IT!

      character*12 limit        !description of what limits time step


C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------

      real aerodens_PSSA
      external aerodens_PSSA

      real gasdiff
      external gasdiff

      real alpha_nh3
      external alpha_nh3

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(R=8.314) !gas constant (J/mol K)
      parameter(Neps=1.0d+10, zeta13=0.98483, cthresh=1.d-16)
                         !Neps is set a little higher than multicoag 
                         !to avoid redundant time step segregations.
cjgj      parameter(corfactor=1.4) ! a correction factor
      parameter(corfactor=1.0) ! a correction factor
      parameter(cvt = 1.0d-15) ! a convert factor from ug/m3 to kg/cm3

C-----CODE--------------------------------------------------------------

cdbg      write(*,*)'xk(j)',(xk(j),j=1,ibins)
cdbg      print*,'boxmass=',boxmass
cdbg      pause
cdbg      write(*,*) '\n Coord.(i,j,k)=',ichm,jchm,kchm

      !Initialization
      k=0
      jj=0
cdbg      limit='null'
cneg_numb_cons      dNerr=0.0
      time=0.0         !subroutine exits when time=dt
      tdt=2.d0/3.d0
cdbg      itr=0
cdbg      itrts=0
cnogas      do j=1,icomp-1
cnogas        Gcf(j)=Gci(j)
cnogas      enddo
      do k=1,ibins
        Nkf(k)=Nki(k)
        do j=1,icomp
          Mkf(k,j)=Mki(k,j)
        enddo
      enddo

      do j=1,icomp
        Gcflag(j)=0    !=0, Do insignificant gas concentration command
      enddo

      !xk is provided
      do k=1,ibins+1
        xk(k)=xkDMAN(k)
      enddo

cdbg        do k=1, ibins+1
cdbg          print*,xk(k)
cdbg        enddo

      !Converting gas unit from ppt to kg
cdbg      if (icond_test .eq. 1) then
cdbg        ygas(mgsvi)=1.0E+3 ! [=] ppt
cdbg        gmw(srtso4)=100. ! [=]g mol-1
cdbg        Gc(srtso4)=boxmass*ygas(mgsvi)*1.0e-12*gmw(srtso4)/28.9
cdbg      else
cdbg        Gc(srtso4)=boxmass*ygas(mgsvi)*1.0e-12*gmw(srtso4)/28.9
cdbg        Gc(srtnh3)=boxmass*ygas(mgnh3)*1.0e-12*gmw(srtnh3)/28.9
cdbg      endif

      !Skip negative gas concentration if command
cnogas      if (Gci(srtso4) .lt. cthresh*boxmass) Gcflag(srtso4)=1 ! =1, Skip 
cnogas      if (Gci(srtnh3) .lt. cthresh*boxmass) Gcflag(srtnh3)=1 ! =1, Skip

      !If PSSA is on, turn on Gcflag(srtso4) for ezcond does H2SO4 condensatoin.
cinclude_ammonia      Gcflag(srtnh3)=1 !Skip an ammonia condensation

C Repeat from this point if multiple internal time steps are needed
 10   continue

      !Set dp equals to zero
cnogas      do k=1,ibins
cnogas         do j=1,icomp-1
cnogas            dp(k,j)=0.0
cnogas         enddo
cnogas      enddo

C Set dp for nonvolatile species
cnogas      do k=1,ibins
cnogas        if (icond_test .eq. 1) then
cnogas          dp(k,srtso4)=(Gcf(srtso4)/100.)/(boxmass/28.9)*pres
cnogas        else
cnogas          dp(k,srtso4)=(Gcf(srtso4)/98.)/(boxmass/28.9)*pres
cnogas        endif
cnogas      enddo

C Ammonia Condensation Strategy by jgj
      !If ammonium concentration is less than 0.375 times of sulfate,
      !then pressure difference, dp is calculated. Where 0.375 is equal 
      !to mass ratio of ammonium and sulfate of ammonium sulfate.
cnonh3      do k=1,ibins
cnonh3        if (Mkf(k,srtso4) .gt. 0.) then
cnonh3          ratio=Mkf(k,srtnh3)/Mkf(k,srtso4)
cnonh3        else
cnonh3          ratio=0.375 ! skip calculating of dp(k,srtnh3)
cnonh3        endif
cnonh3        if (ratio .lt. 0.375) then ! 0.375 is max ratio
cnonh3          dp(k,srtnh3)=(Gcf(srtnh3)/molwt(srtnh3))/(boxmass/28.9)*pres
cnonh3        endif
cnonh3      enddo

C Calculate tj and tk factors needed to calculate tau values
cnogas      mu=2.5277e-7*temp**0.75302
cnogas      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temp)))  !S&P eqn 8.6
cnogas      do j=1,icomp-1
cnogas        if (icond_test .eq. 1) then
cnogas          Di=0.1*1.0d-4 ! m^2/s ! cond_test 6/12/04 jgj
cnogas        else
cnogas          Di=gasdiff(temp,pres,gmw(j),Sv(j))
cnogas        endif
cnogas        tj(j)=2.*pi*Di*molwt(j)*1.0d-3/(R*temp)
cnogas      enddo

cnogas      sK(srtso4)=0.0d0
cnogas      sK(srtnh3)=0.0d0

cnogas      do k=1,ibins
cnogas        if (Nkf(k) .gt. Neps) then
cnogas          if (icond_test .eq. 1) then
cnogas            density=1400.0 ! Mk's are calculated based on this density
                           ! in initbounds. Dpk can be derived from this
                           ! density.
cnogas          else
cnogas            density=1400.0 ! [=] kg/m3
cnogas          endif
cnogas          mp=(Mkf(k,srtso4)+Mkf(k,srtnh3))/Nkf(k)
cnogas        else
          !nothing in this bin - set to "typical value"
cnogas          if (icond_test .eq. 1) then !cond test 6/24/04 jgj
cdbg            density = 1000.
cnogas            density=1400.0 ! Mk's are calculated based on this density
                           ! in initbounds. Dpk can be derived from this
                           ! density.
cnogas          else
cnogas            density=1400. ! was 1500. jgj 6/28/04
cnogas          endif
cnogas          mp=1.414*xk(k)
cnogas        endif
cnogas        Dpk(k)=((mp/density)*(6./pi))**(0.333)
cnogas        Dpk(k)=h2ogrowth*Dpk(k)
cnogas        if (icond_test .eq. 1) then
cnogas          beta(srtso4)=1.
cnogas        else
cnogas          Kn=2.0*mfp/Dpk(k)                             !S&P eqn 11.35 (text)
cnogas          beta(srtso4)=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha(srtso4)) 
                                ! S&P eqn 11.35, Table 11.1 Dahneke
cnogas          beta(srtnh3)=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha(srtnh3))
cnogas        endif
cnogas        tk(k,srtso4)=(6./(pi*density))**(1./3.)*beta(srtso4)
cnogas        tk(k,srtnh3)=(6./(pi*density))**(1./3.)*beta(srtnh3)
cnogas        if (Nkf(k) .gt. 0.0) then
cnogas          sK(srtso4)=sK(srtso4)+tk(k,srtso4)*Nkf(k)*(Mkf(k,srtso4)
cnogas     &             /Nkf(k))**(1.d0/3.d0)
cnogas          sK(srtnh3)=sK(srtnh3)+tk(k,srtnh3)*Nkf(k)*(Mkf(k,srtnh3)
cnogas     &             /Nkf(k))**(1.d0/3.d0)
cnogas        endif
cnogas      enddo

cnogas      sK(srtso4)=sK(srtso4)*zeta13*tj(srtso4)*R*temp/(molwt(srtso4)
cnogas     &      *1.d-3)/(boxvol*1.d-6)
cnogas      sK(srtnh3)=sK(srtnh3)*zeta13*tj(srtnh3)*R*temp/(molwt(srtso4)
cnogas     &      *1.d-3)/(boxvol*1.d-6)

C Choose appropriate time step

      !Try to take full time step
      cdt=dt-time
cdbg      limit='complete'

      !Not more than 15 minutes
      if (cdt .gt. 900.) then
        cdt=900.
cdbg      limit='max'
      endif

 20   continue   !if time step is shortened, repeat from here

C-----Calculate tau values for all species/bins

      !Set atauc and atau equal to zero
      do k=1,ibins
         do j=1,icomp-1
            atauc(k,j)=0.0  
            atau(k,j)=0.0
         enddo

        !atauc is tau that is a parameter to describe condensation
        !driving force when pressure difference is equal to constant.
        !atau is same to atauc except having a exponential decaying
        !pressure term. The current ammonia condensation uses atauc as atau.
        !In case atauc is bigger than taumax that atau is forced to be taumax.
 
        !Calculate a driving force for ammonia condensation
cnogas        atauc(k,srtso4)=tj(srtso4)*tk(k,srtso4)*dp(k,srtso4)*cdt
cnogas        atauc(k,srtnh3)=tj(srtnh3)*tk(k,srtnh3)*dp(k,srtnh3)*cdt
cnogas        Mknh3max=0.375*Mkf(k,srtso4) ! maximally allowable NH4 mass
cdbg        if (Nkf(k) .gt. 0.) then
cnogas        if (Nkf(k) .gt. Neps) then
cnogas          taumax=1.5*((Mknh3max**tdt)-(Mkf(k,srtnh3)**tdt))/(Nkf(k)
cnogas     &     **tdt) 
            ! See Eq.(9) in Adams and Seinfeld (2002, JGR)
cnogas        else
cnogas          taumax=0.
cnogas        endif
cnogas        if (atauc(k,srtnh3) .gt. taumax) then
cnogas          if (taumax .ge. 0.) then
cnogas            atau(k,srtnh3)=taumax
cnogas          else
cnogas            atau(k,srtnh3)=0.  
cnogas          endif
cnogas        else
cnogas          atau(k,srtnh3)=atauc(k,srtnh3)
cnogas        endif

        !Calculate a driving force for sulfuric acid condensation
cnogas        if (sK(srtso4) .gt. 0.0) then
cnogas          atau(k,srtso4)=tj(srtso4)*R*temp/(molwt(srtso4)*1.d-3)
cnogas     &      /(boxvol*1.d-6)*tk(k,srtso4)*Gcf(srtso4)/sK(srtso4)
cnogas     &      *(1.d0-exp(-1.d0*sK(srtso4)*cdt))
cnogas        else
cnogas          atau(k,srtso4)=0.0  !nothing to condense onto
cnogas        endif
      enddo

      !Update atau considering the amount added by aqueous chemistry
      do k=iact,ibins
         do j=1,icomp-1
            Mko(k,j)=Mkf(k,j)+(moxid(k,j)*(1./Nkf(k)))
            atau(k,j)=1.5*((Mko(k,j)**tdt)-(Mkf(k,j)**tdt)) 
         enddo
      enddo
      do k=1,ibins
         do j=1,icomp-1
            Mko(k,j)=Mkf(k,j)
         enddo
      enddo


C-----Adjust a time step 

      tr=1.0 !The following sections limit the condensation time step
             !when necessary.  tr is a factor that describes by
             !how much to reduce the time step.
cdbg      nsteps = 1 ! number of steps for condensation

      !Make sure masses of individual species don't change too much
cdbg      do j=1,icomp-1
cdbg        if(Gcflag(j).eq.1) goto 30 ! If gas concentration is tiny, then skip 
                                   ! condensation. 12/06/07 jgj
c_include_inert        if(j.eq.srtorg) goto 30 ! OM does not have condensation process.
cdbg        do k=1,ibins
cdbg          if (Nkf(k) .gt. Neps) then
cdbg            mc=0.0
cdbg            do jj=1,icomp-1
cdbg              mc=mc+Mkf(k,jj)/Nkf(k)
cdbg            enddo
cdbg            if (mc/xk(k) .gt. 1.0d-3) then
cdbg              !species has significant mass in particle - limit change
c
c In aqueous chemistry, it is frequent that mass increases more than 10%.
c This is treated as nsteps as ezcond, so it moved to the loop included
c tmcond. jgj 2/20/08
c
cdbg              if (abs(atau(k,j))/(mc**tdt) .gt. 0.1) then
cdbg                floor=int(abs(atau(k,j))/(mc**tdt)/0.05)
cdbg                ceil=floor + 1
cdbg                nsteps=max(nsteps,ceil)
cmass                ttr=abs(atau(k,j))/(mc**tdt)/0.05
cmass                if (ttr. gt. tr) then 
cmass                  limit='amass'
cmass                  write(limit(7:11),'(I2,X,I2)') k,j
cmass                  tr=ttr
cmass                endif
cdbg              endif
cdbg            else
cdbg              !species is new to particle - set max time step
cdbg              if ((cdt/tr .gt. 0.1) .and. (atau(k,j).gt. 0.0)) then 
cdbg                tr=cdt/0.1
cdbg                  limit='nspec'
cdbg                  write(limit(7:11),'(I2,X,I2)') k,j
cdbg              endif
cdbg            endif
cdbg          endif
cdbg        enddo
cdbg        !Make sure gas phase concentrations don't change too much
cdbg        !control gas reduction rate on 12/15/07, jgj
cdbg        if (dt.le.120.) then
cdbg           gasfrac=0.25 ! not allow less than 75% of reduction theoretically
cdbg        else
cdbg           if (dt.le.600.) then
cdbg             gasfrac=0.989 ! not allow less than 1.1% of reduction exoponentally
cdbg           elseif (dt.le.700.) then
cdbg             gasfrac=0.991 ! not allow less than 0.9% of reduction exponentially
cdbg           elseif (dt.le.800.) then
cdbg             gasfrac=0.993 ! not allow less than 0.7% of reduction exponentially
cdbg           elseif (dt.le.850.) then
cdbg             gasfrac=0.995 ! not allow less than 0.5% of reduction exponentially
cdbg           elseif (dt.lt.900.) then
cdbg             gasfrac=0.998 ! not allow less than 0.2% of reduction exponentially
cdbg           else
cdbg             gasfrac=0.999 ! not allow less than 0.1% of reduction theoretically
cdbg           endif
cbdg        endif
cdbg        if (j .eq. srtso4) then ! Deal only sulfuric acid by jgj
cdbg          if (exp(-1.d0*sK(srtso4)*cdt) .lt. 0.25) then
cdbg            ttr=-2.d0*cdt*sK(j)/log(0.25)
cnogas        if (exp(-1.d0*sK(j)*cdt) .lt. gasfrac) then
cnogas          ttr=-2.d0*cdt*sK(j)/log(gasfrac)
cnogas          if (ttr .gt. tr) then 
cnogas            tr=ttr
cdbg              limit='gphas'
cdbg              write(limit(7:8),'(I2)') j
cnogas          endif
cnogas        endif
cdbg        endif
cdbg 30   continue
cdbg      enddo

      !Never shorten timestep by less than half
cdbg      if (tr .gt. 1.0) tr=max(tr,2.0)

      !Repeat for shorter time step if necessary
cdbg      if (tr .gt. 1.0) then
cdbg         cdt=cdt/tr
cdbg         !Iteration timestep
cdbg         itrts=itrts+1
cdbg         if (itrts.gt.5000) then
cdbg            write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
cdbg            write(*,*) 'Iterations of timestep in so4cond exceed 5000'
cdbg            write(*,*) 'dt=',dt,'time=',time,'cdt=',cdt
cnogas            write(*,*) 'exponential decaying frac=',
cnogas     &               exp(-1.d0*sK(srtnh4)*cdt)
cdbg            write(*,*) 'limit=',limit
cdbg            write(*,*) 'moxd sulfate produced from aqueous chemistry'
cdbg            do j=1, icomp-1
cdbg              write(*,*)'j=',j
cdbg              do k=1, ibins
cdbg                write(*,*) moxid(k,j),moxid(k,j)/cvt/boxvol 
cdbg              enddo
cdbg            enddo
cdbg            write(*,*)'atau='
cdbg            do j=1, icomp-1
cdbg              write(*,*)'j=',j
cdbg              do k=1, ibins
cdbg                write(*,*) atau(k,j) 
cdbg              enddo
cdbg            enddo
cdbg            write(*,*)'Nk='
cdbg            do k=1,ibins
cdbg               write(*,*)Nk(k)
cdbg            enddo
cdbg            write(*,*)'Mk='
cdbg            do j=1,icomp
cdbg               write(*,*)'j=',j
cdbg               do k=1,ibins
cdbg                  write(*,*)Mk(k,j)
cdbg               enddo
cdbg            enddo
cdbg            STOP
cdbg         endif
cdbg         goto 20
cdbg      endif
      
C Call condensation subroutine to do mass transfer

      do j=1,icomp-1  !Loop over all aerosol components
c        if(Gcflag(j).eq.1) goto 40 ! If gas concentration is tiny, then skip 
                                   ! condensation. 12/06/07 jgj
c_include_inert        if(j.eq.srtorg) goto 40 ! OM does not have condensation process.
        
        !Swap tau values for this species into array for cond
        do k=1,ibins
          if (icond_test .eq. 1) then
            tau(k)=atauc(k,j) ! for cond_test 6/24/04 jgj
          else
            tau(k)=atau(k,j)
          endif
          tau(k)=corfactor*tau(k) ! A correction factor is applied.
        enddo

cdbg        do insteps=1,nsteps ! do steps of condensation

        call mnfix_PSSA(Nkf,Mkf,ichm,jchm,kchm)
            ! adjust average mass in each size bin within boundaries

        !Call condensation routine
cneg_numb_cons        Ntotf=0.0
cneg_numb_cons        do k=1,ibins
cneg_numb_cons          Ntotf=Ntotf+Nkf(k)
cneg_numb_cons        enddo

        !oxidated mass calculation, dummy, not using it, jgj
        do k=1, ibins
          moxd(k)=0.0 !oxidated mass
        enddo
 
cdbg        print*,'j=',j
cdbg        print*,'tau=',tau
cdbg        print*,'Gc=',Gcf
cdbg        print*,'tk=',tk
cdbg        print*,'tj=',tj
cdbg        print*,'sK=',sK
cdbg        print*,'temp=',temp
cdbg        print*,'boxvol=',boxvol
cdbg        print*,'cdt=',cdt
cdbg        print*,'xk=',xk
cdbg        print*,'Nk='
cdbg        pause

        call tmcond(tau,xk,Mkf,Nkf,Mko,Nko,j,moxd)

cdbg        if ((ichm.eq.53).and.(jchm.eq.19).and.(kchm.eq.1)) then
cdbg           write(*,*)'Coord.=',ichm,jchm,kchm
cdbg           write(*,*)'addtau, tau='
cdbg           do k=1,ibins
cdbg              write(*,*)addtau(k),tau(k)
cdbg           enddo
cdbg           write(*,*)'Nko='
cdbg           do k=1,ibins
cdbg              write(*,*)Nko(k)
cdbg           enddo
cdbg           write(*,*)'j=',j
cdbg           do k=1,ibins
cdbg              write(*,*)Mko(k,j)
cdbg           enddo
cdbg        endif

        !Check for number conservation
cneg_numb        Ntoto=0.0
cneg_numb        do k=1,ibins
cneg_numb          Ntoto=Ntoto+Nko(k)
cneg_numb        enddo
cdbg         write(*,*) 'Time=', time
cdbg         write(*,*) 'Ntoto=', Ntoto
cdbg         write(*,*) 'Ntotf=', Ntotf
cdbg         dNerr=dNerr+Ntotf-Ntoto
cneg_numb        dNerr=Ntotf-Ntoto
cneg_numb        if (abs(dNerr/Ntoto) .gt. 1.d-4) then
cneg_numb          write(*,*)'ERROR in so4cond: Number not conserved'
cneg_numb          write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
cneg_numb          write(*,*)'Ntoto, Ntotf, dNerr/Ntoto'
cneg_numb     &             ,Ntoto, Ntotf, dNerr/Ntoto
cneg_numb          if (abs(dNerr/Ntoto) .gt. 1.d-2) then
cneg_numb            write(*,*)'Serious Error in so4cond: Number not conserved 
cneg_numb&less than 1 %'
cneg_numb            write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
cneg_numb            write(*,*) 'moxd sulfate produced from aqueous chemistry'
cneg_numb            do k=1, ibins
cneg_numb              write(*,*) moxid(k),moxid(k)/cvt/boxvol 
cneg_numb            enddo
cneg_numb            write(*,*) 'Nko='
cneg_numb            do k=1, ibins
cneg_numb              write(*,*) Nko(k)
cneg_numb            enddo
cneg_numb            write(*,*) 'Mko='
cneg_numb            do jj=1,icomp
cneg_numb              write(*,*) 'jj=',jj
cneg_numb              do k=1, ibins
cneg_numb                write(*,*) Mko(k,jj)
cneg_numb              enddo
cneg_numb            enddo
cneg_numb            STOP
cneg_numb          endif
cneg_numb        endif

        !Update gas phase concentration
cnogas        mi(j)=0.0
cnogas        mf(j)=0.0
cnogas        do k=1,ibins
cnogas          mi(j)=mi(j)+Mkf(k,j)
cnogas          mf(j)=mf(j)+Mko(k,j)
cnogas        enddo
cnogas        Gcf(j)=Gcf(j)+(mi(j)-mf(j))*gmw(j)/molwt(j)

        !Swap into Nk, Mk
        do k=1,ibins
          Nkf(k)=Nko(k)
          do jj=1,icomp-1
            Mkf(k,jj)=Mko(k,jj)
          enddo
        enddo

        !Update water concentrations
c        call ezwatereqm(Mkf)
cdbg      enddo !nsteps
 40   continue
      enddo

C Update time
      time=time+cdt
cdbg      write(*,*) 'so4cond - time = ', time, ' ',limit
cdbg      write(*,*) 'H2SO4(g)= ', Gcf(srtso4)
 
      !Check sulfuric acid
cnogas      if (Gcflag(srtso4) .eq. 0) then
cnogas        if (Gcf(srtso4) .lt. 0.0) then
cnogas          if (abs(Gcf(srtso4)) .gt. 1.d-5) then
            !Gcf is substantially less than zero - this is a problem
cnogas            write(*,*) 'ERROR in so4cond: H2SO4(g) < 0'
cnogas            write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
cnogas            write(*,*) 'time=',time
cnogas            write(*,*) 'Gcf()=',Gcf(srtso4)
cnogas            write(*,*) 'H2SO4 [=] ppt',Gcf(srtso4)/boxmass/1.0e-12
cnogas     &               /gmw(srtso4)*28.9
cnogas            write(*,*) 'Gass consumed in kg',(mi(srtso4)-mf(srtso4))
cnogas     &               *gmw(srtso4)/molwt(srtso4)
cnogas            write(*,*) 'mi,mf,gmw(j),molwt(j)=',mi(srtso4),mf(srtso4)
cnogas     &               ,gmw(srtso4),molwt(srtso4)
cnogas            write(*,*) 'time=',time,'cdt=',cdt
cnogas            write(*,*) 'exponential decaying frac=',
cnogas     &                exp(-1.d0*sK(srtso4)*cdt)
cnogas            STOP
cnogas          else
            !Gcf is negligibly less than zero - probably roundoff error
cnogas            Gcf(srtso4)=0.0
cnogas          endif
cnogas        endif
cnogas      endif

      !Check ammonia gas
cnogas      if (Gcflag(srtnh4) .eq. 0) then
cnogas        if (Gcf(srtnh4) .lt. 0.0) then
cnogas          if (abs(Gcf(srtnh4)) .gt. 1.d-5) then
            !Gcf is substantially less than zero - this is a problem
cnogas            write(*,*) 'ERROR in so4cond: NH3(g) < 0'
cnogas            write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
cnogas            write(*,*) 'time=',time
cnogas            write(*,*) 'Gcf()=',Gcf(srtnh4)
cnogas            write(*,*) 'NH3 [=] ppt',Gcf(srtnh4)/boxmass/1.0d-12
cnogas     &               /gmw(srtnh4)*28.9
cnogas            write(*,*) 'Gass consumed in kg',(mi(srtnh4)-mf(srtnh4))
cnogas     &               *gmw(srtnh4)/molwt(srtnh4)
cnogas            write(*,*) 'time=',time,'cdt=',cdt
cnogas            write(*,*) 'exponential decaying frac=',
cnogas     &               exp(-1.d0*sK(srtnh4)*cdt)
cnogas            STOP
cnogas          else
            !Gcf is negligibly less than zero - probably roundoff error
cnogas            Gcf(srtnh4)=0.0
cnogas          endif
cnogas        endif
cnogas      endif

C Repeat process if necessary
cdbg      if (time .lt. dt) then
        !Iteration
cdbg        itr=itr+1
cdbg        if (itr.gt.5000) then
cdbg           write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
cdbg           write(*,*) 'An iteration in so4cond exceeds 5000'
cdbg           write(*,*) 'dt=',dt,'time=',time,'cdt=',cdt
cdbg           write(*,*) 'exponential decaying frac=',
cdbg     &               exp(-1.d0*sK(srtnh4)*cdt)
cdbg           STOP
cdbg        endif
cdbg        goto 10
cdbg      endif

cdbg      write(*,*) 'Number cons. error was ', dNerr

 100  continue   !skip to here if there is no gas phase to condense
c
cdbg      if (icond_test .eq. 1) then
cdbg          ygas(mgsvi)=Gcf(srtso4)*1.0e+12/boxmass*28.9/100.
                 !const H2SO4 in case of condensation tests
cdbg      else
cdbg          ygas(mgsvi)=Gcf(srtso4)*1.0e+12/boxmass*28.9/gmw(srtso4)
cdbg          ygas(mgnh3)=Gcf(srtnh3)*1.0e+12/boxmass*28.9/gmw(srtnh3) 
                                   !eznh3eqm in PSSA could substitue
cdbg      endif

      ! Adjust finally before return
      call mnfix_PSSA(Nkf,Mkf,ichm,jchm,kchm)

      RETURN
      END

