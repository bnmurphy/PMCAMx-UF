
C     **************************************************
C     *  so4cond                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, June 2000
C     Based on thermocond.f

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
      SUBROUTINE so4cond(Nki,Mki,Gci,Nkf,Mkf,Gcf,dt,ichm,jchm,kchm,
     &           iflagez)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'
      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

cdbg      real ygas(ngas) ! gas species mixing ratio
      double precision Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      double precision Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      real dt ! timestep
cdbg      double precision xkDMAN(ibins+1) ! xk as inputs
      integer ichm ! i coordinate in PMCAMx
      integer jchm ! j coordinate in PMCAMx
      integer kchm ! k coordinate in PMCAMx
      integer iflagez ! If iteration exceeds more than thousand
                      ! , return with iflagez = 1.

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer j,k,jj,kk        !counters
      integer Gcflag(icomp)    !flags for checking insignificant Gc
      integer itr, itrts       !iterations for master loop and timestep

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
      double precision Ntotf, Ntoto, dNerr  !used to track number cons.

      double precision Mknh3max !Maximum allowable NH3
      double precision taumax   !Maxmum tau for ammonia condensation
      double precision moxd(ibins) !moxid/Nact
      double precision Mktot(icomp) ! Total mass of each species
      double precision Mtot     !Total mass of each size bin
      double precision Gcknh3(ibins) ! Fractional ammonia gas to be condensed
                                     ! when ammonia is limited.
      double precision sumataunh3    ! sum of atauc

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
      parameter(corfactor=1.4) ! a correction factor

C-----CODE--------------------------------------------------------------

cdbg      write(*,*)'xk(j)',(xk(j),j=1,ibins)
cdbg      print*,'boxmass=',boxmass
cdbg      pause
cdbg      write(*,*) '\n Coord.(i,j,k)=',ichm,jchm,kchm

      !Initialization
      k=0
      jj=0
      iflagez=0
cdbg      limit='null'
      dNerr=0.0
      time=0.0         !subroutine exits when time=dt
      tdt=2.d0/3.d0
      itr=0
      itrts=0
      do j=1,icomp-1
        Gcf(j)=Gci(j)
      enddo
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
cdbg      do k=1,ibins+1
cdbg        xk(k)=xkDMAN(k)
cdbg      enddo

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
      if (Gci(srtso4) .lt. cthresh*boxmass) Gcflag(srtso4)=1 ! =1, Skip 
      if (Gci(srtnh3) .lt. cthresh*boxmass) Gcflag(srtnh3)=1 ! =1, Skip

      !If PSSA is on, turn on Gcflag(srtso4) for ezcond does H2SO4 condensatoin.
      Gcflag(srtso4)=1 !PSSA

C Repeat from this point if multiple internal time steps are needed
 10   continue

      !Set dp equals to zero
      do k=1,ibins
         do j=1,icomp-1
            dp(k,j)=0.0
         enddo
      enddo

C Set dp for nonvolatile species
      do k=1,ibins
        if (icond_test .eq. 1) then
          dp(k,srtso4)=(Gcf(srtso4)/100.)/(boxmass/28.9)*pres
        else
          dp(k,srtso4)=(Gcf(srtso4)/98.)/(boxmass/28.9)*pres
        endif
      enddo

C Ammonia Condensation Strategy by jgj
      !If ammonium concentration is less than 0.375 times of sulfate,
      !then pressure difference, dp is calculated. Where 0.375 is equal 
      !to mass ratio of ammonium and sulfate of ammonium sulfate.
      do k=1,ibins
        if (Mkf(k,srtso4) .gt. 0.) then
          ratio=Mkf(k,srtnh3)/Mkf(k,srtso4)
        else
          ratio=0.375 ! skip calculating of dp(k,srtnh3)
        endif
        if (ratio .lt. 0.375) then ! 0.375 is max ratio
          dp(k,srtnh3)=(Gcf(srtnh3)/molwt(srtnh3))/(boxmass/28.9)*pres
        endif
      enddo

C Calculate tj and tk factors needed to calculate tau values
      mu=2.5277e-7*temp**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temp)))  !S&P eqn 8.6
      do j=1,icomp-1
        if (icond_test .eq. 1) then
          Di=0.1*1.0d-4 ! m^2/s ! cond_test 6/12/04 jgj
        else
          Di=gasdiff(temp,pres,gmw(j),Sv(j))
        endif
        tj(j)=2.*pi*Di*molwt(j)*1.0d-3/(R*temp)
      enddo

      sK(srtso4)=0.0d0
      sK(srtnh3)=0.0d0

      do k=1,ibins
        if (Nkf(k) .gt. Neps) then
          if (icond_test .eq. 1) then
cdbg            density=1.0e+3 ! cond_test 6/24/04 jgj
            density=1400.0 ! Mk's are calculated based on this density
                           ! in initbounds. Dpk can be derived from this
                           ! density.
          else
cdbg             density=aerodens_PSSA(Mkf(k,srtso4),0.0,Mkf(k,srtso4),
cdbg     &           0.0,Mkf(k,srtnh3)) 
            density=1400.0 ! [=] kg/m3
          endif
c          mp=(Mkf(k,srtso4)+Mkf(k,srtnh3))/Nkf(k)
          ! JJ bugfix: Nkf refers to total particle number so must
          ! include mass of all particulate species (except water) or we end up with mp below
          ! the mass bin lower bound. Although why not take the geometric mean
          ! here as well?
          mp=(Mkf(k,srtso4)+Mkf(k,srtsoa1)+Mkf(k,srtsoa2)+
     &     +Mkf(k,srtsoa3)+Mkf(k,srtsoa4)
     &     +Mkf(k,srtnh3)+Mkf(k,srtinrt))/Nkf(k)
        else
          !nothing in this bin - set to "typical value"
          if (icond_test .eq. 1) then !cond test 6/24/04 jgj
cdbg            density = 1000.
            density=1400.0 ! Mk's are calculated based on this density
                           ! in initbounds. Dpk can be derived from this
                           ! density.
          else
            density=1400. ! was 1500. jgj 6/28/04
          endif
          mp=1.414*xk(k)
        endif
        Dpk(k)=((mp/density)*(6./pi))**(1./3.)
        Dpk(k)=h2ogrowth*Dpk(k)
        if (icond_test .eq. 1) then
          beta(srtso4)=1.
        else
          Kn=2.0*mfp/Dpk(k)                             !S&P eqn 11.35 (text)
          beta(srtso4)=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha(srtso4)) 
                                ! S&P eqn 11.35, Table 11.1 Dahneke
cdbg           alphanh3=alpha_nh3(Mkf(srtnh3),Mkf(srtso4),rh) ! Pathak et al.
cdbg           beta(srtnh3)=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alphanh3) ! method
          beta(srtnh3)=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha(srtnh3))
        endif
        tk(k,srtso4)=(6./(pi*density))**(1./3.)*beta(srtso4)
        tk(k,srtnh3)=(6./(pi*density))**(1./3.)*beta(srtnh3)
        if (Nkf(k) .gt. 0.0) then
          Mtot=0.0
          do jj=1, icomp
            Mtot=Mtot+Mkf(k,jj)
          enddo
          sK(srtso4)=sK(srtso4)+tk(k,srtso4)*Nkf(k)*(Mtot
     &             /Nkf(k))**(1.d0/3.d0)
          sK(srtnh3)=sK(srtnh3)+tk(k,srtnh3)*Nkf(k)*(Mtot
     &             /Nkf(k))**(1.d0/3.d0)
        endif
      enddo
      sK(srtso4)=sK(srtso4)*zeta13*tj(srtso4)*R*temp/(molwt(srtso4)
     &      *1.d-3)/(boxvol*1.d-6)
      sK(srtnh3)=sK(srtnh3)*zeta13*tj(srtnh3)*R*temp/(molwt(srtnh3)
     &      *1.d-3)/(boxvol*1.d-6)

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
        atauc(k,srtso4)=tj(srtso4)*tk(k,srtso4)*dp(k,srtso4)*cdt
        atauc(k,srtnh3)=tj(srtnh3)*tk(k,srtnh3)*dp(k,srtnh3)*cdt

        !Calculate a driving force for sulfuric acid condensation
        if (sK(srtso4) .gt. 0.0) then
          atau(k,srtso4)=tj(srtso4)*R*temp/(molwt(srtso4)*1.d-3)
     &      /(boxvol*1.d-6)*tk(k,srtso4)*Gcf(srtso4)/sK(srtso4)
     &      *(1.d0-exp(-1.d0*sK(srtso4)*cdt))
        else
          atau(k,srtso4)=0.0  !nothing to condense onto
        endif

        ! JJ bugfix: since atau is the one that is being used for the 
        ! driving force we need to set that to zero when there is already
        ! enough nh4 to neutralize the so4. This check was done above with
        ! dp, but that is used to calculate atauc, not atau!

        !Calculate a driving force for ammonia condensation
        if (sK(srtnh3) .gt. 0.0 .and. atauc(k,srtnh3).gt.0.D0) then
          atau(k,srtnh3)=tj(srtnh3)*R*temp/(molwt(srtnh3)*1.d-3)
     &      /(boxvol*1.d-6)*tk(k,srtnh3)*Gcf(srtnh3)/sK(srtnh3)
     &      *(1.d0-exp(-1.d0*sK(srtnh3)*cdt))
        else
          atau(k,srtnh3)=0.0  !nothing to condense onto or already enough ammonia
        endif

        ! JJ bugfix: if we do not end the do loop over the bins here, the following 
        ! calculation with atau for sumtaunh3 is not using correct ataus for k>current step

      end do ! added by JJ

      ! it is needless to recalculate these Mktot:s within an outer loop so start next
      ! loop over k after these / JJ
      Mktot(srtso4)=0.0
      do kk=1,ibins
         Mktot(srtso4)=Mktot(srtso4)+Mkf(kk,srtso4)
      enddo
      Mktot(srtnh4)=0.0
      do kk=1,ibins
         Mktot(srtnh4)=Mktot(srtnh4)+Mkf(kk,srtnh4)
      enddo

      do k=1,ibins ! added by JJ: now we can loop over the bins and do the next part
         ! JJ change: instead of comparing total so4 to total ammonia in order to get
         ! taumax for the current bin, check the so4 and ammonia in the current bin.
         ! For this purpose calculate Gcknh3 already outside of the if statement
         sumataunh3=0.0
         do kk=1,ibins
            sumataunh3=sumataunh3+atau(kk,srtnh4)
         enddo
         Gcknh3(k)=Gcf(srtnh4)*atau(k,srtnh4)/sumataunh3
         
        !Separate the cases of total ammonia is greater than existing sulfate
        ! or not.
c        if (0.375*Mktot(srtso4).gt.(Mktot(srtnh4)+Gcf(srtnh4))) then
         if (0.375*Mkf(k,srtso4).gt.Mkf(k,srtnh4)+Gcknh3(k)) then !JJ change
        !Sulfate rich condition
           
            Mknh3max=Mkf(k,srtnh3)+Gcknh3(k)
                                         ! maximally allowable NH4 mass
         else                   !Total ammonia rich condition
            Mknh3max=0.375*Mkf(k,srtso4) ! maximally allowable NH4 mass
         endif
cdbg        if (Nkf(k) .gt. 0.) then
         if (Nkf(k) .gt. Neps) then
            taumax=1.5*((Mknh3max**tdt)-(Mkf(k,srtnh3)**tdt))/(Nkf(k)
     &           **tdt) 
            ! See Eq.(9) in Adams and Seinfeld (2002, JGR)
         else
            taumax=0.           !For safety
         endif
         ! JJ change: since atau is the one that is used for tmcond, there seems
         ! to be no reason to compare taumax to atauc here
c         if (atauc(k,srtnh3) .gt. taumax) then
         if (atau(k,srtnh3) .gt. taumax) then
            if (taumax .ge. 0.) then
               atau(k,srtnh3)=taumax
            else
               atau(k,srtnh3)=0. !For safety 
            endif
         endif
      enddo

C-----Adjust a time step 

      tr=1.0 !The following sections limit the condensation time step
             !when necessary.  tr is a factor that describes by
             !how much to reduce the time step.

      !Make sure masses of individual species don't change too much
      do j=1,icomp-1
        if(Gcflag(j).eq.1) goto 30 ! If gas concentration is tiny, then skip 
                                   ! condensation. 12/06/07 jgj
        if(j.eq.srtorg) goto 30 ! OM does not have condensation process.
        do k=1,ibins
          if (Nkf(k) .gt. Neps) then
            mc=0.0
            do jj=1,icomp-1
              mc=mc+Mkf(k,jj)/Nkf(k)
            enddo
            if (mc/xk(k) .gt. 1.0d-3) then
              !species has significant mass in particle - limit change
              if (abs(atau(k,j))/(mc**tdt) .gt. 0.1) then
                ttr=abs(atau(k,j))/(mc**tdt)/0.05
                if (ttr. gt. tr) then 
cdbg                  limit='amass'
cdbg                  write(limit(7:11),'(I2,X,I2)') k,j
                  tr=ttr
                endif
              endif
            else
              !species is new to particle - set max time step
              if ((cdt/tr .gt. 0.1) .and. (atau(k,j).gt. 0.0)) then 
                tr=cdt/0.1
cdbg                  limit='nspec'
cdbg                  write(limit(7:11),'(I2,X,I2)') k,j
              endif
            endif
          endif
        enddo
        !Make sure gas phase concentrations don't change too much
        !control gas reduction rate on 12/15/07, jgj
        if (dt.le.120.) then
           gasfrac=0.25 ! not allow less than 75% of reduction theoretically
        else
cdbg           if (dt.le.600.) then
cdbg             gasfrac=0.989 ! not allow less than 1.1% of reduction exoponentally
cdbg           elseif (dt.le.700.) then
cdbg             gasfrac=0.991 ! not allow less than 0.9% of reduction exponentially
cdbg           elseif (dt.le.800.) then
cdbg             gasfrac=0.993 ! not allow less than 0.7% of reduction exponentially
cdbg           elseif (dt.le.850.) then
cdbg             gasfrac=0.995 ! not allow less than 0.5% of reduction exponentially
cdbg           elseif (dt.lt.900.) then
           if (dt.lt.900.) then
cdbg             gasfrac=0.998 ! not allow less than 0.2% of reduction exponentially
             gasfrac=0.50 ! not allow less than 10% of reduction exponentially
           else
cdbg             gasfrac=0.999 ! not allow less than 0.1% of reduction theoretically
             gasfrac=0.55 ! not allow less than 5% of reduction theoretically
           endif
        endif
cdbg        if (j .eq. srtso4) then ! Deal only sulfuric acid by jgj
cdbg          if (exp(-1.d0*sK(srtso4)*cdt) .lt. 0.25) then
cdbg            ttr=-2.d0*cdt*sK(j)/log(0.25)
        if (exp(-1.d0*sK(j)*cdt) .lt. gasfrac) then
          ttr=-2.d0*cdt*sK(j)/log(gasfrac)
          if (ttr .gt. tr) then 
            tr=ttr
cdbg              limit='gphas'
cdbg              write(limit(7:8),'(I2)') j
          endif
        endif
cdbg        endif
 30   continue
      enddo

! JJ comment: this is not actually working now, 2 has to be double because tr is double
! but left it for now because there might be good reason for this?
      !Never shorten timestep by less than half
c     changed by LA
      if (tr .gt. 1.0) tr=max(tr,2.0)
c      if (tr .gt. 1.0) tr=max(tr,2.D0)
c     end changed by LA

      !Repeat for shorter time step if necessary
      if (tr .gt. 1.0) then
         cdt=cdt/tr
         !Iteration timestep
         itrts=itrts+1
         if (itrts.gt.5000) then
            write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
            write(*,*) 'Iterations of timestep in so4cond exceed 5000'
            write(*,*) 'dt=',dt,'time=',time,'cdt=',cdt
            write(*,*) 'exponential decaying frac=',
     &               exp(-1.d0*sK(srtnh4)*cdt)
            STOP
         endif
         goto 20
      endif
      
C Call condensation subroutine to do mass transfer

      do j=1,icomp-1  !Loop over all aerosol components
        if(Gcflag(j).eq.1) goto 40 ! If gas concentration is tiny, then skip 
                                   ! condensation. 12/06/07 jgj
        if(j.eq.srtorg) goto 40 ! OM does not have condensation process.

        !Swap tau values for this species into array for cond
        do k=1,ibins
          if (icond_test .eq. 1) then
            tau(k)=atauc(k,j) ! for cond_test 6/24/04 jgj
          else
            tau(k)=atau(k,j)
          endif
          ! JJ comment: if atau and consequently tau got the taumax value earlier,
          ! a corfactor > 1 will try to condense more than should be allowed.
          ! Leave it for now
          tau(k)=corfactor*tau(k) ! A correction factor is applied.
        enddo

        call mnfix_PSSA(Nkf,Mkf,ichm,jchm,kchm)
            ! adjust average mass in each size bin within boundaries

        !Call condensation routine
        Ntotf=0.0
        do k=1,ibins
          Ntotf=Ntotf+Nkf(k)
        enddo

        !oxidated mass calculation
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
cdbg        do k=1, ibins
cdbg          print*,Nkf(k)
cdbg        enddo
cdbg        print*,'Mk='
cdbg        do k=1, ibins
cdbg          print*,Mkf(k,srtso4)
cdbg        enddo
cdbg        pause
        call tmcond(tau,xk,Mkf,Nkf,Mko,Nko,j,moxd)

        !Check for number conservation
        Ntoto=0.0
        do k=1,ibins
          Ntoto=Ntoto+Nko(k)
        enddo
cdbg         write(*,*) 'Time=', time
cdbg         write(*,*) 'Ntoto=', Ntoto
cdbg         write(*,*) 'Ntotf=', Ntotf
cdbg         dNerr=dNerr+Ntotf-Ntoto
        dNerr=Ntotf-Ntoto
        if (abs(dNerr/Ntoto) .gt. 1.d-4) then
          write(*,*)'ERROR in so4cond: Number not conserved'
          write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
          write(*,*)'Ntoto, Ntotf, dNerr/Ntoto'
     &             ,Ntoto, Ntotf, dNerr/Ntoto
          if (abs(dNerr/Ntoto) .gt. 1.d-2) then
            write(*,*)'Serious Error in so4cond: Number not conserved 
     &less than 1 %'
            write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
            STOP
          endif
        endif

        !Update gas phase concentration
        mi(j)=0.0
        mf(j)=0.0
        do k=1,ibins
          mi(j)=mi(j)+Mkf(k,j)
          mf(j)=mf(j)+Mko(k,j)
        enddo
        Gcf(j)=Gcf(j)+(mi(j)-mf(j))*gmw(j)/molwt(j)

        !Swap into Nk, Mk
        do k=1,ibins
          Nkf(k)=Nko(k)
          do jj=1,icomp-1
            Mkf(k,jj)=Mko(k,jj)
          enddo
        enddo

        !Update water concentrations
c        call ezwatereqm(Mkf)
 40   continue
      enddo

C Update time
      time=time+cdt
cdbg      write(*,*) 'so4cond - time = ', time, ' ',limit
cdbg      write(*,*) 'H2SO4(g)= ', Gcf(srtso4)
 
      !Check sulfuric acid
      if (Gcflag(srtso4) .eq. 0) then
        if (Gcf(srtso4) .lt. 0.0) then
          if (abs(Gcf(srtso4)) .gt. 1.d-5) then
            !Gcf is substantially less than zero - this is a problem
            write(*,*) 'ERROR in so4cond: H2SO4(g) < 0'
            write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
            write(*,*) 'time=',time
            write(*,*) 'Gcf()=',Gcf(srtso4)
            write(*,*) 'H2SO4 [=] ppt',Gcf(srtso4)/boxmass/1.0e-12
     &               /gmw(srtso4)*28.9
            write(*,*) 'Gass consumed in kg',(mi(srtso4)-mf(srtso4))
     &               *gmw(srtso4)/molwt(srtso4)
            write(*,*) 'mi,mf,gmw(j),molwt(j)=',mi(srtso4),mf(srtso4)
     &               ,gmw(srtso4),molwt(srtso4)
            write(*,*) 'time=',time,'cdt=',cdt
            write(*,*) 'exponential decaying frac=',
     &                exp(-1.d0*sK(srtso4)*cdt)
            STOP
          else
            !Gcf is negligibly less than zero - probably roundoff error
            Gcf(srtso4)=0.0
          endif
        endif
      endif

      !Check ammonia gas
      if (Gcflag(srtnh4) .eq. 0) then
        if (Gcf(srtnh4) .lt. 0.0) then
          if (abs(Gcf(srtnh4)) .gt. 1.d-5) then
            !Gcf is substantially less than zero - this is a problem
            write(*,*) 'ERROR in so4cond: NH3(g) < 0'
            write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
            write(*,*) 'time=',time
            write(*,*) 'Gcf()=',Gcf(srtnh4)
            write(*,*) 'NH3 [=] ppt',Gcf(srtnh4)/boxmass/1.0d-12
     &               /gmw(srtnh4)*28.9
            write(*,*) 'Gass consumed in kg',(mi(srtnh4)-mf(srtnh4))
     &               *gmw(srtnh4)/molwt(srtnh4)
            write(*,*) 'time=',time,'cdt=',cdt
            write(*,*) 'exponential decaying frac=',
     &               exp(-1.d0*sK(srtnh4)*cdt)
            STOP
          else
            !Gcf is negligibly less than zero - probably roundoff error
            Gcf(srtnh4)=0.0
          endif
        endif
      endif

C Repeat process if necessary
      if (time .lt. dt) then
        !Iteration
        itr=itr+1
c        if (itr.gt.5000) then
        if (itr.gt.500) then
           write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
           write(*,*) 'An iteration in so4cond exceeds 500'
           write(*,*) 'dt=',dt,'time=',time,'cdt=',cdt
           write(*,*) 'exponential decaying frac=',
     &               exp(-1.d0*sK(srtnh4)*cdt)
c$$$           iflagez=1
c$$$           return
           ! JJ change : if there is a problem just stop the program and
           ! fix this routine. Condensing just with whatever routine does not 
           ! fail is not a good design.
           STOP
        endif
        goto 10
      endif

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


      RETURN
      END

