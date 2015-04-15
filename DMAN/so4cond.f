
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

C     Added amine condensation,gas phase amine (passed in as dmappt) will 
C     condense to particle phase ammonium JJ 2015/02

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
     &           iflagez,dmappt)

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

      double precision dmappt ! dma concentration in ppt
      double precision dmaGc  ! dma gas conc in kg/gridcell
      real dmaMw  ! molecular weight of dma (real because of gasdiff subr.)
      real molwtdma ! mol. weight dma ion form
      double precision dpdma(ibins) ! pressure difference for dma (assume 0 sat.vap.pres)
      double precision tjdma ! tj factor for dma
      double precision tkdma(ibins) ! tk factor for dma
      double precision sKdma ! sk decay constant for dma
      double precision ataudma(ibins) !cond. parameter for dma
      double precision ataucdma(ibins) !another cond parameter for dma
      real betadma ! trans. regime correction for dma
      real alphadma ! accommodation coefficient for dma
      integer Gcflagdma ! flag to check insignificant gas conc.
      double precision Gckdma(ibins) ! Fractional dma gas to be condensed
                                     ! when ammonia+dma is limited.
      double precision sumataudma    ! sum of atau for dma
      double precision ataufrac ! nh3 atau divided by dma+nh3 ataus in a bin
      double precision midma,mfdma ! condensed phase nh4 mass before and after
                                   ! dma condensation (used to update dmaGc)

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

      !initialize dma related constants JJ 150228
      dmaMw=45.0
      molwtdma=46.0
      alphadma=1.0
      betadma=0.0
      dpdma=0.d0
      Gcflagdma=0
      tkdma=0.0d0

!dma converting from ppt to kg/gridcell
      dmaGc=dmappt*pres*boxvol*dmaMw*1.0d-21/(R*temp)

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

! Skip also dma condensation if too small gas concentration
      if (dmaGc .lt. cthresh*boxmass) Gcflagdma=1 ! =1, Skip
!!!!!!!!!!!!!!!!!

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
      
      !dp is only used for calculation of "atauc", which is not used for
      !the condensation, only for checking if taumax should be used
      do k=1,ibins
        if (Mkf(k,srtso4) .gt. 0.) then
          ratio=Mkf(k,srtnh3)/Mkf(k,srtso4)
        else
          ratio=0.375 ! skip calculating of dp(k,srtnh3)
        endif
        if (ratio .lt. 0.375) then ! 0.375 is max ratio
          dp(k,srtnh3)=(Gcf(srtnh3)/molwt(srtnh3))/(boxmass/28.9)*pres
        endif

       !for dma
        if (ratio .lt. 0.375) then ! 0.375 is max ratio
          dpdma(k)=(dmaGc/molwtdma)/(boxmass/28.9)*pres
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

C Calculate diffusivity and tj for dma, Sv is the same as for all species currently
      Di=gasdiff(temp,pres,dmaMw,42.88)
      tjdma=2.d0*pi*Di*molwtdma*1.0d-3/(R*temp)
      sKdma=0.0d0
CCCCCCCCCCCCC

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
          mp=(Mkf(k,srtso4)+Mkf(k,srtnh3))/Nkf(k)
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
        Dpk(k)=((mp/density)*(6./pi))**(0.333)
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
          ! DMA
          betadma=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alphadma)
          !
        endif
        tk(k,srtso4)=(6./(pi*density))**(1./3.)*beta(srtso4)
        tk(k,srtnh3)=(6./(pi*density))**(1./3.)*beta(srtnh3)
        ! DMA
        tkdma(k)=(6./(pi*density))**(1./3.)*betadma
        !
        if (Nkf(k) .gt. 0.0) then
          Mtot=0.0
          do jj=1, icomp
            Mtot=Mtot+Mkf(k,jj)
          enddo
          sK(srtso4)=sK(srtso4)+tk(k,srtso4)*Nkf(k)*(Mtot
     &             /Nkf(k))**(1.d0/3.d0)
          sK(srtnh3)=sK(srtnh3)+tk(k,srtnh3)*Nkf(k)*(Mtot
     &             /Nkf(k))**(1.d0/3.d0)
          ! DMA
          sKdma=sKdma+tkdma(k)*Nkf(k)*(Mtot
     &             /Nkf(k))**(1.d0/3.d0)
          !
        endif
      enddo
      sK(srtso4)=sK(srtso4)*zeta13*tj(srtso4)*R*temp/(molwt(srtso4)
     &      *1.d-3)/(boxvol*1.d-6)
      sK(srtnh3)=sK(srtnh3)*zeta13*tj(srtnh3)*R*temp/(molwt(srtnh3)
     &      *1.d-3)/(boxvol*1.d-6)
      ! DMA
      sKdma=sKdma*zeta13*tjdma*R*temp/(molwtdma
     &      *1.d-3)/(boxvol*1.d-6)
      !

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
         !DMA 
         ataucdma(k)=0.0d0
         ataudma(k)=0.d0
         ! 

        !atauc is tau that is a parameter to describe condensation
        !driving force when pressure difference is equal to constant.
        !atau is same to atauc except having a exponential decaying
        !pressure term. The current ammonia condensation uses atauc as atau.
        !In case atauc is bigger than taumax that atau is forced to be taumax.
 
        !Calculate a driving force for ammonia condensation
        atauc(k,srtso4)=tj(srtso4)*tk(k,srtso4)*dp(k,srtso4)*cdt
        atauc(k,srtnh3)=tj(srtnh3)*tk(k,srtnh3)*dp(k,srtnh3)*cdt
        !DMA
        ataucdma(k)=tjdma*tkdma(k)*dpdma(k)*cdt
        !since dma is condensing to nh4 scaling is necessary
        ataucdma(k)=molwt(srtnh3)/molwtdma*ataucdma(k)
        !

        !Calculate a driving force for sulfuric acid condensation
        if (sK(srtso4) .gt. 0.0) then
          atau(k,srtso4)=tj(srtso4)*R*temp/(molwt(srtso4)*1.d-3)
     &      /(boxvol*1.d-6)*tk(k,srtso4)*Gcf(srtso4)/sK(srtso4)
     &      *(1.d0-exp(-1.d0*sK(srtso4)*cdt))
        else
          atau(k,srtso4)=0.0  !nothing to condense onto
        endif

        !Calculate a driving force for ammonia condensation
        if (sK(srtnh3) .gt. 0.0) then
          atau(k,srtnh3)=tj(srtnh3)*R*temp/(molwt(srtnh3)*1.d-3)
     &      /(boxvol*1.d-6)*tk(k,srtnh3)*Gcf(srtnh3)/sK(srtnh3)
     &      *(1.d0-exp(-1.d0*sK(srtnh3)*cdt))
        else
          atau(k,srtnh3)=0.0  !nothing to condense onto
        endif

        !Calculate a driving force for amine condensation
        if (sKdma .gt. 0.0 .and. ataucdma(k).gt.0.0) then
          ataudma(k)=tjdma*R*temp/(molwtdma*1.d-3)
     &      /(boxvol*1.d-6)*tkdma(k)*dmaGc/sKdma
     &      *(1.d0-exp(-1.d0*sKdma*cdt))
          ! again scale because of dma->nh4
          ataudma(k)=molwt(srtnh3)/molwtdma*ataudma(k)
        else
          ataudma(k)=0.0  !nothing to condense onto or too much nh3 already
        endif

        Mktot(srtso4)=0.0
        do kk=1,ibins
          Mktot(srtso4)=Mktot(srtso4)+Mkf(kk,srtso4)
        enddo
        Mktot(srtnh4)=0.0
        do kk=1,ibins
          Mktot(srtnh4)=Mktot(srtnh4)+Mkf(kk,srtnh4)
        enddo
        

c$$$        !DMA
c$$$        sumataudma=0.0
c$$$        do kk=1,ibins
c$$$           sumataudma=sumataudma+ataudma(kk)
c$$$        enddo
c$$$        Gckdma(k)=dmaGc*ataudma(k)/sumataudma
c$$$        !

        !Separate the cases of total ammonia and amine (!) is greater than existing sulfate
        ! or not.
        !Note that the gas phase mass of dma is multiplied by the ammonia/dma mass ratio
        !since the current treatment has the amine transforming into ammonium in the particle
        !phase
        !The following procedure has ammonia condensing before dma: the order is otherwise 
        !irrelevant, but if we have ammonia rich conditions, and if the individual taus are 
        !defined through the final and initial masses as for taumax below they will depend on the order of condensation

        !150305 lets do this the old fashioned way, if there is no space for amine to condense then too bad...

        if (0.375*Mktot(srtso4).gt.(Mktot(srtnh4)+Gcf(srtnh4))) then
c        if (0.375*Mktot(srtso4).gt.(Mktot(srtnh4)+Gcf(srtnh4)+dmaGc*17.d0/45.d0)) then
        !Sulfate rich condition
           sumataunh3=0.0
           do kk=1,ibins
              sumataunh3=sumataunh3+atau(kk,srtnh4)
           enddo
           Gcknh3(k)=Gcf(srtnh4)*atau(k,srtnh4)/sumataunh3
           Mknh3max=Mkf(k,srtnh3)+Gcknh3(k)
c           Mknh3max=Mkf(k,srtnh3)+Gcknh3(k)+Gckdma(k)*17.d0/45.d0
                                         ! maximally allowable NH4 mass
c        else if (Mkf(k,srtnh3).gt.0.375*Mkf(k,srtso4)) then !Total ammonia rich condition
        else
           Mknh3max=0.375*Mkf(k,srtso4) ! maximally allowable NH4 mass
c$$$        else
c$$$           ! this is the tricky part, current solution: let ammonia and dma condense in the 
c$$$           ! same ratio that they would if the condition for sulfate neutralization is not in effect
c$$$
c$$$           Mknh3max=Mkf(k,srtnh3)+Gcknh3(k)/(Gcknh3(k)+Gckdma(k))*
c$$$     &      (0.375*Mkf(k,srtso4)-Mkf(k,srtnh3)) ! maximally allowable NH4 mass
        endif
cdbg        if (Nkf(k) .gt. 0.) then
        if (Nkf(k) .gt. Neps) then
          taumax=1.5*((Mknh3max**tdt)-(Mkf(k,srtnh3)**tdt))/(Nkf(k)
     &     **tdt) 
            ! See Eq.(9) in Adams and Seinfeld (2002, JGR)
        else
          taumax=0. !For safety
        endif
c        if (atauc(k,srtnh3) .gt. taumax) then
!     taumax here refers to nh3
        if ((atauc(k,srtnh3)) .gt. taumax) then
          if (taumax .ge. 0.) then
c            ataufrac=atau(k,srtnh3)/(atau(k,srtnh3)+ataudma(k)) 
            atau(k,srtnh3)=taumax
c            atau(k,srtnh3)=ataufrac*taumax
c            ataudma(k)=(1-ataufrac)*taumax
          else
            atau(k,srtnh3)=0. !For safety
c            ataudma(k)=0. !for safety
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
              ! check if species ammonia
               if (j.eq.srtnh4) then
                  !first check with ammonia
                  if (abs(atau(k,j))/(mc**tdt) .gt. 0.1) then
                     ttr=abs(atau(k,j))/(mc**tdt)/0.05
                     if (ttr. gt. tr) then 
                        tr=ttr
                     endif
                  endif
                  ! check also with dma if there is dma to condense
                  !for now lets just assume dma will change particle mass minimally
                  ! because the timestep iteration seems to have trouble
c$$$                  if (Gcflagdma.eq.0) then
c$$$                     if (abs(ataudma(k))/(mc**tdt) .gt. 0.1) then
c$$$                        ttr=abs(ataudma(k))/(mc**tdt)/0.05
c$$$                        if (ttr. gt. tr) then 
c$$$                           tr=ttr
c$$$                        endif
c$$$                     endif
c$$$                  endif               
               else ! other species
                  if (abs(atau(k,j))/(mc**tdt) .gt. 0.1) then
                     ttr=abs(atau(k,j))/(mc**tdt)/0.05
                     if (ttr. gt. tr) then 
cdbg                  limit='amass'
cdbg                  write(limit(7:11),'(I2,X,I2)') k,j
                        tr=ttr
                     endif
                  endif
               endif
            else
              !maybe should do something here to account for dma but this else and tr
              !does not seem to make much sense anyway, and the apparently huge tr
              !will be limited below to max value of 2 / JJ 150228

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

      !DMA
      !check the change in dma gas phase concentration and adjust tr if necessary
      if (Gcflagdma.eq.0) then
         if (exp(-1.d0*sKdma*cdt) .lt. gasfrac) then
            ttr=-2.d0*cdt*sKdma/log(gasfrac)
            if (ttr .gt. tr) then 
               tr=ttr
            endif
         endif
      endif
      !

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

C Finally condense also dma

      if(Gcflagdma.eq.0) then ! If there is gas to condense
         !first check the taumax for dma
         Mktot(srtso4)=0.0
         do kk=1,ibins
            Mktot(srtso4)=Mktot(srtso4)+Mkf(kk,srtso4)
         enddo
         Mktot(srtnh4)=0.0
         do kk=1,ibins
            Mktot(srtnh4)=Mktot(srtnh4)+Mkf(kk,srtnh4)
         enddo

         sumataudma=0.0
         do kk=1,ibins
            sumataudma=sumataudma+ataudma(kk)
         enddo
         do k=1,ibins
            Gckdma(k)=dmaGc*ataudma(k)/sumataudma

            !we dont need to care about gas phase ammonia anymore for this time step
            if (0.375*Mktot(srtso4).gt.(Mktot(srtnh4)+dmaGc*17.d0/45.d0)) then
            !Sulfate rich condition

               Mknh3max=Mkf(k,srtnh3)+Gckdma(k)*17.d0/45.d0
                                         ! maximally allowable NH4 mass
            else                !Total ammonia rich condition
               !unlike nh3 earlier, we can now just let enough dma condense to finish
               !neutralizing so4
               Mknh3max=0.375*Mkf(k,srtso4) ! maximally allowable NH4 mass
            endif
cdbg        if (Nkf(k) .gt. 0.) then
            if (Nkf(k) .gt. Neps) then
               taumax=1.5*((Mknh3max**tdt)-(Mkf(k,srtnh3)**tdt))/(Nkf(k)
     &              **tdt) 
            ! See Eq.(9) in Adams and Seinfeld (2002, JGR)
            else
               taumax=0.        !For safety
            endif
!     taumax here refers to dma
            if ((ataudma(k)) .gt. taumax) then !changed atauc -> atau
               if (taumax .ge. 0.) then 
                  ataudma(k)=taumax
               else
                  ataudma(k)=0. !for safety
               endif
            endif
         enddo

        !we use atau for tau variable for dma as well (and we ignore icond_test which should
        !always be zero
        do k=1,ibins
           !make sure we are not condensing to bins where we have enough already to neutralize so4
           if (ataucdma(k).gt.0.0d0) then
              tau(k)=ataudma(k)
           else
              tau(k)=0
           endif
c          tau(k)=corfactor*tau(k) ! A correction factor is applied. (use same as for ammonia/other species)
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
           moxd(k)=0.0          !oxidated mass
        enddo

        !we call tmcond with tau for dma, and species index for nh4
        !this way we condense gas phase dma into particle phase nh4
        call tmcond(tau,xk,Mkf,Nkf,Mko,Nko,srtnh4,moxd)

        !Check for number conservation
        Ntoto=0.0
        do k=1,ibins
          Ntoto=Ntoto+Nko(k)
        enddo

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

        !Update gas phase concentration for dma
        !change in particle phase is in nh4 mass
        midma=0.0
        mfdma=0.0
        do k=1,ibins
          midma=midma+Mkf(k,srtnh4)
          mfdma=mfdma+Mko(k,srtnh4)
        enddo
        dmaGc=dmaGc+(midma-mfdma)*dmaMw/molwt(srtnh4)

        !Swap into Nk, Mk
        do k=1,ibins
          Nkf(k)=Nko(k)
          do jj=1,icomp-1
            Mkf(k,jj)=Mko(k,jj)
          enddo
        enddo

      endif

C End of DMA condensing

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

      !Check dma gas
      if (Gcflagdma .eq. 0) then
        if (dmaGc .lt. 0.0) then
          if (abs(dmaGc) .gt. 1.d-5) then
            !Gcf is substantially less than zero - this is a problem
            write(*,*) 'ERROR in so4cond: DMA(g) < 0'
            write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
            write(*,*) 'time=',time
            write(*,*) 'Gcf()=',dmaGc
            write(*,*) 'DMA [=] ppt',dmaGc/boxmass/1.0d-12
     &               /dmaMw*28.9
            write(*,*) 'Gass consumed in kg',(midma-mfdma)
     &               *dmaMw/molwt(srtnh4)
            write(*,*) 'time=',time,'cdt=',cdt
            write(*,*) 'exponential decaying frac=',
     &               exp(-1.d0*sKdma*cdt)
            STOP
          else
            !Gcf is negligibly less than zero - probably roundoff error
            dmaGc=0.0
          endif
        endif
      endif

C Repeat process if necessary
      if (time .lt. dt) then
        !Iteration
        itr=itr+1
c     Commenting this out: iteration should never take this long, but if it fails
c     the subroutine would return with "halfway" updated values for the various
c     concentrations (and the try eznh3), which sounds bad / JJ 150228
c        if (itr.gt.5000) then
c$$$        if (itr.gt.500) then
c$$$           write(*,*) 'Coord.(i,j,k)=',ichm,jchm,kchm
c$$$           write(*,*) 'An iteration in so4cond exceeds 5000'
c$$$           write(*,*) 'dt=',dt,'time=',time,'cdt=',cdt
c$$$           write(*,*) 'exponential decaying frac=',
c$$$     &               exp(-1.d0*sK(srtnh4)*cdt)
c$$$           iflagez=1
c$$$           return
c$$$c           STOP
c$$$        endif
        goto 10
      endif

      ! change dma from kg/grid cell to ppt so that dman gets the updated value
      dmappt=1.0d21*R*temp*dmaGc/(pres*boxvol*dmaMw)
      

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

