
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

      SUBROUTINE so4cond_oxd(Nki,Mki,Nkf,Mkf,dt,moxid,iact,ichm,jchm    ! cf
     &                      ,kchm,xkDMAN)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'
      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

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
      parameter(corfactor=1.0) ! a correction factor
      parameter(cvt = 1.0d-15) ! a convert factor from ug/m3 to kg/cm3

C-----CODE--------------------------------------------------------------

      !Initialization
      k=0
      jj=0

      time=0.0         !subroutine exits when time=dt
      tdt=2.d0/3.d0

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

C Repeat from this point if multiple internal time steps are needed
 10   continue


C Choose appropriate time step

      !Try to take full time step
      cdt=dt-time

      !Not more than 15 minutes
      if (cdt .gt. 900.) then
        cdt=900.

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

      enddo

      !Update atau considering the amount added by aqueous chemistry
      do k=iact,ibins
         do j=1,icomp-1
            ! JJ bugfix: both the initial and final masses need to be 
            ! average mass of single particle (Adams & Seinfeld 2002, Eq. 9&10)
            Mko(k,j)=Mkf(k,j)+moxid(k,j)!*(1./Nkf(k)))
            atau(k,j)=1.5*((Mko(k,j)**tdt)-(Mkf(k,j)**tdt))/(Nkf(k)**tdt) 
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
      
C Call condensation subroutine to do mass transfer

      do j=1,icomp-1  !Loop over all aerosol components

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

        !oxidated mass calculation, dummy, not using it, jgj
        do k=1, ibins
          moxd(k)=0.0 !oxidated mass
        enddo

        call tmcond(tau,xk,Mkf,Nkf,Mko,Nko,j,moxd)

        !Swap into Nk, Mk
        do k=1,ibins
          Nkf(k)=Nko(k)
          do jj=1,icomp-1
            Mkf(k,jj)=Mko(k,jj)
          enddo
        enddo

 40   continue
      enddo

C Update time
      time=time+cdt

 100  continue   !skip to here if there is no gas phase to condense

      ! Adjust finally before return
      call mnfix_PSSA(Nkf,Mkf,ichm,jchm,kchm)

      RETURN
      END

