
C     **************************************************
C     *  cond_nuc                                      *
C     **************************************************

C     WRITTEN BY Jeff Pierce, May 2007

C     This subroutine calculates the change in the aerosol size distribution
C     due to so4 condensation and binary/ternary nucleation during the
C     overal microphysics timestep.

C     ADD MORE HERE!

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Nnuci - number of nucleation size particles per size bin in grid cell
C     Mnuci - mass of given species in nucleation pseudo-bin (kg/grid cell)
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell
C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water
C     H2SO4rate - rate of H2SO4 chemical production [kg s^-1]
C     dt - total model time step to be taken (s)

C-----OUTPUTS-----------------------------------------------------------

C     Nkf, Mkf, Gcf - same as above, but final values

      SUBROUTINE cond_nuc(Nki,Mki,Gci,Nkf,Mkf,Gcf,H2SO4rate,dmappt,dt,ichm,jchm
     & ,kchm)             
      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      double precision Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      double precision H2SO4rate

      real dt
c      double precision dt

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      integer num_iter !number of iteration
      integer nuc_bin ! the nucleation bin
      integer iter ! number of iteration

      integer ichm ! i coordinate in PMCAMx
      integer jchm ! j coordinate in PMCAMx
      integer kchm ! k coordinate in PMCAMx by jgj 01/27/08

      double precision pi, R    ! pi and gas constant (J/mol K)
      double precision CSi,CSa   ! intial and average condensation sinks
      double precision CS1,CS2       ! guesses for condensation sink [s^-1]
      double precision CStest   !guess for condensation sink
      double precision Nk1(ibins), Mk1(ibins, icomp), Gc1(icomp-1)
      double precision Nk2(ibins), Mk2(ibins, icomp), Gc2(icomp-1)
      double precision Nk3(ibins), Mk3(ibins, icomp), Gc3(icomp-1)
      double precision mcond,mcond1    !mass to condense [kg]
      double precision tol      !tolerance
      double precision eps      !small number
      double precision sinkfrac(ibins) !fraction of condensation sink coming from bin k
      double precision totmass  !the total mass of H2SO4 generated during the timestep
      double precision tmass
      double precision CSch     !fractional change in condensation sink
      double precision CSch_tol !tolerance in change in condensation sink
      double precision addt     !adaptive timestep time
      double precision time_rem !time remaining
      double precision sumH2SO4 !used for finding average H2SO4 conc over timestep
      double precision fn, rnuc ! nucleation rate [# cm-3 s-1] and critical radius [nm]
      double precision gasConc  !gas concentration [kg]
      double precision mass_change !change in mass during nucleation.f
      double precision total_nh4_1,total_nh4_2
      double precision min_tstep !minimum timestep [s]
      double precision dmappt   !mixing ratio of dimethyl amine (for nucleation)

      logical nflg ! returned from nucleation, says whether nucleation occurred or not
 
C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------


C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(eps=1E-20)
      parameter(CSch_tol=0.01)
      parameter(min_tstep=1.0d0)

C-----CODE--------------------------------------------------------------

C Initialize values of Nkf, Mkf, Gcf, and time
      do j=1,icomp-1
         Gc1(j)=Gci(j)
      enddo
      do k=1,ibins
         Nk1(k)=Nki(k)
         do j=1,icomp
            Mk1(k,j)=Mki(k,j)
         enddo
      enddo

C     Get initial condensation sink
      CS1 = 0.d0
      call getCondSink(Nk1,Mk1,srtso4,CS1,sinkfrac)

C     Get initial H2SO4 concentration guess (assuming no nucleation)
C     Make sure that H2SO4 concentration doesn't exceed the amount generated
C     during that timestep (this will happen when the condensation sink is very low)

      !Get the steady state H2SO4 concentration
      call getH2SO4conc(H2SO4rate,CS1,Gc1(srtnh4),gasConc,dmappt)
      Gc1(srtso4) = gasConc
      addt = min_tstep
c      addt = 3600.d0
      totmass = H2SO4rate*addt*96.d0/98.d0

      !Get change size distribution due to nucleation with initial guess
      call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,nuc_bin,addt,CS1,dmappt)

      mass_change = 0.d0
      do k=1,ibins
         mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
      enddo
      mcond = totmass-mass_change ! mass of h2so4 to condense

      if (mcond.lt.0.d0)then
         tmass = 0.d0
         do k=1,ibins
            do j=1,icomp-idiag
               tmass = tmass + Mk2(k,j)
            enddo
         enddo
         if (abs(mcond).gt.totmass*1.0d-8) then
            if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
               tmass = 0.d0
               do j=1,icomp-idiag
                  tmass = tmass + Mk2(nuc_bin,j)
               enddo
               Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
               Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
               mcond = 0.d0
            else
            print*,'mcond < 0 in cond_nuc', mcond, totmass
            stop
            endif
         else
            mcond = 0.d0
         endif
      endif

      tmass = 0.d0
      do k=1,ibins
         do j=1,icomp-idiag
            tmass = tmass + Mk2(k,j)
         enddo
      enddo

C     Get guess for condensation
      call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3,ichm,jchm,kchm)

      Gc3(srtnh4) = Gc1(srtnh4)   

cdynamic      call eznh3eqm(Gc3,Mk3) ! NH3 condensation is calculated dynamically.-on
      call ezwatereqm(Mk3) !Do not count water amount for a while.

      ! check to see how much condensation sink changed
      call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac)
      CSch = abs(CS2 - CS1)/CS1    
       
c      if (CSch.gt.CSch_tol) then ! condensation sink didn't change much use whole timesteps
         ! get starting adaptive timestep to not allow condensationk sink
         ! to change that much
         addt = addt*CSch_tol/CSch/2
         addt = min(addt,dt)
         addt = max(addt,min_tstep)
         time_rem = dt ! time remaining
         
         num_iter = 0
         sumH2SO4=0.d0


         !DO ADAPTIVE TIMESTEPS
         do while (time_rem .gt. 0.d0)
            num_iter = num_iter + 1

C     Get the steady state H2SO4 concentration
            if (num_iter.gt.1)then ! no need to recalculate for first step
               call getH2SO4conc(H2SO4rate,CS1,Gc1(srtnh4),gasConc,dmappt)
               Gc1(srtso4) = gasConc
            endif

            sumH2SO4 = sumH2SO4 + Gc1(srtso4)*addt
            totmass = H2SO4rate*addt*96.d0/98.d0
            call nucleation(Nk1,Mk1,Gc1,Nk2,Mk2,Gc2,nuc_bin,addt,CS1,dmappt) 
            
            mass_change = 0.d0
            do k=1,ibins
               mass_change = mass_change + (Mk2(k,srtso4)-Mk1(k,srtso4))
            enddo
            mcond = totmass-mass_change ! mass of h2so4 to condense

            if (mcond.lt.0.d0)then
               tmass = 0.d0
               do k=1,ibins
                  do j=1,icomp-idiag
                     tmass = tmass + Mk2(k,j)
                  enddo
               enddo
               if (abs(mcond).gt.totmass*1.0d-8) then
                  if (-mcond.lt.Mk2(nuc_bin,srtso4)) then
                     tmass = 0.d0
                     do j=1,icomp-idiag
                        tmass = tmass + Mk2(nuc_bin,j)
                     enddo
                     Nk2(nuc_bin) = Nk2(nuc_bin)*(tmass+mcond)/tmass
                     Mk2(nuc_bin,srtso4) = Mk2(nuc_bin,srtso4) + mcond
                     mcond = 0.d0
                  else
                  print*,'mcond < 0 in cond_nuc', mcond, totmass
                  stop
                  endif
               else
                  mcond = 0.d0
               endif
            endif
            
c            Gc2(srtnh4) = Gc1(srtnh4)
c            call eznh3eqm(Gc2,Mk2,Mnuc2)
c            call ezwatereqm(Mk2,Mnuc2)

c            call getCondSink(Nk2,Mk2,Nnuc2,Mnuc2,srtso4,CStest,sinkfrac)  

            call ezcond(Nk2,Mk2,mcond,srtso4,Nk3,Mk3,ichm,jchm,kchm)
            Gc3(srtnh4) = Gc1(srtnh4)
cdynamic            call eznh3eqm(Gc3,Mk3) ! NH3 is calculated dynamically.-on
            call ezwatereqm(Mk3)
            
            ! check to see how much condensation sink changed
            call getCondSink(Nk3,Mk3,srtso4,CS2,sinkfrac)  

            time_rem = time_rem - addt
            if (time_rem .gt. 0.d0) then
               CSch = abs(CS2 - CS1)/CS1
Cjrp               if (CSch.lt.0.d0) then
Cjrp                  print*,''
Cjrp                  print*,'CSch LESS THAN ZERO!!!!!', CS1,CStest,CS2
Cjrp                  print*,'Nnuc',Nnuc1,Nnuc2
Cjrp                  print*,''
Cjrp
Cjrp                  addt = min(addt,time_rem)
Cjrp               else
c               addt = min(addt*CSch_tol/CSch,addt*1.5d0) ! allow adaptive timestep to change  commented out by LA
               if(CSch.le.eps) then   ! added by LA
                  addt=addt*1.5d0
               else
                  addt=min(addt*CSch_tol/CSch,addt*1.5d0)
               end if
                     
               addt = min(addt,time_rem) ! allow adaptive timestep to change
               addt = max(addt,min_tstep)
Cjrp               endif
               CS1 = CS2
               Gc1(srtnh4)=Gc3(srtnh4)
               do k=1,ibins
                  Nk1(k)=Nk3(k)
                  do j=1,icomp
                     Mk1(k,j)=Mk3(k,j)
                  enddo
               enddo         
            endif
         enddo
         Gcf(srtso4)=sumH2SO4/dt
Cjrp      else
Cjrp         num_iter = 1
Cjrp         Gcf(srtso4)=Gc1(srtso4)
Cjrp      endif
      
      do k=1,ibins
         Nkf(k)=Nk3(k)
         do j=1,icomp
            Mkf(k,j)=Mk3(k,j)
         enddo
      enddo      
      Gcf(srtnh4)=Gc3(srtnh4)


      RETURN
      END
