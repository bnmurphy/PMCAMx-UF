
C     **************************************************
C     * org_cond                                       *
C     **************************************************

C     WRITTEN BY Peter Adams, June 2000
C     Based on thermocond.f

c     This routine determines the condensational driving force for
C     mass transfer of organic species between gas and aerosol phases.  It then calls
C     a mass- and number-conserving algorithm for condensation (or
C     evaporation) of aerosol.

C     An adaptive time step is used to prevent numerical difficulties.
C     To account for the changing gas phase concentration of organic
C     species, its decrease during a condensational time step is well-
C     approximated by an exponential decay with a constant, sK (Hz).
C     sK is calculated from the mass and number distribution of the
C     aerosol.  Not only does this approach accurately take into account
C     the changing organic species concentration, it is also used to
C     predict (and limit) the final organic species concentration.

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
      SUBROUTINE org_cond(Nki,Mki,Gci,Nkf,Mkf,Gcf,dt,xkDMAN,
     &       ichm, jchm, kchm)


      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------
      include 'sizecode.COM'
      include 'aervaria.inc'


C-----ARGUMENT DECLARATIONS---------------------------------------------

      real dt ! timestep
cdbg      real ygas(ngas) ! gas species mixing ratio
      double precision Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      double precision Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      double precision Nka(ibins), Mka(ibins, icomp), Gca(icomp-1)
      double precision xkDMAN(ibins+1), Mkmin ! xk is 
      double precision  cond_step2
C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,jj,kk,fff        !counters
      integer Gcflag(icomp)    !flags for checking insignificant Gc
      integer evaflag(ibins,icomp)
      integer  mg_srt,dii   

      integer ichm ! i coordinate in PMCAMx
      integer jchm ! j coordinate in PMCAMx
      integer kchm ! k coordinate in PMCAMx

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
      real tunning
      real temp_298k


      double precision csat1,csat2,csat3,psat1,psat2,psat3  !david 
      double precision psat4,psat5,psat6,psat7,psat8,cond_step
      double precision psat9,psat10,psat11,psat12,ssss
      double precision psat(icomp), Gcsat(icomp)
      double precision Gcsat1,Gcsat2,Gcsat3,Gcsat4      !david
      double precision Gcsat5,Gcsat6,Gcsat7,Gcsat8
      double precision Gcsat9,Gcsat10,Gcsat11,Gcsat12
      double precision kelvin(ibins),DH_elthap(icomp)             !david
      double precision sum_organic1,sum_organic2,sum_organic3
      double precision sum_tot_organic(ibins)
cd      double precision x_mol(ibins,13)    
      double precision x_mol_1(ibins),x_mol_2(ibins),x_mol_3(ibins)  
      double precision x_mol_4(ibins),x_mol_5(ibins),x_mol_6(ibins)
      double precision x_mol_7(ibins),x_mol_8(ibins),x_mol_9(ibins)
      double precision x_mol_10(ibins),x_mol_11(ibins),x_mol_12(ibins)
      double precision x_mol(ibins,icomp)
      double precision A_kelvin,dav1,dav2
      double precision evap(ibins,icomp)
 
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
      double precision mc, ttr,max_ttr
      double precision Neps     !value below which Nk is insignificant
      double precision cthresh  !determines minimum gas conc. for cond.
      double precision Ntotf, Ntoto, dNerr  !used to track number cons.

      double precision Mknh3max !Maximum allowable NH3
      double precision taumax   !Maxmum tau for ammonia condensation
      double precision moxd(ibins) !moxid/Nact
      double precision Mktot(icomp) ! Total mass of each species
      double precision Mktot1(icomp) ! Total mass of each species
      double precision Mktot2(icomp) ! Total mass of each species
      double precision Mtot ! Total mass of each section
      double precision Gcknh3(ibins) ! Fractional ammonia gas to be condensed
                                     ! when ammonia is limited.
      double precision sumataunh3    ! sum of atuc(ibins,srtnh3)
      double precision diameter_atau(41)
      double precision tau_max_evap(ibins,icomp)
    

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
      parameter(Neps=1.0e+10, zeta13=0.98483, cthresh=1.d-16)
                         !Neps is set a little higher than multicoag 
                         !to avoid redundant time step segregations.
      parameter(corfactor=1.0) ! a correction factor

C-----CODE---------------------i-----------------------------------------
      ssss= 0.025
      cond_step= 0.01      

      cond_step2=cond_step/2.d0

      do dii=srtsoa1,srtsoa5
         DH_elthap(dii)= 30.    !(KJ/mol )
      end do

      temp_298k=298.15     

      A_kelvin=4.*1.e-3*ssss*mol_weight/R/temp/liquid_density


       psat(srtsoa1) = 1.238786d-5  !Pa =>1 ug/m3
       psat(srtsoa2) = 1.238786d-4  !Pa =>10 ug/m3
       psat(srtsoa3) = 1.238786d-3  !Pa =>100 ug/m3
       psat(srtsoa4) = 1.238786d-2  !Pa =>1000 ug/m3
       psat(srtsoa5) = 1.238786e-8  !Pa =>10^-3 ug/m3  !EXLVOCS

        
      !new Psat    
      do dii=srtsoa1,srtsoa5
          psat(dii)= psat(dii)*exp(DH_elthap(dii)*1000./R
     &       *(1./temp_298k -1./temp))
      end do 

      !Initialization
      k=0
      jj=0

      dNerr=0.0
      time=0.0         !subroutine exits when time=dt
      tdt=2.d0/3.d0
      tunning=10.0     !multiplication for Mktot
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
      do k=1,ibins+1
        xk(k)=xkDMAN(k)
      enddo


      !Skip negative gas concentration if command
      if (Gci(srtso4) .lt. cthresh*boxmass) Gcflag(srtso4)=1 ! =1, Skip 
      if (Gci(srtnh3) .lt. cthresh*boxmass) Gcflag(srtnh3)=1 ! =1, Skip

      if (Gci(srtsoa1) .lt. cthresh*boxmass) Gcflag(srtsoa1)=1 ! =1, Skip
      if (Gci(srtsoa2) .lt. cthresh*boxmass) Gcflag(srtsoa2)=1 ! =1, Skip
      if (Gci(srtsoa3) .lt. cthresh*boxmass) Gcflag(srtsoa3)=1 ! =1, Skip
      if (Gci(srtsoa4) .lt. cthresh*boxmass) Gcflag(srtsoa4)=1 ! =1, Skip      
      if (Gci(srtsoa5) .lt. cthresh*boxmass) Gcflag(srtsoa5)=1 ! EXLVOCs 

      !If PSSA is on, turn on Gcflag(srtso4) for ezcond does H2SO4 condensatoin.
       Gcflag(srtso4)=0 !PSSA          

C Repeat from this point if multiple internal time steps are needed
 10   continue

      !Set dp equals to zero
      do k=1,ibins
         do j=1,icomp-1
            dp(k,j)=0.0
         enddo
      enddo


C Calculate tj and tk factors needed to calculate tau values
      mu=2.5277e-7*temp**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temp)))  !S&P eqn 8.6
      do j=1,icomp-1
        if (icond_test .eq. 1) then
          Di=0.1*1.0e-4 ! m^2/s ! cond_test 6/12/04 jgj
        else
          Di=gasdiff(temp,pres,gmw(j),Sv(j))
        endif
        tj(j)=2.*pi*Di*molwt(j)*1.0d-3/(R*temp)
      enddo
      sK(srtso4)=0.0d0
      sK(srtnh3)=0.0d0

      do dii=srtsoa1,srtsoa5     ! EXLVOCs  
         sK(dii)=0.0d0
      end do

 
      do k=1,ibins
         if (Nkf(k) .gt. Neps) then

            density=1400.0 ! [=] kg/m3

c            mp=(Mkf(k,srtso4)+Mkf(k,srtnh3))/Nkf(k)
          mp=(Mkf(k,srtso4) +Mkf(k,srtsoa1)+Mkf(k,srtsoa2)+
     &       +Mkf(k,srtsoa3)+Mkf(k,srtsoa4)+Mkf(k,srtsoa5)+      !  EXLVOCs
     &       +Mkf(k,srtinrt)+Mkf(k,srtnh3))/Nkf(k)

         else
           !nothing in this bin - set to "typical value"
             density=1400.0 ! Mk's are calculated based on this density
                           ! in initbounds. Dpk can be derived from this
                           ! density.
           mp=1.414*xk(k)
         endif

         Dpk(k)=((mp/density)*(6./pi))**(1.d0/3.d0)
         diameter_atau(k)=((mp/density)*(6./pi))**(0.333)
         Dpk(k)=h2ogrowth*Dpk(k)
         if (icond_test .eq. 1) then
           beta(srtso4)=1.
         else
           Kn=2.0*mfp/Dpk(k)                             !S&P eqn 11.35 (text)

           do dii=srtsoa1,srtsoa5    !!WITH EXLVOCs
              beta(dii)=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha(dii)) 
           end do            
         endif

         do dii=srtsoa1,srtsoa5      !!WITH EXLVOCS
             tk(k,dii)=(6./(pi*density))**(1./3.)*beta(dii)
         end do

         if (Nkf(k) .gt. 0.0) then
           Mtot=0.0
           do jj=1, icomp
             Mtot=Mtot+Mkf(k,jj)
           enddo

           do dii=srtsoa1,srtsoa5    !!WITH EXLVOCS
              sK(dii)=sK(dii)+tk(k,dii)*Nkf(k)*(Mtot
     &             /Nkf(k))**(1.d0/3.d0)
           end do
         endif
      enddo

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

       do k=1,ibins
         kelvin(k)=exp(A_kelvin/dpk(k))

         sum_tot_organic(k)=Mkf(k,srtsoa1) +Mkf(k,srtsoa2)+  !
     &    +Mkf(k,srtsoa3)+Mkf(k,srtsoa4)+ Mkf(k,srtsoa5)     ! EXLVOCs
     &    +Mkf(k,srtso4)+Mkf(k,srtnh3)                       ! 

          do dii=srtsoa1,srtsoa5   !!WITH EXLVOCS
            x_mol(k,dii)=Mkf(k,dii)/sum_tot_organic(k)
            if (x_mol(k,dii).ne.x_mol(k,dii)) then
               x_mol(k,dii)= 1.d-36
            endif
          end do
       enddo

       do dii=srtsoa1,srtsoa5  !!WITH EXLVOCS
            Gcsat(dii)=(psat(dii)*200.*boxmass)/(28.9*pres)
       end do 

       do k=1,ibins

          do dii=srtsoa1,srtsoa5   !!WITH EXLVOCS
             dp(k,dii)=(Gcf(dii)/200.)/(boxmass/28.9)*pres
     &                -(psat(dii)*kelvin(k))*x_mol(k,dii) !!WITH EXLVOCS
          end do 

       enddo

cdav   organic
       do dii=srtsoa1,srtsoa5    ! !WITH EXLVOCS
               sK(dii)=sK(dii)*zeta13*tj(dii)*R*temp/
     &           (molwt(dii) *1.d-3)/(boxvol*1.d-6)
       end do
          
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
            atauc(k,j)=0.0d0  
            atau(k,j)=0.0d0
         enddo

        !atauc is tau that is a parameter to describe condensation
        !driving force when pressure difference is equal to constant.
        !atau is same to atauc except having a exponential decaying
        !pressure term. The current ammonia condensation uses atauc as atau.
        !In case atauc is bigger than taumax that atau is forced to be taumax.
 

c        organic
         do dii=srtsoa1,srtsoa5   !!WITH EXLVOCS
             atauc(k,dii)=tj(dii)*tk(k,dii)*dp(k,dii)*cdt
         end do

c---------------------------------------------------------------------------
         do dii=srtsoa1,srtsoa5    !!WITH EXLVOCS
          if (sK(dii) .gt. 0.0) then
            evap(k,dii)=Gcf(dii)-(kelvin(k)*Gcsat(dii))*x_mol(k,dii)

            atau(k,dii)=tj(dii)*R*temp/(molwt(dii)*1.d-3)
     &      /(boxvol*1.d-6)*tk(k,dii)
     &      *(Gcf(dii)-(kelvin(k)*Gcsat(dii))*x_mol(k,dii))/sK(dii)
     &      *(1.d0-exp(-1.d0*sK(dii)*cdt))

            tau_max_evap(k,dii) = tj(dii)*R*temp/(molwt(dii)*1.d-3)
     &      /(boxvol*1.d-6)*tk(k,dii)
     &      *(-0.01d0*Mkf(k,dii))/sK(dii)
     &      *(1.d0-exp(-1.d0*sK(dii)*cdt))

            if (atau(k,dii).lt.0.0) then
             if (abs(atau(k,dii)).ge.abs(tau_max_evap(k,dii))) then
                 atau(k,dii)=tau_max_evap(k,dii)
             end if
            end if

          else
           atau(k,dii)=0.0d0  !nothing to condense onto
          endif
         end do 
      enddo
        
C-----Adjust a time step 

      tr=1.0 !The following sections limit the condensation time step
             !when necessary.  tr is a factor that describes by
             !how much to reduce the time step.

      !Make sure masses of individual species don't change too much
      do j= srtsoa1,srtsoa5   !!WITH EXLVOCS
        do k=1,ibins
          if (Nkf(k) .gt. Neps) then
            mc=0.0
            do jj=1,icomp-1
              mc=mc+Mkf(k,jj)/Nkf(k)
            enddo
            if (mc/xk(k) .gt. 1.0d-3) then
                 if (ttr.gt.max_ttr)then
                  max_ttr=ttr
                 end if      

cd              !species has significant mass in particle - limit change
              if (abs(atau(k,j))/(mc**tdt) .gt.0.1 ) then
                ttr=abs(atau(k,j))/(mc**tdt)/0.05

                if (ttr. gt. tr) then 
                  tr=ttr
                endif
              endif
            else
              !species is new to particle - set max time step
              if ((cdt/tr .gt. 0.1) .and. (atau(k,j).gt. 0.0)) then 
                tr=cdt/0.1
              endif
            endif
          endif
        enddo
        !Make sure gas phase concentrations don't change too much
        !control gas reduction rate on 12/15/07, jgj
        if (dt.gt.120.) then
           gasfrac=0.50 ! not allow less than 50% of reduction theoretically
        else
           gasfrac=0.25 ! not allow less than 75% of reduction theoretically
        endif
        if (exp(-1.d0*sK(j)*cdt) .lt. gasfrac) then
          ttr=-2.d0*cdt*sK(j)/log(gasfrac)
          if (ttr .gt. tr) then 
            tr=ttr
          endif
           if (ttr .gt.100000.)then
                   ttr = tr
                 end if

        endif

 30   continue
      enddo

      !Never shorten timestep by less than half
      if (tr .gt. 1.0) tr=max(tr,2.0)

      !Repeat for shorter time step if necessary
      if (tr .gt. 1.0) then
         cdt=cdt/tr
         goto 20
      endif

C Call condensation subroutine to do mass transfer

      do j= srtsoa1,srtsoa5  !!WITH EXLVOCS    !Loop over all aerosol components

        !Swap tau values for this species into array for cond
        do k=1,ibins
          tau(k)=atau(k,j)

          tau(k)=corfactor*tau(k) ! A correction factor is applied.
        enddo
        
        call mnfix_PSSA(Nkf,Mkf,ichm,jchm,kchm,6)
            ! adjust average mass in each size bin within boundaries

        !Call condensation routine
        Ntotf=0.0d0
        do k=1,ibins
          Ntotf=Ntotf+Nkf(k)
        enddo

        !oxidated mass calculation
        do k=1, ibins
          moxd(k)=0.0d0 !oxidated mass
        enddo

        call tmcond(tau,xk,Mkf,Nkf,Mko,Nko,j,moxd)

        !Check for number conservation
        Ntoto=0.0d0
        do k=1,ibins
          Ntoto=Ntoto+Nko(k)
        enddo

        dNerr=Ntotf-Ntoto

        if (abs(dNerr/Ntoto) .gt. 1.e-4) then
          write(*,*)'ERROR in so4cond: Number not conserved'
          write(*,*)'Ntoto, Ntotf, dNerr/Ntoto'
     &             ,Ntoto, Ntotf, dNerr/Ntoto
          if (abs(dNerr/Ntoto) .gt. 1e-2) then
            write(*,*)'Serious Error in so4cond: Number not conserved 
     &       less than 1 %'
            STOP
          endif
        endif

        !Update gas phase concentration
        mi(j)=0.0d0
        mf(j)=0.0d0
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

c=======================================================================    
       !Check organic gas
      do fff=srtsoa1,srtsoa5    !!WITH EXLVOCS
      if (Gcflag(fff) .eq. 0.d0) then
        if (Gcf(fff) .lt. 0.0d0) then
          if (abs(Gcf(fff)) .gt. 1.d-5) then
            !Gcf is substantially less than zero - this is a problem
            write(*,*) 'ERROR in organis cond: organic(g) < 0'
            write(*,*) 'time=',time
            write(*,*) 'Gcf()=',Gcf(fff),fff
            write(*,*) 'organic [=] ppt',Gcf(fff)/boxmass/1.0e-12
     &               /gmw(fff)*28.9
            write(*,*) 'Gass consumed in kg',(mi(fff)-mf(fff))
     &               *gmw(fff)/molwt(fff)
            write(*,*) 'time=',time,'cdt=',cdt
            write(*,*) 'exponential decaying frac=',
     &               exp(-1.d0*sK(fff)*cdt)
            Mktot(fff)=0.0d0
            do kk=1,ibins
               Mktot(fff)=Mktot(fff)+Mk(kk,fff)
            enddo
            write(*,*)'Max. Possible Total Mkorgan =',
     &         0.375d0*Mktot(fff)
            write(*,*)'Initial organic mass=',Gci(fff)
            STOP
          else
            !Gcf is negligibly less than zero - probably roundoff error
            Gcf(fff)=0.0d0
          endif
        endif
      endif
      end do
   
C Repeat process if necessary
      if (time .lt. dt) goto 10
     
 100  continue   !skip to here if there is no gas phase to condense

      RETURN
      END

