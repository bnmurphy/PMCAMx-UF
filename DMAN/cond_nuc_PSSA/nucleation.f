
C     **************************************************
C     *  nucleation                                    *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     UPDATED BY Ben Murphy & Jan Julin, October-November 2014

C     This subroutine calls the Vehkamaki 2002 nucleation
C     parameterizations and the NH3 and amine nucleation routines based on the 
C     Atmospheric Cluster Dynamic Code (ACDC)and gets the nucleation rates.

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell
C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water
C     dt - total model time step to be taken (s)

C-----OUTPUTS-----------------------------------------------------------

C     Nkf, Mkf, Gcf - same as above, but final values

      SUBROUTINE nucleation(Nki,Mki,Gci,Nkf,Mkf,Gcf,nuc_bin,dt,fn_all,cs)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer j,i,k
      double precision Nki(ibins), Mki(ibins, icomp), Gci(icomp-1)
      double precision Nkf(ibins), Mkf(ibins, icomp), Gcf(icomp-1)
      integer nuc_bin
      double precision dt

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision nh3ppt   ! gas phase ammonia in pptv
      double precision nh3_molec! gas phase ammonia in molec cm-3
      double precision dmappt   ! gas phase dimethyl amine in pptv
      double precision dma_molec! gas phase dimethyl amine in molec cm-3
      double precision h2so4    ! gas phase h2so4 in molec cc-1
      double precision fn       ! nucleation rate cm-3 s-1
      double precision rnuc     ! critical nucleation radius [nm]
      double precision gtime    ! time to grow to first size bin [s]
      double precision ltc, ltc1, ltc2 ! coagulation loss rates [s-1]
      double precision Mktot    ! total mass in bin
      double precision neps
      double precision meps
      double precision density  ! density of particle [kg/m3]
      double precision pi
      double precision frac     ! fraction of particles growing into first size bin
      double precision d1,d2    ! diameters of particles [m]
      double precision mp       ! mass of particle [kg]
      double precision mold     ! saved mass in first bin
      double precision mnuc     !mass of nucleation
      double precision cs       ! condensation sink [s-1]
      double precision d_dma    ! consumed dimethyl amine by the nucleation
                                ! process (kg m-3)
      double precision mfrac(icomp) ! fraction of the nucleated mass that will be assigned
                                    ! to each icomp species (passed to nucMassUpdate)
      double precision fn_all(nJnuc) !Magnitude of nucleation rates resolved by pathway
      double precision N_Vehk   !  critical cluster from Vehkamäki parameterization
      double precision amu      !atomic mass unit [kg]
      double precision max_addedH2SO4 ! the largest amount of SA mass added in this call by a pathway, used for nuc_bin
      double precision addedH2SO4 ! amount of SA added by a NPF pathway
      integer bin_addmnuc  ! internal variable used to determine to which bin the mass is added
                           ! for an individual pathway. The output variable nuc_bin will currently
                           ! contain the bin corresponding to the pathway that adds most mass /JJ 23/05/2016


C-----EXTERNAL FUNCTIONS------------------------------------------------

      double precision aerodens_PSSA
      external aerodens_PSSA

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (neps=1E8, meps=1E-8)
      parameter (pi=3.14159d0)
      parameter (amu=1.660539040d-27)

C-----CODE--------------------------------------------------------------

      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23
      ! ACDC Lookup table ternary nuclation does not need this in ppt
c      nh3ppt= (1.0e+21*8.314)*Gci(srtnh4)*temp/(pres*boxvol*gmw(srtnh4))
      nh3_molec = Gci(srtnh4)/boxvol*1000.d0/17.d0*6.022d23 
      dma_molec = Gci(srtdma)/boxvol*1000.d0/45.d0*6.022d23 !molec cm-3

      fn = 0.d0
      rnuc = 0.d0
      fn_all = 0.d0

      ! if no nucleation occurs the final arrays will be same as the initial arrays
      Mkf=Mki
      Nkf=Nki
      Gcf=Gci

      !nuc_bin tells cond_nuc which bin to remove sulfuric acid from if the mass to condense becomes negative
      nuc_bin=1
      max_addedH2SO4=0.d0 !keeping track which pathway resulted in most SA mass (for cond_nuc.f & nuc_bin) 

C     if requirements for nucleation are met, call nucleation subroutines

      if (amine_nuc.eq.1 .and.h2so4.gt.minval(amine_nuc_tbl_H2SO4).and.
     &     dma_molec.gt.minval(amine_nuc_tbl_DMA)) then
         call amine_nucl(temp,cs,h2so4,dma_molec,fn) !amine nuc

         if (fn.gt.0.d0) then
            !ACDC formation rate (the 5Feb2014 lookup table) corresponds to 
            !clusters with 4 acids+5 DMA
            mnuc=(4.d0*gmw(srtso4)+5.d0*gmw(srtdma))*amu !mass of 4A+5DMA particle
            !Determine to which bin the mass is added
            bin_addmnuc = 1
            do while (mnuc .gt. xk(bin_addmnuc+1))
               bin_addmnuc = bin_addmnuc + 1
            enddo

            !update mass and number
            addedH2SO4=4.d0*gmw(srtso4)*amu*fn*boxvol*dt
            Mkf(bin_addmnuc,srtso4) = Mkf(bin_addmnuc,srtso4)+
     &           +addedH2SO4
            Mkf(bin_addmnuc,srtdma) = Mkf(bin_addmnuc,srtdma)+
     &           +5.d0*gmw(srtdma)*amu*fn*boxvol*dt
            Nkf(bin_addmnuc) = Nkf(bin_addmnuc)+fn*boxvol*dt
            
            !checking if nuc_bin needs to be updated
            if (addedH2SO4.gt.max_addedH2SO4) then
               max_addedH2SO4=addedH2SO4
               nuc_bin=bin_addmnuc
            end if
            
            !Update DMA gas phase concentration
            Gcf(srtdma)=Gcf(srtdma)-Mkf(bin_addmnuc,srtdma) ! kg/grid cell
            !nucleation could use all of the DMA so check that we do not get negative
            if (Gcf(srtdma).lt.0.d0) then
c               print*, 'Neg DMA in nucleation', Gcf(srtdma)
               Gcf(srtdma) = 0.d0
            end if

            !we set amine nucleation rate as number 3 (it is
            !the latest addition to fn_all
            fn_all(3)=fn
         endif
      endif       

      if (h2so4.gt.minval(tern_nuc_tbl_H2SO4).and.
     &     nh3_molec.gt.minval(tern_nuc_tbl_NH3).and.tern_nuc.eq.1) then
c$$$         if (nh3ppt.gt.0.1.and.tern_nuc.eq.1) then
c$$$            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc
         call tern_nucl_acdc(temp,rh,cs,h2so4,nh3_molec,fn)

         if (fn.gt.0.d0) then
            !ACDC formation rate (the 2014-12-04 lookup table) corresponds to 
            !clusters with 4 acids+3 NH3
            mnuc=(4.d0*gmw(srtso4)+3.d0*gmw(srtnh4))*amu !mass of 4A+3NH3 particle
            !Determine to which bin the mass is added
            bin_addmnuc = 1
            do while (mnuc .gt. xk(bin_addmnuc+1))
               bin_addmnuc = bin_addmnuc + 1
            enddo            

            !update mass and number
            addedH2SO4=4.d0*gmw(srtso4)*amu*fn*boxvol*dt
            Mkf(bin_addmnuc,srtso4) = Mkf(bin_addmnuc,srtso4)+
     &           +addedH2SO4
            Mkf(bin_addmnuc,srtnh4) = Mkf(bin_addmnuc,srtnh4)+
     &           +3.d0*gmw(srtnh4)*amu*fn*boxvol*dt
            Nkf(bin_addmnuc) = Nkf(bin_addmnuc)+fn*boxvol*dt

            !checking if nuc_bin needs to be updated
            if (addedH2SO4.gt.max_addedH2SO4) then
               max_addedH2SO4=addedH2SO4
               nuc_bin=bin_addmnuc
            end if

            fn_all(1) = fn
         endif
      endif

      if (h2so4.gt.1.d4.and.bin_nuc.eq.1) then
         call vehk_nucl(temp,rh,h2so4,fn,rnuc,N_Vehk) !binary nuc

         if (fn.gt.0.d0 .and. N_Vehk.ge.4.d0) then !Vehkamäki param. only valid for N>4
            !update mass and number
            !nuclei are assumed to be sulfuric acid
            mnuc=N_Vehk*gmw(srtso4)*amu
            !Determine to which bin the mass is added
            bin_addmnuc = 1
            do while (mnuc .gt. xk(bin_addmnuc+1))
               bin_addmnuc = bin_addmnuc + 1
            enddo

            addedH2SO4=mnuc*fn*boxvol*dt
            Mkf(bin_addmnuc,srtso4) = Mkf(bin_addmnuc,srtso4)+addedH2SO4
            Nkf(bin_addmnuc) = Nkf(bin_addmnuc)+fn*boxvol*dt

            !checking if nuc_bin needs to be updated
            if (addedH2SO4.gt.max_addedH2SO4) then
               max_addedH2SO4=addedH2SO4
               nuc_bin=bin_addmnuc
            end if

            fn_all(2) = fn
         endif
      endif


      RETURN
      END

