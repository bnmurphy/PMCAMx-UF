
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


C-----EXTERNAL FUNCTIONS------------------------------------------------

      double precision aerodens_PSSA
      external aerodens_PSSA

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (neps=1E8, meps=1E-8)
      parameter (pi=3.14159d0)

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

cdbg      print*,'h2so4',h2so4,'nh3ppt',nh3ppt

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (h2so4.gt.1.d4) then
         if (amine_nuc.eq.1.and.dma_molec.gt.1.d4) then
            call amine_nucl(temp,cs,h2so4,dma_molec,fn,rnuc) !amine nuc

c$$$            !Update DMA Concentration
c$$$	    !The 0.31 factor should be the mass fraction of DMA in
c$$$	    !the nucleated particles.
c$$$            d_dma = 0.31d0 * (4.d0/3.d0*pi*(rnuc*1D-9)**3)*1350.d0*fn*1.d6*dt !kg m-3
c$$$            dmappt = dmappt - d_dma/0.045d0 / (pres/(8.314d0*temp)) * 1.d12 !pptv
c$$$
c$$$            !nucleation could use all of the DMA so check that we do not get negative
c$$$            if (dmappt.lt.0.d0) then
c$$$               print*, 'Neg DMA in nucleation', dmappt
c$$$               dmappt = 0.d0
c$$$            end if

            if (fn.gt.0.d0) then
               !update mass and number
               !nuclei consist of sulfate and amine compounds
	       !assume 1-to-1 ratio of DMA and sulfuric acid
               mfrac = (/0.69, 0.0, 0.0, 0.31, 0.0/)
               call nuc_massupd(Nkf,Mkf,Gcf,nuc_bin,dt,fn,rnuc,mfrac)
               
               !Update DMA gas phase concentration
               Gcf(srtdma)=Gcf(srtdma)-Mkf(nuc_bin,srtdma)   ! kg/grid cell
               !nucleation could use all of the DMA so check that we do not get negative
               if (Gcf(srtdma).lt.0.d0) then
                  print*, 'Neg DMA in nucleation', Gcf(srtdma)
                  Gcf(srtdma) = 0.d0
               end if

               !set amine nucleation rate as number 3 for now (it is
               !the latest addition to fn_all
               fn_all(3)=fn
            endif
         endif       

         if (nh3_molec.gt.1.d6.and.tern_nuc.eq.1) then
c$$$         if (nh3ppt.gt.0.1.and.tern_nuc.eq.1) then
c$$$            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc
            call tern_nucl_acdc(temp,rh,cs,h2so4,nh3_molec,fn,rnuc)

            if (fn.gt.0.d0) then
               !update mass and number
               !nuclei are assumed as ammonium bisulfte
               mfrac = (/0.8144, 0.0, 0.1856, 0.0, 0.0/)
               call nuc_massupd(Nkf,Mkf,Gcf,nuc_bin,dt,fn,rnuc,mfrac)
               fn_all(1) = fn
            endif
         endif

         if (bin_nuc.eq.1) then
            call vehk_nucl(temp,rh,h2so4,fn,rnuc) !binary nuc

            if (fn.gt.0.d0) then
               !update mass and number
               !nuclei are assumed to be sulfuric acid
               mfrac = (/1.0, 0.0, 0.0, 0.0, 0.0/)
               call nuc_massupd(Nkf,Mkf,Gcf,nuc_bin,dt,fn,rnuc,mfrac)
               fn_all(2) = fn
            endif
         endif    
      endif


      RETURN
      END

