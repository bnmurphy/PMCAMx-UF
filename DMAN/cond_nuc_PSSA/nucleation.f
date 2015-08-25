
C     **************************************************
C     *  nucleation                                    *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     UPDATED BY Ben Murphy & Jan Julin, October-November 2014

C     This subroutine calls the Vehkamaki 2002 and Napari 2002 nucleation
C     parameterizations and the amine nucleation routine based on the 
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
      
      double precision mfrac(icomp) ! fraction of the nucleated mass that will be assigned
      double precision fn_all(2) !Magnitude of nucleation rates resolved by pathway
                                    ! to each icomp species (passed to nucMassUpdate)

C-----EXTERNAL FUNCTIONS------------------------------------------------

      double precision aerodens_PSSA
      external aerodens_PSSA

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (neps=1E8, meps=1E-8)
      parameter (pi=3.14159)

C-----CODE--------------------------------------------------------------

      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23
      ! ACDC Lookup table ternary nuclation does not need this in ppt
c      nh3ppt= (1.0e+21*8.314)*Gci(srtnh4)*temp/(pres*boxvol*gmw(srtnh4))
      nh3_molec = Gci(srtnh4)/boxvol*1000.d0/17.d0*6.022d23 

      fn = 0.d0
      rnuc = 0.d0
      gtime = 0.d0
      fn_all = 0.

      ! if no nucleation occurs the final arrays will be same as the initial arrays
      Mkf=Mki
      Nkf=Nki
      Gcf=Gci

cdbg      print*,'h2so4',h2so4,'nh3ppt',nh3ppt

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (h2so4.gt.1.d4) then         

         if (nh3_molec.gt.1.d6.and.tern_nuc.eq.1) then
c$$$         if (nh3ppt.gt.0.1.and.tern_nuc.eq.1) then
c$$$            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc
            call tern_nucl_acdc(temp,rh,cs,h2so4,nh3_molec,fn,rnuc)

            if (fn.gt.0.d0) then
               !update mass and number
               !nuclei are assumed as ammonium bisulfate
               mfrac = (/0.8144, 0.0, 0.1856, 0.0/)
               call nuc_massupd(Nkf,Mkf,Gcf,nuc_bin,dt,fn,rnuc,mfrac)
               fn_all(1) = fn
            endif
         endif

         if (bin_nuc.eq.1) then
            call vehk_nucl(temp,rh,h2so4,fn,rnuc) !binary nuc

            if (fn.gt.0.d0) then
               !update mass and number
               !nuclei are assumed to be sulfuric acid
               mfrac = (/1.0, 0.0, 0.0, 0.0/)
               call nuc_massupd(Nkf,Mkf,Gcf,nuc_bin,dt,fn,rnuc,mfrac)
               fn_all(2) = fn
            endif
         endif    
      endif

C     if nucleation occured, see how many particles grow to join the first size
C     section
c$$$      if (fn.gt.0.d0) then
cdbg         print*, 'entered fn'
cdbg         print*, 'fn', fn, 'rnuc',rnuc
         ! get the time it takes to grow to first size bin
Cjrp         d1 = rnuc*2.d0*1E-9
Cjrp         d2 = (6.d0*sqrt(xk(1)*xk(2))/1800.d0/pi)**(1.d0/3.d0)
Cjrp         print*, 'd1', d1, 'd2', d2
Cjrp         call getGrowthTime(d1,d2,Gc(srtso4),temp,
Cjrp     &        boxvol,gtime)
Cjrp         print*, 'gtime', gtime
Cjrp         ! get the first order loss rate of these particles due to coagulation
Cjrp         ! get loss rate for cluster size
Cjrp         d1 = rnuc*2.d0*1E-9
Cjrp         call getCoagLoss(d1,ltc1)
Cjrp         print*, 'd1', d1, 'ltc1', ltc1
Cjrp         ! we want the loss rate for particles at the size where they are added
Cjrp         ! to bin 1
Cjrp         Mktot = 0.d0
Cjrp         do j=1,icomp
Cjrp            Mktot=Mktot+Mk(1,j)
Cjrp         enddo
Cjrp         Mktot=Mktot+2.d0*Mk(1,srtso4)/96.d0-Mk(1,srtnh4)/18.d0 ! account for h+
Cjrp
Cjrp         if (Mktot.gt.meps)then
Cjrp            density=aerodens_PSSA(Mk(1,srtso4),0.d0,Mk(1,srtnh4),
Cjrp     &           0.d0,Mk(1,srth2o)) !assume bisulfate
Cjrp         else
Cjrp            density = 1400.
Cjrp         endif
Cjrp
Cjrp         if(Nk(1).gt.neps .and. Mktot.gt.meps)then
Cjrp            mp=Mktot/Nk(1)
Cjrp         else
Cjrp            mp=sqrt(xk(1)*xk(2))
Cjrp         endif
Cjrp         
Cjrp         d1 = (6.d0*mp/density/pi)**(1.d0/3.d0)
Cjrp         call getCoagLoss(d1,ltc2)
Cjrp         print*, 'd1', d1, 'ltc2', ltc2 
Cjrp
Cjrp         ! take average of loss rate constants
Cjrp         ltc = (ltc1 + ltc2)/2.d0
Cjrp         print*, 'ltc', ltc
Cjrp
Cjrp	 if (gtime.lt.0.d0)then
Cjrp	    gtime=0.d0
Cjrp         endif
Cjrp
Cjrp         frac = exp(-ltc*gtime)
Cjrp         print*,'frac', frac
Cjrpc         frac = 1.d0

cdbg         print*,'Nk_NUC',Nki(nuc_bin),Nkf(nuc_bin)
c$$$         Gcf(srtso4) = Gci(srtso4)! - (Mkf(nuc_bin,srtso4)-mold)
c$$$                                  !PSSA decide Gc as a diagnostic way


C there is a chance that Gcf will go less than zero because we are artificially growing
C particles into the first size bin.  don''t let it go less than zero.
c         if (Gcf(srtso4).lt.0.d0)then
c            Mkf(nuc_bin,srtso4) = Mki(nuc_bin,srtso4) + 
c     &           Gci(srtso4)*96./98.
c            Nkf(nuc_bin) = Nki(1) + Gci(srtso4)*96./98./mnuc
c            Gcf(srtso4) = 0.d0
c            print*,'used up h2so4 in nuc subroutine'
c         endif

c$$$      endif

      RETURN
      END

