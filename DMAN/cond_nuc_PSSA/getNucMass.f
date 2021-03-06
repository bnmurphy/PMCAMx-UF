

C     **************************************************
C     *  getNucMass                                    *
C     **************************************************

C     WRITTEN BY Jan Julin, December 2014
C     on top of the getNucRate subroutine by Jeff Pierce, April 2007

C     This subroutine calls the Vehkamaki 2002 nucleation
C     parameterization and the ternary and amine nucleation routines based
C     on ACDC and returns the total nucleated mass.

C     Called by getH2SO4conc


C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water

C-----OUTPUTS-----------------------------------------------------------

C     massnuc - nucleated mass [kg s-1 box-1]
C     nflg - says if nucleation happened

      SUBROUTINE getNucMass(Gci,massnuc,nflg,cs)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer j,i,k
      double precision Gci(icomp-1)
      double precision massnuc ! mass being removed by nucleation [kg s-1 box-1]
      logical nflg

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision nh3ppt   ! gas phase ammonia in pptv
      double precision nh3_molec! gas phase ammonia in molec cm-3
      double precision h2so4    ! gas phase h2so4 in molec cc-1
      double precision cs       ! condensation sink [s-1]
      double precision fn       ! nucleation rate cm-3 s-1
      double precision rnuc     ! critical nucleation radius [nm]
      double precision massnuc2 ! variable to store mass from individual
                                ! nucleation pathways [kg s-1 box-1]
      double precision pi
      double precision N_Vehk   ! critical cluster from Vehkamäki parameterization
      double precision amu      !atomic mass unit [kg]

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654)
      parameter (amu=1.660539040d-27)

C-----CODE--------------------------------------------------------------

      !Convert Gas concentrations from kg/(grid cell) to...
      h2so4  = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23 ![molec cm-3]
c      nh3ppt= (1.0e+21*8.314)*Gci(srtnh4)*temp/(pres*boxvol*gmw(srtnh4)) ![ppt]      
      nh3_molec = Gci(srtnh4)/boxvol*1000.d0/17.d0*6.022d23  ![molec cm-3]

      fn = 0.d0
      rnuc = 0.d0
      massnuc = 0.d0

      nflg=.false.

C     if requirements for nucleation are met, call nucleation subroutines
C     and calculate nucleated mass
C     Total mass is sum from all the nucleation pathways

      if (h2so4.gt.minval(tern_nuc_tbl_H2SO4).and.
     &     nh3_molec.gt.minval(tern_nuc_tbl_NH3).and.tern_nuc.eq.1) then
c$$$         if ((h2so4.gt.1.d4).and.(nh3ppt.gt.0.1).and.(tern_nuc.eq.1)) then
c$$$            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc
         call tern_nucl_acdc(temp,rh,cs,h2so4,nh3_molec,fn) !ternary nuc
         nflg=.true.
         !ACDC formation rate (the 2014-12-04 lookup table) corresponds to 
         !clusters with 4 acids+3 NH3
         massnuc2 = 4.d0*gmw(srtso4)*amu*fn*boxvol
         massnuc = massnuc + massnuc2
      endif

      if (h2so4.gt.1.d4.and.bin_nuc.eq.1) then
         call vehk_nucl(temp,rh,h2so4,fn,rnuc,N_Vehk) !binary nuc
         if (fn.gt.0.d0 .and. N_Vehk.ge.4.d0) then 
            nflg=.true.
            !nuclei are assumed to be sulfuric acid as in nucleation.f
            massnuc2 = N_Vehk*gmw(srtso4)*amu*fn*boxvol
            massnuc = massnuc + massnuc2
         endif
      endif

      return
      end
