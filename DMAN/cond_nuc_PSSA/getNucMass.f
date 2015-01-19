

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
C     nflg - says if nucleation happend

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

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654)

C-----CODE--------------------------------------------------------------

      !Convert Gas concentrations from kg/(grid cell) to...
      h2so4  = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23 ![molec cm-3]
      nh3ppt= (1.0e+21*8.314)*Gci(srtnh4)*temp/(pres*boxvol*gmw(srtnh4)) ![ppt]      
c      nh3_molec = Gci(srtnh4)/boxvol*1000.d0/17.d0*6.022d23  ![molec cm-3]

      fn = 0.d0
      rnuc = 0.d0
      massnuc = 0.d0

      nflg=.false.

C     if requirements for nucleation are met, call nucleation subroutines
C     and calculate nucleated mass
C     Total mass is sum from all the nucleation pathways
      if (h2so4.gt.1.d4) then

c	 if ((nh3_molec.gt.1.d6).and.(tern_nuc.eq.1)) then
         if ((nh3ppt.gt.0.1).and.(tern_nuc.eq.1)) then
            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc
            nflg=.true.
            massnuc2 = 4.d0/3.d0*pi*(rnuc*1.d-9)**3*1350.*fn*boxvol*
     &        98.d0/96.d0
            massnuc = massnuc + massnuc2
         endif

         if (bin_nuc.eq.1) then
            call vehk_nucl(temp,rh,h2so4,fn,rnuc) !binary nuc
            if (fn.gt.1.0d-7)then ! JJ changed this from 1d-6 to the actual lower limit
                                    ! of the Vehkamäki param. This should actually be also
                                    ! checked in the vehk_nucl routine?
               nflg=.true.
               massnuc2 = 4.d0/3.d0*pi*(rnuc*1.d-9)**3*1350.*fn*boxvol*
     &           98.d0/96.d0
               massnuc = massnuc + massnuc2
            endif
         endif
      endif

      return
      end
