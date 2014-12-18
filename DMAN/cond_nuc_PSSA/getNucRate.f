

C     **************************************************
C     *  getNucRate                                    *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine calls the Vehkamaki 2002 and Napari 2002 nucleation
C     parameterizations and gets the binary and ternary nucleation rates.
C     November 2014: also amine nucleation
C     December 2014: calls the ternary nucleation rate routine based on ACDC
C     instead of the Napari param.
C     JJ December 2014: this subroutine replaced by getNucMass.f: since 
C     several nucleation pathways can in the current version occur at the same time
C     a subroutine that returns total nucleated mass is more meaningful for the
C     purposes of getH2SO4conc.f

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water

C-----OUTPUTS-----------------------------------------------------------

C     fn - nucleation rate [# cm-3 s-1]
C     rnuc - radius of nuclei [nm]
C     nflg - says if nucleation happend

      SUBROUTINE getNucRate(Gci,fn,rnuc,nflg,cs,dmappt)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer j,i,k
      double precision Gci(icomp-1)
      double precision fn       ! nucleation rate cm-3 s-1
      double precision rnuc     ! critical nucleation radius [nm]
      logical nflg

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision nh3ppt   ! gas phase ammonia in pptv
      double precision nh3_molec! gas phase ammonia in molec cm-3
      double precision dmappt   ! gas phase dimethyl amine in pptv
      double precision dma_molec! gas phase dimethyl amine in molec cm-3
      double precision h2so4    ! gas phase h2so4 in molec cc-1
      double precision cs       ! condensation sink [s-1]
      double precision fn2      ! temporary nucleation rate storage cm-3 s-1

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      !Convert Gas concentrations from kg/(grid cell) to...
      h2so4  = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23 ![molec cm-3]
      !tern_nuc_acdc.f needs nh3 in molec cm-3
c      nh3ppt = (1.0e+21*8.314)*Gci(srtnh4)*temp/(pres*boxvol*gmw(srtnh4)) ![ppt]
      nh3_molec = Gci(srtnh4)/boxvol*1000.d0/17.d0*6.022d23  ![molec cm-3]
      dma_molec = dmappt*1.e-18 * pres/8.314/temp * 6.022d23 !molec cm-3

      fn = 0.d0
      rnuc = 0.d0
c      fn2 = 0.d0
      nflg=.false.

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
C     Total nucleation rate is sum from all the nucleation pathways
      if (h2so4.gt.1.d4) then
         if (amine_nuc.eq.1.and.dma_molec.gt.1.e4) then
            call amine_nucl(temp,cs,h2so4,dma_molec,fn,rnuc) !amine nuc
            nflg=.true.
c            fn=fn+fn2
c         endif

	 else if ((nh3_molec.gt.1.d6).and.(tern_nuc.eq.1)) then
c            print*, 'napari'
c            call napa_nucl(temp,rh,h2so4,nh3ppt,fn2,rnuc) !ternary nuc
            call tern_nucl_acdc(temp,rh,cs,h2so4,nh3_molec,fn,rnuc)
            nflg=.true.
c            fn=fn+fn2
c         endif

         else if (bin_nuc.eq.1) then
c            print*, 'vehk'
            call vehk_nucl(temp,rh,h2so4,fn,rnuc) !binary nuc
            if (fn2.gt.1.0d-7)then  ! JJ changed this from 1d-6 to the actual lower limit
                                    ! of the Vehkamäki param. Would make more sense to have 
                                    ! vehk_nucl do this               
               nflg=.true.
c               fn=fn+fn2
            endif
         endif
      endif

      return
      end
