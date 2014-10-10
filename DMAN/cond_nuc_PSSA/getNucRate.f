

C     **************************************************
C     *  getNucRate                                    *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine calls the Vehkamaki 2002 and Napari 2002 nucleation
C     parameterizations and gets the binary and ternary nucleation rates.

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Gci(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water

C-----OUTPUTS-----------------------------------------------------------

C     fn - nucleation rate [# cm-3 s-1]
C     rnuc - radius of nuclei [nm]
C     nflg - says if nucleation happend

      SUBROUTINE getNucRate(Gci,fn,rnuc,nflg)

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
      double precision h2so4    ! gas phase h2so4 in molec cc-1

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      h2so4 = Gci(srtso4)/boxvol*1000.d0/98.d0*6.022d23
      nh3ppt= (1.0e+21*8.314)*Gci(srtnh4)*temp/(pres*boxvol*gmw(srtnh4))
cjgj	nh3ppt = NH3ppt_o
c      nh3ppt = Gci(srtnh4)/17.d0/(boxmass/29.d0)*1d12

      fn = 0.d0
      rnuc = 0.d0

c      print*,'h2so4',h2so4,'nh3ppt',nh3ppt

C     if requirements for nucleation are met, call nucleation subroutines
C     and get the nucleation rate and critical cluster size
      if (h2so4.gt.1.d4) then
         if (amine_nuc.eq.1.and.amineppt.gt.0.001) then
            call amine_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !amine nuc
            nflg=.true.
	 elseif ((nh3ppt.gt.0.1).and.(tern_nuc.eq.1)) then
c            print*, 'napari'
            call napa_nucl(temp,rh,h2so4,nh3ppt,fn,rnuc) !ternary nuc
            nflg=.true.
         elseif (bin_nuc.eq.1) then
c            print*, 'vehk'
            call vehk_nucl(temp,rh,h2so4,fn,rnuc) !binary nuc
            if (fn.gt.1.0d-6)then
               nflg=.true.
            else
               fn = 0.d0
               nflg=.false.
            endif

         else
            nflg=.false.
         endif
      else
         nflg=.false.
      endif

      return
      end
