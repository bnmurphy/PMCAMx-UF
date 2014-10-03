
C     *************************************************
C     * Mtot_sub                                      *
C     *************************************************

C     WRITTEN BY JaeGun Jung, November 2007

C     This subroutine calculates total mass concentrations of all size 
C     bins with respect to each speices.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

C     Mtot - total mass concentration of all size bins for each speices

      SUBROUTINE Mtot_sub(Mtot)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'
      
C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Mtot(icomp)

C-----VARIABLE DECLARATIONS---------------------------------------------

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      integer i ! a counter for size bins
      integer j ! a counter for species

      real cvt ! a unit converter from kg/cm3 to ug/m3

      cvt = 1./boxvol*1.0e+15 ! conversion from kg/cm3 to ug/m3

      do j=1, icomp
        Mtot(j)=0.
      enddo
      do i=1,ibins
        do j=1, icomp
          Mtot(j)=Mtot(j)+Mk(i,j)*cvt !ug/m3
        enddo
      enddo


      RETURN
      END
