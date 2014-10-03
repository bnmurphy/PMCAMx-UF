
C     *************************************************
C     * Nk2pn                                         *
C     *************************************************

C     WRITTEN BY JaeGun Jung, November 2007

C     This subroutine converts number of particles in a grid cell, Nk
C     to particle number concentration, pn, and mass of particles in 
C     a grid cell, Mk to mass concentration, pm. Then, calculate cn 
C     from pn.

C-----INPUTS------------------------------------------------------------

C     Nk(ibins) - number of particles per size bin in a grid cell
C     Mk(ibins,icomp) - mass of a given species per size bin/grid cell

C-----OUTPUTS-----------------------------------------------------------

C     pn(nsect) - number concentrations per size bin as particles per cm3
C     pm(nsect,icomp) - mass concentrations per size bin as ug per m3

      SUBROUTINE Nk2pn(pn,cn)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real pn(nsect)
      real cn

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k
      real cvt

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------
  
C-----CODE--------------------------------------------------------------

cdbg      print*,'In Nk2pn'
cdbg      print*,'xko=',xko
cdbg      print*,'xk=',xk
cdbg      print*,'Nk=',Nk
cdbg      print*,'Mk=',Mk
cdbg      pause

      ! Initialization
      cn = 0.0

      ! Set a conversion factor
      cvt = 1./boxvol *1.0e+15 !conversion from kg/cm3 to ug/m3

      ! Do conversion using a box interpolation
      j=1
      do i=1,nsect
         if (xko(i+1).lt.xk(j+1)) then
            pn(i)=Nk(j)/boxvol*(log(xko(i+1))
     &       -log(xko(i)))/(log(xk(j+1))-log(xk(j)))
            do k=1,icomp-1
               pm(i,k)=Mk(j,k)*(log(xko(i+1))
     &          -log(xko(i)))/(log(xk(j+1))-log(xk(j)))
     &          *cvt ! ug/m3
            enddo
         else
            pn(i)=Nk(j)/boxvol*(log(xk(j+1))-log(xko(i)))
     &       /(log(xk(j+1))-log(xk(j)))+Nk(j+1)
     &       /boxvol*(log(xko(i+1))-log(xk(j+1)))
     &       /(log(xk(j+2))-log(xk(j+1)))
            do k=1,icomp-1
               pm(i,k)=(Mk(j,k)*(log(xk(j+1))-log(xko(i)))
     &          /(log(xk(j+1))-log(xk(j)))+Mk(j+1,k)
     &          *(log(xko(i+1))-log(xk(j+1)))
     &          /(log(xk(j+2))-log(xk(j+1))))*cvt !ug/m3
            enddo
            j=j+1
         endif
      enddo

      ! Get cn from pn
      do i=1,nsect
        if (dpmean(i) .lt. 0.1) then
          cn = cn + pn(i)
        endif
      enddo

      RETURN
      END
