
C     ************************************************
C     * pn2Nk                                        *
C     ************************************************

C     WRITTEN BY JaeGun Jung, November 2007

C     This subroutine adds number concentrations increased by nucleation 
C     to processing number of particles, Nk using by a box interpolation.
C     For the box interpolation, see figure 1 in Koo et al. (2003)*.

C     * Koo, B., Gaydos, M. T., and Pandis, N. S., Evaluation of the
C     Equilibrium, Dynamic, and Hybrid Aerosol Modeling Approaches,
C     Aero. Sci. Tech., 37, 53-64, 2003.

C-----INPUTS------------------------------------------------------------

C     pn_p(nsect) - number concentration increased by nucleation
C                 [=] particles cm-3
c
      SUBROUTINE pn2Nk(pn_p)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real pn_p(nsect) ! a number concentrations of particles

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k
      integer icsect

      real rso4, rnh4
      real density

      real neps, cvt

c-----EXTERNAL FUCTIONS-------------------------------------------------

      double precision aerodens_PSSA
      external aerodens_PSSA

C-----PARAMETERS--------------------------------------------------------

      parameter(icsect = 1) ! a currently-treated number of cells
      parameter(cvt = 1.0d-18) ! a convert factor (um3/m3 to m3/um3)
      parameter(rso4 = 0.8144) ! a sulfate mass ratio of nuclei 
      parameter(rnh4 = 0.1875) ! a ammonium mass ratio of nuclei 
                               ! I assumed that nuclei is ammonium 
                               ! bisulfate for rso4 and rnh4 from
                               !(Napari et al. 2002).
      
      !ADJUSTABLE PARAMETERS
      parameter(neps = 1.0e-20)

C-----CODE--------------------------------------------------------------

cdbg      density=aerodens_PSSA(Mk(srtso4),0.0,Mk(srtnh3),0.0,Mk(srth2o)) ! [=] kg/m3
      density = 1400 ! density of H2SO4-H2O droplets (20% solution)

      !First check pn_p whether it is negative
      do i = 1, icsect
         if (pn_p(i) .lt. 0.) then
            if (pn_p(i) .gt. -1.0d5) then
               pn_p(i) = neps
            else
               print*,'A pn_p is significantly negative.',pn_p(i)
               STOP
            endif
         endif       
      enddo

C-----Convert number from nucleation to processing variables
      j=1
      do i=1,icsect
        if (xko(i).lt.xk(j+1)) then
          Nk(j)=Nk(j)+pn_p(i)*boxvol
          Mk(j,srtso4)=Mk(j,srtso4)+pn_p(i)*density
     &                *boxvol*(dpmean(i)**3)*pi6*cvt*rso4*96./97.
          Mk(j,srtnh3)=Mk(j,srtso4)*rnh4
        else
          Nk(j)=Nk(j)+(pn_p(i)*(xk(j+1)
     &         -xko(i-1))/(xko(i)-xko(i-1)))*boxvol
          Nk(j+1)=Nk(j+1)+(pn_p(i)*(xko(i)-xk(j+1))
     &           /(xko(i)-xko(i-1)))*boxvol
          Mk(j,srtso4)=Mk(j,srtso4)+(pn_p(i)*(xk(j+1)
     &                -xko(i-1))/(xko(i)-xko(i-1)))*boxvol
     &                *density*(dpmean(i)**3)*pi6*cvt*rso4*96./97.
          Mk(j,srtnh3)=Mk(j,srtso4)*rnh4
          Mk(j+1,srtso4)=Mk(j+1,srtso4)+(pn_p(i)
     &                  *(xko(i)-xk(j+1))/(xko(i)-xko(i-1)))*boxvol
     &                  *density*(dpmean(i)**3)*pi6*cvt*rso4*96./97.
          Mk(j+1,srtnh3)=Mk(j+1,srtso4)*rnh4
          j=j+1
        endif
      enddo



      RETURN
      END
