
C     *************************************************
C     * initconv                                      *
C     *************************************************

C     WRITTEN BY JaeGun Jung, November 2007

C     This subroutine converts number concentration variables read from 
C     input files to TOMAS number and mass varibles. The subroutine is
C     implemented at the beginning of simulation for initialization. 
C     A box interpolation is used for conversion between the variables.
C     For the box interpolation, see figure 1 in Koo et al. (2003)*.

C     * Koo, B., Gaydos, M. T., and Pandis, N. S., Evaluation of the
C     Equilibrium, Dynamic, and Hybrid Aerosol Modeling Approaches,
C     Aero. Sci. Tech., 37, 53-64, 2003.

C-----INPUTS------------------------------------------------------------

C     Initial variables of
C     ====================

C     pn(nsect) - number concentrations per size bin as particles per cm3

C-----OUTPUTS-----------------------------------------------------------

C     Nk(ibins) - number of particles per size bin in a grid cell
C     Mk(ibins,icomp) - mass of a given species per size bin/grid cell
 
      SUBROUTINE initconv(pn,ygas)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'
      include 'IO.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real pn(nsect)
      real ygas(ngas)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k

      real cvt ! a unit converter from kg/cm3 to ug/m3
      real ratiobm
      real R   ! gas constant J/mol-K

      double precision density
      double precision Neps

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS--------------------------------------------

      parameter(R=8.314) !J/mol-K
      parameter(Neps=1.0e-20)

C-----CODE-------------------------------------------------------------

      density=1400.0  ! density of H2SO4-H2O droplets (20% solution)

      !Set a conversion factor
      cvt = 1./boxvol*1.0e15 ! conversion from kg/cm3 to ug/m3

C-----Mass of particles from input files
      do i=1,nsect+1
        xko(i)=pi6*dpbound(i)**3.*density*1.0e-18 ! kg
      enddo

C-----Output header for coagulation and condensation tests
      if (icoag_test .eq. 1) then
        do i=1, ibins
          read(idat+7,*)Nk(i) 
        enddo
      elseif (icond_test .eq. 1) then
        do i=1, ibins
          read(idat+8,*)Nk(i) !for boxvol=3.0e+20 cm3
        enddo
      else

C-----Convert the input number concentrations, pn(nsect) to the processing 
C    number concentrations, Nk(ibins)
        j=1
        do i=1,nsect
          if (xko(i).lt.xk(j+1)) then
            Nk(j)=Nk(j)+pn(i)*boxvol
          else
            Nk(j)=Nk(j)+(pn(i)*(xk(j+1)-xko(i-1))/(xko(i)-xko(i-1)))
     &           *boxvol
            Nk(j+1)=Nk(j+1)+(pn(i)*(xko(i)-xk(j+1))/(xko(i)-xko(i-1)))
     &           *boxvol
            j=j+1
          endif
        enddo
        !Put Neps if the input variables are less than processing ones
        if (j.lt.ibins) then 
          do i=(j+1),ibins
            Nk(i)=Neps 
          enddo   
        endif
      endif

C-----Initialize mass concentrations
      !Account only sulfate when a test are on
      if ((icoag_test .eq. 1) .or. (icond_test .eq. 1))then
        do k=1,ibins
          Mk(k,srtso4)=1.414*xk(k)*Nk(k)
          Mk(k,srtorg)=0.
          Mk(k,srtnh3)=0.
          Mk(k,srth2o)=0.
        enddo
      else !Assume half of mass consist of ammonium sulfate 
        do k=1,ibins
          !0.72723 is a fraction of sulfate in ammonium sulfate
          Mk(k,srtso4)=0.5*1.414*xk(k)*Nk(k)*0.727273  ! sulfate 
          Mk(k,srtorg)=0.5*1.414*xk(k)*Nk(k)  ! organic matter
          !0.272727 is a fraction of ammonium in ammonium sulfate
          Mk(k,srtnh3)=0.5*1.414*xk(k)*Nk(k)*0.272727 ! ammonium
          Mk(k,srth2o)=0.0 ! water
        enddo
      endif

C-----Copy gas unit from ppt to kg
      Gc(srtso4)=1.0e-21*pres*boxvol*gmw(srtso4)/(R*temp)*ygas(mgsvi) !H2SO4
      Gc(srtnh3)=1.0e-21*pres*boxvol*gmw(srtnh3)/(R*temp)*ygas(mgnh3) !NH3


      RETURN
      END
