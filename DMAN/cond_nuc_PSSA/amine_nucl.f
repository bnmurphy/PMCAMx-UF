
C     **************************************************
C     *  amine_nucl                                     *
C     **************************************************

C     WRITTEN BY Jan and Ben, October 2014

C     This subroutine calculates the ternary nucleation rate and radius of the 
C     critical nucleation cluster using the parameterization of...

c     Napari, I., M. Noppel, H. Vehkamaki, and M. Kulmala. "Parametrization of 
c     Ternary Nucleation Rates for H2so4-Nh3-H2o Vapors." Journal of Geophysical 
c     Research-Atmospheres 107, no. D19 (2002).

      SUBROUTINE amine_nucl(tempi,rhi,cnai,dma_i,fn,rnuc)

      IMPLICIT NONE

C-----INPUTS------------------------------------------------------------

      double precision tempi                ! temperature of air [K]
      double precision rhi                  ! relative humidity of air as a fraction
      double precision cnai                 ! concentration of gas phase sulfuric acid [molec cm-3]
      double precision dma_i                ! concentration of dimethyl amine [ppt]

C-----OUTPUTS-----------------------------------------------------------

      double precision fn                   ! nucleation rate [cm-3 s-1]
      double precision rnuc                 ! critical cluster radius [nm]
c      double precision tonset

C-----INCLUDE FILES-----------------------------------------------------

C-----ARGUMENT DECLARATIONS---------------------------------------------

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision fnl                  ! natural log of nucleation rate
      double precision temp                 ! temperature of air [K]
      double precision rh                   ! relative humidity of air as a fraction
      double precision cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      double precision dma                  ! concentration of gas phase dimethyl amine [ppt]
c
      integer ii                 ! counter

C
      temp = tempi
      rh   = rhi
      cna  = cnai
      dma  = dma_i
       
c     Check for Lower Bounds on LookUp Table
      if ((cna .lt. 1.d4).or.(dma.lt.0.1)) then ! limit sulf acid and nh3 conc
c      if ((cna .lt. 5.d4).or.(nh3ppt.lt.0.1)) then ! limit sulf acid and nh3 conc

         fn = 0.
         rnuc = 1
         goto 10
      endif  
      
      if (cna .gt. 1.0d9) cna=1.0d9 ! limit sulfuric acid conc
      if (nh3ppt .gt. 100.) nh3ppt=100. ! limit nh3ppt
c      if (nh3ppt .gt. 1000.) nh3ppt=1000. ! limit nh3ppt
      if (temp .lt. 240.) temp=240. ! limit temp
c      if (temp .lt. 235.) temp=235. ! limit temp
      if (temp .gt. 300.) temp=300. ! limit temp
c      if (temp .gt. 295.) temp=295. ! limit temp
      if (rh .lt. 0.05) rh=0.05 ! limit rh 
      if (rh .gt. 0.95) rh=0.95 ! limit rh

      do i=1,20
         fa(i)=aa0(i)+a1(i)*temp+a2(i)*temp**2.+a3(i)*temp**3.
      enddo

      fnl=-84.7551+fa(1)/log(cna)+fa(2)*log(cna)+fa(3)*(log(cna))**2.
     &  +fa(4)*log(nh3ppt)+fa(5)*(log(nh3ppt))**2.+fa(6)*rh
     &  +fa(7)*log(rh)+fa(8)*log(nh3ppt)/log(cna)+fa(9)*log(nh3ppt)
     &  *log(cna)+fa(10)*rh*log(cna)+fa(11)*rh/log(cna)
     &  +fa(12)*rh
     &  *log(nh3ppt)+fa(13)*log(rh)/log(cna)+fa(14)*log(rh)
     &  *log(nh3ppt)+fa(15)*(log(nh3ppt))**2./log(cna)+fa(16)*log(cna)
     &  *(log(nh3ppt))**2.+fa(17)*(log(cna))**2.*log(nh3ppt)
     &  +fa(18)*rh
     &  *(log(nh3ppt))**2.+fa(19)*rh*log(nh3ppt)/log(cna)+fa(20)
     &  *(log(cna))**2.*(log(nh3ppt))**2.
c    c
      fn=exp(fnl)
      fn=fn*1.0e-6 ! For adjusting

cc   Cap at 10^6 particles/cm3-s, limit for parameterization
       if (fn.gt.1.0d6) then
        fn=1.0d6
        fnl=log(fn)
       endif

      rnuc=0.141027-0.00122625*fnl-7.82211d-6*fnl**2.
     &     -0.00156727*temp-0.00003076*temp*fnl
     &     +0.0000108375*temp**2.
     
 

      return

      end

