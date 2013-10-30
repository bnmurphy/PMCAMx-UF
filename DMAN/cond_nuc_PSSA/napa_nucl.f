
C     **************************************************
C     *  napa_nucl                                     *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine calculates the ternary nucleation rate and radius of the 
C     critical nucleation cluster using the parameterization of...

c     Napari, I., M. Noppel, H. Vehkamaki, and M. Kulmala. "Parametrization of 
c     Ternary Nucleation Rates for H2so4-Nh3-H2o Vapors." Journal of Geophysical 
c     Research-Atmospheres 107, no. D19 (2002).

      SUBROUTINE napa_nucl(tempi,rhi,cnai,nh3ppti,fn,rnuc)

      IMPLICIT NONE

C-----INPUTS------------------------------------------------------------

      double precision tempi                ! temperature of air [K]
      double precision rhi                  ! relative humidity of air as a fraction
      double precision cnai                 ! concentration of gas phase sulfuric acid [molec cm-3]
      double precision nh3ppti              ! concentration of gas phase ammonia

C-----OUTPUTS-----------------------------------------------------------

      double precision fn                   ! nucleation rate [cm-3 s-1]
      double precision rnuc                 ! critical cluster radius [nm]

C-----INCLUDE FILES-----------------------------------------------------

C-----ARGUMENT DECLARATIONS---------------------------------------------

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision aa0(20),a1(20),a2(20),a3(20),fa(20) ! set parameters
      double precision fnl                  ! natural log of nucleation rate
      double precision temp                 ! temperature of air [K]
      double precision rh                   ! relative humidity of air as a fraction
      double precision cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      double precision nh3ppt               ! concentration of gas phase ammonia
c
      integer i                 ! counter

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      data aa0 /-0.355297, 3.13735, 19.0359, 1.07605, 6.0916,
     $         0.31176, -0.0200738, 0.165536,
     $         6.52645, 3.68024, -0.066514, 0.65874,
     $         0.0599321, -0.732731, 0.728429, 41.3016,
     $         -0.160336, 8.57868, 0.0530167, -2.32736        /

      data a1 /-33.8449, -0.772861, -0.170957, 1.48932, -1.25378,
     $         1.64009, -0.752115, 3.26623, -0.258002, -0.204098,
     $         -7.82382, 0.190542, 5.96475, -0.0184179, 3.64736,
     $         -0.35752, 0.00889881, -0.112358, -1.98815, 0.0234646/
     
      data a2 /0.34536, 0.00561204, 0.000479808, -0.00796052,
     $         0.00939836, -0.00343852, 0.00525813, -0.0489703,
     $         0.00143456, 0.00106259, 0.0122938, -0.00165718,
     $         -0.0362432, 0.000147186, -0.027422, 0.000904383,
     $         -5.39514d-05, 0.000472626, 0.0157827, -0.000076519/
     
      data a3 /-0.000824007, -9.74576d-06, -4.14699d-07, 7.61229d-06,
     $         -1.74927d-05, -1.09753d-05, -8.98038d-06, 0.000146967,
     $         -2.02036d-06, -1.2656d-06, 6.18554d-05, 3.41744d-06,
     $         4.93337d-05, -2.37711d-07, 4.93478d-05, -5.73788d-07, 
     $         8.39522d-08, -6.48365d-07, -2.93564d-05, 8.0459d-08   /

C-----CODE--------------------------------------------------------------

cdbg      print*,'temp=',tempi,' rh=',rhi,' cna=',cnai,' nh3ppt=',nh3ppti
      temp=tempi
      rh=rhi
      cna=cnai
      nh3ppt=nh3ppti

c     Napari's parameterization is only valid within limited area
      if ((cna .lt. 1.d4).or.(nh3ppt.lt.0.1)) then ! limit sulf acid and nh3 conc
         fn = 0.
         rnuc = 1
         goto 10
      endif  
      if (cna .gt. 1.0d9) cna=1.0d9 ! limit sulfuric acid conc
      if (nh3ppt .gt. 100.) nh3ppt=100. ! limit nh3ppt  
      if (temp .lt. 240.) temp=240. ! limit temp
      if (temp .gt. 300.) temp=300. ! limit temp
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
c
c
      fn=exp(fnl)
      fn=fn*1.0e-5 ! For adjusting
c   Cap at 10^6 particles/cm3-s, limit for parameterization
      if (fn.gt.1.0d6) then
        fn=1.0d6
        fnl=log(fn)
      endif

      rnuc=0.141027-0.00122625*fnl-7.82211d-6*fnl**2.
     &     -0.00156727*temp-0.00003076*temp*fnl
     &     +0.0000108375*temp**2.

 10   return
      end
