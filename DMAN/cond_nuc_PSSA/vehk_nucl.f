
C     **************************************************
C     *  vehk_nucl                                     *
C     **************************************************

C     WRITTEN BY Jeff Pierce, April 2007

C     This subroutine calculates the binary nucleation rate and radius of the 
C     critical nucleation cluster using the parameterization of...

c     Vehkamaki, H., M. Kulmala, I. Napari, K. E. J. Lehtinen, C. Timmreck, 
C     M. Noppel, and A. Laaksonen. "An Improved Parameterization for Sulfuric 
C     Acid-Water Nucleation Rates for Tropospheric and Stratospheric Conditions." 
C     Journal of Geophysical Research-Atmospheres 107, no. D22 (2002).

      SUBROUTINE vehk_nucl(tempi,rhi,cnai,fn,rnuc)

      IMPLICIT NONE

C-----INPUTS------------------------------------------------------------

      double precision tempi                ! temperature of air [K]
      double precision rhi                  ! relative humidity of air as a fraction
      double precision cnai                 ! concentration of gas phase sulfuric acid [molec cm-3]

cjgj      real tempi                ! temperature of air [K]
cjgj      real rhi                  ! relative humidity of air as a fraction
cjgj      real cnai                 ! concentration of gas phase sulfuric acid [molec cm-3]
C-----OUTPUTS-----------------------------------------------------------

      double precision fn                   ! nucleation rate [cm-3 s-1]
      double precision rnuc                 ! critical cluster radius [nm]

C-----INCLUDE FILES-----------------------------------------------------

C-----ARGUMENT DECLARATIONS---------------------------------------------

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision fb0(10),fb1(10),fb2(10),fb3(10),fb4(10),fb(10)
      double precision gb0(10),gb1(10),gb2(10),gb3(10),gb4(10),gb(10) ! set parameters
      double precision temp                 ! temperature of air [K]
      double precision rh                   ! relative humidity of air as a fraction
      double precision cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
cdbg      real temp                 ! temperature of air [K]
cdbg      real rh                   ! relative humidity of air as a fraction
cdbg      real cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      double precision xstar                ! mole fraction sulfuric acid in cluster
      double precision ntot                 ! total number of molecules in cluster
      integer i                 ! counter

C-----ADJUSTABLE PARAMETERS---------------------------------------------

c     Nucleation Rate Coefficients
c
      data fb0 /0.14309, 0.117489, -0.215554, -3.58856, 1.14598,
     $          2.15855, 1.6241, 9.71682, -1.05611, -0.148712        /
      data fb1 /2.21956, 0.462532, -0.0810269, 0.049508, -0.600796,
     $       0.0808121, -0.0160106, -0.115048, 0.00903378, 0.00283508/
      data fb2 /-0.0273911, -0.0118059, 0.00143581, -0.00021382, 
     $       0.00864245, -0.000407382, 0.0000377124, 0.000157098,
     $       -0.0000198417, -9.24619d-6 /
      data fb3 /0.0000722811, 0.0000404196, -4.7758d-6, 3.10801d-7,
     $       -0.0000228947, -4.01957d-7, 3.21794d-8, 4.00914d-7,
     $       2.46048d-8, 5.00427d-9 /
      data fb4 /5.91822, 15.7963, -2.91297, -0.0293333, -8.44985,
     $       0.721326, -0.0113255, 0.71186, -0.0579087, -0.0127081  / 

c     Coefficients of total number of molecules in cluster     
c
      data gb0 /-0.00295413, -0.00205064, 0.00322308, 0.0474323,
     $         -0.0125211, -0.038546, -0.0183749, -0.0619974,
     $         0.0121827, 0.000320184 /
      data gb1 /-0.0976834, -0.00758504, 0.000852637, -0.000625104,
     $         0.00580655, -0.000672316, 0.000172072, 0.000906958,
     $         -0.00010665, -0.0000174762 /      
      data gb2 /0.00102485, 0.000192654, -0.0000154757, 2.65066d-6,
     $         -0.000101674, 2.60288d-6, -3.71766d-7, -9.11728d-7,
     $         2.5346d-7, 6.06504d-8 /
      data gb3 /-2.18646d-6, -6.7043d-7, 5.66661d-8, -3.67471d-9,
     $         2.88195d-7, 1.19416d-8, -5.14875d-10, -5.36796d-9,
     $         -3.63519d-10, -1.42177d-11 /
      data gb4 /-0.101717, -0.255774, 0.0338444, -0.000267251,
     $         0.0942243, -0.00851515, 0.00026866, -0.00774234,
     $         0.000610065, 0.000135751 /

C-----CODE--------------------------------------------------------------

      temp=tempi
      rh=rhi
      cna=cnai

c     Respect the limits of the parameterization
      if (cna .lt. 1.d4) then ! limit sulf acid conc
         fn = 0.
         rnuc = 1
         goto 10
      endif
      if (cna .gt. 1.0d11) cna=1.0e11 ! limit sulfuric acid conc  
      if (temp .lt. 230.15) temp=230.15 ! limit temp
      if (temp .gt. 305.15) temp=305.15 ! limit temp
      if (rh .lt. 1d-4) rh=1d-4 ! limit rh
      if (rh .gt. 1.) rh=1. ! limit rh
c
c     Mole fraction of sulfuric acid
      xstar=0.740997-0.00266379*temp-0.00349998*log(cna)
     &   +0.0000504022*temp*log(cna)+0.00201048*log(rh)
     &   -0.000183289*temp*log(rh)+0.00157407*(log(rh))**2.
     &   -0.0000179059*temp*(log(rh))**2.
     &   +0.000184403*(log(rh))**3.
     &   -1.50345d-6*temp*(log(rh))**3.
c 
c     Nucleation rate coefficients 
      do i=1, 10
         fb(i) = fb0(i)+fb1(i)*temp+fb2(i)*temp**2.
     &        +fb3(i)*temp**3.+fb4(i)/xstar
      enddo
c
c     Nucleation rate (1/cm3-s)
      fn = exp(fb(1)+fb(2)*log(rh)+fb(3)*(log(rh))**2.
     &    +fb(4)*(log(rh))**3.+fb(5)*log(cna)
     &    +fb(6)*log(rh)*log(cna)+fb(7)*(log(rh))**2.*log(cna)
     &    +fb(8)*(log(cna))**2.+fb(9)*log(rh)*(log(cna))**2.
     &    +fb(10)*(log(cna))**3.)
c
c   Cap at 10^6 particles/s, limit for parameterization
      if (fn.gt.1.0d6) then
         fn=1.0d6
      endif
c
c     Coefficients of total number of molecules in cluster 
      do i=1, 10
         gb(i) = gb0(i)+gb1(i)*temp+gb2(i)*temp**2.
     &        +gb3(i)*temp**3.+gb4(i)/xstar
      enddo
c     Total number of molecules in cluster
      ntot=exp(gb(1)+gb(2)*log(rh)+gb(3)*(log(rh))**2.
     &    +gb(4)*(log(rh))**3.+gb(5)*log(cna)
     &    +gb(6)*log(rh)*log(cna)+gb(7)*log(rh)**2.*log(cna)
     &    +gb(8)*(log(cna))**2.+gb(9)*log(rh)*(log(cna))**2.
     &    +gb(10)*(log(cna))**3.)

c     cluster radius
      rnuc=exp(-1.6524245+0.42316402*xstar+0.3346648*log(ntot)) ! [nm]

 10   return
      end
