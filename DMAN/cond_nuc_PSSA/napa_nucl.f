
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
c      double precision tonset

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
c
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
cc
c      data aa0 /-358.2337705052991,-980.923146020468,1200.472096232311,
c     &         -14.833042158178936,-4.39129415725234d6,4.905527742256349
c     $         ,-231375.56676032578,75061.15281456841,-3180.5610833308,
c     $         -100.21645273730675,5599.912337254629,2.360931724951942d6
c     $         ,16597.75554295064,-89.38961120336789,-629.7882041830943,
c     $         -732006.8180571689,40751.075322248245,-1911.0303773001353
c     $         ,2.792313345723013, 3.1712136610383244     /
c
c      data a1 /4.8630382337426985,10.054155220444462,-17.37107890065621,
c     $        0.2932631303555295,56383.93843154586,-0.05463019231872484
c     $         ,2919.2852552424706,-931.8802278173565,39.08268568672095,
c     $         0.977886555834732,-70.70896612937771,-29752.130254319443,
c     $         -175.2365504237746,1.153344219304926,7.772806552631709,
c     $         9100.06398573816,-501.66977622013934,23.6903969622286,
c     $         -0.03422552111802899,-0.037822330602328806      /
c
c      data a2 /-0.02175548069741675,-0.03306644502023841,
c     $     0.08170681335921742,-0.0016497524241142845,-239.835990963361
c     $     ,0.00020258394697064567,-12.286497122264588,3.863266220840964
c     $  ,-0.16048521066690752,-0.0030511783284506377,0.29788016132269466
c     $    ,125.04965118142027,0.6033215603167458,-0.004954549700267233,
c     $    -0.031974053936299256,-37.771091915932004,2.063469732254135,
c     $ -0.09807872005428583,0.00014019195277521142,0.0001500555743561457
c     $      /
c
c      data a3 /0.00003212869941055865,0.000034274041225891804,
c     $         -0.00012534476159729881,2.844074805239367d-6,
c     $         0.33765136625580167,-2.502406532869512d-7,
c     $         0.017249301826661612,-0.005349472062284983,
c     $         0.00022031380023793877,2.967320346100855d-6,
c     $         -0.00041866525019504,-0.1752996881934318,
c     $         -0.0006731787599587544,7.096309866238719d-6,
c     $         0.00004383764128775082,0.05235455395566905,
c     $         -0.002836873785758324,0.00013564560238552576,
c     $         -1.9201227328396297d-7,-1.9828365865570703d-7   /
c
C-----CODE--------------------------------------------------------------

cdbg      print*,'temp=',tempi,' rh=',rhi,' cna=',cnai,' nh3ppt=',nh3ppti
      temp=tempi
      rh=rhi
      cna=cnai
      nh3ppt=nh3ppti

c     Napari's parameterization is only valid within limited area
      if ((cna .lt. 1.d4).or.(nh3ppt.lt.0.1)) then ! limit sulf acid and nh3 conc
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
c
c      fnl=-12.861848898625231+fa(1)*rh+fa(2)*log(rh)+fa(3)*(log(cna))
c     &  +fa(4)*(log(cna))**2.+fa(5)/((log(cna))**2.)+fa(6)*(nh3ppt)
c     &  +fa(7)*log(nh3ppt)+fa(8)*(log(nh3ppt))**2.
c     &  +fa(9)*(log(nh3ppt))**3.+fa(10)*rh*log(nh3ppt)
c     &  +fa(11)*log(nh3ppt)*log(cna)+fa(12)*log(nh3ppt)/log(cna)
c     &  +fa(13)*log(rh)/log(cna)+fa(14)*log(rh)*log(nh3ppt)
c     &  +fa(15)*rh/(((nh3ppt)**3.)*log(cna))
c     &  +fa(16)*(log(nh3ppt))**2./log(cna)
c     &  +fa(17)*(log(nh3ppt))**3./log(cna)
c     &  +fa(18)*log(cna)*(log(nh3ppt))**2.
c     &  +fa(19)*(log(cna))**2.*(log(nh3ppt))**3.
c     &  +fa(20)*log(rh)*(log(nh3ppt))**3.
c
c

c      tonset = 143.6002929064716 + 1.0178856665693992*rh +
c     &  10.196398812974294*log(cna) -
c     &  0.1849879416839113*log(cna)**2 - 17.161783213150173*log(nh3ppt)
c     &  + (109.92469248546053*log(nh3ppt))/log(cna) +
c     &  0.7734119613144357*log(cna)*log(nh3ppt)
c     &  - 0.15576469879527022*log(nh3ppt)**2
c
c      if(tonset.gt.temp) then
c
c      fnl=-12.861848898625231 + 4.905527742256349*nh3ppt-
c     &  358.2337705052991*rh -
c     &  0.05463019231872484*nh3ppt*temp + 4.8630382337426985*rh*temp +
c     &  0.00020258394697064567*nh3ppt*temp**2 -
c     &  0.02175548069741675*rh*temp**2 -
c     &  2.502406532869512e-7*nh3ppt*temp**3 +
c     &  0.00003212869941055865*rh*temp**3 -
c     &  4.39129415725234e6/log(cna)**2 +
c     &  (56383.93843154586*temp)/log(cna)**2 -
c     &  (239.835990963361*temp**2)/log(cna)**2 +
c     &  (0.33765136625580167*temp**3)/log(cna)**2 -
c     &  (629.7882041830943*rh)/(nh3ppt**3*log(cna)) +
c     &  (7.772806552631709*rh*temp)/(nh3ppt**3*log(cna)) -
c     &  (0.031974053936299256*rh*temp**2)/(nh3ppt**3*log(cna)) +
c     &  (0.00004383764128775082*rh*temp**3)/(nh3ppt**3*log(cna)) +
c     &  1200.472096232311*log(cna) - 17.37107890065621*temp*log(cna) +
c     &  0.08170681335921742*temp**2*log(cna) -
c     &  0.00012534476159729881*temp**3*log(cna) -
c     &  14.833042158178936*log(cna)**2 +
c     &  0.2932631303555295*temp*log(cna)**2 -
c     &  0.0016497524241142845*temp**2*log(cna)**2 +
c     &  2.844074805239367e-6*temp**3*log(cna)**2 -
c     &  231375.56676032578*log(nh3ppt) -
c     &  100.21645273730675*rh*log(nh3ppt) +
c     &  2919.2852552424706*temp*log(nh3ppt) +
c     &  0.977886555834732*rh*temp*log(nh3ppt) -
c     &  12.286497122264588*temp**2*log(nh3ppt) -
c     &  0.0030511783284506377*rh*temp**2*log(nh3ppt) +
c     &  0.017249301826661612*temp**3*log(nh3ppt) +
c     &  2.967320346100855e-6*rh*temp**3*log(nh3ppt) +
c     &  (2.360931724951942e6*log(nh3ppt))/log(cna) -
c     &  (29752.130254319443*temp*log(nh3ppt))/log(cna) +
c     &  (125.04965118142027*temp**2*log(nh3ppt))/log(cna) -
c     &  (0.1752996881934318*temp**3*log(nh3ppt))/log(cna) +
c     &  5599.912337254629*log(cna)*log(nh3ppt) -
c     &  70.70896612937771*temp*log(cna)*log(nh3ppt) +
c     &  0.2978801613269466*temp**2*log(cna)*log(nh3ppt) -
c     &  0.00041866525019504*temp**3*log(cna)*log(nh3ppt) +
c     &  75061.15281456841*log(nh3ppt)**2 -
c     &  931.8802278173565*temp*log(nh3ppt)**2 +
c     &  3.863266220840964*temp**2*log(nh3ppt)**2 -
c     &  0.005349472062284983*temp**3*log(nh3ppt)**2 -
c     &  (732006.8180571689*log(nh3ppt)**2)/log(cna) +
c     &  (9100.06398573816*temp*log(nh3ppt)**2)/log(cna) -
c     &  (37.771091915932004*temp**2*log(nh3ppt)**2)/log(cna) +
c     &  (0.05235455395566905*temp**3*log(nh3ppt)**2)/log(cna) -
c     &  1911.0303773001353*log(cna)*log(nh3ppt)**2 +
c     &  23.6903969622286*temp*log(cna)*log(nh3ppt)**2 -
c     &  0.09807872005428583*temp**2*log(cna)*log(nh3ppt)**2 +
c     &  0.00013564560238552576*temp**3*log(cna)*log(nh3ppt)**2 -
c     &  3180.5610833308*log(nh3ppt)**3 +
c     &  39.08268568672095*temp*log(nh3ppt)**3 -
c     &  0.16048521066690752*temp**2*log(nh3ppt)**3 +
c     &  0.00022031380023793877*temp**3*log(nh3ppt)**3 +
c     &  (40751.075322248245*log(nh3ppt)**3)/log(cna) -
c     &  (501.66977622013934*temp*log(nh3ppt)**3)/log(cna) +
c     &  (2.063469732254135*temp**2*log(nh3ppt)**3)/log(cna) -
c     &  (0.002836873785758324*temp**3*log(nh3ppt)**3)/log(cna) +
c     &  2.792313345723013*log(cna)**2*log(nh3ppt)**3 -
c     &  0.03422552111802899*temp*log(cna)**2*log(nh3ppt)**3 +
c     &  0.00014019195277521142*temp**2*log(cna)**2*log(nh3ppt)**3 -
c     &  1.9201227328396297e-7*temp**3*log(cna)**2*log(nh3ppt)**3 -
c     &  980.923146020468*log(rh) + 10.054155220444462*temp*log(rh) -
c     &  0.03306644502023841*temp**2*log(rh) +
c     &  0.000034274041225891804*temp**3*log(rh) +
c     &  (16597.75554295064*log(rh))/log(cna) -
c     &  (175.2365504237746*temp*log(rh))/log(cna) +
c     &  (0.6033215603167458*temp**2*log(rh))/log(cna) -
c     &  (0.0006731787599587544*temp**3*log(rh))/log(cna) -
c     &  89.38961120336789*log(nh3ppt)*log(rh) +
c     &  1.153344219304926*temp*log(nh3ppt)*log(rh) -
c     &  0.004954549700267233*temp**2*log(nh3ppt)*log(rh) +
c     &  7.096309866238719e-6*temp**3*log(nh3ppt)*log(rh) +
c     &  3.1712136610383244*log(nh3ppt)**3*log(rh) -
c     &  0.037822330602328806*temp*log(nh3ppt)**3*log(rh) +
c     &  0.0001500555743561457*temp**2*log(nh3ppt)**3*log(rh) -
c     &  1.9828365865570703e-7*temp**3*log(nh3ppt)**3*log(rh)
c
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
     

c      rnuc=0.328886-0.00337417*temp+0.0000183474*temp**2.
c     &     +0.00254198*log(cna)-0.0000949811*temp*log(cna)
c     &     +0.000744627*(log(cna))**2.+0.0243034*log(nh3ppt)
c     &     +0.0000158932*temp*log(nh3ppt)
c     &     -0.00203460*log(cna)*log(nh3ppt)
c     &     -0.000559304*(log(nh3ppt))**2.
c     &     -4.88951d-7*temp*(log(nh3ppt))**2.
c     &     +0.000138470*(log(nh3ppt))**3.+4.14108d-6*fnl
c     &     -0.0000268131*temp*fnl+0.00128791*log(nh3ppt)*fnl
c     &     -3.80352d-6*temp*log(nh3ppt)*fnl-0.0000187902*(fnl)**2.
c


c       rnuc=3.2888553966535506e-10 - 3.374171768439839e-12*temp +
c     &  1.8347359507774313e-14*temp**2+2.5419844298881856e-12*log(cna)-
c     &  9.498107643050827e-14*temp*log(cna) +
c     &  7.446266520834559e-13*log(cna)**2 +
c     &  2.4303397746137294e-11*log(nh3ppt) +
c     &  1.589324325956633e-14*temp*log(nh3ppt) -
c     &  2.034596219775266e-12*log(cna)*log(nh3ppt) -
c     &  5.59303954457172e-13*log(nh3ppt)**2 -
c     &  4.889507104645867e-16*temp*log(nh3ppt)**2 +
c     &  1.3847024107506764e-13*log(nh3ppt)**3 +
c     &  4.141077193427042e-15*log(fn) -
c     &  2.6813110884009767e-14*temp*log(fn) +
c     &  1.2879071621313094e-12*log(nh3ppt)*log(fn) -
c     &  3.80352446061867e-15*temp*log(nh3ppt)*log(fn) -
c     &   1.8790172502456827e-14*log(fn)**2
c
c      else
c  ! Nucleation rate less that 5E-6, setting fnl arbitrary small
cc       fnl=-300.
c       fn = 0.
c       rnuc = 1
c      end if

 10   return
      end
