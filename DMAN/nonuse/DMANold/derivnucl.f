
C     **************************************************
C     * derivnucl                                         *
C     **************************************************

C     RE-WRITTEN BY JaeGun Jung, Novemeber 2007

C     This subroutine calculates Right Hand Side of ODE's.

      SUBROUTINE derivnucl(neq,time,y,yprime,pn,ygas,inucl,deltat)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'
      include 'nucleation.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer inucl ! a nucleation flag
      integer neq 
      real time ! time of a main routine
      real y(neq)
      real yprime(neq) !It is same to dfdx(n)
      real pn(nsect)
      real ygas(ngas)
      real deltat

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i ! a counter
      integer indx1 ! indx1 =first bin of particles
      integer indx2 ! indx2=sulfuric acid
      integer indx3 ! indx3=ammonia gas

      real alpha1 ! a reaction constant for Clement and Ford theory
      real cna ! sulfuric acid concentration as molec cm-3
      real cnasav ! sulfuric acid concentration before nucleation
      real cvt1 ! a conversion factor from kg/m3 to ug/um3
      real cvt2 ! a conversion factor from 1/cm3 to 1/m3
      real cvt3 ! a conversion factor from ppm to ppt
      real density ! density of fresh particles
      real dmassnucl ! mass of all nuclei particles
      real dh2so4 ! sulfuric acid consumed by nucleation
      real dnh3 ! ammonia consumed by nucleation
      real eps  ! a general epslion value
      real fn   ! nucleation rate as particles cm-3 hr-1
      real fnl  ! log10(fn) for JV&M, log(fn) for Napari et al.
      real k ! a rate constant in Spracklen et al. (2006) 
      real Neps ! an epsilon value of pn
      real nh2so4 ! number of sulfuric acid molecules in nuclei
      real nh3pptsav ! ammonia ppt value before nucleation
      real ntot ! total number of molecules in cluster
      real ntotl! log(ntot)
      real oneppt ! molec cm-3 of 1 ppt sulfuric acid at T=298.15K and P=1atm
      real R ! gas constant
      real rnuc ! radius of nuclei
      real rhsav ! rh saved before nucleation
      real tuner ! nucleation tuner
      real vol  ! volume of one nuclei
      real volnucl ! volume of all nuclei 
      real xstar ! mole fraction of sulfuric acid

      ! Ion induced nucleation of Modgil et al.(2005)
      real q ! ion-pair source rate, ion-pair cm-3 s-1
      real sa ! surface area, um3 cm-3
      real h1,h2,h3,h4,h5,h6 ! outputs of the ion induced nucleation

C-----EXTERNAL FUNCTIONS------------------------------------------------

      double precision aerodens_PSSA
      external aerodens_PSSA

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (cvt1 = 1.0e-9) ! a conversion factor from kg/m3 to ug/um3
      parameter (cvt2 = 1.0e+6) ! a conversion factor from 1/cm3 to 1/m3
      parameter (cvt3 = 1.0e+6) ! a conversion factor from ppm to ppt
      parameter (eps = 1.0e-10) ! a general epsilon value
      parameter (Neps = 1.0e-20)! an epsilon value of pn
      parameter (oneppt = 2.463e+7 ) ! molec cm3 corresponding to 1 ppt
                                     ! at T = 298 K and p =  1 atm
      parameter (R = 8.314) ! gas constant, J/mol-K
      parameter (tuner=1.0e+7) ! nucleation tuner

C-----CODE--------------------------------------------------------------

C-----SET ALL THE DERIVATIVES EQUAL TO ZERO

      do i=1, neq
         if(y(i).lt.0.)y(i)=0
         yprime(i)=0.0
      enddo

      ! Setting indx values
      indx1=1  ! first bin of particle distributions
      indx2=2  ! sulfuric acid
      indx3=3  ! ammonia gas.

      ! Initialization
      vol = 0.0

C-----CALCULATE THE RATE OF CHANGE OF H2SO4(g) DUE TO LINEARIZED GAS-PHASE
C     CHEMISTRY (OPERATOR SPLITTING SCHEME)

      yprime(indx2)=yprime(indx2)+dsulfdt

C-----CALCULATE THE RATE OF CHANGE OF NH3(g) FROM EMISSION

      if (nh3flag.eq.1) then
        yprime(indx3)=yprime(indx3)+nh3flux
      endif

      cna = y(indx2)*oneppt*(temp0/temp)*(pres/pres0) 
                               ! H2SO4 concentration in molec/cm3
      if ((iin_flag.eq.1).or.(CF_flag.eq.1)) goto 10 
                               ! If the theories, goto 10.

C-----CALCULATE NULEATION RATE FOR VARIOUS NUCLEATION THEORIES

C-----Binary nucleation theory

C-----Jaeker-Voirol and Mirabel (1989) Atmos. Environ. 23(9):2054-2057

c      if (cna. le. 0) goto 30
c      fnl=-(64.24+4.7*rh)+(6.13+1.95*rh)*
c     &      (LOG10(cna)+(298.-temp)/25.)
c      fn = tuner*(10.**fnl)*3600. ! nucleation rate in p/(cm3 hr)
csensitivity      fn = 1.0e10*(10.**fnl)*3600. ! A tuner value is changed to 1.e10
c      goto 20

C-----Spracklen et al. (2006)  Atmos. Chem. Phys. 6:5631-5648

c       k = 2.0e-6
csensitivity      k = 0.5e-6
csensitivity_not_used      k = 0.5e-7
csensitivity      k = 0.1e-7 
c       fn = k * cna * 3600 ! p/(cm3 hr)
c       goto 20

C-----IF [NH3] is below threshold of ternary nucleation, 
C     then do the binary nucleation of Vehkamaki et al.(2002).

      IF (y(indx3) .lt. 0.1) THEN ! Comment out when do binary of Vehkamaki.
cdbg      IF (ygas(mgnh3) .lt. 0.1) THEN ! Comment out when do binary of Vehkamaki.

C-----Vehkamaki et al. (2002) J. Geophys. Res. 107(D22):4622-4632

      !Save original values
      cnasav=cna
      rhsav =rh

      !Valid Boundary of each parameters
      if (cna .lt. 1.e-4) goto 30
      cna = min(1.0e+11,cna)
      if (temp .gt. 305.15) then
         temp=305.15
      elseif (temp .lt. 230.15) then
         temp=230.15
      endif
      if (rh .lt. 1e-4) rh=1e-4

      !Mole fraction of sulfuric acid
      xstar=0.740997-0.00266379*temp-0.00349998*log(cna)
     &   +0.0000504022*temp*log(cna)+0.00201048*log(rh)
     &   -0.000183289*temp*log(rh)+0.00157407*(log(rh))**2.
     &   -0.0000179059*temp*(log(rh))**2.
     &   +0.000184403*(log(rh))**3.
     &   -1.50345E-6*temp*(log(rh))**3.

      !Nucleation rate coefficients 
      do i=1, 10
        fb(i) = fb0(i)+fb1(i)*temp+fb2(i)*temp**2.
     &       +fb3(i)*temp**3.+fb4(i)/xstar
      enddo

      !Log nucleation rate
      fnl = fb(1)+fb(2)*log(rh)+fb(3)*(log(rh))**2.
     &    +fb(4)*(log(rh))**3.+fb(5)*log(cna)
     &    +fb(6)*log(rh)*log(cna)+fb(7)*(log(rh))**2.*log(cna)
     &    +fb(8)*(log(cna))**2.+fb(9)*log(rh)*(log(cna))**2.
     &    +fb(10)*(log(cna))**3.

      !Nucleation rate (1/cm3-s)
      fn = exp(fnl)

      !Cap at 10^6 particles/cm3-s, limit for parameterization
      if (fn.gt.1.0d6) then
        fn=1.0d6
      endif

      !Coefficients of total number of molecules in cluster 
      do i=1, 10
        gb(i) = gb0(i)+gb1(i)*temp+gb2(i)*temp**2.
     &       +gb3(i)*temp**3.+gb4(i)/xstar
      enddo

      !log total number of molecules in cluster
      ntotl=gb(1)+gb(2)*log(rh)+gb(3)*(log(rh))**2.
     &    +gb(4)*(log(rh))**3.+gb(5)*log(cna)
     &    +gb(6)*log(rh)*log(cna)+gb(7)*log(rh)**2.*log(cna)
     &    +gb(8)*(log(cna))**2.+gb(9)*log(rh)*(log(cna))**2.
     &    +gb(10)*(log(cna))**3.

      !Total number of molecules in cluster
      ntot=exp(ntotl)

      !If H2SO4 concentration is small, then nucleation rate equal to zero.
      !Otherwise, set nuclei radius based on parameterization.
      if (cna.lt.1.0d4) then
        fn=0.0d0
        rnuc=0.
      else
        fn=3600.*fn  ! [1/cm3-hr]
        rnuc=exp(-1.6524245+0.42316402*xstar+0.3346648*log(ntot)) ! [nm]
      endif
csensitivity       fn=1.0e10*fn ! tuner = 1.0e10
csensitivity       fn=1.0e15*fn ! tuner = 1.0e15

      !idnuc=0 set nuclei to size of smallest section
      !idnuc=1 allows size to change
      !currently idnuc=0 kept as 0 due to numerical issues

      if (idnuc .eq. 0) then
        yprime(indx1)=yprime(indx1)+fn*FLOAT(inucl)
cdbg      else
cdbg        i=1
cdbg        do while(dpmean(i) .le. 2.*rnuc*1.0d-3)
cdbg          i=i+1
cdbg        enddo
cdbg        yprime(i)=yprime(i) + fn*FLOAT(inucl)
      endif

      !In binary nucleation case, there is a lot of
      !water molecule inside of nuclei. So, the dh2so4 is
      !calculated from nh2so4 rather than a total mass nucleated
      !such as a way of ternary nucleation.

      nh2so4=xstar*ntot !total sulfuric acid molecules in nuclei cluster
      dh2so4= fn*nh2so4*(1.0/oneppt)*(temp/temp0)*(pres0/pres) ![ppt/hr]
      yprime(indx2)=yprime(indx2)-dh2so4*FLOAT(inucl)

      !Recover original values      
      cna=cnasav
      rh =rhsav
      ELSE ! Comment out when do binary of Vehkamaki.
cternary      goto 30 ! COMMENT OUT WHEN TERNARY NUCLEATION IS IMPLEMENTED!!

C-----Ternary nucleation theory
C     Napari et al (2002)  J. Geophys. Res. 107(D19):4381-4386

      !Save original values of parameters

      nh3pptsav=ygas(mgnh3) ! save original NH3 mixing ratio.
      cnasav=cna
      rhsav =rh

      !Napari's parameterization is only valid within limited area. 
      !So variables changes to fix maximum, minimum H2SO4 concentration
      if (cna.ge.1.0d9) cna=1.0d9 !limit sulfuric acid conc
      if (rh.lt.0.05) rh = 5.0d-2 !limit rh
      if (rh.gt.0.95) rh = 9.5d-1 !limit rh
      ygas(mgnh3)=min(1.0e+2,ygas(mgnh3)) ! limit NH3 mixing ratio

      !Calculates parameters
      do i=1,20
       fa(i)=aa0(i)+a1(i)*temp+a2(i)*temp**2.+a3(i)*temp**3.
      enddo

      !Derive log nucleation rate
      fnl=-84.7551+fa(1)/log(cna)+fa(2)*log(cna)+fa(3)*(log(cna))**2.
     &  +fa(4)*log(y(indx3))+fa(5)*(log(y(indx3)))**2.+fa(6)*rh
     &  +fa(7)*log(rh)+fa(8)*log(y(indx3))/log(cna)+fa(9)
     &  *log(y(indx3))*log(cna)+fa(10)*rh*log(cna)+fa(11)*rh/log(cna)
     &  +fa(12)*rh*log(y(indx3))+fa(13)*log(rh)/log(cna)+fa(14)
     &  *log(rh)*log(y(indx3))+fa(15)*(log(y(indx3)))**2./log(cna)
     &  +fa(16)*log(cna)*(log(y(indx3)))**2.+fa(17)*(log(cna))**2.
     &  *log(y(indx3))+fa(18)*rh*(log(y(indx3)))**2.+fa(19)*rh
     &  *log(y(indx3))/log(cna)+fa(20)*(log(cna))**2.
     &  *(log(y(indx3)))**2.

      !Get rid of log value
      fn=exp(fnl)
 
      !Cap at 10^6 particles/cm3-s, limit for parameterization
      if (fn.gt.1.0d6) then
        fn=3600.*1.0d6
        fnl=log(fn/3600.)
      else
        fn=3600.*fn
      endif

ctmg  If H2SO4 concentration is small set nucleation rate equal to zero
ctmg  Otherwise, set nuclei radius based on parameterization.
      if (cna.lt.1.0d4) then
        fn=0.0d0
        fnl=0.0d0
      else 
        rnuc=0.141027-0.00122625*fnl-7.82211d-6*fnl**2.
     &     -0.00156727*temp-0.00003076*temp*fnl
     &     +0.0000108375*temp**2.
      endif

      !idnuc=0 sets nuclei to size of smallest section
      !idnuc=1 allows size to change
      !currently idnuc=0 kept as 0 due to numerical issues
      if(idnuc .eq. 0) then
        yprime(indx1)=yprime(indx1)+fn*FLOAT(inucl)
cdbg       else
cdbg         i=1
cdbg         do while(dpmean(i).le.2.*rnuc*1.0d-3)
cdbg           i=i+1
cdbg         enddo
cdbg         yprime(i)=yprime(i)+fn*float(inucl)
      endif

      !number of sulfuric acid molecules in a nucleus
cdbg      nh2so4=38.1645+0.774106*fnl+0.00298879*fnl**2.-0.357605
cdbg     &     *temp-0.00366358*temp*fnl+0.0008553*temp**2.
      ! nh2so4 is not used. Instead volume of nuclei is used to 
      ! calculate loss of H2SO4. It is because the uncertainty
      ! of parameterization increases more when you use these 
      ! twice than once.

      if (idnuc .eq. 0) then
        vol = (1./6.)*pi*dpmean(1)**3           ! in um3 / particle
      else                 
        vol = (1./6.)*pi*dpmean(i)**3           ! in um3 / particle
      endif
      volnucl = fn*FLOAT(inucl)*vol     ! um3 / cm3-hr 
cdbg     density = aerodens_PSSA(Mk(srtso4),0.0,Mk(srtnh3),0.0,Mk(srth2o))
      density = 1400 ! [=] kg/m3 
      dmassnucl = density * volnucl * cvt1 * cvt2 ! ug / m3-hr

      !nuclei are assumed as ammonium bisulfate
      dh2so4 = 0.8144*dmassnucl*(R*temp)/(pres*97.0)*cvt3 ! in ppt/hr
      dnh3 = 0.1856*dmassnucl*(R*temp)/(pres*18.0)*cvt3 ! in ppt/hr
      yprime(indx2) = yprime(indx2) - dh2so4*float(inucl)
cdbg      if (time.ge.6.338) then
cdbg        print*,'time=',time, ' deltat=',deltat
cdbg        print*,'y=',y
cdbg        print*,'yprime=',yprime
cdbg        pause
cdbg      endif
      if (nh3flag.eq.1) then
         yprime(indx3) = yprime(indx3) - dnh3*float(inucl) !NH3 mass balance
      endif

      !Recover original values
      ygas(mgnh3)=nh3pptsav
      cna=cnasav
      rh=rhsav
      goto 30

      ENDIF ! Comment out when do binary of Vehkamaki.
c
 10   continue

C-----OTHER NUCLEATION THEORIES

C-----Ion-induced of Modgil et al. (2005) J. Geophys. Res. 110:D19205

      !If the ion induced nucleation is requested, then excute

      if (iin_flag .eq. 1) then
        q = 2! ion pairs cm-3 s-1 at ground level from Yu and Turco (2001)
             ! The range of parameterization of q = 2 - 50 cm-3 s-1
csensitivity        q = 50 ! for the maximum value
cnon_used        q = 1.0e3 ! for the maximum value
        sa = 0
        sa = sa+y(1)*3.14*dpmean(i)**2 ! because y(1) is kept changing during RKQC.
        do i=2,nsect
          if (pn(i) .ge. Neps) then
            sa=sa+pn(i)*3.14*dpmean(i)**2 !um2/cm3
          endif
        enddo

        !The boundaries of parameterization
        if (cna .lt. 1e5) goto 30
        if (cna .gt. 1e8) cna = 1e8
        if (sa .lt. 2) sa = 2
        if (sa .gt. 100) sa = 100
        if (temp .lt. 190) temp = 190
        if (temp .gt. 300) temp = 300
        if (rh .lt. 0.05) rh = 0.05
        if (rh .gt. 0.95) rh = 0.95

        call ion_nucl(cna,sa,temp,q,rh,h1,h2,h3,h4,h5,h6) 
        if (h1.le.1.0e-6) goto 30 ! Auxiliary material in Modgil et al.(2005)
        if (h1.ge.eps) then
          if (h1.gt.1.0e6) then
              h1 = 1.0e6 ! A maximum boundary is not suggested in Modgil et al.
                         ! (2005). Although the maximum value shown in the 
                         ! paper is less than 100, the value of 1.0e6 cm-3s-1 
                         ! is sellect following a cap value in Vehkamaki et al.
                         ! (2002) and Napari et al. (2002).
          endif
          fn = h1 * 3600.0 ! particle cm-3 hr-1
        else
          fn = 0
        endif
        yprime(1) = yprime(1) + fn*FLOAT(inucl)
        if (h2.ge.eps) then
          dh2so4=h2*(1.0/oneppt)*(temp/temp0)*(pres0/pres)*3600.0 ![ppt/hr]
cdbg          dh2so4=h1*h3*(1/oneppt)*(temp/temp0)*(pres0/pres)*3600.0 ![ppt/hr]
        else
          dh2so4 = 0
        endif
        yprime(indx2) = yprime(indx2) - dh2so4*float(inucl)

        goto 30
      endif

C-----Clement and Ford (1999) Atmos. Environ. 33:489-499

      if (CF_flag.eq.1) then
        if (ygas(mgnh3) .lt. 0.1) then
          alpha1=4.276e-10*sqrt(temp/293.15) ! For sulfuric acid
        else
          alpha1=3.684e-10*sqrt(temp/293.15) ! For ammonium sulfate
        endif
        fn = alpha1 * cna**2. *3600.
csensitivity       fn = 1.e-3 * fn ! 10^-3 tuner
        if (fn.gt.1.0e9) fn=1.0e9 ! For numerical conversion
      endif

      !Calculate derivatives 

 20   yprime(indx1) = yprime(indx1) + fn*FLOAT(inucl)
      vol = (1./6.) *pi*dpmean(1)**3.
      volnucl = fn * FLOAT(inucl) * vol
cdbg      density = aerodens_PSSA(Mk(srtso4),0.0,Mk(srtnh3),0.0,Mk(srth2o))
      density = 1400 ! [=] kg/m3
      dmassnucl = density * 0.001 * volnucl
      dh2so4 = dmassnucl*(R*temp)/(pres*97.0)*cvt3 ! in ppt/hr
      yprime(indx2)=yprime(indx2)-dh2so4*float(inucl)


 30   continue
      RETURN
      END
