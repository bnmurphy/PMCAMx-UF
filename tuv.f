      SUBROUTINE tuv(nz, z, airlev, alb, coszen, 
     $     odcld, omcld, gcld, 
     $     odaer1, odaer2, omaer1, omaer2, gaer,
     $     rafcld, actflx)

      IMPLICIT NONE

**** Inputs:

* altitude grid

      INTEGER nz
      REAL z(nz)

* solar zenith angle 

      REAL coszen

* air profile

      REAL airlev(nz)

* cloud and aerosol 

      REAL odcld(nz), omcld(nz), gcld(nz)
      REAL odaer1(nz), omaer1(nz), gaer(nz)
      REAL odaer2(nz), omaer2(nz)

* surface albedo

      REAL alb

**** Output: Actinic flux cloudy/clear ratio as a funct of altitude

      REAL rafcld(nz)

**** Internal:

* internally computed optical depths (to avoid overwriting odcld, odaer)

      REAL dtrl(nz), dtcld(nz), dtaer(nz)

* Rayleigh scattering cross section from WMO 1985 (originally from
* Nicolet, M., On the molecular scattering in the terrestrial atmosphere:
* An empirical formula for its calculation in the homoshpere, Planet.
* Space Sci., 32, 1467-1468, 1984.
* computed at the selected wavelenth wc(nm) using the equations:
c         wmicrn =  wc/1.E3
c         IF( wmicrn .LE. 0.55) THEN
c            xx = 3.6772 + 0.389*wmicrn + 0.09426/wmicrn
c         ELSE
c            xx = 4. + 0.04
c         ENDIF
c         srayl(iw) = 4.02e-28/(wmicrn)**xx

      REAL srayl
      DATA srayl /2.922E-26/

* scale heights for air density

      REAL hscale
      DATA hscale /8.05e5/

* actinic flux (relative units) 

      REAL fdir(nz), fdn(nz), fup(nz)
      REAL actflx(2,nz)

* loop indices

      INTEGER i
* ___

* compute Rayleigh optical depths:
* include exponential tail integral at top (iz = nz)

      DO i = 1, nz-1
         dtrl(i) = (airlev(i+1)+airlev(i)) * 0.5 *
     $        (z(i+1)-z(i)) * 1.E5 * srayl
      ENDDO
      dtrl(nz) = hscale * airlev(nz) * 1.E5

**** clear sky calculation
** monochromatic radiative transfer:
* call rtlink for plane-parallel 2-stream
      
      DO i = 1, nz-1
         dtcld(i) = 0.
         dtaer(i) = odaer1(i)
      ENDDO
      CALL rtlpp2(nz, alb, coszen,
     $     dtrl,
     $     dtcld, omcld, gcld,
     $     dtaer,omaer1,gaer,
     $     fdir, fdn, fup)
      DO i = 1, nz
         actflx(1,i) = fdir(i) + fdn(i) + fup(i)
      ENDDO

**** cloudy sky calculation:

      DO i = 1, nz-1
         dtcld(i) = odcld(i)
         dtaer(i) = odaer2(i)
      ENDDO
      CALL rtlpp2(nz, alb, coszen,
     $     dtrl,
     $     dtcld, omcld, gcld,
     $     dtaer,omaer2,gaer,
     $     fdir, fdn, fup)
      DO i = 1, nz
         actflx(2,i) = fdir(i) + fdn(i) + fup(i)
      ENDDO

      DO i = 1, nz
         rafcld(i) = actflx(2,i)/actflx(1,i)
c      print *, rafcld(i) !added by Elham
      ENDDO
      
      
      RETURN
      END

*==============================================================================

      SUBROUTINE rtlpp2(nz, ag, m0, 
     $     dtrl, 
     $     dtcld, omcld, gcld,
     $     dtaer,omaer,gaer,
     $     fdir, fdn, fup)
*_______________________________________________________________________

      IMPLICIT NONE

* input:

      INTEGER nz
      REAL ag
      REAL m0
      REAL dtrl(nz)
      REAL dtcld(nz), omcld(nz), gcld(nz)
      REAL dtaer(nz), omaer(nz), gaer(nz)

* output

      REAL fdir(nz), fdn(nz), fup(nz)

* local:
      REAL dt(nz), om(nz), g(nz)
      REAL fdiri(nz), fdni(nz), fupi(nz)
      REAL daaer, dtsct, dtabs, dsaer, dscld, dacld
      INTEGER i, ii

      LOGICAL delta
      DATA delta /.true./

* largest number of the machine:
      REAL largest
      PARAMETER(largest=1.E+36)

*_______________________________________________________________________

* initialize:

      DO 5 i = 1, nz
         fdir(i) = 0.
         fup(i) = 0.
         fdn(i) = 0.
 5    CONTINUE

* ----

*  set here any coefficients specific to rt scheme, 
* ----

      DO 10, i = 1, nz - 1

         dscld = dtcld(i)*omcld(i)
         dacld = dtcld(i)*(1.-omcld(i))

         dsaer = dtaer(i)*omaer(i)
         daaer = dtaer(i)*(1.-omaer(i))

         dtsct = dtrl(i) + dscld + dsaer
         dtabs = dacld + daaer

         dtabs = AMAX1(dtabs,1./largest)
         dtsct = AMAX1(dtsct,1./largest)

* invert z-coordinate:

         ii = nz - i
         dt(ii) = dtsct + dtabs
         om(ii) = dtsct/(dtsct + dtabs)
           IF(dtsct .EQ. 1./largest) om(ii) = 1./largest
         g(ii) = (gcld(i)*dscld + gaer(i)*dsaer)/dtsct

   10 CONTINUE

*  call rt routine:

      CALL pp2str(nz,m0,ag,dt,om,g,
     $        fdiri, fupi, fdni)

* put on upright z-coordinate

      DO 20, i = 1, nz
         ii = nz - i + 1
         fdir(i) = fdiri(ii)
         fup(i) = fupi(ii)
         fdn(i) = fdni(ii)
 20   CONTINUE
*_______________________________________________________________________

      RETURN
      END   

*==============================================================================

      SUBROUTINE pp2str(nz,mu,rsfc,tauu,omu,gu,
     $     fdr, fup, fdn)
*_______________________________________________________________________
* TWO-STREAM EQUATIONS FOR MULTIPLE LAYERS
* based on equations from Toon et al.,  Journal of Geophysical Research
* Volume 94, #D13  Nov. 20, 1989 Issue 
* programmed by:  Kathleen G. Mosher
* for Sasha Madronich, N.C.A.R. A.C.D.
*_______________________________________________________________________
* Now it contains 9 two-stream methods to choose from.
* programmed  on 05.26.94 by:  Irina V.Petropavlovskikh
* for Sasha Madronich, N.C.A.R. A.C.D.
*_______________________________________________________________________

      IMPLICIT NONE

*******
* input:
*******
      INTEGER nz
      REAL mu, rsfc
      REAL tauu(nz), omu(nz), gu(nz)

*******
* output:
*******
      REAL fup(nz),fdn(nz),fdr(nz)
      REAL eup(nz),edn(nz),edr(nz)

*******
* local:
*******

* internal coefficients and matrix
      REAL lam(nz),taun(nz),bgam(nz)
      REAL e1(nz),e2(nz),e3(nz),e4(nz)
      REAL cup(nz),cdn(nz),cuptn(nz),cdntn(nz)
      REAL tauc, mu1(nz)
      INTEGER row
      REAL a(2*nz),b(2*nz),d(2*nz),e(2*nz),y(2*nz)

*******
* other:
*******
      REAL pifs, fdn0
      REAL f, g, tau, om, tempg
      REAL gam1, gam2, gam3, gam4

* For calculations of Associated Legendre Polynomials for GAMA1,2,3,4
* in delta-function, modified quadrature, hemispheric constant,
* Hybrid modified Eddington-delta function metods, p633,Table1.
* W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
* W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440, 
* uncomment the following two lines and the appropriate statements further
* down.
C     REAL YLM0, YLM2, YLM4, YLM6, YLM8, YLM10, YLM12, YLMS, BETA0,
C    >     BETA1, BETAn, amu1, subd

      REAL expon, expon0, expon1, divisr, temp, up, dn
      REAL ssfc
      REAL taug
      INTEGER nlayer, mrows, lev

      INTEGER i, j

* Some additional program constants:
      REAL eps, precis
      PARAMETER (eps = 1.E-3, precis = 1.E-7) 

*_______________________________________________________________________

* MU = cosine of solar zenith angle
* RSFC = surface albedo
* TAUU =  unscaled optical depth of each layer
* OMU  =  unscaled single scattering albedo
* GU   =  unscaled asymmetry factor
* KLEV = max dimension of number of layers in atmosphere
* NLAYER = number of layers in the atmosphere
* NLEVEL = nlayer + 1 = number of levels

* initial conditions:  pi*solar flux = 1;  diffuse incidence = 0

      pifs = 1.      
      fdn0 = 0.

      nlayer = nz - 1

************** compute coefficients for each layer:
* GAM1 - GAM4 = 2-stream coefficients, different for different approximations
* EXPON0 = calculation of e when TAU is zero
* EXPON1 = calculation of e when TAU is TAUN
* CUP and CDN = calculation when TAU is zero
* CUPTN and CDNTN = calc. when TAU is TAUN
* DIVISR = prevents division by zero

      tauc = 0.
      DO 10, i = 1, nlayer

         g = gu(i)
         tau = tauu(i)
         om = omu(i)

* stay away from 1 by precision.  For g, also stay away from -1

         tempg = AMIN1(abs(g),1. - precis)
         g = SIGN(tempg,g)
         om = AMIN1(om,1-precis)

* delta-scaling. Have to be done for delta-Eddington approximation, 
* delta discrete ordinate, Practical Improved Flux Method, delta function,
* and Hybrid modified Eddington-delta function methods approximations

         f = g*g
         g = (g - f)/(1 - f)
         taun(i) = (1 - om*f)*tau
         om = (1 - f)*om/(1 - om*f)       

*** the following gamma equations are from pg 16,289, Table 1

* Eddington approximation(Joseph et al., 1976, JAS, 33, 2452):

         gam1 = (7. - om*(4. + 3.*g))/4.
         gam2 = -(1. - om*(4. - 3.*g))/4.
         gam3 = (2. - 3.*g*mu)/4.
         gam4 = 1. - gam3
         mu1(i) = 0.5

* quadrature (Liou, 1973, JAS, 30, 1303-1326; 1974, JAS, 31, 1473-1475):

c         gam1 = 1.7320508*(2. - om*(1. + g))/.2
c         gam2 = 1.7320508*om*(1. - g)/2.
c         gam3 = (1. - 1.7320508*g*mu)/2
c         gam4 = 1. - gam3
         
* hemispheric mean (Toon et al., 1089, JGR, 94, 16287):

c          gam1 = 2. - om*(1. + g)
c          gam2 = om*(1. - g)
c          gam3 = (2. - g*mu)/4.
c          gam4 = 1. - gam3

* PIFM  (Zdunkovski et al.,1980, Conrib.Atmos.Phys., 53, 147-166):
c         GAM1 = 0.25*(8. - OM*(5. + 3.*G))
c         GAM2 = 0.75*OM*(1.-G)
c         GAM3 = 0.25*(2.-3.*G*MU)
c         GAM4 = 1. - GAM3

* delta discrete ordinates  (Schaller, 1979, Contrib.Atmos.Phys, 52, 17-26):
c         GAM1 = 0.5*1.7320508*(2. - OM*(1. + G))
c         GAM2 = 0.5*1.7320508*OM*(1.-G)
c         GAM3 = 0.5*(1.-1.7320508*G*MU)
c         GAM4 = 1. - GAM3

* Calculations of Associated Legendre Polynomials for GAMA1,2,3,4
* in delta-function, modified quadrature, hemispheric constant,
* Hybrid modified Eddington-delta function metods, p633,Table1.
* W.E.Meador and W.R.Weaver, GAS,1980,v37,p.630
* W.J.Wiscombe and G.W. Grams, GAS,1976,v33,p2440
c      YLM0 = 2.
c      YLM2 = -3.*G*MU
c      YLM4 = 0.875*G**3*MU*(5.*MU**2-3.)
c      YLM6=-0.171875*G**5*MU*(15.-70.*MU**2+63.*MU**4)
c     YLM8=+0.073242*G**7*MU*(-35.+315.*MU**2-693.*MU**4
c    *+429.*MU**6)
c     YLM10=-0.008118*G**9*MU*(315.-4620.*MU**2+18018.*MU**4
c    *-25740.*MU**6+12155.*MU**8)
c     YLM12=0.003685*G**11*MU*(-693.+15015.*MU**2-90090.*MU**4
c    *+218790.*MU**6-230945.*MU**8+88179.*MU**10)
c      YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
c      YLMS=0.25*YLMS
c      BETA0 = YLMS
c
c         amu1=1./1.7320508
c      YLM0 = 2.
c      YLM2 = -3.*G*amu1
c      YLM4 = 0.875*G**3*amu1*(5.*amu1**2-3.)
c      YLM6=-0.171875*G**5*amu1*(15.-70.*amu1**2+63.*amu1**4)
c     YLM8=+0.073242*G**7*amu1*(-35.+315.*amu1**2-693.*amu1**4
c    *+429.*amu1**6)
c     YLM10=-0.008118*G**9*amu1*(315.-4620.*amu1**2+18018.*amu1**4
c    *-25740.*amu1**6+12155.*amu1**8)
c     YLM12=0.003685*G**11*amu1*(-693.+15015.*amu1**2-90090.*amu1**4
c    *+218790.*amu1**6-230945.*amu1**8+88179.*amu1**10)
c      YLMS=YLM0+YLM2+YLM4+YLM6+YLM8+YLM10+YLM12
c      YLMS=0.25*YLMS
c      BETA1 = YLMS
c
c         BETAn = 0.25*(2. - 1.5*G-0.21875*G**3-0.085938*G**5
c    *-0.045776*G**7)


* Hybrid modified Eddington-delta function(Meador and Weaver,1980,JAS,37,630):
c         subd=4.*(1.-G*G*(1.-MU))
c         GAM1 = (7.-3.*G*G-OM*(4.+3.*G)+OM*G*G*(4.*BETA0+3.*G))/subd
c         GAM2 =-(1.-G*G-OM*(4.-3.*G)-OM*G*G*(4.*BETA0+3.*G-4.))/subd
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
*****
* delta function  (Meador, and Weaver, 1980, JAS, 37, 630):
c         GAM1 = (1. - OM*(1. - beta0))/MU
c         GAM2 = OM*BETA0/MU
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3
*****
* modified quadrature (Meador, and Weaver, 1980, JAS, 37, 630):
c         GAM1 = 1.7320508*(1. - OM*(1. - beta1))
c         GAM2 = 1.7320508*OM*beta1
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3

* hemispheric constant (Toon et al., 1989, JGR, 94, 16287):
c         GAM1 = 2.*(1. - OM*(1. - betan))
c         GAM2 = 2.*OM*BETAn
c         GAM3 = BETA0
c         GAM4 = 1. - GAM3

*****

*  save mu1 for use in converting irradiance to actinic flux

c         mu1(i) = (1-om)/(gam1 - gam2)

* lambda = pg 16,290 equation 21
* big gamma = pg 16,290 equation 22
 
         lam(i) = sqrt(gam1*gam1 - gam2*gam2)

         IF( gam2 .NE. 0.) THEN
            bgam(i) = (gam1 - lam(i))/gam2
         ELSE
            bgam(i) = 0.
         ENDIF

         expon = EXP(-lam(i)*taun(i))

* e1 - e4 = pg 16,292 equation 44
         
         e1(i) = 1. + bgam(i)*expon
         e2(i) = 1. - bgam(i)*expon
         e3(i) = bgam(i) + expon
         e4(i) = bgam(i) - expon

* the following sets up for the C equations 23, and 24
* found on page 16,290
* prevent division by zero (if LAMBDA=1/MU, shift 1/MU^2 by EPS = 1.E-3
* which is approx equiv to shifting MU by 0.5*EPS* (MU)**3

         expon0 = EXP(-(tauc          )/mu)
         expon1 = EXP(-(tauc + taun(i))/mu)

         divisr = lam(i)*lam(i) - 1./(mu*mu)
         temp = AMAX1(eps,abs(divisr))
         divisr = SIGN(temp,divisr)

         up = om*pifs*((gam1 - 1./mu)*gam3 + gam4*gam2)/divisr
         dn = om*pifs*((gam1 + 1./mu)*gam4 + gam2*gam3)/divisr
         
* cup and cdn are when tau is equal to zero
* cuptn and cdntn are when tau is equal to taun

         cup(i) = up*expon0
         cdn(i) = dn*expon0
         cuptn(i) = up*expon1
         cdntn(i) = dn*expon1

         tauc = tauc + taun(i)

   10 CONTINUE

***************** set up matrix ******
* ssfc = pg 16,292 equation 37  where pi Fs is one (unity).

      ssfc = rsfc*mu*EXP(-tauc/mu)*pifs

* MROWS = the number of rows in the matrix

      mrows = 2*nlayer     
      
* the following are from pg 16,292  equations 39 - 43.
* set up first row of matrix:

      i = 1
      a(1) = 0.
      b(1) = e1(i)
      d(1) = -e2(i)
      e(1) = fdn0 - cdn(i)

      row=1

* set up odd rows 3 thru (MROWS - 1):

      i = 0
      DO 20, row = 3, mrows - 1, 2
         i = i + 1
         a(row) = e2(i)*e3(i) - e4(i)*e1(i)
         b(row) = e1(i)*e1(i + 1) - e3(i)*e3(i + 1)
         d(row) = e3(i)*e4(i + 1) - e1(i)*e2(i + 1)
         e(row) = e3(i)*(cup(i + 1) - cuptn(i)) + 
     $        e1(i)*(cdntn(i) - cdn(i + 1))
   20 CONTINUE

* set up even rows 2 thru (MROWS - 2): 

      i = 0
      DO 30, row = 2, mrows - 2, 2
         i = i + 1
         a(row) = e2(i + 1)*e1(i) - e3(i)*e4(i + 1)
         b(row) = e2(i)*e2(i + 1) - e4(i)*e4(i + 1)
         d(row) = e1(i + 1)*e4(i + 1) - e2(i + 1)*e3(i + 1)
         e(row) = (cup(i + 1) - cuptn(i))*e2(i + 1) - 
     $        (cdn(i + 1) - cdntn(i))*e4(i + 1)
   30 CONTINUE

* set up last row of matrix at MROWS:

      row = mrows
      i = nlayer
      
      a(row) = e1(i) - rsfc*e3(i)
      b(row) = e2(i) - rsfc*e4(i)
      d(row) = 0.
      e(row) = ssfc - cuptn(i) + rsfc*cdntn(i)

* solve tri-diagonal matrix:

      CALL tridag(a, b, d, e, y, mrows)

**** unfold solution of matrix, compute output fluxes:

      row = 1 
      lev = 1
      j = 1
      taug = 0.
      
* the following equations are from pg 16,291  equations 31 & 32

      fdr(lev) = 1.
      edr(lev) = mu
      edn(lev) = fdn0
      eup(lev) =  y(row)*e3(j) - y(row + 1)*e4(j) + cup(j)
      fdn(lev) = edn(lev)/mu1(lev)
      fup(lev) = eup(lev)/mu1(lev)

      DO 60, lev = 2, nlayer + 1
         taug = taug + taun(j)
         fdr(lev) = EXP(-taug/mu)
         edr(lev) = mu*fdr(lev)
         edn(lev) =  y(row)*e3(j) + y(row + 1)*e4(j) + cdntn(j)
         eup(lev) =  y(row)*e1(j) + y(row + 1)*e2(j) + cuptn(j)
         fdn(lev) = edn(lev)/mu1(j)
         fup(lev) = eup(lev)/mu1(j)

         row = row + 2
         j = j + 1
   60 CONTINUE
*_______________________________________________________________________

      RETURN
      END

*==============================================================================

      SUBROUTINE tridag(a,b,c,r,u,n)
*_______________________________________________________________________
* solves tridiagonal system.  From Numerical Recipies, p. 40
*_______________________________________________________________________

      IMPLICIT NONE

* input:
      INTEGER n
      REAL a, b, c, r
      DIMENSION a(n),b(n),c(n),r(n)

* output:
      REAL u
      DIMENSION u(n)

* local:
      INTEGER j
      REAL bet, gam
      DIMENSION gam(n)
*_______________________________________________________________________

      IF (b(1) .EQ. 0.) STOP 1001
      bet   = b(1)
      u(1) = r(1)/bet
      DO 11, j = 2, n   
         gam(j) = c(j - 1)/bet
         bet = b(j) - a(j)*gam(j)
         IF (bet .EQ. 0.) STOP 2002 
         u(j) = (r(j) - a(j)*u(j - 1))/bet
   11 CONTINUE
      DO 12, j = n - 1, 1, -1  
         u(j) = u(j) - gam(j + 1)*u(j + 1)
   12 CONTINUE
*_______________________________________________________________________

      RETURN
      END
