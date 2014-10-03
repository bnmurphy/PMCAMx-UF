
C     **************************************************
C     * alpha_nh3                                      *
C     **************************************************

C     WRITTEN BY JaeGun Jung, November 2007

C     This function calculates accomodation coefficient of ammonia from
C     ammonium and sulfate mass ratio.

C-----Literature cited--------------------------------------------------
C     Pathak et al., Characteristics of aerosol acidity in Hong Kong,
C       Atmos. Environ., 38, 2965-2974,2004 
C     Swartz et al., Uptake of Gas-Phase Ammonia. 2. Uptake by Sulfuric
C       Acid Surfaces, J. Phys. Chem. A, 103, 8824-8833

C-----INPUTS------------------------------------------------------------

C     mnh4 - mass of ammonium [=] kg
C     mso4 - mass of sulfate  [=] kg
C     rh   - relative humidity

C-----OUTPUTS-----------------------------------------------------------

      real FUNCTION alpha_nh3(mnh4,mso4,rh)

      IMPLICIT NONE

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real mnh4
      real mso4
      real rh   !relative humidity

C-----VARIABLE DECLARATIONS---------------------------------------------

      real AS   !ammonium molar ratio to sulfate
      real pH   !pH of particles
      real log_alpha_nh3

C-----CODE--------------------------------------------------------------

      !ammonium molar ratio to sulfate
      AS=mnh4/mso4*2.667

      !Derive pH of particle. 
      !The detail equation is received by email from Dr. Pathak
      pH=-3.50E-3*rh**2.-7.5E-5*rh+0.5*AS**2-0.25
     &       *AS-6.0E-3*(rh*AS)-1+0.5

      !The following polynormial extracted from figure 3 of the second
      !literature by hand. Temperature range is from 285 - 291 K.
      log_alpha_nh3=-4.0536E-6*pH**6.+1.177E-4*pH**5.-0.0012*pH**4.
     &       +0.0042*pH**3.+0.0219*pH**2.-0.3838*pH-0.7223
      alpha_nh3=10.**log_alpha_nh3




      RETURN
      END
