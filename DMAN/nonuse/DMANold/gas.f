
C     *************************************************
C     * diff                                          *
C     *************************************************

C      RE-WRITTEN BY JaeGun Jung, November 2007

C      This subroutine is differential equation for gas concentration

C-----INPUTS------------------------------------------------------------      

C     time - time in a main program
C     ygas - gas species mixing ratio as ppt
C     avgtemp - average temperature during gas phase reaction

C-----OUTPUTS-----------------------------------------------------------

C     fgas - change rate of gas species per hour

      SUBROUTINE diff(time,ygas,fgas,avgtemp)
 
      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real time
      real ygas(ngas)
      real avgtemp
      real fgas(ngas) ! [=] ppt hr-1

C-----VARIABLES DECLARATIONS--------------------------------------------

      real coh ! concentration of OH radical [=] molec cm-3
      real cno3 ! concentration of NO3 radical [=] molec cm-3
      real cgas(ngas) ! gas species mixing ratio as ppt
      real lr(ngas) ! loss rate of gas species per unit ppt of the gas
                    ! [=] hr-1
      real fr(ngas) ! production rate of gas species [=] ppt hr-1

      integer i ! a counter for gas speices

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      ! Get OH and NO3 radical concentration
      call OH(time,coh,cno3)

      ! Calculate the contribution due to chemical reaction
      do i=1,ngas
        cgas(i)=ygas(i)
      enddo
      call gaschem(coh,cno3,cgas,fr,lr,avgtemp)

      ! Generate Right Hand Side of differential equations for solver
      do i=1, ngas
        fgas(I) = fr(I) - lr(I)*ygas(I)
      enddo


      RETURN
      END

C-----------------------------------------------------------------------

C     *************************************************
C     * gaschem                                       *
C     *************************************************

C      RE-WRITTEN BY JaeGun Jung, November 2007

C      This subroutine determines the formation and loss rates of 
C      gas-phase species due to gas-phase reactions. It could be 
C      replaced by difun and bldup if mechanisms of Carter's form are 
C      used to include more gas phase chemical reactions.

C-----INPUTS------------------------------------------------------------      

C     coh - OH radical concentration, molec cm-3
C     cno3 - NO3 radical concentration, molec cm-3
C     cgas - gas species mixing ratio as ppt
C     avgtemp - average temperature during gas phase reaction

C-----OUTPUTS-----------------------------------------------------------

C     fr - gas production rate
C     lr - gas loss rate

      SUBROUTINE gaschem(coh,cno3,cgas,fr,lr,avgtemp)
 
      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real coh         ![=] molec cm-3
      real cno3        ![=] molec cm-3
      real cgas(ngas)  !array of active and buildup gas-phase species
      real fr(ngas)    !formation rates of gas-phase species
      real lr(ngas)    !loss rates of gas-phase species
      real avgtemp     !average temperture during chemical reactions

C-----VARIABLES DECLARATIONS--------------------------------------------

      real rkso2       !an oxidation  reaction constant from SO2 to sulfate
      real a0          !a parameter for rkso2
      real ainf        !a parameter for rkso2
      real m           !a parameter for rkso2

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      !SO2 gas phase oxidation to sulfate from IMAGES
      a0=3e-31*(300/avgtemp)**3.3
      ainf=1.5e-12
      m=6.02e23*22.4/1000*298/avgtemp   !in molecules/cm3
      rkso2=a0*m/(1+a0*m/ainf)*0.6**
     &    (1/(1+log10((a0*m/ainf)**2)))*3600. ! in cm3/(molec hr)

      !SO2 chemical reaction rate
      lr(mgso2) = rkso2*coh

      !H2SO4 chemical reaction rate
      fr(mgsvi) = rkso2*coh*cgas(mgso2)
      lr(mgsvi) = 0.0

      RETURN
      END

C-----------------------------------------------------------------------

C     *************************************************
C     * OH                                            *
C     *************************************************

C      RE-WRITTEN BY JaeGun Jung, November 2007

C      This subroutine calculates variation of OH. See Gaydos et al.
C      (2005) in J. Geophys. Res. for details. For NO3 mixing ratio,
C      it seems consulting Lin et al, JGR. 
C                           - I could not find exact source 11/24/07

C-----INPUTS------------------------------------------------------------      

C     time - time of a main program

C-----OUTPUTS-----------------------------------------------------------

C     coh - OH radical concentration, molec cm-3
C     cno3 - NO3 radical concentration, molec cm-3

      SUBROUTINE OH(time,coh,cno3)
 
      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real time
      real coh
      real cno3

C-----VARIABLES DECLARATIONS--------------------------------------------

      integer iday ! an order of days from the beginning of the simulation
      integer ihour ! an order of hour during a day

      real t ! time during a day

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      ! Calculate day and hour
      iday = INT(time/24) + 1
      ihour = INT(time-(iday-1)*24)

C-----Set radical concentration during night time
      if (ihour .ge. 0 .AND. ihour .lt. 6) then
        coh = 0.0
        cno3 = 0.25/3.72e-8          !in molec/cm3
        go to 100
      endif

      if (ihour .ge. 18 .AND. ihour .le. 24) then
        coh=0.0
        cno3 = 0.25/3.72e-8          !in molec/cm3
        go to 100
      endif

      ! calculate time during a day
      t = time - (iday-1)*24.

C-----Set radical concentration during day time
      ! Calculate coh. ADJUST MAXIMUM DEPENDING ON SEASON!!!!
      coh = 5.E6*rad/3.815d1 ! The maximum of coh is 5.0E6 and 1.0E6 
                             ! for summer and winter, respectively.
                             ! The average coh is 2.0E6 molec cm-3.
ctmg      coh = 5.E6 * SIN(pi*(t-6.)/12.)    !in molec/cm3

      cno3 = 0.0             ! NO3 is turned off during day time.


 100  continue

      RETURN
      END
