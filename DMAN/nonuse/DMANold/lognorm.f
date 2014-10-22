
C     *************************************************
C     * logdist                                       *
C     *************************************************

C     RE-WRITTEN BY JaeGun Jung, November 2007

C     This subroutine loads a lognormal distribution into matirx pn(nsect)

C-----INPUTS------------------------------------------------------------

C     numbr - a total number concentration of the log-normal distribution 
C     dpm - a mean diameter of the log-normal distribution 
C     sigma - a standard deviation of the distribution

C-----OUTPUTS-----------------------------------------------------------

C     pn - each portion of the log-normal distribution corresponding to 
C         specific size sections     

      SUBROUTINE logdist(numbr,dpm,sigma,pn)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real numbr
      real dpm
      real sigma
      real pn(nsect)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i ! a counter for size section

      real cntot ! total number of particles form the distribution
      real dp ! mean diameters of each size section
      real f ! a kind of dN/dlog10(Dp), whose sum equal to 1

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      ! If there are no particles then return
      if (numbr .lt. 0.1) then
        do i = 1, nsect
          pn(i) =0.0
        enddo
        return
      endif

      ! Call lognorm, then get the portion from the distribution
      do i=1, nsect
        dp = dpmean(i)
        call lognorm(dp, dpm, sigma, f)
        pn(i) = f*(log10(dpbound(i+1))-log10(dpbound(i)))
      enddo

      ! Calculate total CN fraction from pn, which is still fraction of 1
      cntot =0.0
      do i=1, nsect
        cntot = cntot+pn(i) ! cntot would be close to 1.
      enddo

      ! Calculate the final number concentration from the fraction
      do i=1, nsect
      pn(i)=numbr*pn(i)/cntot
      enddo



      RETURN
      END


C-----------------------------------------------------------------------


C     *************************************************     
C     * lognorm                                       *
C     *************************************************

C     RE-WRITTEN BY JaeGun Jung, November 2007

C     This subroutine gives a kind of dN/dlog10(Dp) value corresponding 
C     to each size bins, whose summation is 1. It uses a mean diameter of 
C     the log-normal distribution, dp and a standard deviation of sigma. 
C     dpm is a mean diameter of each size section.

C-----INPUTS------------------------------------------------------------

C     dp - a mean diameter of log-normal distribution
C     dpm - mean diameters of each size bin
C     sigma - a standard deviation

C-----OUTPUTS-----------------------------------------------------------

C     f - dN/dlog10(Dp) value for the size bin

      SUBROUTINE lognorm(dp,dpm,sigma,f)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real dp  ! a mean diameter of log-normal distribution 
      real dpm ! mean diameters of each size bin
      real sigma ! a standard deviation
      real f !dN/dlog10(Dp) value for the size bin

C-----VARIABLE DECLARATIONS---------------------------------------------

      real a
      real b
      real c
      real pi

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (pi = 3.14)

C-----CODE--------------------------------------------------------------

      a = SQRT(2.*pi)*log10(sigma)
      b = log10(dp) - log10(dpm)
      c = 2.*log10(sigma)*log10(sigma)
      f=  exp(-b*b/c)/a


      RETURN
      END
