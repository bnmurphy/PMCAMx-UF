c----------------------------------------------------------------------
c Function DMDT_INT:  Here we apply the analytic solution to the
c droplet growth equation in mass space for a given scale length which
c mimics the inclusion of gas kinetic effects
c
c M0 ......... initial mass
c L0 ......... length scale
c Tau ........ forcing from vapor field
c

Cpja I have changed the length scale.  Non-continuum effects are
Cpja assumed to be taken into account in choice of tau (in so4cond
Cpja subroutine).

Cpja I have also added another argument to the function call, WR.  This
Cpja is the ratio of wet mass to dry mass of the particle.  I use this
Cpja information to calculate the amount of growth of the wet particle,
Cpja but then return the resulting dry mass.  This is the appropriate
Cpja way to implement the condensation algorithm in a moving sectional
Cpja framework.

      double precision FUNCTION DMDT_INT(M0,TAU,WR)
 
      IMPLICIT NONE
      double precision M0,TAU,X,L0,C,ZERO,WR,MH2O
      PARAMETER (C=2.d0/3.d0,L0=0.0d0,ZERO=0.0d0)
 
      MH2O=(WR-1.d0)*M0
      X=((M0+MH2O)**C+L0)
c      X=MAX(ZERO,SQRT(MAX(ZERO,C*TAU+X))-L0)
c      X=MAX(ZERO,SQRT(C*TAU+X))    ! commented out by LA
      if((C*TAU+X).gt.ZERO) then
         X=SQRT(C*TAU+X)
      else
         X=ZERO
      end if
      DMDT_INT=X*X*X-MH2O

Cpja Perform some numerical checks on dmdt_int
      if ((tau .gt. 0.0) .and. (dmdt_int .lt. m0)) dmdt_int=m0
      if ((tau .lt. 0.0) .and. (dmdt_int .gt. m0)) dmdt_int=m0

      RETURN
      END
