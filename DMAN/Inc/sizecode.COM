
C	*************************************************************
C	* Size-Resolved Aerosol Microphysics header file            *
C	*************************************************************

C	This header file includes all the variables used by the
C	size-resolved aerosol microphysics code incorporated into
C	the GISS GCM II' by Peter Adams.  The microphysics algorithm
C	conserves aerosol number and mass using the schemes developed
C	by Graham Feingold and others.

C	Tzivion, S., Feingold, G., and Levin, Z., An Efficient
C	   Numerical Solution to the Stochastic Collection Equation,
C	   J. Atmos. Sci., 44, 3139-3149, 1987.
C	Feingold, G., Tzivion, S., and Levin, Z., Evolution of
C	   Raindrop Spectra. Part I: Solution to the Stochastic
C	   Collection/Breakup Equation Using the Method of Moments,
C	   J. Atmos. Sci., 45, 3387-3399, 1988.
C	Tzivion, S., Feingold, G., and Levin, Z., The Evolution of
C	   Raindrop Spectra. Part II: Collisional Collection/Breakup
C	   and Evaporation in a Rainshaft, J. Atmos. Sci., 46, 3312-
C	   3327, 1989.
C	Feingold, G., Levin, Z., and Tzivion, S., The Evolution of
C	   Raindrop Spectra. Part III: Downdraft Generation in an
C	   Axisymmetrical Rainshaft Model, J. Atmos. Sci., 48, 315-
C	   330, 1991.

C	The algorithms described in these papers have been extended
C	to include multicomponent aerosols and modified for a moving
C	sectional approach.  Using this approach, the boundaries
C	between size bins are defined in terms of dry aerosol mass
C	such that the actual sizes of the sections move as water
C	is added to or lost from the aerosol.

C	All of the subroutines needed for this aerosol microphysics
C	algorithm use only their own internal variables or the ones
C	listed here.  GISS GCM II' variables are not used (a driver
C	subroutine performs the necessary swapping between the GCM
C	and the microphysics code).  The microphysics code is,
C	therefore, completely modular.

C--------------------------------------------------------------------

C ibins is the number of size bins and icomp the number of chemical
C species.  There are 4 aerosol species (SO4--, NO3-, NH4+, and H2O).
C Water is assumed to always be the last component and is mostly
C ignored (see above) because this algorithm uses a moving sectional
C approach.  Most loops are from 1 to icomp-idiag as a result.

      integer ibins, icomp
      integer idiag ! number of diagnostic aerosol species
cjgj      parameter(ibins=43, icomp=4, idiag=2)
      parameter(ibins=41, icomp=4, idiag=1)

      integer srtso4, srtna, srtnh4, srth2o
      parameter (srtso4=1,
     &           srtna=2,
     &           srtnh4=3,
     &           srth2o=4)

C Nk and Mk contain the number and mass size distributions of the
C aerosol.  Units are #/grid cell or kg/grid cell, respectively.
C Gc are gas phase concentrations (kg/grid cell) of species
C corresponding to the aerosol species (e.g. H2SO4 for sulfate).
C Nkd and Mkd store values of Nk and Mk for diagnostic purposes.

      double precision Nk(ibins), Mk(ibins,icomp), Gc(icomp-1)
      double precision Nkd(ibins), Mkd(ibins,icomp), Gcd(icomp-1) 

C    Nnuc is the number of particles in the nucleation pseudo-bin
C    #/grid cell 
c      double precision Nnuc
C    Mnuc is the mass of the particles in the nucleation pseudo-bin
C    kg/grid cell
c      double precision Mnuc(icomp)

C    The following array contains the mass boundaries of the size
C    bins.  xk(1) is lower boundary of the first size bins in terms
C    of mass (kg).  xk(2) is the boundary between the first and
C    second bins, etc...  The scheme used here assumes mass doubling:

C	xk(k+1) = 2 * xk(k)

      double precision xk(ibins+1)

C The following variables describe the grid cell in which the
C microphysics is operating.

      double precision boxvol     !volume of grid cell (cm3)
      double precision boxmass    !volume of grid cell (kg)
cdbg      double precision temp       !temperature (K) of grid cell
cdbg      double precision pres       !air pressure (Pa) of grid cell
cdbg      double precision rh         !relative humidity (0-1)
c
      real temp       !temperature (K) of grid cell
      real pres       !air pressure (Pa) of grid cell
      real rh         !relative humidity (0-1)

C Physical properties of aerosol components

      real molwt(icomp)
cjgj      data molwt/96.,58.45,18.,18./ ! The 2nd is NaCl
      data molwt/96.,50.,18.,18./ ! The 2nd is organic matter. NEED TO BE UPDATED. 
      real gmw(icomp)
      data gmw/98.,50.,17.,18./
      real alpha(icomp) ! accomodation coef  
      data alpha/0.65,0.,0.08,0./
           ! alpha(H2SO4) is from U. Poschl et al., J. Phys. Chem. A, 102,
                                                        !10082-10089, 1998
           ! alpha(NH3) is from E. Swartz et al., J. Phys. Chem. A, 103,
                                                        !8824-8833, 1999
      real Sv(icomp)         !parameter used for estimating diffusivity
      data Sv /42.88,42.88,42.88,42.88/

C Flag for which nucleation parameterizations to use
      integer bin_nuc, tern_nuc, amine_nuc   !flags for binary and ternary nuc
      parameter(bin_nuc=1, tern_nuc=1, amine_nuc=1) ! 1 = on

cdbg	  double precision NH3ppt_o
cdbg	  double precision h2so4_print

C Make variables common

      common /sizecode/ Nk, Mk, Gc, Nkd, Mkd, Gcd, xk, 
     &                  boxvol, boxmass, temp, pres, rh
cdbg     &                  h2so4_print
cdbg     &                  nh3ppt_o, h2so4_print




