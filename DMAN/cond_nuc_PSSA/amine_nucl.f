
C     **************************************************
C     *  amine_nucl                                     *
C     **************************************************

C     WRITTEN BY Jan and Ben, October 2014

C     This subroutine calculates the ternary nucleation rate and radius of the 
C     critical nucleation cluster using the parameterization of...

c     Napari, I., M. Noppel, H. Vehkamaki, and M. Kulmala. "Parametrization of 
c     Ternary Nucleation Rates for H2so4-Nh3-H2o Vapors." Journal of Geophysical 
c     Research-Atmospheres 107, no. D19 (2002).

      SUBROUTINE amine_nucl(tempi,csi,cnai,dma_i,fn,rnuc)

      IMPLICIT NONE
 
      include 'sizecode.COM'

C-----INPUTS------------------------------------------------------------

      double precision tempi                ! temperature of air [K]
      double precision csi                  ! condensation sink [s-1]
      double precision cnai                 ! concentration of gas phase sulfuric acid [molec cm-3]
      double precision dma_i                ! concentration of dimethyl amine [ppt]

C-----OUTPUTS-----------------------------------------------------------

      double precision fn                   ! nucleation rate [cm-3 s-1]
      double precision rnuc                 ! critical cluster radius [nm]
c      double precision tonset

C-----INCLUDE FILES-----------------------------------------------------

C-----ARGUMENT DECLARATIONS---------------------------------------------

C-----VARIABLE DECLARATIONS---------------------------------------------

      double precision fnl                  ! natural log of nucleation rate
      double precision temp                 ! temperature of air [K]
      double precision cs                   ! condensation sink [s-1]
      double precision cna                  ! concentration of gas phase sulfuric acid [molec cm-3]
      double precision dma                  ! concentration of gas phase dimethyl amine [ppt]
c
      integer ii                 ! counter

C
      temp = tempi  !K
      cs   = csi    !s-1
      cna  = cnai   !molec cm-3
      dma  = dma_i  !ppt
       
      !Limit All Values to Upper Bound on Lookup Table
      temp = max(min(temp,320.0),180.0)
      cs   = max(min(cs,2.0e-1),1.0e-5)
      cna  = max(min(cna,3.16e9),1.e4)
      dma  = max(min(dma,1.0e9),1.e4)

      !Locate the indices of all the independent variables
      itemp = locate(amine_nuc_tbl_TEMP, temp)

      !Define Nearest Neighbors for Each Index
      itemp0 = floor(itemp)
      itemp1 = ceiling(itemp)

      !Use Multilinear Interpolation even though it is a rather coarse,
      !inaccurate method
      ii1  =  (temp - amine_nuc_tbl(itemp0)) / (amine_nuc_tbl(itemp1) - amine_nuc_tbl(itemp0) )
      

      return

      end

