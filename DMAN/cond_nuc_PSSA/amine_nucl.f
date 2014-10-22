
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

      end subroutine

!=============================================================
!
      subroutine read_amine_nuc_table
c
c
c     read_amine_nuc_table opens the lookup table for amine-sulfuric
c     acid nucleation and reads in the table values
c
c     Written: Ben Murphy and Jan Julin 10/17/14
c
c     Input arguments: 
c        none 
c 
c     Output arguments: 
c        All captured variables go to common block in sizecode.COM
c        amine_nuc_tbl_H2SO4 - sulfuric acid conc. [molec cm-3]
c        amine_nuc_tbl_DMA   - dimethyl amine conc. [ppt]
c        amine_nuc_tbl_CS    - condensaiton sink [s-1]
c        amine_nuc_tbl_TEMP  - temperature [K]
c        amine_nuc_tbl_J     - Nucleation Rate [Particles cm-3 s-1]
c            
c     Called by:
c        CHMPREP
c
      include 'sizecode.COM'
c
      integer iH2SO4, iDMA, iCS, iTEMP

      open(unit=98,file='DMAN/cond_nuc_PSSA/ACDC_H2SO4_DMA_05Feb2014.txt')

      !First read header and toss it
      read (98)

      !Now Start Reading in the Table
      do iTEMP = 1,amine_nuc_nTEMP
        do iCS = 1,amine_nuc_nCS
          do iH2SO4 = 1,amine_nuc_nH2SO4
            do iDMA = 1,amine_nuc_nDMA
	      read (98,*) amine_nuc_H2SO4(iH2SO4), amine_nuc_DMA(iDMA),
     &                    amine_nuc_TEMP(iTEMP),   amine_nuc_CS(iCS),
     &                    amine_nuc_J(iDMA, iH2SO4, iCS, iTEMP)
            enddo
          enddo
        enddo
      enddo

      close(98)

      return

      end subroutine
        
