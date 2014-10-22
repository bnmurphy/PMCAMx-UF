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
        

c 
