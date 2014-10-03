      subroutine soap(ntot,caer,cgas,tempk,convfac,
     &                idiag,iout,igrdchm,ichm,jchm,kchm,lppm,
     &                cpre,mwpre,csatT)
      implicit none
c
c**********************************************************************
c                                                                     * 
c                                                                     * 
c  SOAP - A SUBROUTINE TO CALCULATE GAS-AEROSOL PARTITIONING OF       *
c                     SECONDARY ORGANIC SPECIES                       *
c                                                                     *
c                           DEVELOPED BY:                             *
c                                                                     * 
c                           ROSS STRADER                              *
c            DEPARTMENT OF CIVIL/ENVIRONMENTAL ENGINEERING            *
c             CARNEGIE MELLON UNIVERSITY, 5000 FORBES AVE             *
c                   PITTSBURGH, PENNSYLVANIA, 15213                   *
c                                                                     * 
c                               AND                                   *
c                                                                     * 
c                         SPYROS N. PANDIS                            *
c       DEPARTMENTS OF CHEMICAL ENGINEERING AND ENGINEERING AND       *
c                          PUBLIC POLICY                              *
c             CARNEGIE MELLON UNIVERSITY, 5000 FORBES AVE             *
c                   PITTSBURGH, PENNSYLVANIA, 15213                   *
c                                                                     * 
c                                                                     * 
c                                                                     * 
c                                                                     * 
c      ASSOCIATED SUBROUTINES:                                        *
c      CALC       - CALCULATES THE MOLE FRACTION OF EACH SOLUTION-    *
c                   FORMING COMPOUND IN THE PSEUDO-IDEAL SOLUTION     * 
c                   (AEROSOL PHASE)                                   *
c      FCN        - NONLINEAR SYSTEM TO BE SOLVED BY HYBRD            *
c      RNDNUM     - GENERATES RANDOM NUMBERS                          *
c      HYBRD      - SOLVES A SYSTEM OF NONLINEAR EQUATIONS            *
c                                                                     * 
c                                                                     * 
c**********************************************************************
c                                                                     * 
c                                                                     * 
c     THIS MODEL USES A PSEUDO-IDEAL SOLUTION THEORY TO PARTITION     *
c     SECONDARY ORGANIC SPECIES BETWEEN THE AEROSOL AND GAS PHASE.    *
c     THE THEORY IS OUTLINED IN:                                      * 
c                                                                     * 
c     STRADER, R. (1998), 'EVALUATION OF SECONDARY ORGANIC AEROSOL    *
c          FORMATION IN WINTER,' M.S. THESIS, DEPT. OF CIVIL/         *
c          ENVIRONMENTAL ENGINEERING, CARNEGIE MELLON UNIVERSITY      *
c                                                                     * 
c                                                                     * 
c**********************************************************************
c
c
c  MODIFICATIONS
c
c     Revised for new SOAP driver by bkoo (03/09/03)
c     Modified by bkoo (01/20/03)
c     Revised for CAMx by Greg Yarwood, January 2003
c
c DEFINITION OF VARIABLES:
c
c  INPUTS
c
c     ntot    - total number of CG/SOA species pairs
c     caer    - aerosol-phase concentrations of SOA species (ug/m3)
c     cgas    - gas-phase concentrations of CG species (ppm or ug/m3)
c     tempk   - temperature (K)
c     convfac - conversion factor: umol/m3 = ppm * convfac 
c     idiag   - diagnostic output file unit
c     iout    - standard output file unit
c     igrdchm - index for grid containing the grid cell
c     ichm    - i index of grid cell
c     jchm    - j index of grid cell
c     kchm    - k index of grid cell
c     lppm    - gases in ppm if true, ugm3 if false
c
c   OUTPUTS
c
c     caer    - aerosol-phase concentrations of SOA species (ug/m3)
c     cgas    - gas-phase concentrations of CG species (ppm or ug/m3)
c
c   VARIABLES USED WITHIN SOAP
c
c     i       - counter
c     icont   - counter
c     sum     - counter
c     nsol    - total number of solution-forming SOA species
c     cstemp  - temperatures corresponding to saturation concentrations
c               of CG/SOA species (K)
c     csat    - saturation concentrations of CG/SOA species (ug/m3)
c     deltah  - enthalpy of vaporization of CG/SOA species (kJ/mol)
c     flagsoap- set to 1 if CG/SOA species forms solutions; 0 if not
c     mwsoap  - molecular weights of CG/SOA species (g/mol)
c     scaer   - aerosol-phase concentrations of solution-forming 
c               SOA species (ug/m3)
c     scgas   - gas-phase concentrations of solution-forming 
c               SOA species (ug/m3)
c     scsat   - saturation concentrations of solution-forming 
c               SOA species (ug/m3)
c     sctot   - total concentrations of solution-forming SOA species
c               (ug/m3)
c     smw     - molecular weights of solution-forming SOA species 
c               (g/mol)
c     sx      - mole fraction of solution-forming SOA species in the
c               pseudo-ideal solution
c     znum    - number of iterations HYBRD goes through before finding
c               solution
c     cpre    - concentration of pre-existing organic aerosol (ug/m3)
c     mwpre   - molecular weight of pre-existing organic aerosol (g/mol)
c     caerin  - input aerosol-phase concentrations of SOA species (ug/m3)
c     conmin  - use simple solution for species below this level(ug/m3)
c     lae3    - true to emulate CMAQ AE3 algorithm (no evporation)
c
c
c***********************************************************************
c
c
c VARIABLE DECLARATION
c
      include 'soap.com'
c
      integer     ntot
      real        caer(ntot), cgas(ntot)
      real        ctot(NSOAP), caerin(NSOAP)
      real        smw(NSOAP),scsat(NSOAP)
      real        sctot(NSOAP),scaer(NSOAP)
      real        scgas(NSOAP),sx(NSOAP)
      real        csatT(NSOAP)
      integer     idx(NSOAP)
c
      real        mwpre,cpre,tempk,sum,convfac,conmin
      integer     iout, idiag, igrdchm, ichm, jchm, kchm
      integer     i,icont,nsol,znum
      logical     lppm
c
c***********************************************************************
c
c Entry point
c
      if(ntot.ne.NSOAP) then 
        write(iout,'(//,A)') 'ERROR in SOAP:'
        write(iout,*) 'ERROR: ntot.ne.NSOAP in subroutine soap'
        call camxerr()
      endif
      do i=1,ntot
        if (lppm) cgas(i) = cgas(i)*convfac*mwsoap(i)
        ctot(i) = caer(i) + cgas(i)
        caerin(i) = caer(i)
      enddo
c
c Assume no pre-existing aerosol phase, such as primary
c organic aerosol.
c now check pflag to set cpre - bkoo (11/13/03)
      if(pflag.eq.0) cpre = 0.0

      znum = 0
c
c CHANGE SATURATION CONCENTRATIONS ACCORDING TO CURRENT TEMPERATURE
c
      do i=1,ntot
         csatT(i)=csat(i)*(cstemp(i)/tempk)*exp((deltah(i)/8.314)
     &                *(1/cstemp(i)-1/tempk))
      enddo     
c
c MAP COMPOUNDS THAT FORM SOLUTIONS ONTO ARRAYS TO PASS TO SUBROUTINE 
c CALC.  COMPOUNDS THAT HAVE A CONCENTRATION OF LESS THAN 1E-04 ARE 
c IGN0RED.
c
c CALCULATE AEROSOL-PHASE CONCENTRATION (CAER) AND GAS-PHASE 
c CONCENTRATION (CGAS) FOR NON-SOLUTION-FORMING COMPOUNDS
c
      icont=0
      conmin=0.0001
      do i=1,ntot
         if (flagsoap(i).eq.0) then
            cgas(i) = amin1(ctot(i), csatT(i))
            caer(i) = ctot(i) - cgas(i)
         elseif (ctot(i).lt.conmin) then
            caer(i) = caerin(i)
            cgas(i) = ctot(i) - caer(i)
         else
            icont=icont+1
            idx(icont) = i
            smw(icont)=mwsoap(i)
            scsat(icont)=csatT(i)
            sctot(icont)=ctot(i)
            scaer(icont)=caer(i)
         endif
      enddo
      nsol=icont
c
c Check for a trivial solution
c
      if (nsol.eq.0) goto 1000
      if (nsol.eq.1 .and. cpre.eq.0.0) then
         scgas(1) = amin1(sctot(1), scsat(1))
         scaer(1) = sctot(1) - scgas(1)
         goto 900
      endif
      sum=0.0
      do i=1,nsol
         sum = sum + sctot(i)/scsat(i)
      enddo
      if (cpre.eq.0.0 .and. sum.lt.1.0) then
         do i=1,nsol
            scgas(i)=sctot(i)
            scaer(i)=0.0
         enddo
         goto 900
      endif
c
c Initial guess for the aerosol organic phase mole fractions
c
c Good guess, sometimes converges poorly from here
c
c      sum = 0.0
c      do i=1,nsol
c         sx(i) = amax1(0.01,scaer(i)/smw(i))
c         sum = sum + sx(i)
c      enddo
c      do i=1,nsol
c         sx(i) = sx(i)/sum
c      enddo
c
c Simplistic guess is more robust
c
      if (pflag.eq.1) then ! absorbing POA
        do i=1,nsol
          sx(i) = 0.1/nsol
        enddo
      else
        do i=1,nsol
          sx(i) = 1.0/nsol
        enddo
      endif
c      
c CALL SUBROUTINE CALC TO CALCULATE THE MOLE FRACTION (SX) OF EACH
c SOLUTION-FORMING COMPOUND IN THE PSEUDO-IDEAL SOLUTION (AEROSOL 
c PHASE)
c
      call calc(nsol,sx,smw,scsat,sctot,mwpre,cpre,znum,
     &          idiag,iout,igrdchm,ichm,jchm,kchm,tempk)
c
c CALCULATE AEROSOL-PHASE CONCENTRATION (SCAER) AND GAS-PHASE 
c CONCENTRATION (SCGAS) FOR SOLUTION-FORMING COMPOUNDS
c
      do i=1,nsol
         if (sctot(i).le.sx(i)*scsat(i)) then
            scaer(i) = 0.0
            scgas(i) = sctot(i)
         else
            scgas(i) = sx(i)*scsat(i)
            scaer(i) = sctot(i) - scgas(i)
         endif
      enddo
c
c REMAP COMPOUNDS THAT FORM SOLUTIONS BACK ONTO ORIGINAL ARRAYS
c      
 900  continue
      do i=1,nsol
            caer(idx(i))=scaer(i)
            cgas(idx(i))=scgas(i)
      enddo
c
c --- For emulation of the CMAQ AE3 algorithm, do not allow
c     any SOA to evaporate back to the gas phase.  The flag
c     that controls this (lae3) is set in soapdat.f
c
 1000 continue
      if (lae3) then
         do i=1,ntot
            caer(i) = amax1(caer(i), caerin(i))
            cgas(i) = ctot(i) - caer(i)
         enddo
      endif
c
      if (lppm) then
         do i=1,ntot
            cgas(i) = cgas(i)/(convfac*mwsoap(i))
         enddo
      endif
c
      return
      end

