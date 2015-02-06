
c***********************************************************************
c
c            'INCLUDE' FILE WITH AEROSOL PARAMETERS 
c
c                                            Last change:  JGJ 6/10/05
c                                                       :  TMG 5/27/03
c                                                       :  COS 4/24/02 
c***********************************************************************
c
c Include TOMAS library file
c
cPMCAMx      include 'sizecode.COM'
c
c Setting size sections from input files
c
      integer dres    ! a resolution below measurement
      integer nsectbm ! a number of size sections below measurements
      integer nsectm  ! a number of measured sections with measurements
                      ! It was 186 until 10 um.
      integer nsecto  ! a number of total size section with dres of 1
      integer nsect   ! the number of aerosol section with variable dres
      integer inmax   ! a maximum number of data for input variables
      integer intime  ! a number of row for hourly changable input variables

      parameter (dres=1)
      parameter (nsectbm = 35*dres) ! number of sections below measurements
      parameter (nsectm = 195) ! number of measured sections with measurements
                               ! for initial distribution. It was 186 by 10 um.
      parameter (nsecto = nsectbm+nsectm)
      parameter (nsect  = nsectbm+nsectm*dres)
      parameter (inmax = 150) ! It was 145.
      parameter (intime = 25)
c
c USEFUL CONSTANTS
c
      real pi
      real pi6
      real pres0 ! 1.0133E+5 Pa = 1atm
      real temp0 ! 298.15 K
      real dpmin ! (um) lower boundary of smallest section
      
      parameter (pi  = 3.14159)
      parameter (pi6 = pi/6.0)
      parameter (pres0 = 1.0133E+5) ! [=] Pa
      parameter (temp0 = 298.15) ! [=] K
      parameter (dpmin = 0.0008)
c
c Environmental setting and speices flags
c
      integer ngas  ! number of active and buildup species in gas-phase array
      integer mgso2 ! sulfur dioxide
      integer mgsvi ! sulfuric acid
      integer mgnh3 ! ammonia
      integer mgdma ! dimethyl amine
      integer srtorg ! a species flag of organic matter in aerosol phase
      integer srtnh3 ! a species flag of ammonium in aerosol phase

      parameter (ngas = 4)   
      parameter (mgso2 = 1)  ! sulfur dioxide
      parameter (mgsvi = 2)  ! sulfuric acid
      parameter (mgnh3 = 3)  ! ammonia
      parameter (mgdma = 4)  ! dimethyl amine
      parameter (srtorg=2, srtnh3=3)
cdbg      parameter (srtso4=1, srth2o=4) ! moved to sizecode.COM
c
c NUMBER OF EQUATIONS FOR THE AEROSOL OPERATOR
c
      integer keqna

      parameter (keqna = 3)         ! number of variables dealt in nucl_sub
                                    ! pn(1), H2SO4, and NH3 are currently dealt
                                    ! there. It was 2 pn(1) and H2SO4.
c
c THE FLAGS BELOW ARE USED TO TURN ON AND OFF PROCESSES IN THE MODEL
c
      integer idnuc
      integer iwknd
      integer kprint ! a interval for priniting output files
      integer icoag_test
      integer icond_test
      integer iin_flag
      integer ievaluation ! swith to decide whether it makes output file 
      integer iexptest    ! experimental test
      integer CF_flag ! a flag for nucleation of Clement and Ford (1999)
      integer nh3flag ! ammonia emission flag

      parameter (idnuc = 0)  ! =1 calculate nuclei diam using parameterization 
                             ! =0 set fresh nuclei diam =dpmean(1)
      parameter (iwknd = 0)       ! =0 weekday
      parameter (kprint = 60)            ! PRINT FREQUENCY
      parameter (icoag_test = 0)  ! =0 disable coag test 6/24/04 jgj
      parameter (icond_test = 0)  ! =0 disable cond test jgj
      parameter (ievaluation = 1) ! evaluation for print output files
      parameter (iexptest = 0)    ! =1  Do experimental test
      parameter (iin_flag = 0)    ! =1  Do ion induced nucleation
      parameter (CF_flag  = 0)    ! =1  Do Clement and Ford (1999)
      parameter (nh3flag = 1)    ! =1  Do ammonia mass balance
c
c Variables to interface with the TOMAS algorithm     
c
      real dp4xk(ibins+1) ! boundary diameters of TOMAS size section
      real dpmean4xk(ibins) ! mean diameter of TOMAS size section
      double precision xko(nsect) ! mass boundaries of an input file
c
c Variables in common statements
c
      integer kcount ! a counter for printing output files
      real dpmean(nsect) ! a mean diameter of both sectional boundaries, um
      real dpbound(nsect+1) ! sectional boundaries, um
      real dpboundo(nsecto)
      real h2ogrowth
      real itcommon
      real rad
      real pm(nsect,icomp)
      real nh3flux
cPMCAMx      real dsulfdt !sulfuric acid production rate, moved to dman.f
cdbg      real nh3ppt
cdbg      real temp, rh ! moved to sizecode.COM
c
c Diagnostic variables in TOMAS algorithm
c
      double precision mnfixdiagn(ibins) 
                       ! number concentration change by mnfix subroutine
      double precision mnfixdiagm(ibins, icomp-1)
                       ! mass concentration change by mnfix subroutine
      double precision conddiagn(ibins)
                       ! number concentration change by condensation, tmconds
      double precision conddiagm(ibins, icomp-1)
                       ! mass concentration change by condesation, tmconds
      double precision coagdiagn(ibins)
                       ! number concentration change by coagulation, multicoag
      double precision coagdiagm(ibins, icomp-1)
                       ! mass concentration change by coagulation, multicoag
c
c COMMON STATEMENTS
c
      common /kcount/ kcount
      common /diametros/ dpmean,dpbound, dpboundo
      common /wetdiameter/ h2ogrowth
      common /xronos/ itcommon
      common /meteorology/ rad
      common /convert/ xko
      common /pm/ pm
      common /nh3flux/ nh3flux
cPMCAMx      common /sulfuric/ dsulfdt
cPMCAMx      common /xk/ dp4xk, dpmean4xk
      common /tomasdiag/ mnfixdiagn, mnfixdiagm, conddiagn, conddiagm,
     &		      coagdiagn, coagdiagm
cdbg      common /nh3ppt/ nh3ppt
