
C     **************************************************      
C     * the Dynamic Model for Aerosol Nucleation       *
C     *                                   (DMAN)       *
C     **************************************************      
 
C     Written by JaeGun Jung, November 2007

C     This program is a box model which is used in PMCAMx for
C     aerosol microphysical process module. This program is 
C     developed taking advantages of two moment aerosol sectional
C     algorithm based on Gaydos et al. (2005)* model. 

C     *Gaydos et al., 2005. Modeling of in-situ ultrafine atmospheric
C     particle formation in the eastern United State. J. Geophys. Res.
C     110, D07S12.

C     For the details of development history, see history.help file.

      PROGRAM DMAN

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'
      include 'IO.inc'

C-----VARIABLE DECLARATIONS---------------------------------------------

C     Definition of Variables
C     =======================

C     deltat : operator time step (hr)
C     pn(i) : aerosol number concentration (particles cm-3) 
C     tstart : simulation start time from midnight (hr)
C     tend : simulation end time from midnight (hr)

      real cn
      real deltat
      real dtgas
      real dndlog(nsect)
      real dndlogo(nsecto)
      real dt ! [=] sec
      real fgas(ngas) ! [=] ppt
      real h2ostd ! a standard h2o growth rate
      real mettime(inmax)
      real Neps
      real nh3time(intime)
      real nh3dat(intime)
      real pn(nsect) ! [=] particles cm-3
      real pn_p(nsect) ! [=] cm-3 change by nucleation or deposition
      real R ! gas constant
      real raddat(inmax)
      real ratio
      real rhdat(inmax)
      real so2dat(inmax)
      real so2time(inmax)
      real sulf0
      real sulf1
      real t1
      real t2
      real tempdat(inmax)
      real tend
      real time
      real tgas
      real tstart
      real winddat(inmax)
      real ygas(ngas) ! [=] ppt

      ! for PSSA_cond_nuc
      double precision Nkout(ibins) ! Number output
      double precision Mkout(ibins,icomp) ! Mass output
      double precision Gcout(icomp-1) ! Gas output
      double precision H2SO4rate ! sulfuric acid production rate [kg box-1 s-1]
      double precision NH3rate ! ammonia emission rate [kg hr-1]

      integer i ! a counter for size section
      integer icoag ! a coagulation flag
      integer icond ! a condensation flag
      integer inucl ! a nucleation flag
      integer indx ! index of high size resolved particle bins
      integer iterso4 ! How many times iterated in so4cond.
                      ! It is printed in so4cond.itr.
      integer j ! a counter for species
      integer metcount ! a meteorology counter for input
      integer metterms ! a number of lines in meteorology input files
      integer so2count ! a SO2 counter for input
      integer so2terms ! a number of lines in SO2 input files
      integer nh3count ! a NH3 counter for input
      integer nh3terms ! a number of lines in NH3 input files
      integer PM25cutpoint ! PM2.5 cutpoint in terms of TOMAS

      character*12 limit ! What kind of limitation is applied in so4ocnd.
                         ! It is printed in so4cond.itr.
      character*6 date ! simulation date

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (R = 8.314) ! gas constant, J/mol-K
      parameter (Neps = 1.0e-20) ! an epsilon of number concentration

C-----CODE--------------------------------------------------------------

C-----Initialize variables
    
      !Initialize flag variables
      icond = 1
      icoag = 1
      inucl = 1

      !Initialize count varables
      metcount = 0
      so2count = 0
      nh3count = 0

      !Initialize Nk, Mk, and Gc
      do i=1,ibins
        Nk(i)=0.0
        do j=1,icomp
          Mk(i,j)=0.0
        enddo
      enddo
      do j=1,icomp
        Gc(j)=0.0
      enddo

      !Initialize ygas
      do j=1,3
        ygas(j)=0.0
      enddo

      !Initialize parameters
      iterso4 = 0
      limit = ''
      metterms= 0 !meteorology total number of terms
      so2terms= 0 !so2 total number of terms
      nh3terms= 0 !nh3 total number of terms

      !Set environmental variables
      date = '091202' ! set date
      dt = 15. ! a time step [sec]
      deltat = dt/3600. ! a time step [hr]
      pres=1.01325e5 ! initial pressure [=] Pa , set as 1 atm
      temp=298.15    ! initial temperature [=] K
      boxvol = 3.0e+20 ! a volume of a grid cell, cm3
      PM25cutpoint = 35 ! It corresponds to 2.6 um.

      !Ammonia emission rate in case of a calculation of ammonia mass balance
      nh3flux = 3.0 ! [=]ppt/hr
      ygas(mgnh3) = 0.0
cdbg      ygas(mgnh3) = 10000.0

      !Set initial sulfuric acid conditions
      if (icond_test .eq. 1) then  ! cond_test 06/24/04 jgj
        ygas(mgsvi) = 1.0e+3     ! in ppt 
      else
        ygas(mgsvi) = 0.
      endif

      !Set a print counter
      kcount = kprint

C-----OPEN I/O UNITS 
      call startup(date)

C-----Get initial aerosol distribution from file, then calculate sectional
C    boundaries and mean diameters for aerosol distribution by tmg (05/27/03)

      ratio = (3.0d-3/dpmin)**(1./FLOAT(nsectbm))
      dpbound(1)=dpmin

      do i=2,nsectbm 
        dpbound(i)=dpbound(i-1)*ratio
      enddo

      do i=nsectbm+1,nsecto+1
        if (i.lt.nsecto+1) then
          read(idat,*)dpboundo(i),dndlogo(i)
          if (dpboundo(i) .ge. 2.)dndlogo(i)=0.0 ! to fit with measurement
        else
          read(idat,*)dpboundo(i)
        endif
      enddo

      do i=1, nsectbm
          dndlog(i)=Neps ! not to make nan by jgj
      enddo

      do i=nsectbm+1,nsecto
        ratio=(dpboundo(i+1)/dpboundo(i))**(1./FLOAT(dres))
        indx=(i-nsectbm-1)*dres+nsectbm
        dpbound(indx+1)=dpboundo(i)
        if (dndlogo(i) .eq. 0.) dndlogo(i)=Neps ! not to make nan
        dndlog(indx+1)=dndlogo(i)
      enddo
      dpbound(nsect+1)=dpboundo(nsecto+1)

      do i=1,nsect
        pn(i)=dndlog(i)*log10(dpbound(i+1)/dpbound(i))
      enddo

      do i=1, nsect
        dpmean(i) = SQRT(dpbound(i)*dpbound(i+1))
      enddo

C-----Read and initialize

      ! Set start and end times for simulation
      if ((icond_test .eq. 1) .or.(icoag_test .eq. 1)) then
        tstart = 0.0
        tend = 2.
      else
        read(idat,*)tstart
cchg        tstart = 7.0
        tend = 24.
      endif

      time = tstart

      !Convert pn to Nk
      Call initbounds()
      Call initconv(pn,ygas)
      
      !Read meteorology, SO2, and NH3 data from input files
      Call initinput(mettime,tempdat,rhdat,raddat,winddat,so2time,
     & so2dat,nh3time,nh3dat,metterms,so2terms,nh3terms,date,tend)
 
      !Print screen header
      write(6,*) '   time           T          RH             CN      
&H2SO4 [ppt]   NH3 [ppt]'

C
C-----DO A MAIN LOOP----------------------------------------------------  
C

      do while (time.le.tend)        !Outer Loop: time

      itcommon = int(time+deltat/2.) ! Only hour value goes to itcommon

C-----Experimental Evaluation Requested by Dr.Spyros Pandis 6/2/05
C     Only do when experimental test swith is turned on.

      if (iexptest .eq. 1) then  
        !Initialize distribution at 12:00 P.M.
        if ((time.gt.(12.-deltat)) .and. (time.le.(12.+deltat))) then 
          call experiment(ygas,icoag,icond,inucl)
        endif
      endif

C-----Read input files

      Call readinput(time,mettime,tempdat,rhdat,raddat,winddat,so2time,
     & so2dat,nh3time,nh3dat,metterms,so2terms,nh3terms,metcount,
     & so2count,nh3count,ygas,date)
 
      ! Update boxmass from a new temperature
      boxmass=0.0289*pres*boxvol*1.0e-6*(1./(R*temp))
                                             ! a mass of grid cell(kg)

C-----Calculate water growth for a calculation of a wet diameter

      !Changed standard rh to 0.5 (indicating dry diameter given)
      !by tmg (5/16/03)

      h2ostd = (17.8*(0.5-0.5)**3.5 + 1.28)
      if (rh.gt.0.5) then
        h2ogrowth = (17.8*(rh-0.5)**3.5 + 1.28)/h2ostd
      else
        h2ogrowth =1.0d0
      endif

C-----Print output files

      Call printf(deltat, pn, tend, time, tstart, ygas)

C-----SO2 oxidation by OH - preprocessing

      !Step 1 - Save H2SO4(g) before doing oxidation, initialize variables
      sulf0= ygas(mgsvi)     ! [ppt]
      dtgas= deltat / 8.    ! in hours must be an integer increment of deltat
      tgas = time

      !Step 2 - Call oxidation subroutine
      do while(tgas.le.time+deltat)
        Call diff(tgas,ygas,fgas,temp)
        do i=1,ngas-1 ! exclude ammonia
          ygas(i)=ygas(i)+fgas(i)*dtgas
        enddo
        tgas = tgas + dtgas
      enddo

      !Step 3 - Get oxidation rate, dsulfdt
      sulf1 = ygas(mgsvi)
      dsulfdt = (sulf1-sulf0)/deltat    ! [ppt/hr]
      if (inucl.eq.1) ygas(mgsvi) = sulf0
        ! If inucl = 0, use a sulf1 value as ygas
      t1=time
      t2=time+deltat

C-----Ammonia concentration increase by a emission - preprocessing
      !If PSSA, turn on, otherwise nh3flux will be considered in derivnucl.
      NH3rate=(1.0e-12/28.9)*nh3flux*boxmass*gmw(srtnh3) ![kg/hr]
      Gc(srtnh3)=Gc(srtnh3)+NH3rate*deltat ![kg]

C-----Micro physical processes 6/9/04

      !Call condensation with dt/2
ceznh3eqm      if (icond .eq. 1) then
ceznh3eqm        Call so4cond(dt/2.,ygas,iterso4) !Also ammonia
ceznh3eqm        call neutral(1)
cdbg        write(*,*)'1'
ceznh3eqm      endif

      !Call coagulation with dt/2
      if (icoag .eq. 1) then
        call multicoag(dt/2.,time) ! coagulation
        call neutral(2)
cdbg        write(*,*)'2'
      endif
      
      !Convert sulfuric acid production rate from ppt hr-1
      H2SO4rate = dsulfdt*1.0e-12*boxmass/28.9*gmw(srtso4)/3600. !PSSA
                                              ! [=] kg box-1 s-1

cdbg       print*,'Before cond_nuc'
cdbg       do i=1,ibins
cdbg          write(*,*)Nk(i)
cdbg       enddo
cdbg       do j=1,icomp
cdbg          write(*,*)'species=',j          
cdbg          do i=1,ibins
cdbg             write(*,*)Mk(i,j)
cdbg          enddo
cdbg       enddo
cdbg       pause

      !Call Pseudo Steady State condensation and nucleation
      if (inucl .eq. 1) then
        call cond_nuc(Nk,Mk,Gc,Nkout,Mkout,Gcout,H2SO4rate,dt,xk,
     &   temp,pres,rh,boxvol)
      endif

      do i=1,ibins
        do j=1,icomp
          Mk(i,j)=Mkout(i,j)
        enddo
        Nk(i)=Nkout(i)
      enddo
      Gc(srtnh3) = Gcout(srtnh3)
      Gc(srtso4) = Gcout(srtso4)

cdbg       print*,'After cond_nuc'
cdbg       print*,'Nk='
cdbg       do i=1,ibins
cdbg          write(*,*)Nk(i)
cdbg       enddo
cdbg       do j=1,icomp
cdbg         write(*,*)'species=',j          
cdbg         do i=1,ibins
cdbg           write(*,*)Mk(i,j)
cdbg         enddo
cdbg       enddo
cdbg       pause

      !CALL NUCLEATION AND GAS RXN WITH dt

cPSSA      if (inucl .eq. 1) then
cPSSA        do i=1,nsect
cPSSA          pn_p(i)=0.
cPSSA        enddo
cdbg        Call Nk2pn(pn,cn) !Convert Nk to pn. Turn on when doing IIN.
cPSSA        Call nucl_sub(ygas,pn_p,t1,t2,pn,inucl)
cPSSA        call pn2Nk(pn_p)
cPSSA        do i=1,nsect
cPSSA          pn_p(i)=0.
cPSSA        enddo
cPSSA      endif
cPSSA      call neutral(3)
cdbg      write(*,*)'3'

      !Call coagulation with dt/2
      if (icoag .eq. 1) then
        call multicoag(dt/2,time) !coagulation
        call neutral(4) 
cdbg        write(*,*)'4'
      endif

      !Call condensation with dt/2
ceznh3eqm      if (icond .eq. 1) then
ceznh3eqm        Call so4cond(dt/2.,ygas,iterso4) !Also ammonia
ceznh3eqm        call neutral(5) 
cdbg        write(*,*)'5'
ceznh3eqm      endif

      if (ygas(mgsvi) .lt. 0.) ygas(mgsvi)=0. 
                     ! This is hard-wired stuff. I need to improve this.

C-----Print cn, cn20, ccn and H2SO4 [kg] on screen
      if (time.le.(itcommon+deltat/2)) then
        call Nk2pn(pn,cn) ! converting Nk to pn
        ygas(mgsvi)=Gc(srtso4)*1.0e+12/boxmass*28.9/gmw(srtso4)
        ygas(mgnh3)=Gc(srtnh3)*1.0e+12/boxmass*28.9/gmw(srtnh3)
        write(6,101)time,temp,rh,cn,ygas(mgsvi),ygas(mgnh3)
      endif
 101  format(6g13.5)
   
      !Re-initialize pn
      do i=1, nsect
        if (pn(i) .lt. 0.) pn(i)=0.
      enddo

      time = time +deltat
      kcount = kcount +1

      enddo    ! main loop



      STOP
      END
