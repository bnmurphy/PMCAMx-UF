
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

cPMCAMx      PROGRAM DMAN
      SUBROUTINE dman(tstart,tend,Nki,Mki,h2so4,nh3ppt,rhi,tempi,presi,
     &                fsurfi,dsulfdt,ichm, jchm, kchm, fndt)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'
      include 'aervaria.inc'
cPMCAMx      include 'IO.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real tstart ! start time
      real tend   ! end time
      double precision Nki(ibins)       !Number in a box 
      double precision Mki(ibins,icomp) !Mass in a box as kg
      real h2so4   ! h2so4 from PMCAMx [=] ppt
      real nh3ppt  ! nh3ppt from PMCAMx
      real rhi     ! relative humidity from PMCAMx
      real tempi   ! temperature from PMCAMx
      real presi   ! pressure from PMCAMx
      real fsurfi(11)  ! LandUse Category from PMCAMx
      real dsulfdt ! sulfuric acid production rate
      integer ichm ! i coordinate in PMCAMx
      integer jchm ! j coordinate in PMCAMx
      integer kchm ! k coordinate in PMCAMx
      double precision fndt(nJnuc)  !nucleation diagnostic

C-----VARIABLE DECLARATIONS---------------------------------------------

C     Definition of Variables
C     =======================

C     deltat : operator time step (hr)
C     pn(i) : aerosol number concentration (particles cm-3) 
C     tstart : simulation start time from midnight (hr)
C     tend : simulation end time from midnight (hr)

ctest      real cn
      double precision cn
      real deltat
cPMCAMx      real dtgas
cPMCAMx      real dndlog(nsect)
cPMCAMx      real dndlogo(nsecto)
      real dt ! [=] sec
      real fgas(ngas) ! [=] ppt
      real h2ostd ! a standard h2o growth rate
cPMCAMx      real mettime(inmax)
      real Neps
cPMCAMx      real nh3time(intime)
cPMCAMx      real nh3dat(intime)
      real pn(nsect) ! [=] particles cm-3
cPSSA      real pn_p(nsect) ! [=] cm-3 change by nucleation or deposition
      real R ! gas constant
cPMCAMx      real raddat(inmax)
cPMCAMx      real ratio
cPMCAMx      real rhdat(inmax)
      real remainder ! for print step size
cPMCAMx      real so2dat(inmax)
cPMCAMx      real so2time(inmax)
      real sulf0
      real sulf1
      real t1
      real t2
cPMCAMx      real tempdat(inmax)
cPMCAMx      real tend !moved up as an argument variable
      real time
      real tgas
cPMCAMx      real tstart !moved up as argument variable
cPMCAMx      real winddat(inmax)
      real ygas(ngas) ! [=] ppt

      ! for PSSA_cond_nuc
      double precision Nkout(ibins) ! Number output
      double precision Mkout(ibins,icomp) ! Mass output
      double precision Gcout(icomp-1) ! Gas output
      double precision H2SO4rate ! sulfuric acid production rate [kg box-1 s-1]
cPMCAMx      double precision NH3rate ! ammonia emission rate [kg hr-1]

      integer i ! a counter for size section
      integer icoag ! a coagulation flag
      integer icond ! a condensation flag
      integer inucl ! a nucleation flag
cPMCAMx      integer indx ! index of high size resolved particle bins
      integer j ! a counter for species
      integer k ! a counter for size section
      integer iflagez ! if =1, call eznh3eqm
cPMCAMx      integer metcount ! a meteorology counter for input
cPMCAMx      integer metterms ! a number of lines in meteorology input files
cPMCAMx      integer so2count ! a SO2 counter for input
cPMCAMx      integer so2terms ! a number of lines in SO2 input files
cPMCAMx      integer nh3count ! a NH3 counter for input
cPMCAMx      integer nh3terms ! a number of lines in NH3 input files
cPMCAMx      integer PM25cutpoint ! PM2.5 cutpoint in terms of TOMAS

cPMCAMx      character*6 date ! simulation date

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (R = 8.314) ! gas constant, J/mol-K
      parameter (Neps = 1.0d-10) ! an epsilon of number concentration

C-----CODE--------------------------------------------------------------

C-----Initialize variables
    
      !Initialize flag variables
      icond = 1
      icoag = 1
      inucl = 1

cdbg      print*,'In dman, check arguments'
cdbg      print*,'tstart=', tstart ! start time
cdbg      print*,'tend=', tend   ! end time
cdbg      print*,'h2so4=', h2so4   ! h2so4 from PMCAMx [=] ppt
cdbg      print*,'nh3ppt=', nh3ppt  ! nh3ppt from PMCAMx
cdbg      print*,'rhi=', rhi     ! relative humidity from PMCAMx
cdbg      print*,'tempi=', tempi   ! temperature from PMCAMx
cdbg      print*,'presi=',presi   ! pressure from PMCAMx
cbdg      print*,'dsulfdt=', dsulfdt ! sulfuric acid production rate
cdbg      write(*,*)'coordinate (i,j,k)=', ichm, jchm, kchm ! i coordinate in PMCAMx
      !Debug for dsulfdt
cdbg      if ((ichm.eq.65).and.(jchm.eq.51).and.(kchm.eq.1)) then
cdbg        write(*,*)'dsulfdt in Pittsburgh=',dsulfdt
cdbg        write(*,*)'ygas(mgsvi)=',ygas(mgsvi)
cdbg      endif

      do i=1,ibins
        Nk(i)=Nki(i)
        do j=1,icomp
          Mk(i,j)=Mki(i,j)
        enddo
      enddo

      !For a debuging purpose
cdbg      if ((ichm.eq.65).and.(jchm.eq.51).and.(kchm.eq.1)) then
cdbg         write(*,*)'In dman at the very beginning' 
cdbg         write(*,*)'coordinate of (i,j,k)',ichm, jchm, kchm
cdbg         write(*,*)'tempi,presi,dsulfdt=',tempi,presi,dsulfdt
cdbg         write(*,*)'h2so4=',h2so4,'nh3ppt=',nh3ppt
cdbg         write(*,*)'Nk='
cdbg         do i=1, ibins
cdbg           write(*,*)Nk(i)
cdbg         enddo
cdbg         write(*,*)'Mk='
cdbg         do j=1, icomp
cdbg           write(*,*)'j=',j
cdbg           do i=1, ibins
cdbg             write(*,*)Mk(i,j)
cdbg           enddo
cdbg         enddo
cdbg       endif
    
      !Initialize ygas
      ygas(mgsvi)=h2so4
      ygas(mgnh3)=nh3ppt
      rh=rhi
      temp=tempi
      pres=presi
      fsurf=fsurfi

      !Set environmental variables
      deltat=tend-tstart
      dt=deltat*3600.
      boxvol = 3.0d+20 ! a volume of a grid cell, cm3

      ! Update boxmass from a new temperature
      boxmass=0.0289*pres*boxvol*1.0d-6*(1./(R*temp))
                                             ! a mass of grid cell(kg)
C-----Read and initialize

      ! Set start and end times for simulation
      time = tstart

      !Convert pn to Nk
      Call initbounds()

      !Check pt 0 for debug
cdbg      if ((ichm.eq.79).and.(jchm.eq.54).and.(kchm.eq.1).and.
cdbg     &  (time.ge.7.)) then
cdbg      if ((ichm.eq.82).and.(jchm.eq.32).and.(kchm.eq.4).and.
cdbg     &  (time.gt.2.25)) then
cdbg        write(*,*)'dt=',dt,'Coord. (i,j,k)=',ichm,jchm,kchm
cdbg        write(*,*)'dsulfdt=',dsulfdt
cdbg        write(*,*)'rh=',rh,'temp=',temp,'pres=',pres
cdbg      endif

C
C-----DO A MAIN LOOP----------------------------------------------------  
C
      do while (time.lt.tend)        !Outer Loop: time

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

cPMCAMx      Call printf(deltat, pn, tend, time, tstart, ygas)

      t1=time
      t2=time+deltat

C-----Micro physical processes 6/9/04

      !Check pt 1 for debug
cdbg      if ((time.gt.10.75).and.(time.lt.11.25)) then
cdbg        if ((ichm.eq.76).and.(jchm.eq.44).and.(kchm.eq.1)) then
cdbg          write(*,*)'In dman at the beginning' 
cdbg          write(*,*)'coordinate of (i,j,k)',ichm, jchm, kchm
cdbg          write(*,*)'tempi,presi,dsulfdt=',tempi,presi,dsulfdt
cdbg          write(*,*)'h2so4=',h2so4,'nh3ppt=',nh3ppt
cdbg          write(*,*)'deltat=',deltat
cdbg          write(*,*)'Nk='
cdbg          do i=1, ibins
cdbg            write(*,*)Nk(i)
cdbg          enddo
cdbg            write(*,*)'Mk='
cdbg          do j=1, icomp
cdbg            write(*,*)'j=',j
cdbg            do i=1, ibins
cdbg              write(*,*)Mk(i,j)
cdbg            enddo
cdbg          enddo
cdbg        endif
cdbg      endif
      !moved these two lines here from inside the icond if statement
      !otherwise icond=0 means that there is junk values in Gc
      Gc(srtso4)=boxmass*ygas(mgsvi)*1.0d-12*gmw(srtso4)/28.9
      Gc(srtnh3)=boxmass*ygas(mgnh3)*1.0d-12*gmw(srtnh3)/28.9

      !Call condensation with "dt"
      if (icond .eq. 1) then
        !Converting gas unit from ppt to kg
cPMCAMx        if (icond_test.eq.1) then
cPMCAMx          Gc(srtso4)=boxmass*ygas(mgsvi)*1.0e-12*100.0/28.9
cPMCAMx        else
        
cPMCAMx        endif

        Call so4cond(Nk,Mk,Gc,Nkout,Mkout,Gcout,dt,ichm,jchm,kchm,
     &              iflagez) 
                                                           !Only ammonia
        if (iflagez.eq.1) then
          call eznh3(Gc,Mk,Nk,Gcout,Mkout,Nkout,ichm,jchm,kchm)
        endif

        do i=1,ibins
          do j=1,icomp
            Mk(i,j)=Mkout(i,j)
          enddo
          Nk(i)=Nkout(i)
        enddo
        Gc(srtnh3) = Gcout(srtnh3)
        Gc(srtso4) = Gcout(srtso4)

        !Converting Nk, Mk, and Gc
        do i=1,ibins
          do j=1, icomp
            Mk(i,j)=Mkout(i,j)   
          enddo
          Nk(i)=Nkout(i)
        enddo
        if (icond_test .ne. 1) then !Notice that not equal,"ne"
          Gc(srtso4)=Gcout(srtso4)
          Gc(srtnh3)=Gcout(srtnh3)
          ygas(mgsvi)=Gc(srtso4)*1.0d+12/boxmass*28.9/gmw(srtso4)
          ygas(mgnh3)=Gc(srtnh3)*1.0d+12/boxmass*28.9/gmw(srtnh3)
        endif
      endif

      !Check pt 2 for debug
cdbg      if ((time.gt.10.75).and.(time.lt.11.25)) then
cdbg        if ((ichm.eq.76).and.(jchm.eq.44).and.(kchm.eq.1)) then
cdbg          write(*,*)'In dman after so4cond' 
cdbg          write(*,*)'coordinate of (i,j,k)',ichm, jchm, kchm
cdbg          write(*,*)'tempi,presi,dsulfdt=',tempi,presi,dsulfdt
cdbg          write(*,*)'h2so4=',h2so4,'nh3ppt=',nh3ppt
cdbg          write(*,*)'Nk='
cdbg          do i=1, ibins
cdbg            write(*,*)Nk(i)
cdbg          enddo
cdbg            write(*,*)'Mk='
cdbg          do j=1, icomp
cdbg            write(*,*)'j=',j
cdbg            do i=1, ibins
cdbg              write(*,*)Mk(i,j)
cdbg            enddo
cdbg          enddo
cdbg        endif
cdbg      endif

      !Call coagulation with dt/2
      if (icoag .eq. 1) then
        call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm)
        call multicoag(dt/2.,time,Nk,Mk,Nkout,Mkout,ichm,jchm,kchm) !coagulation
c
        do i=1,ibins
          do j=1, icomp
            Mk(i,j)=Mkout(i,j)   
          enddo
          Nk(i)=Nkout(i)
        enddo
cPMCAMx        call neutral(2)
cdbg        write(*,*)'2'
      endif
      
      !Check pt 3 for debug
cdbg      if ((time.gt.10.75).and.(time.lt.11.25)) then
cdbg        if ((ichm.eq.76).and.(jchm.eq.44).and.(kchm.eq.1)) then
cdbg          write(*,*)'In dman after 1st multicoag' 
cdbg          write(*,*)'coordinate of (i,j,k)',ichm, jchm, kchm
cdbg          write(*,*)'tempi,presi,dsulfdt=',tempi,presi,dsulfdt
cdbg          write(*,*)'h2so4=',h2so4,'nh3ppt=',nh3ppt
cdbg          write(*,*)'Nk='
cdbg          do i=1, ibins
cdbg            write(*,*)Nk(i)
cdbg          enddo
cdbg            write(*,*)'Mk='
cdbg          do j=1, icomp
cdbg            write(*,*)'j=',j
cdbg            do i=1, ibins
cdbg              write(*,*)Mk(i,j)
cdbg            enddo
cdbg          enddo
cdbg        endif
cdbg      endif

      !Convert sulfuric acid production rate from ppt hr-1
      H2SO4rate = dsulfdt*1.0d-12*boxmass/28.9*gmw(srtso4)/3600. !PSSA
                                              ! [=] kg box-1 s-1

      !Call Pseudo Steady State condensation and nucleation
      if (inucl .eq. 1) then
        call cond_nuc(Nk,Mk,Gc,Nkout,Mkout,Gcout,H2SO4rate,dt,
     &   ichm,jchm,kchm,fndt)
      endif

      do i=1,ibins
        do j=1,icomp
          Mk(i,j)=Mkout(i,j)
        enddo
        Nk(i)=Nkout(i)
      enddo
      Gc(srtnh3) = Gcout(srtnh3)
      Gc(srtso4) = Gcout(srtso4)

      ygas(mgsvi)=Gc(srtso4)*1.0d+12/boxmass*28.9/gmw(srtso4)
      ygas(mgnh3)=Gc(srtnh3)*1.0d+12/boxmass*28.9/gmw(srtnh3)

      !Check pt 4 for debug
cdbg      if ((time.gt.10.75).and.(time.lt.11.25)) then
cdbg        if ((ichm.eq.76).and.(jchm.eq.44).and.(kchm.eq.1)) then
cdbg          write(*,*)'In dman after cond_nuc' 
cdbg          write(*,*)'coordinate of (i,j,k)',ichm, jchm, kchm
cdbg          write(*,*)'tempi,presi,dsulfdt=',tempi,presi,dsulfdt
cdbg          write(*,*)'h2so4=',h2so4,'nh3ppt=',nh3ppt
cdbg          write(*,*)'Nk='
cdbg          do i=1, ibins
cdbg            write(*,*)Nk(i)
cdbg          enddo
cdbg            write(*,*)'Mk='
cdbg          do j=1, icomp
cdbg            write(*,*)'j=',j
cdbg            do i=1, ibins
cdbg              write(*,*)Mk(i,j)
cdbg            enddo
cdbg          enddo
cdbg        endif
cdbg      endif

      !Call coagulation with dt/2
      if (icoag .eq. 1) then
        call multicoag(dt/2.,time,Nk,Mk,Nkout,Mkout,ichm,jchm,kchm) !coagulation
c
        do i=1,ibins
          do j=1, icomp
            Mk(i,j)=Mkout(i,j)   
          enddo
          Nk(i)=Nkout(i)
        enddo

      !Check pt 4.5 for debug
cdbg      if ((time.gt.10.75).and.(time.lt.11.25)) then
cdbg        if ((ichm.eq.76).and.(jchm.eq.44).and.(kchm.eq.1)) then
cdbg          write(*,*)'In dman after 2nd multicoag' 
cdbg          write(*,*)'coordinate of (i,j,k)',ichm, jchm, kchm
cdbg          write(*,*)'tempi,presi,dsulfdt=',tempi,presi,dsulfdt
cdbg          write(*,*)'h2so4=',h2so4,'nh3ppt=',nh3ppt
cdbg          write(*,*)'Nk='
cdbg          do i=1, ibins
cdbg            write(*,*)Nk(i)
cdbg          enddo
cdbg            write(*,*)'Mk='
cdbg          do j=1, icomp
cdbg            write(*,*)'j=',j
cdbg            do i=1, ibins
cdbg              write(*,*)Mk(i,j)
cdbg            enddo
cdbg          enddo
cdbg        endif
cdbg      endif
        call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm)
cPMCAMx        call neutral(4) 
      endif

C-----Print cn, cn20, ccn and H2SO4 [kg] on screen
cPMCAMx      if (icond_test.eq.1) then
cPMCAMx        remainder=mod(time,1.667d-2) ! every 1 min
cPMCAMx        if ((remainder.le.deltat/2.).or.((1.667d-2-remainder)
cPMCAMx     &    .le.deltat/2.)) then
cPMCAMx          call Nk2pn(pn,cn) ! converting Nk to pn
cPMCAMx          write(6,101)time,temp,rh,cn,ygas(mgsvi),ygas(mgnh3)
cPMCAMx        endif
cPMCAMx      else  
cPMCAMx        if (time.le.(itcommon+deltat/2)) then
cPMCAMx          call Nk2pn(pn,cn) ! converting Nk to pn
cPMCAMx          cn=0.0
cPMCAMx          do i=1,22 ! about 0.1 um
cPMCAMx            cn=cn+Nk(i)/boxvol
cPMCAMx          enddo
cPMCAMx          print*,'boxvol=',boxvol
cPMCAMx          pause
cPMCAMx          ygas(mgsvi)=Gc(srtso4)*1.0d+12/boxmass*28.9/gmw(srtso4)
cPMCAMx          ygas(mgnh3)=Gc(srtnh3)*1.0d+12/boxmass*28.9/gmw(srtnh3)
cdbg          write(6,101)time,temp,rh,cn,ygas(mgsvi),ygas(mgnh3)
cPMCAMx        endif
cPMCAMx      endif
cPMCAMx 101  format(6g13.5)
   
      !Re-initialize pn
cPMCAMx      do i=1, nsect
cPMCAMx        if (pn(i) .lt. 0.) pn(i)=0.
cPMCAMx      enddo

      !Check pt 5 for debug
cdbg      if ((time.gt.10.75).and.(time.lt.11.25)) then
cdbg        if ((ichm.eq.76).and.(jchm.eq.44).and.(kchm.eq.1)) then
cdbg          write(*,*)'In dman after 2nd multicoag and mnfix' 
cdbg          write(*,*)'coordinate of (i,j,k)',ichm, jchm, kchm
cdbg          write(*,*)'tempi,presi,dsulfdt=',tempi,presi,dsulfdt
cdbg          write(*,*)'h2so4=',h2so4,'nh3ppt=',nh3ppt
cdbg          write(*,*)'Nk='
cdbg          do i=1, ibins
cdbg            write(*,*)Nk(i)
cdbg          enddo
cdbg            write(*,*)'Mk='
cdbg          do j=1, icomp
cdbg            write(*,*)'j=',j
cdbg            do i=1, ibins
cdbg              write(*,*)Mk(i,j)
cdbg            enddo
cdbg          enddo
cdbg        endif
cdbg      endif

      time = time +deltat
cPMCAMx      kcount = kcount +1

      enddo    ! main loop

      !Return gas concentrations
      ygas(mgsvi)=Gc(srtso4)*1.0d+12/boxmass*28.9/gmw(srtso4)
      ygas(mgnh3)=Gc(srtnh3)*1.0d+12/boxmass*28.9/gmw(srtnh3)

      do i=1,ibins
        Nki(i)=Nk(i)
        do j=1,icomp
          Mki(i,j)=Mk(i,j)
        enddo
      enddo

      h2so4 = ygas(mgsvi)
      nh3ppt = ygas(mgnh3)

cdbg      write(*,*)'h2so4=',h2so4,' nh3ppt=',nh3ppt

    


      RETURN
      END
