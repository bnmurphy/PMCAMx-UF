
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
      SUBROUTINE dman(tstart,tend,Nki,Mki,h2so4,nh3ppt,rhi,tempi,presi
     &               ,dsulfdt,organic_ppt,ichm, jchm, kchm, fndt)

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
      real dsulfdt ! sulfuric acid production rate
      integer ichm ! i coordinate in PMCAMx
      integer jchm ! j coordinate in PMCAMx
      integer kchm ! k coordinate in PMCAMx
cd    david new variable
      integer icond_soa  

      double precision fndt(njnuc) !nucleation diagnostic

C-----VARIABLE DECLARATIONS---------------------------------------------

C     Definition of Variables
C     =======================

C     deltat : operator time step (hr)
C     pn(i) : aerosol number concentration (particles cm-3) 
C     tstart : simulation start time from midnight (hr)
C     tend : simulation end time from midnight (hr)

ctest      real cn
      double precision cn
      double precision xkDMAN(ibins+1),Moo
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
      real organic_ppt(4)
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
      double precision sum_mass_dman,sum_num_dman
      double precision sum_so4_dman,sum_init_dman,sum_nh4_dman

      double precision sum_soa1_dman,sum_soa2_dman,sum_soa3_dman
      double precision sum_soa4_dman 
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
      icond_soa = 1 
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
     
       Moo = 3.7531528494783419e-25 ! 0.8 nm diameter particle assuming 1400 kg m-3 density

      do k=1,ibins+1
         xkDMAN(k)=Moo*2.d0**(k-1)
      enddo
        


      do i=1,ibins
        Nk(i)=Nki(i)
        do j=1,icomp
          Mk(i,j)=Mki(i,j)
        enddo
      enddo

    
      !Initialize ygas
      ygas(mgsvi)=h2so4
      ygas(mgnh3)=nh3ppt
      rh=rhi
      temp=tempi
      pres=presi

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


       Gc(srtsoa1)=boxmass*organic_ppt(1)*1.0d-12*200.0/28.9
       Gc(srtsoa2)=boxmass*organic_ppt(2)*1.0d-12*200.0/28.9
       Gc(srtsoa3)=boxmass*organic_ppt(3)*1.0d-12*200.0/28.9
       Gc(srtsoa4)=boxmass*organic_ppt(4)*1.0d-12*200.0/28.9
          

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


      !Call condensation with "dt"
      if (icond .eq. 1) then
c
        Gc(srtso4)=boxmass*ygas(mgsvi)*1.0d-12*gmw(srtso4)/28.9
        Gc(srtnh3)=boxmass*ygas(mgnh3)*1.0d-12*gmw(srtnh3)/28.9


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
      

      !Convert sulfuric acid production rate from ppt hr-1
      H2SO4rate = dsulfdt*1.0d-12*boxmass/28.9*gmw(srtso4)/3600. !PSSA
                                              ! [=] kg box-1 s-1

      !Call Pseudo Steady State condensation and nucleation
      if (inucl .eq. 1) then
        call cond_nuc(Nk,Mk,Gc,Nkout,Mkout,Gcout,H2SO4rate,dt,
     &   ichm,jchm,kchm, fndt )
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

        call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm)


c---------------------------------------------------------------------------
c     !Call condensation with "dt" of organics 
      if (icond_soa.eq.1) then

        !Converting gas unit from ppt to kg
       Gc(srtsoa1)=boxmass*organic_ppt(1)*1.0d-12*200.0d0/28.9d0
       Gc(srtsoa2)=boxmass*organic_ppt(2)*1.0d-12*200.0d0/28.9d0
       Gc(srtsoa3)=boxmass*organic_ppt(3)*1.0d-12*200.0d0/28.9d0
       Gc(srtsoa4)=boxmass*organic_ppt(4)*1.0d-12*200.0d0/28.9d0


        call org_cond(Nk,Mk,Gc,Nkout,Mkout,Gcout,dt,xkDMAN,
     &     ichm,jchm,kchm)

        do i=1,ibins
          do j=1,icomp
            Mk(i,j)=Mkout(i,j)
          enddo
          Nk(i)=Nkout(i)
        enddo
        Gc(srtsoa1) = Gcout(srtsoa1)
        Gc(srtsoa2) = Gcout(srtsoa2)
        Gc(srtsoa3) = Gcout(srtsoa3)
        Gc(srtsoa4) = Gcout(srtsoa4)

       organic_ppt(1)=Gc(srtsoa1)*1.0d+12/boxmass*28.9d0/200.0d0 !gmw(srtsoa1)
       organic_ppt(2)=Gc(srtsoa2)*1.0d+12/boxmass*28.9d0/200.0d0 !gmw(srtsoa2)
       organic_ppt(3)=Gc(srtsoa3)*1.0d+12/boxmass*28.9d0/200.0d0 !gmw(srtsoa3)
       organic_ppt(4)=Gc(srtsoa4)*1.0d+12/boxmass*28.9d0/200.0d0 !gmw(srtsoa4)

       call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm)
      end if
c---------------------------------------------------------------------------

cd      call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm)

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

        call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm)
cPMCAMx        call neutral(4) 
      endif

C-----Print cn, cn20, ccn and H2SO4 [kg] on screen

      time = time +deltat

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

      
       organic_ppt(1)=Gc(srtsoa1)*1.0d+12/boxmass*28.9d0/200.0d0 !gmw(srtsoa1)
       organic_ppt(2)=Gc(srtsoa2)*1.0d+12/boxmass*28.9d0/200.0d0 !gmw(srtsoa2)
       organic_ppt(3)=Gc(srtsoa3)*1.0d+12/boxmass*28.9d0/200.0d0 !gmw(srtsoa3)
       organic_ppt(4)=Gc(srtsoa4)*1.0d+12/boxmass*28.9d0/200.0d0 !gmw(srtsoa4)



      !For a debuging purpose
      if ((ichm.eq.90).and.(jchm.eq.97).and.(kchm.eq.1)) then
         write(*,*)'-------------------------------------------------'
         write(*,*)'In dman at the finishing dman finokalia'
         write(*,*)'coordinate of (i,j,k)',ichm, jchm, kchm
         write(*,*)'tempi,dsulfdt=',tempi,dsulfdt
         write(*,*)'organic(ppt)=',organic_ppt
         write(*,*)'h2so4=',h2so4,'nh3ppt=',nh3ppt
         write(*,*)'Nk='
         sum_num_dman=0.0
         do i=1, ibins
           sum_num_dman=sum_num_dman+Nki(i)
         enddo
         write(*,*)'sum_num_dman=',sum_num_dman
         write(*,*)'Mk='
         sum_mass_dman=0.0
         sum_soa1_dman=0.0
         sum_soa2_dman=0.0
         sum_soa3_dman=0.0
         sum_soa4_dman=0.0

         sum_so4_dman=0.0
         sum_init_dman=0.0
         sum_nh4_dman=0.0

         do j=1, icomp-1
           do i=1, ibins
             sum_mass_dman=sum_mass_dman+Mki(i,j)
             if (j.eq.srtsoa1) then
             sum_soa1_dman=sum_soa1_dman+Mki(i,j)
             end if
             if (j.eq.srtsoa2) then
             sum_soa2_dman=sum_soa2_dman+Mki(i,j)
             end if
             if (j.eq.srtsoa3) then
             sum_soa3_dman=sum_soa3_dman+Mki(i,j)
             end if
             if (j.eq.srtsoa4) then
             sum_soa4_dman=sum_soa4_dman+Mki(i,j)
             end if

             if (j.eq.srtso4) then
             sum_so4_dman=sum_so4_dman+Mki(i,j)
             end if
             if (j.eq.srtinrt) then
             sum_init_dman=sum_init_dman+Mki(i,j)
             end if
             if (j.eq.srtnh4) then
             sum_nh4_dman=sum_nh4_dman+Mki(i,j)
             end if
           enddo
         enddo

         write(*,*)'sum_so4_dman=',sum_so4_dman,srtso4
         write(*,*)'sum_soa1_dman=',sum_soa1_dman,srtsoa1
         write(*,*)'sum_soa2_dman=',sum_soa2_dman,srtsoa2
         write(*,*)'sum_soa3_dman=',sum_soa3_dman,srtsoa3
         write(*,*)'sum_soa4_dman=',sum_soa4_dman,srtsoa4

         write(*,*)'sum_init_dman=',sum_init_dman
         write(*,*)'sum_nh4_dman=',sum_nh4_dman
         write(*,*)'========================================'
      endif    


      RETURN
      END
