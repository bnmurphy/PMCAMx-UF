
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

      SUBROUTINE dman(tstart,tend,Nki,Mki,h2so4,nh3ppt,dmappt,rhi,tempi,
     &               presi,dsulfdt,organic_ppt, ichm, jchm, kchm, fndt)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'
      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real tstart ! start time
      real tend   ! end time
      double precision Nki(ibins)       !Number in a box 
      double precision Mki(ibins,icomp) !Mass in a box as kg
      real h2so4   ! h2so4 from PMCAMx [=] ppt
      real nh3ppt  ! nh3ppt from PMCAMx
      double precision dmappt  ! dmappt from PMCAMx
      real rhi     ! relative humidity from PMCAMx
      real tempi   ! temperature from PMCAMx
      real presi   ! pressure from PMCAMx
      real dsulfdt ! sulfuric acid production rate
      integer ichm ! i coordinate in PMCAMx
      integer jchm ! j coordinate in PMCAMx
      integer kchm ! k coordinate in PMCAMx
cd    david new variable
      integer icond_soa  

      double precision fndt(nJnuc) !nucleation diagnostic

C-----VARIABLE DECLARATIONS---------------------------------------------

C     Definition of Variables
C     =======================

C     deltat : operator time step (hr)
C     pn(i) : aerosol number concentration (particles cm-3) 
C     tstart : simulation start time from midnight (hr)
C     tend : simulation end time from midnight (hr)

      double precision cn
      double precision xkDMAN(ibins+1),Moo
      real deltat
      real dt ! [=] sec
      real fgas(ngas) ! [=] ppt
      real h2ostd ! a standard h2o growth rate
      real Neps
      real pn(nsect) ! [=] particles cm-3
      real R ! gas constant
      real remainder ! for print step size
      real sulf0
      real sulf1
      real t1
      real t2
      real time
      real tgas
      real ygas(ngas) ! [=] ppt
      real organic_ppt(5)  !4 semivol + 1 extr.-low.-volat.
      ! for PSSA_cond_nuc
      double precision Nkout(ibins) ! Number output
      double precision Mkout(ibins,icomp) ! Mass output
      double precision Gcout(icomp-1) ! Gas output - DONT INCLUDE WATER
      double precision H2SO4rate ! sulfuric acid production rate [kg box-1 s-1]

      integer i ! a counter for size section
      integer icoag ! a coagulation flag
      integer icond ! a condensation flag
      integer inucl ! a nucleation flag
      integer j ! a counter for species
      integer k ! a counter for size section
      integer iflagez ! if =1, call eznh3eqm

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
      ygas(mgdma)=dmappt
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

      Gc(srtsoa1)=boxmass*organic_ppt(1)*1.0d-12*gmw(srtsoa1)/28.9
      Gc(srtsoa2)=boxmass*organic_ppt(2)*1.0d-12*gmw(srtsoa2)/28.9
      Gc(srtsoa3)=boxmass*organic_ppt(3)*1.0d-12*gmw(srtsoa3)/28.9
      Gc(srtsoa4)=boxmass*organic_ppt(4)*1.0d-12*gmw(srtsoa4)/28.9
      Gc(srtsoa5)=boxmass*organic_ppt(5)*1.0d-12*gmw(srtsoa5)/28.9 !EXLVOCS


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

      t1=time
      t2=time+deltat

C-----Micro physical processes 6/9/04

 
      !moved these two lines here from inside the icond if statement
      !otherwise icond=0 means that there is junk values in Gc
      Gc(srtso4)=boxmass*ygas(mgsvi)*1.0d-12*gmw(srtso4)/28.9
      Gc(srtnh3)=boxmass*ygas(mgnh3)*1.0d-12*gmw(srtnh3)/28.9
      Gc(srtdma)=boxmass*ygas(mgdma)*1.0d-12*gmw(srtdma)/28.9

      !Call condensation with "dt"
      if (icond .eq. 1) then

        Call so4cond(Nk,Mk,Gc,Nkout,Mkout,Gcout,dt,ichm,jchm,kchm,
     &              iflagez) 
                                                           !Only ammonia and amine
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
        Gc(srtdma) = Gcout(srtdma)

        if (icond_test .ne. 1) then !Notice that not equal,"ne"
          Gc(srtso4)=Gcout(srtso4)
          Gc(srtnh3)=Gcout(srtnh3)
          Gc(srtdma)=Gcout(srtdma)
          ygas(mgsvi)=Gc(srtso4)*1.0d+12/boxmass*28.9/gmw(srtso4)
          ygas(mgnh3)=Gc(srtnh3)*1.0d+12/boxmass*28.9/gmw(srtnh3)
          ygas(mgdma)=Gc(srtdma)*1.0d+12/boxmass*28.9/gmw(srtdma)
        endif
      endif

c---------------------------------------------------------
      !Call coagulation with dt/2
      if (icoag .eq. 1) then
cJJ        call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm,4)
        call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm)
        call multicoag(dt/2.,time,Nk,Mk,Nkout,Mkout,ichm,jchm,kchm) !coagulation
c
        do i=1,ibins
          do j=1, icomp
            Mk(i,j)=Mkout(i,j)   
          enddo
          Nk(i)=Nkout(i)
        enddo
      endif

      !Convert sulfuric acid production rate from ppt hr-1
      H2SO4rate = dsulfdt*1.0d-12*boxmass/28.9*gmw(srtso4)/3600. !PSSA
                                              ! [=] kg box-1 s-1

      !Call Pseudo Steady State condensation and nucleation
      if (inucl .eq. 1) then
        call cond_nuc(Nk,Mk,Gc,Nkout,Mkout,Gcout,H2SO4rate,dt,
     &   ichm,jchm,kchm, fndt)

        do i=1,ibins
           do j=1,icomp
              Mk(i,j)=Mkout(i,j)
           enddo
           Nk(i)=Nkout(i)
        enddo
        Gc(srtnh3) = Gcout(srtnh3)
        Gc(srtso4) = Gcout(srtso4)
        Gc(srtdma) = Gcout(srtdma)

        ygas(mgsvi)=Gc(srtso4)*1.0d+12/boxmass*28.9/gmw(srtso4)
        ygas(mgnh3)=Gc(srtnh3)*1.0d+12/boxmass*28.9/gmw(srtnh3)
        ygas(mgdma)=Gc(srtdma)*1.0d+12/boxmass*28.9/gmw(srtdma)

cJJ Not called here in master either - not necessary (?)
cJJ      call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm,11)

      endif  

c--------------------------------------------------------

      if (icond_soa.eq.1) then

        !Converting gas unit from ppt to kg
         Gc(srtsoa1)=boxmass*organic_ppt(1)*1.0d-12*gmw(srtsoa1)/28.9
         Gc(srtsoa2)=boxmass*organic_ppt(2)*1.0d-12*gmw(srtsoa2)/28.9
         Gc(srtsoa3)=boxmass*organic_ppt(3)*1.0d-12*gmw(srtsoa3)/28.9
         Gc(srtsoa4)=boxmass*organic_ppt(4)*1.0d-12*gmw(srtsoa4)/28.9
         Gc(srtsoa5)=boxmass*organic_ppt(5)*1.0d-12*gmw(srtsoa5)/28.9 !EXLVOCS

         call org_cond(Nk,Mk,Gc,Nkout,Mkout,Gcout,dt,xkDMAN
     &     ,ichm,jchm,kchm)

        !Converting Nk, Mk, and G                                                           !Only ammonia
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
         Gc(srtsoa5) = Gcout(srtsoa5) !EXLVOCS

         organic_ppt(1)=Gc(srtsoa1)*1.0d+12/boxmass*28.9/gmw(srtsoa1)
         organic_ppt(2)=Gc(srtsoa2)*1.0d+12/boxmass*28.9/gmw(srtsoa2)
         organic_ppt(3)=Gc(srtsoa3)*1.0d+12/boxmass*28.9/gmw(srtsoa3)
         organic_ppt(4)=Gc(srtsoa4)*1.0d+12/boxmass*28.9/gmw(srtsoa4)
         organic_ppt(5)=Gc(srtsoa5)*1.0d+12/boxmass*28.9/gmw(srtsoa5) !EXLVOCS

      end if

cJJ       call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm,10)
      call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm)

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

cJJ        call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm,5)
        call mnfix_PSSA(Nk,Mk,ichm,jchm,kchm)

      endif

c--------------------------------------------------------      

      time = time +deltat

      enddo    ! main loop

      !Return gas concentrations
      ygas(mgsvi)=Gc(srtso4)*1.0d+12/boxmass*28.9/gmw(srtso4)
      ygas(mgnh3)=Gc(srtnh3)*1.0d+12/boxmass*28.9/gmw(srtnh3)
      ygas(mgdma)=Gc(srtdma)*1.0d+12/boxmass*28.9/gmw(srtdma)

      do i=1,ibins
        Nki(i)=Nk(i)
        do j=1,icomp
          Mki(i,j)=Mk(i,j)
        enddo
      enddo

      h2so4 = ygas(mgsvi)
      nh3ppt = ygas(mgnh3)
      dmappt = ygas(mgdma)

      
      organic_ppt(1)=Gc(srtsoa1)*1.0d+12/boxmass*28.9/gmw(srtsoa1) 
      organic_ppt(2)=Gc(srtsoa2)*1.0d+12/boxmass*28.9/gmw(srtsoa2) 
      organic_ppt(3)=Gc(srtsoa3)*1.0d+12/boxmass*28.9/gmw(srtsoa3) 
      organic_ppt(4)=Gc(srtsoa4)*1.0d+12/boxmass*28.9/gmw(srtsoa4) 
      organic_ppt(5)=Gc(srtsoa5)*1.0d+12/boxmass*28.9/gmw(srtsoa5) !EXLVOCS
c==============================================================

      RETURN
      END
