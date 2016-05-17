
c     **************************************************
c     *  CAMx2dman                                     *
c     **************************************************

c     WRITTEN BY JaeGun Jung, March 2006
c                                                                   
c     This subroutine is a linker between PMCAMx and DMAN            
c     converting following variables.
c  
c     +-+-------------------------+-----------+----------------+
c     | | name                    | input     | output         |
c     +-+-------------------------+-----------+----------------+
c     |1| pressure                | pressure  | pres           |
c     |2| relative humdity        | rh        | relh           |
c     |3| H2SO4 conc.             | q(+ih2so4)| h2so4          | 
c     |4| NH3 conc.               | q(+inh3)  | nh3ppt         |
c     |5| number conc.            | q(+knum)  | Nk(ibins)      |
c     |6| mass conc.              | q(+k...)  | Mk(ibins,icomp)| 
c     +-+-------------------------+-----------+----------------+
c     

      subroutine CAMx2dman(q,t0,t1,tempK,pressure,dsulfdt,ich,jch,kch,fndt)

c
c-----Include files
c
      include 'dynamic.inc'
c
c-----Argument declarations
c
      double precision q(ntotal) !particles mass [=] ugr/m3, number [=]#/cm3, gas [=]ppm
      double precision t0, t1
      real tempK      ! temperature [=] K
      real pressure   ! atm
      real dsulfdt    !sulfuric acid production rate
      integer ich, jch, kch ! coordiate, x, y, z

c
c-----Variable declarations
c
      integer ibins, icomp
      parameter (ibins=41, icomp=10)


      integer srtso4, srtinrt, srtnh3, srtdma, srth2o !species indicators
      integer srtsoa1, srtsoa2, srtsoa3, srtsoa4 
      integer srtsoa5    !!EXLVOCS
c
      parameter (srtso4=1, srtinrt=2)
      parameter (srtsoa1=3, srtsoa2= 4, srtsoa3=5, srtsoa4=6)
      parameter (srtsoa5=7)  !!EXLVOCS
      parameter (srtnh3=8, srtdma=9, srth2o=10)
      integer i, j ! counter variables
      integer ii, jj ! counter variables
      integer iflag ! a flag for checking a negative mass
c
      double precision Nk(ibins)       ! Number concentrations
      double precision Mk(ibins,icomp) ! Mass concentrations
      double precision xk(ibins+1)     ! bins' boundary mass [=] kg
      double precision Neps, Meps
      double precision totmass, newmass! For updating Nk when Mk is less
                                       ! than zero.
      double precision eps
c
C-----------------------------------------------------------
      real organic_ppt(5) !organic [=] ppt
c-----------------------------------------------------------
 
      real h2so4      ! sulfuric acid gas [=] ppt
      real nh3ppt     ! ammonia gas [=] ppt
      double precision dmappt ! dma gas [=] ppt
      real relh       ! relative humidity
      real tstart     ! starting time and ending time of simulation
      real tend       ! starting time and ending time of simulation
      real boxvol     ! A volume of arbitrary box
      real pres       ! Pa 
      real rt_pom(ibins), rt_ec(ibins), rt_crst(ibins), rt_cl(ibins),
     &    rt_na(ibins), rt_no3(ibins) ! Ratios of each inert component


      real tot_soa1,tot_soa2   
      real rt_soa1(ibins),rt_soa2(ibins), rt_soa3(ibins),rt_soa4(ibins)
      real rt_soa5(ibins) !!EXLVOCS      

      real tot_inert ! total inert mass
      real tot_inert2 ! total inert mass after calling dman
c
      real cvt, cvt2
      double precision fndt(3) !Nucleation diagnostic
c
c-----Adjustable parameters
c
      parameter (Neps=1.0d-10, Meps=1.0d-20)
cdbg      parameter (eps = 1.0d-20)
      parameter (eps = 1.0d-10)
      parameter (boxvol = 3.0d20) ! cm3, arbitrary value
      parameter (cvt = 1.0d-15)   ! a convert factor from ug/m3 to kg/cm3
      parameter (cvt2 = 1.0d+15)  ! a covert factor from kg/cm3 to ug/m3
c
c-----Code--------------------------------------------------------------
c
      ! Set variables
      iflag = 0 
      ! Set xk's
      xk(1) = 3.75315d-25 ! kg
      do i = 1, ibins
         xk(i+1) = 2.0 * xk(i)
      enddo

      ! Adjust unit of timestep 
      tstart = t0/3600. ![=]hr
      tend = t1/3600.   ![=]hr

      ! Converting PMCAMx variables to DMAN variables
      pres = pressure * 1.01325d5 ! Pa
      relh = rh ! Change a relative humidity variable for DMAN

c-----------------------------------------------------------------------
c                                 ORGANIC
c-----------------------------------------------------------------------
        if (q(naer+icg1).ge.0.0) then
          organic_ppt(1) = q(naer+icg1) * 1.0d6   ! organic 1 [=] ppt, q [=] ppm
        else   
          if (q(naer+icg1).gt.(-eps*1.0d-6)) then
            organic_ppt(1) = eps
            q(naer+icg1) = eps * 1.0d-6
          else
            write(*,*)'organic_ppt1 is less than zero'
            write(*,*)'q(naer+icg1) [ppm]',q(naer+icg1)
            write(*,*)'Coordinate =', ich, jch, kch
            write(*,*)'organic_ppt1=',organic_ppt(1)
            STOP
          endif
        endif

c-----------------------------------------------------------------
        if (q(naer+icg2).ge.0.0) then
          organic_ppt(2) = q(naer+icg2) * 1.0d6   ! organic 2 [=] ppt, q [=] ppm
        else
          if (q(naer+icg2).gt.(-eps*1.0d-6)) then
            organic_ppt(2) = eps
            q(naer+icg2) = eps * 1.0d-6
          else
            write(*,*)'organic_ppt2 is less than zero'
            write(*,*)'q(naer+icg2) [ppm]',q(naer+icg2)
            write(*,*)'Coordinate =', ich, jch, kch
            write(*,*)'organic_ppt2=',organic_ppt(2)
            STOP
          endif
        endif

c------------------------------------------------------------------
        if (q(naer+icg3).ge.0.0) then
          organic_ppt(3) = q(naer+icg3) * 1.0d6   ! organic 3 [=] ppt, q [=] ppm
        else
          if (q(naer+icg3).gt.(-eps*1.0d-6)) then
            organic_ppt(3) = eps
            q(naer+icg3) = eps * 1.0d-6
          else
            write(*,*)'organic_ppt3 is less than zero'
            write(*,*)'q(naer+icg3) [ppm]',q(naer+icg3)
            write(*,*)'Coordinate =', ich, jch, kch
            write(*,*)'organic_ppt3=',organic_ppt(3)
            STOP
          endif
        endif

c------------------------------------------------------------------
        if (q(naer+icg4).ge.0.0) then
          organic_ppt(4) = q(naer+icg4) * 1.0d6   ! organic 4 [=] ppt, q [=] ppm
        else
          if (q(naer+icg4).gt.(-eps*1.0d-6)) then
            organic_ppt(4) = eps
            q(naer+icg4) = eps * 1.0d-6
          else
            write(*,*)'organic_ppt4 is less than zero'
            write(*,*)'q(naer+icg4) [ppm]',q(naer+icg4)
            write(*,*)'Coordinate =', ich, jch, kch
            write(*,*)'organic_ppt4=',organic_ppt(4)
            STOP
          endif
        endif

c---------------------------------------------------------------------
        !          !EXLVOCS   
        if (q(naer+icg5).ge.0.0) then
          organic_ppt(5) =q(naer+icg5) * 1.0d6
        else
          if (q(naer+icg5).gt.(-eps*1.0d-6)) then
            organic_ppt(5) = eps
            q(naer+icg5) = eps * 1.0d-6
          else
            write(*,*)'organic_ppt5 is less than zero'
            write(*,*)'q(naer+icg5) [ppm]',q(naer+icg5)
            write(*,*)'Coordinate =', ich, jch, kch
            write(*,*)'organic_ppt5=',organic_ppt(5)
            STOP
          endif
        endif
C========================================================================
      if (q(naer+ih2so4).ge.0.0) then 
        h2so4 = q(naer+ih2so4) * 1.0d6   ! h2so4 [=] ppt, q [=] ppm 
      else
        if (q(naer+ih2so4).gt.(-eps*1.0d-6)) then
          h2so4 = eps
          q(naer+ih2so4) = eps * 1.0d-6
        else
          write(*,*)'H2SO4 is less than zero'
          write(*,*)'q(naer+ih2so4) [ppm]',q(naer+ih2so4)
          write(*,*)'Coordinate =', ich, jch, kch
          write(*,*)'dsulfdt=',dsulfdt
          STOP
        endif
      endif

      !Ammonia
      if (q(naer+inh3).ge.0.d0) then
        nh3ppt = q(naer+inh3) * 1.0d6 
      else
        if (q(naer+inh3).gt.(-eps*1.0d-6)) then
          nh3ppt = eps
          q(naer+inh3) = eps * 1.0d-6
        else
          write(*,*)'NH3 is less than zero'
          write(*,*)'q(naer+inh3) [ppm]=',q(naer+inh3)
          write(*,*)'Coordinate =', ich, jch, kch
          write(*,*)'dsulfdt=',dsulfdt
          STOP
        endif
      endif

      !Dimethyl Amine
      if (q(naer+idma).ge.0.d0) then
        dmappt = q(naer+idma) * 1.0d6 
      else
        if (q(naer+idma).gt.(-eps*1.0d-6)) then
          dmappt = eps
          q(naer+idma) = eps * 1.0d-6
        else
          write(*,*)'DMA is less than zero'
          write(*,*)'q(naer+inh3) [ppm]=',q(naer+idma)
          write(*,*)'Coordinate =', ich, jch, kch
          write(*,*)'dsulfdt=',dsulfdt
          STOP
        endif
      endif
 
      !Load Aerosol Variables
      do i=1, ibins
         ! First check for negative tracers
         if (q((i-1)*nsp+knum).lt.0.0) then
            if (q((i-1)*nsp+knum).gt.(-Neps*(1./boxvol)))then
               q((i-1)*nsp+knum) = Neps*(1./boxvol) 
               do j=2, nsp-1 ! all mass species wo H2O
                  q((i-1)*nsp+j) = Neps*(1./boxvol)*1.4*xk(i)*(1./real(nsp-2))
               enddo
            else
               write(*,*)'Coordinate =', ich, jch, kch
               write(*,*)'Negative tracer in DMAN after dman'
               write(*,*)'sizesection=', i
               write(*,*)'q(+knum)='
               do ii=1,ibins
                  write(*,*)q((ii-1)*nsp+knum)
               enddo
               write(*,*)'q(+k...)='
               do jj=2, nsp-1 ! all mass species wo H2O
                  write(*,*)'species=',jj
                  do ii=1,ibins
                     write(*,*)q((ii-1)*nsp+jj)
                  enddo
               enddo
               STOP
            endif
         endif
c---------------------------------------------------------------------         
         totmass=0.0
         do j=2, nsp-1 ! all mass species wo H2O
            totmass=totmass+q((i-1)*nsp+j)
            if (q((i-1)*nsp+j).lt.0.0) then
               if (q((i-1)*nsp+j).gt.(-Meps*(1./(cvt*boxvol)))) then
                  q((i-1)*nsp+j) = Meps*(1./(cvt*boxvol))
                  iflag=1
               else
                  write(*,*)'Coordinate =', ich, jch, kch
                  write(*,*)'Negative tracer in DMAN after dman'
                  write(*,*)'sizesection=', i
                  write(*,*)'q(+knum)='
                  do ii=1,ibins
                     write(*,*)q((ii-1)*nsp+knum)
                  enddo
                  write(*,*)'q(+k...)='
                  do jj=2, nsp-1 ! all mass species wo H2O
                     write(*,*)'species=',jj
                     do ii=1,ibins
                        write(*,*)q((ii-1)*nsp+jj)
                     enddo
                  enddo
                  STOP
               endif
            endif
         enddo
c--------------------------------------------------------------------------         
         if (iflag.eq.1) then
            newmass=0.0
            do jj=2, nsp-1 ! all mass species wo H2O
               newmass=newmass+q((i-1)*nsp+jj)
            enddo
            q((i-1)*nsp+knum)=q((i-1)*nsp+knum)*newmass/totmass
            iflag=0
         endif
c-----------------------------------------------------------------------
         Nk(i)=q((i-1)*nsp+knum) * boxvol
          ! Nk [=] #, q [=] #/cm3, and boxvol [=] cm3

         Mk(i,srtso4)=q((i-1)*nsp+kso4) * cvt * boxvol
              ! Mk [=] kg, q [=] ug/m3
         tot_inert = (q((i-1)*nsp+kpom)+q((i-1)*nsp+kec)+
                     ! POA                 EC
     &                q((i-1)*nsp+kcrus)+q((i-1)*nsp+kcl)+
                     ! CRST                Cl
     &                q((i-1)*nsp+kna)+
                     ! Na                  
c     &                q((i-1)*nsp+kcl+1)+q((i-1)*nsp+kcl+2)+    !david
c     &                q((i-1)*nsp+kcl+3)+q((i-1)*nsp+kcl+4)+    !david
c                     ! SOA                                      !david
     &                q((i-1)*nsp+kno3))
                     ! Nitrate

         ! Capture ratios before calling dman 
         rt_pom(i) = q((i-1)*nsp+kpom) * (1.0/tot_inert)
         rt_ec(i) = q((i-1)*nsp+kec) * (1.0/tot_inert)
         rt_crst(i) = q((i-1)*nsp+kcrus) * (1.0/tot_inert)
         rt_cl(i) = q((i-1)*nsp+kcl) * (1.0/tot_inert)
         rt_na(i) = q((i-1)*nsp+kna) * (1.0/tot_inert)
         rt_no3(i) = q((i-1)*nsp+kno3) * (1.0/tot_inert)   
c
         Mk(i,srtinrt) = tot_inert * cvt * boxvol

         Mk(i,srtsoa1) = q((i-1)*nsp+kcl+1) * cvt * boxvol
         Mk(i,srtsoa2) = q((i-1)*nsp+kcl+2) * cvt * boxvol
         Mk(i,srtsoa3) = q((i-1)*nsp+kcl+3) * cvt * boxvol
         Mk(i,srtsoa4) = q((i-1)*nsp+kcl+4) * cvt * boxvol
         Mk(i,srtsoa5) = q((i-1)*nsp+kcl+5) * cvt * boxvol !!EXLVOCS

         Mk(i,srtnh3)=q((i-1)*nsp+knh4) * cvt * boxvol
         Mk(i,srtdma)=q((i-1)*nsp+kpami)* cvt * boxvol     !amine /JJ
         Mk(i,srth2o)=q((i-1)*nsp+kh2o) * cvt * boxvol
      enddo



c================ SUBROUTINE DMAN ======================================
   
      call dman(tstart,tend,Nk,Mk,h2so4,nh3ppt,dmappt,relh,tempK,pres,
     &     dsulfdt,organic_ppt,ich,jch,kch,fndt)

C======================================================================


c-----Return DMAN values to the PMCAMx variable, check I can call initbounds
      do i=1, ibins
         ! First check for negative tracers
         if (Nk(i).lt.0.0) then
            if (Nk(i).gt.-Neps)then
               Nk(i) = 0.0
               do j=1, icomp-1
                  Mk(i,j)=Neps*1.4*xk(i)*(1./(icomp-1))
               enddo
            else
               write(*,*)'Coordinate =', ich, jch, kch
               write(*,*)'Negative tracer in DMAN after dman'
               write(*,*)'sizesection=', i
               write(*,*)'Nk='
               do ii=1,ibins
                 write(*,*)Nk(ii)
               enddo
               write(*,*)'Mk='
               do jj=1,icomp
                  write(*,*)'species=',jj
                  do ii=1,ibins
                     write(*,*)Mk(ii,jj)
                  enddo
               enddo
               STOP
            endif
         endif
         totmass=0.0
         do j=1, icomp-1
            totmass=totmass+Mk(i,j)
            if (Mk(i,j).lt.0.0) then
               if (Mk(i,j).gt.-Meps) then
                  Mk(i,j)=0.0
                  iflag=1
               else
                  write(*,*)'Coordinate =', ich, jch, kch
                  write(*,*)'Negative tracer in DMAN after dman'
                  write(*,*)'sizesection=', i
                  write(*,*)'Nk='
                  do ii=1,ibins
                    write(*,*)Nk(ii)
                  enddo
                  write(*,*)'Mk='
                  do jj=1,icomp
                     write(*,*)'species=',jj
                     do ii=1,ibins
                        write(*,*)Mk(ii,jj)
                     enddo
                  enddo
                  STOP
               endif
            endif
         enddo
         if (iflag.eq.1) then
            newmass=0.0
            do jj=1, icomp-1
               newmass=newmass+Mk(i,jj)
            enddo
            Nk(i)=Nk(i)*newmass/totmass
            iflag=0
         endif

        ! Return the values
        q((i-1)*nsp+knum) = Nk(i) * (1.0/boxvol)
        ! Nk [=] #, q [=] #/cm3, and boxvol [=] cm3
        q((i-1)*nsp+kso4) = Mk(i,srtso4) * cvt2 * (1.0/boxvol)

        q((i-1)*nsp+knh4) = Mk(i,srtnh3) * cvt2 * (1.0/boxvol)
        q((i-1)*nsp+kpami)= Mk(i,srtdma) * cvt2 * (1.0/boxvol)
        q((i-1)*nsp+kh2o) = Mk(i,srth2o) * cvt2 * (1.0/boxvol)

        ! Mk [=] kg, q [=] ug/m3
        q((i-1)*nsp+kcl+1)=Mk(i,srtsoa1) * cvt2 * (1.0/boxvol)
        q((i-1)*nsp+kcl+2)=Mk(i,srtsoa2) * cvt2 * (1.0/boxvol)
        q((i-1)*nsp+kcl+3)=Mk(i,srtsoa3) * cvt2 * (1.0/boxvol)
        q((i-1)*nsp+kcl+4)=Mk(i,srtsoa4) * cvt2 * (1.0/boxvol)
        q((i-1)*nsp+kcl+5)=Mk(i,srtsoa5) * cvt2 * (1.0/boxvol)    !david 
c 
c     Only POA has the sum of POA, EC, CRST, Cl, and Na. The rest of
c     species are set to zero.
c
        tot_inert2 = Mk(i,srtinrt) * cvt2 * (1.0/boxvol)

        q((i-1)*nsp+kpom) = tot_inert2 * rt_pom(i)
        q((i-1)*nsp+kec) = tot_inert2 * rt_ec(i)
        q((i-1)*nsp+kcrus) = tot_inert2 * rt_crst(i)
        q((i-1)*nsp+kcl) = tot_inert2 * rt_cl(i)
        q((i-1)*nsp+kna) = tot_inert2 * rt_na(i)
        q((i-1)*nsp+kno3) = tot_inert2 * rt_no3(i)
c
      enddo      

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c--------------           ORGANIC 1       -------------------
      if (organic_ppt(1).ge.0.0) then
        q(naer+icg1) = organic_ppt(1) * 1.0d-6   ! organic 1 [=] ppt, q [=] ppm
      else
        if (organic_ppt(1).gt.-eps) then
          q(naer+icg1) = eps * 1.0d-6  ! organic 1 [=] ppt, q [=] ppm
          organic_ppt(1) = eps
        else
          write(*,*)'organic_ppt(1) is less than zero'
          write(*,*)'organic_ppt(1) [=]ppt',organic_ppt(1)
          write(*,*)'Coordinate =', ich, jch, kch
          STOP
        endif
      endif

C---------------------    ORGANIC 2   -------------------------  
      if (organic_ppt(2).ge.0.0) then
        q(naer+icg2) = organic_ppt(2) * 1.0d-6   ! organic 2 [=] ppt, q [=] ppm
      else
        if (organic_ppt(2).gt.-eps) then
          q(naer+icg2) = eps * 1.0d-6  ! organic 2 [=] ppt, q [=] ppm
          organic_ppt(2) = eps
        else
          write(*,*)'organic_ppt(2) is less than zero'
          write(*,*)'organic_ppt(2) [=]ppt',organic_ppt(2)
          write(*,*)'Coordinate =', ich, jch, kch
          STOP
        endif
      endif

C---------------------    ORGANIC 3   -------------------------
      if (organic_ppt(3).ge.0.0) then
        q(naer+icg3) = organic_ppt(3) * 1.0d-6   ! organic 3 [=] ppt, q [=] ppm
      else
        if (organic_ppt(3).gt.-eps) then
          q(naer+icg3) = eps * 1.0d-6  ! organic 3 [=] ppt, q [=] ppm
          organic_ppt(3) = eps
        else
          write(*,*)'organic_ppt(3) is less than zero'
          write(*,*)'organic_ppt(3) [=]ppt',organic_ppt(3)
          write(*,*)'Coordinate =', ich, jch, kch
          STOP
        endif
      endif

C---------------------    ORGANIC 4   -------------------------
      if (organic_ppt(4).ge.0.0) then
        q(naer+icg4) = organic_ppt(4) * 1.0d-6   ! organic 4 [=] ppt, q [=] ppm
      else
        if (organic_ppt(4).gt.-eps) then
          q(naer+icg4) = eps * 1.0d-6  ! organic 4 [=] ppt, q [=] ppm
          organic_ppt(4) = eps
        else
          write(*,*)'organic_ppt(4) is less than zero'
          write(*,*)'organic_ppt(4) [=]ppt',organic_ppt(4)
          write(*,*)'Coordinate =', ich, jch, kch
          STOP
        endif
      endif

C---------------------    ORGANIC 5 = EXLVOCS   -------------------------
      if (organic_ppt(5).ge.0.0) then
        q(naer+icg5) = organic_ppt(5) * 1.0d-6               
      else
        if (organic_ppt(5).gt.-eps) then
          q(naer+icg5) = eps * 1.0d-6
          organic_ppt(5) = eps
        else
          write(*,*)'organic_ppt(5) is less than zero'
          write(*,*)'organic_ppt(5)[=]ppt',organic_ppt(5)
          write(*,*)'Coordinate =', ich, jch, kch
          STOP
        endif
      endif          

C---------------------   H2SO4   -------------------------
      if (h2so4.ge.0.0) then
        q(naer+ih2so4) = h2so4 * 1.0d-6   ! h2so4 [=] ppt, q [=] ppm 
      else
        if (h2so4.gt.-eps) then
          q(naer+ih2so4) = eps * 1.0d-6  ! h2so4 [=] ppt, q [=] ppm 
          h2so4 = eps
        else
          write(*,*)'H2SO4 is less than zero'
          write(*,*)'h2so4 [=]ppt',h2so4
          write(*,*)'Coordinate =', ich, jch, kch
          write(*,*)'dsulfdt=',dsulfdt
          STOP
        endif

      endif

C---------------------   NH3     -------------------------     
      if (nh3ppt.ge.0.0) then
        q(naer+inh3) = nh3ppt * 1.0d-6 
      else
        if (nh3ppt.gt.-eps) then
          q(naer+inh3) = eps * 1.0d-6
          nh3ppt = eps
        else
          write(*,*)'NH3 is less than zero'
          write(*,*)'nh3ppt=',nh3ppt
          write(*,*)'Coordinate =', ich, jch, kch
          write(*,*)'dsulfdt=',dsulfdt
          STOP
        endif
      endif

      !Dimethyl Amine
      if (dmappt.ge.0.0) then
        q(naer+idma) = dmappt * 1.0d-6 
      else
        if (dmappt.gt.-eps) then
          q(naer+idma) = eps * 1.0d-6
          dmappt = eps
        else
          write(*,*)'DMA is less than zero'
          write(*,*)'dmappt=',dmappt
          write(*,*)'Coordinate =', ich, jch, kch
          write(*,*)'dsulfdt=',dsulfdt
          STOP
        endif
      endif

      RETURN
      END
