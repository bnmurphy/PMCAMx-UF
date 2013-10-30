
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

      subroutine CAMx2dman(q,t0,t1,tempK,pressure,dsulfdt,ich,jch,kch)
c
c-----Include files
c
      include 'dynamic.inc'
cdbg      include 'aervaria.inc'
c
c-----Argument declarations
c
      real*8 q(ntotal) !particles mass [=] ugr/m3, number [=]#/cm3, gas [=]ppm
      real*8 t0, t1
      real tempK      ! temperature [=] K
      real pressure   ! atm
      real dsulfdt    !sulfuric acid production rate
      integer ich, jch, kch ! coordiate, x, y, z
c
c-----Variable declarations
c
      integer ibins, icomp
      parameter (ibins=41, icomp=4)

      integer srtso4, srtorg, srtnh3, srth2o !species indicators
      parameter (srtso4=1, srtorg=2, srtnh3=3, srth2o=4)

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
      real h2so4      ! sulfuric acid gas [=] ppt
      real nh3ppt     ! ammonia gas [=] ppt
      real relh       ! relative humidity
      real tstart     ! starting time and ending time of simulation
      real tend       ! starting time and ending time of simulation
      real boxvol     ! A volume of arbitrary box
      real pres       ! Pa 
      real rt_pom(ibins), rt_ec(ibins), rt_crst(ibins), rt_cl(ibins),
     &    rt_na(ibins), rt_soa1(ibins), rt_soa2(ibins), rt_soa3(ibins),
     &    rt_soa4(ibins), rt_no3(ibins) ! Ratios of each inert component
      real tot_inert ! total inert mass
      real tot_inert2 ! total inert mass after calling dman
cdbg      real eps
      real cvt, cvt2
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
c
cdbg      write(*,*)'CAMx2dman.f - chkpt 1. at the very beginning'
c     added by LA
c      write(*,*)
c      write(*,*)'naer=',naer
c      write(*,*)'ih2so4=',ih2so4
c      write(*,*)'inh3=',inh3
c      write(*,*)'knum=',knum
c     end added by LA
      if (q(naer+ih2so4).ge.0.0) then 
        h2so4 = q(naer+ih2so4) * 1.0d6   ! h2so4 [=] ppt, q [=] ppm 
      else
        if (q(naer+ih2so4).gt.(-eps*1.0d-6)) then
cdbg          h2so4 = 0.0
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
c     added by LA
c      write(*,*)
c      write(*,*)'h2so4=',h2so4
c     end added by LA

      if (q(naer+inh3).ge.0.d0) then
        nh3ppt = q(naer+inh3) * 1.0d6 
      else
        if (q(naer+inh3).gt.(-eps*1.0d-6)) then
cdbg          nh3ppt = 0.d0
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
c     added by LA
c      write(*,*)
c      write(*,*)'nh3ppt=',nh3ppt
c     end added by LA
     
      do i=1, ibins
         ! First check for negative tracers
         if (q((i-1)*nsp+knum).lt.0.0) then
            if (q((i-1)*nsp+knum).gt.(-Neps*(1./boxvol)))then
               q((i-1)*nsp+knum) = Neps*(1./boxvol)
               do j=2, 13 ! from KNa to KEC, all mass species wo H2O 
                  q((i-1)*nsp+j) = Neps*(1./boxvol)*1.4*xk(i)*(1./12.)
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
               do jj=2, 13 ! from KNa to KEC, all mass species wo H2O 
                  write(*,*)'species=',jj
                  do ii=1,ibins
                     write(*,*)q((ii-1)*nsp+jj)
                  enddo
               enddo
               STOP
            endif
         endif
         totmass=0.0
         do j=2, 13 ! from KNa to KEC, all mass species wo H2O 
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
                  do jj=2, 13 ! from KNa to KEC, all mass species wo H2O 
                     write(*,*)'species=',jj
                     do ii=1,ibins
                        write(*,*)q((ii-1)*nsp+jj)
                     enddo
                  enddo
                  STOP
               endif
            endif
         enddo
c     added by LA
c         write(*,*)
c         write(*,*)'iflag=',iflag
c     end added by LA
         if (iflag.eq.1) then
            newmass=0.0
            do jj=2, 13 ! from KNa to KEC, all mass species wo H2O 
               newmass=newmass+q((i-1)*nsp+jj)
            enddo
            q((i-1)*nsp+knum)=q((i-1)*nsp+knum)*newmass/totmass
            iflag=0
         endif
c
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
     &                q((i-1)*nsp+kcl+1)+q((i-1)*nsp+kcl+2)+
     &                q((i-1)*nsp+kcl+3)+q((i-1)*nsp+kcl+4)+
                     ! SOA
     &                q((i-1)*nsp+kno3))
                     ! Nitrate
         ! Capture ratios before calling dman 
         rt_pom(i) = q((i-1)*nsp+kpom) * (1.0/tot_inert)
         rt_ec(i) = q((i-1)*nsp+kec) * (1.0/tot_inert)
         rt_crst(i) = q((i-1)*nsp+kcrus) * (1.0/tot_inert)
         rt_cl(i) = q((i-1)*nsp+kcl) * (1.0/tot_inert)
         rt_na(i) = q((i-1)*nsp+kna) * (1.0/tot_inert)
         rt_soa1(i) = q((i-1)*nsp+kcl+1) * (1.0/tot_inert)
         rt_soa2(i) = q((i-1)*nsp+kcl+2) * (1.0/tot_inert)
         rt_soa3(i) = q((i-1)*nsp+kcl+3) * (1.0/tot_inert)
         rt_soa4(i) = q((i-1)*nsp+kcl+4) * (1.0/tot_inert)
         rt_no3(i) = q((i-1)*nsp+kno3) * (1.0/tot_inert)   
c
         Mk(i,srtorg) = tot_inert * cvt * boxvol
         Mk(i,srtnh3)=q((i-1)*nsp+knh4) * cvt * boxvol
         Mk(i,srth2o)=q((i-1)*nsp+kh2o) * cvt * boxvol
      enddo      
c
cdbg      write(*,*)'CAMx2dman.f - chkpt 2. after converting before dman'
c
      !For a debuging purpose
cdbg      if ((tstart.gt.0.0).and.(tstart.lt.0.5)) then
cdbg        if ((ich.eq.36).and.(jch.eq.29).and.(kch.eq.1)) then
cdbg          write(*,*)'In CAMx2dman before calling dman' 
cdbg          write(*,*)'coordinate of (i,j,k)',ich, jch, kch
cdbg          write(*,*)'tempK,pressure,dsulfdt=',tempK,pressure,dsulfdt
cdbg          write(*,*)'h2so4=',h2so4,'nh3ppt=',nh3ppt
cdbg          write(*,*)'Nk='
cdbg          do i=1, ibins
cdbg            write(*,*)Nk(i)
cdbg          enddo
cdbg          write(*,*)'Mk='
cdbg          do j=1, icomp
cdbg            write(*,*)'j=',j
cdbg            do i=1, ibins
cdbg              write(*,*)Mk(i,j)
cdbg            enddo
cdbg          enddo
cdbg        endif
cdbg      endif
c     added by LA
c      write(*,*)
c      write(*,*)'nh3ppt before dman=',nh3ppt
c      write(*,*)
c     end added by LA

      call dman(tstart,tend,Nk,Mk,h2so4,nh3ppt,relh,tempK,pres,dsulfdt
     & ,ich,jch,kch)

c     added by LA
c      write(*,*)
c      write(*,*)'nh3ppt after dman=',nh3ppt
c      write(*,*)
c     end added by LA
      !For a debuging purpose
cdbg      if ((tstart.gt.0.0).and.(tstart.lt.0.5)) then
cdbg        if ((ich.eq.36).and.(jch.eq.29).and.(kch.eq.1)) then
cdbg          write(*,*)'In CAMx2dman after calling dman' 
cdbg          write(*,*)'coordinate of (i,j,k)',ich, jch, kch
cdbg          write(*,*)'tempK,pressure,dsulfdt=',tempK,pressure,dsulfdt
cdbg          write(*,*)'h2so4=',h2so4,'nh3ppt=',nh3ppt
cdbg          write(*,*)'Nk='
cdbg          do i=1, ibins
cdbg            write(*,*)Nk(i)
cdbg          enddo
cdbg          write(*,*)'Mk='
cdbg          do j=1, icomp
cdbg            write(*,*)'j=',j
cdbg            do i=1, ibins
cdbg              write(*,*)Mk(i,j)
cdbg            enddo
cdbg          enddo
cdbg        endif
cdbg      endif
c
cdbg      write(*,*)'CAMx2dman.f - chkpt 3. after dman'
c
c-----Return DMAN values to the PMCAMx variable, check I can call initbounds
c
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
        ! Mk [=] kg, q [=] ug/m3
c
c     Only POA has the sum of POA, EC, CRST, Cl, and Na. The rest of
c     species are set to zero.
c
        tot_inert2 = Mk(i,srtorg) * cvt2 * (1.0/boxvol)
        q((i-1)*nsp+kpom) = tot_inert2 * rt_pom(i)
        q((i-1)*nsp+kec) = tot_inert2 * rt_ec(i)
        q((i-1)*nsp+kcrus) = tot_inert2 * rt_crst(i)
        q((i-1)*nsp+kcl) = tot_inert2 * rt_cl(i)
        q((i-1)*nsp+kna) = tot_inert2 * rt_na(i)
        q((i-1)*nsp+kcl+1) = tot_inert2 * rt_soa1(i)
        q((i-1)*nsp+kcl+2) = tot_inert2 * rt_soa2(i)
        q((i-1)*nsp+kcl+3) = tot_inert2 * rt_soa3(i)
        q((i-1)*nsp+kcl+4) = tot_inert2 * rt_soa4(i)
        q((i-1)*nsp+kno3) = tot_inert2 * rt_no3(i)
c
        q((i-1)*nsp+knh4) = Mk(i,srtnh3) * cvt2 * (1.0/boxvol)
        q((i-1)*nsp+kh2o) = Mk(i,srth2o) * cvt2 * (1.0/boxvol)
      enddo      

      if (h2so4.ge.0.0) then
        q(naer+ih2so4) = h2so4 * 1.0d-6   ! h2so4 [=] ppt, q [=] ppm 
      else
        if (h2so4.gt.-eps) then
cdbg          q(naer+ih2so4) = 0.0   ! h2so4 [=] ppt, q [=] ppm 
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

      if (nh3ppt.ge.0.0) then
        q(naer+inh3) = nh3ppt * 1.0d-6 
      else
        if (nh3ppt.gt.-eps) then
cdbg          q(naer+inh3) = 0.0
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

      !For a debuging purpose
cdbg      if ((tstart.gt.0.0).and.(tstart.lt.0.5)) then
cdbg      if ((ich.eq.36).and.(jch.eq.29).and.(kch.eq.1)) then
cdbg         write(*,*)'In CAMx2dman after converting Nk and Mk to q' 
cdbg         write(*,*)'coordinate of (i,j,k)',ich, jch, kch
cdbg         write(*,*)'tempK,pressure,dsulfdt=',tempK,pressure,dsulfdt
cdbg         write(*,*)'q(naer+ih2so4)',q(naer+ih2so4)
cdbg         write(*,*)'q(naer+inh3)',q(naer+inh3)
cdbg         i = 4 !4th size section
cdbg         write(*,*)'q(knum)=',q((i-1)*nsp+knum)   !number
cdbg         write(*,*)'q(kso4)=',q((i-1)*nsp+kso4)   !1
cdbg         write(*,*)'q(kpom)=',q((i-1)*nsp+kpom)   !2
cdbg         write(*,*)'q(kec)=',q((i-1)*nsp+kec)     !3
cdbg         write(*,*)'q(kcrus)=',q((i-1)*nsp+kcrus) !4
cdbg         write(*,*)'q(kcl)=',q((i-1)*nsp+kcl)     !5
cdbg         write(*,*)'q(kna)=',q((i-1)*nsp+kna)     !6
cdbg         write(*,*)'q(ksoa1)=',q((i-1)*nsp+kcl+1) !7
cdbg         write(*,*)'q(ksoa2)=',q((i-1)*nsp+kcl+2) !8
cdbg         write(*,*)'q(ksoa3)=',q((i-1)*nsp+kcl+3) !9
cdbg         write(*,*)'q(ksoa4)=',q((i-1)*nsp+kcl+4) !10
cdbg         write(*,*)'q(kno3)=',q((i-1)*nsp+kno3)   !11
cdbg         write(*,*)'q(knh4)=',q((i-1)*nsp+knh4)   !12
cdbg         write(*,*)'q(kh2o)=',q((i-1)*nsp+kh2o)   !13
cdbg       endif
cdbg       endif

cdbg      write(*,*)'CAMx2dman.f - chkpt 4. after converting after dman'



      RETURN
      END
