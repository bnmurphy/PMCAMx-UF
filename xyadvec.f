      subroutine xyadvec(igrd,xyordr,ncol,nrow,nlay,nspc,nsen,nadv,
     &                   deltat,dx,dy,windu,windv,depth,mapscl,conc,
cjgj     &                   fluxes,sens,tarray2,isaptr,ipa_cel)
     &                   fluxes,sens,tarray2,isaptr,ipa_cel,species)
c
c-----CAMx v4.02 030709
c
c     XYADVEC drives 2-D advection of concentrations.  This version also
c     performs the 2-D advection on DDM sensitivies, if DDM is implemented
c                          
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        4/23/99   Coarse grid outflow BC's set equal to concentration in
c                  the first inner computational cells
c        4/26/99   Added Piecewise Parabolic Method for horizontal advection
c       10/30/01   Revised map scale factor application to OSAT fluxes to be
c                  more consistent with how the fluxes are used.
c       12/07/01   added instructions for OMP and rearranged some loops
c                  to facilitate parallel processing
c       01/22/02   now only calls OSAT routines if not doing RTRAC
c       01/30/02   Added code for RTRAC probing tool
c        4/10/03   X/Y advection now uses layer-dependent timestep
c
c     Input arguments:
c        igrd              grid index
c        xyordr            order of x & y advection
c        ncol              number of columns
c        nrow              number of rows
c        nlay              number of layers
c        nspc              number of species
c        nsen              number of species times number of DDM parameters
c        nadv              number of sub-steps per timestep
c        deltat            time step (s)
c        dx                cell size in x-direction (m)
c        dy                cell size in y-direction (m)
c        windu             wind speed in x-direction (m/s)
c        windv             wind speed in y-direction (m/s)
c        depth             layer depth (m)
c        mapscl            map-scale factor at cell centroids
c        conc              species concentrations (umol/m3)
c        sens              sensitivity coefficient (umol/m3/parameter unit)
c        tarray2           CPU timing arguments (s)
c        isaptr            pointer into tracer conc array for this grid
c        ipa_cel           gridded array to identify if cell is
c                          in a IPRM sub-domain
c
c     Output arguments:
c        conc              species concentrations (umol/m3)
c        fluxes            fluxes across the boundaries (umol)
c        sens              sensitivity coefficient (umol/m3/parameter unit)
c
c     Routines Called:
c        HADVBOT
c        HADVPPM
c        BOTTDDM
c
c     Called by:
c        EMISTRNS
c
      include "camx.prm"
      include "bndary.com"
      include "chmstry.com"
      include "filunit.com"
      include "flags.com"
c
c======================== Source Apportion Begin =======================
c
      include "tracer.com"
c
c======================== Source Apportion End =======================
c
c
c======================== Process Analysis Begin ====================================
c
      include "procan.com"
      integer ipa_cel(ncol,nrow,nlay)
c
c========================= Process Analysis End =====================================
c
      integer xyordr
      integer nadv(nlay)
      dimension conc(ncol,nrow,nlay,nspc)
      dimension sens(ncol,nrow,nlay,nsen)
      real windu(ncol,nrow,nlay),windv(ncol,nrow,nlay),
     &     depth(ncol,nrow,nlay),
     &     mapscl(ncol,nrow),
     &     dx(nrow)
      real c1d(MX1D),v1d(MX1D),flxarr(MX1D),m1d(MX1D),saflux(MX1D)
      real fluxo3(MX1D), fluxvoc(MX1D), fluxnox(MX1D)
      real cnco3(MX1D), cncvoc(MX1D), cncnox(MX1D)
      real c1d0(MX1D),fpc(MX1D),fmc(MX1D)
      real sen1d(MX1D,MXTRSP)
      real*8 fluxes(nspc,11),flux1,flux2
      real*8 fluxtmp(MXSPEC,8,MXLAYA)
      real c2d(MX1D),mn(MX1D) ! by jgj 8/13/06
      dimension tarray2(2)
      character*10 species(nspc) ! by jgj 7/18/06
c
      integer ii,jj,kk ! for generating the order of an advection by jgj 8/8/06
      integer order(nspc) ! the order of species by jgj 8/8/06
      integer nsect ! number of size section by jgj 8/8/06
      integer spec ! variable for calculating species by jgj 8/14/06
      integer lmod
      real nflxarr(MX1D) ! flxarr for number species by jgj 8/14/06
c
c======================== Process Analysis Begin ====================================
c
      real fc1(MX1D), fc2(MX1D)
c
c========================= Process Analysis End =====================================
c
c     For setting the order of advection by jgj 8/8/06.
c     The order(kk) has number of gas concentration first.
c     Then it saves number concentration at the first size bin.
c     Thereafter it saves mass concentrations at the size bin.
c     The order(kk) saves the order of species in the same way 
c     as increasing the order of size bins.
c
      nsect = 43
      do kk=1,ngas
         order(kk) = kk
      enddo

      kk=ngas

      do ii=1,nsect
         do jj=1,naero
            kk=kk+1
cdbg            lmod = ngas + (jj-1)*nsect + ii
cdbg            write(*,*)'spname(lmod)=',spname(lmod)
            if (jj.eq.1) then
               order(kk) = ngas + (naero-1)*nsect + ii
            else
               order(kk) = ngas + (jj-2)*nsect + ii 
            endif
         enddo
cdbg         pause
      enddo
c
c-----Entry point
c
c  ---- intialize the local flux aray to zero ---
c
      do i=1,nspc
        do j=1,8
          do k=1,nlay
            fluxtmp(i,j,k) = 0.
          enddo
        enddo
      enddo
c
      if( xyordr .eq. 0 ) goto 200
c
c-----Advection in x-direction
c
 100  continue
c
      write(*,'(a20,$)') 'x advection ......'
      write(iout,'(a20,$)') 'x advection ......'
      tcpu = dtime(tarray2)
c
cdbg      call checkconc(conc,ncol,nrow,nlay,nspc,species,2)
c
      do 20 k = 1,nlay
        dtuse = deltat/nadv(k)
        do 21 istep = 1,nadv(k)
c
        do 10 j = 2,nrow-1
          i1 = 1
          i2 = ncol
          if (igrd.eq.1) then
            if (ibeg(j).eq.-999) goto 10
            i1 = ibeg(j) - 1
            i2 = iend(j) + 1
          endif
          l = 0
          do i = i1,i2-1
            l = l + 1
            v1d(l) = 2.*windu(i,j,k)/(mapscl(i+1,j) + mapscl(i,j))
            m1d(l) = mapscl(i,j)*mapscl(i,j)
          enddo
          l = l + 1
          v1d(l) = windu(i2,j,k)
          m1d(l) = mapscl(i2,j)*mapscl(i2,j)
          nn = i2 - i1 + 1
c
c======================== Source Apportion Begin =======================
c
c
c  --- initialize the tracer fluxes/concs to zero ----
c
           do i=1,MX1D
             fluxo3(i) = 0.
             fluxvoc(i) = 0.
             fluxnox(i) = 0.
             cnco3(i) = 0.
             cncvoc(i) = 0.
             cncnox(i) = 0.
           enddo
c
c========================= Source Apportion End ========================
c
c      
c          do 22 ispc = 1,nspc
c
          do 22 ii = 1,nspc
            ispc = order(ii)
c
            l = 0
            do i = i1,i2
              l = l + 1
              c1d(l) = conc(i,j,k,ispc)*dy*depth(i,j,k)
            enddo
            if (igrd.eq.1 .and. v1d(1).lt.0.) c1d(1) = c1d(2)
            if (igrd.eq.1 .and. v1d(nn-1).gt.0.) c1d(nn) = c1d(nn-1)
c
c======================== DDM Begin =======================
c
            if( lddm ) then
               do iddm=1,nddmsp
                 l = 0
                 do i = i1,i2
                   l = l + 1
                   c1d0(l) = c1d(l)
                   ioff = iptddm(ispc)+iddm-1
                   sen1d(l,iddm) = sens(i,j,k,ioff)*dy*depth(i,j,k)
                 enddo
                 if (igrd.eq.1 .and. v1d(1).lt.0.) 
     &                                  sen1d(1,iddm) = sen1d(2,iddm)
                 if (igrd.eq.1 .and. v1d(nn-1).gt.0.) 
     &                              sen1d(nn,iddm) = sen1d(nn-1,iddm)
               enddo
            endif
c
c======================== DDM End =======================
c
cdbg            if ((ispc.eq.230).and.(j.eq.2)) then
cdbg               write(*,*)'iadvct=',iadvct
cdbg               write(*,*)'Before call hadvbot'
cdbg	       write(*,*)'c1d(2)=',c1d(2)
cdbg               pause
cdbg            endif
c
cdbg            call mnratio(ii,ngas,naero,c1d,c2d,i1,i2,MX1D,mn,j,k,1) ! by jgj 8/13/06
            call mnratio(ii,c1d,c2d,i1,i2,mn,order,nspc,j,k,1) ! by jgj 8/13/06
c
c           In case of mass species, use the flxarr of the number species.
c           by jgj 8/13/06
c
c           The logic below is little bit hard to follow. The idea is as follows.
c           When ii hits the first number concentration, number concentration flux,
c           nflxarr set to equals to zero. Then the flxarr, the same flux but different variable
c           name, is updated by advection module. After that, nflxarr is updated by the flxarr 
c           updated.
c
c           When ii hits the first mass species concentration, hadvppm2 allocates the flxarr of
c           the mass concentration. Then the flow of program goes to line 400. For y-advection, 
c           the same rule has been applied.
c
            if (ii.gt.ngas) then
               spec=mod(ii-ngas,naero)
               if (spec.eq.1) then ! An number concentration
                  do jj = 1,nn-1
                     nflxarr(jj) = 0
                  enddo
               else
                  call hadvppm2(ii,ngas,naero,MX1D,mn,nflxarr,
     &             nn,dtuse,dx(j),c1d,m1d,flxarr,flux1,flux2,
     &             saflux,fc1,fc2)
                  goto 400
               endif
            endif
c
            if (iadvct.eq.2) then
              call hadvbot(nn,dtuse,dx(j),c1d,v1d,m1d,flxarr,flux1,
     &                                   flux2,saflux,fpc,fmc,fc1,fc2)
            elseif( iadvct .eq. 3) then
              call hadvppm(nn,dtuse,dx(j),c1d,v1d,m1d,flxarr,flux1,
     &                                           flux2,saflux,fc1,fc2)
            endif
c
c           Save flxarr of number species by jgj 8/14/06
c
            if (ii.gt.ngas) then
               spec=mod(ii-ngas,naero)
               if (spec.eq.1) then
                  do jj = 1,nn-1
                     nflxarr(jj) = flxarr(jj)
                  enddo
               endif
            endif
c
 400        continue
c
cdbg            if ((ispc.eq.230).and.(j.eq.2)) then
cdbg               write(*,*)'After call hadvbot'
cdbg	       write(*,*)'c1d(2)=',c1d(2)
cdbg            endif
c======================== DDM Begin =======================
c
c
            if( lddm .AND. iadvct .eq. 2 ) then
              call bottddm(nn,dtuse,dx(j),sen1d,nddmsp,c1d0,c1d,
     &                                            fpc,fmc,v1d,m1d)
            endif
c
c======================== DDM End =======================
c
c
c======================== Process Analysis Begin ====================================
c
            if( lipr ) then
c
c-----Change from X-direction horizontal advection
c 
                l = 1
                do i=i1+1,i2-1
                   l = l+1
c
                   if( ipa_cel(i,j,k) .GT. 0 ) then
                      ipa_idx = ipa_cel(i,j,k)
c
c-----Flux at west boundary
c
                      cipr(IPR_WADV, ipa_idx, ispc) =
     &                   cipr(IPR_WADV, ipa_idx, ispc) + 
     &                                          fc1(l)/dy/depth(i,j,k)
c
c-----Flux at east boundary
c
                      cipr(IPR_EADV, ipa_idx, ispc) =
     &                   cipr(IPR_EADV, ipa_idx, ispc) +
     &                                         fc2(l)/dy/depth(i,j,k)
c
c-----Average volume
                      cipr(IPR_VOL, ipa_idx, ispc) =
     &                   cipr(IPR_VOL, ipa_idx, ispc) +
     &                                         dx(j)*dy*depth(i,j,k)
                      npastep(ipa_idx,ispc) = 
     &                                   npastep(ipa_idx,ispc) + 1
                    endif
                 enddo
            endif
c
c========================= Process Analysis End =====================================
c
c
c======================== Source Apportion Begin =======================
c
            if( ltrace .AND. tectyp .NE. RTRAC ) then
               l = 0
               do i = i1,i2
                 l = l + 1
                 if( lo3sp(ispc) ) then
                     fluxo3(l) = fluxo3(l) + saflux(l)
                     cnco3(l) = cnco3(l) + conc(i,j,k,ispc)
                 else if( lnoxsp(ispc) ) then
                     fluxnox(l) = fluxnox(l) + saflux(l)
                     cncnox(l) = cncnox(l) + conc(i,j,k,ispc)
                 else if( lvocsp(ispc) ) then
                     fluxvoc(l) = fluxvoc(l) + saflux(l) * crbnum(ispc)
                     cncvoc(l) = cncvoc(l) + conc(i,j,k,ispc) *
     &                                                   crbnum(ispc)
                 endif
               enddo
            endif
c
c======================== Source Apportion End =======================
c
c
cdbg            if ((ispc.eq.230).and.(j.eq.2)) then
cdbg               write(*,*)'Before x-advection'
cdbg	       do m=1,14
cdbg                  write(*,*)conc(2,2,m,230)
cdbg               enddo
cdbg               write(*,*)'dy=',dy
cdbg               write(*,*)'depth(2,2,k)=',(depth(2,2,m),m=1,14)
cdbg               write(*,*)'c1d(2)=',c1d(2)
cdbg               pause
cdbg            endif
c
            l = 1
            do i = i1+1,i2-1
              l = l + 1
              conc(i,j,k,ispc) = c1d(l)/dy/depth(i,j,k)
            enddo
c
cjgj For debuging of the negative concentration
c
cdbg              if (conc(i,j,k,ispc).lt.0) then
cdbg                 write(*,*)'Negative concentration occurs during 
cdbg     &		x-advection'
cdbg                 write(*,*)'Species=',spname(ispc),ispc
cdbg                 write(*,*)'Location i, j, k=',i ,j ,k
cdbg                 write(*,*)'Concentration',conc(i,j,k,ispc)
cdbg                 write(*,*)'dy=',dy
cdbg                 write(*,*)'depth=',depth(i,j,k)
cdbg                 pause
cdbg              endif
c
cdbg            if ((ispc.eq.230).and.(j.eq.2)) then
cdbg               write(*,*)'After x-advection'
cdbg	       do m=1,14
cdbg                  write(*,*)conc(2,2,m,230)
cdbg               enddo
cdbg               write(*,*)'dy=',dy
cdbg               write(*,*)'depth(2,2,k)=',(depth(2,2,m),m=1,14)
cdbg               write(*,*)'c1d(2)=',c1d(2)
cdbg               pause
cdbg            endif
c
c======================== DDM Begin =======================
c
            if( lddm ) then
               do iddm=1,nddmsp
                 l = 1
                 do i = i1+1,i2-1
                   l = l + 1
                   ioff = iptddm(ispc)+iddm-1
                   sens(i,j,k,ioff) = sen1d(l,iddm)/dy/depth(i,j,k)
                 enddo
               enddo
             endif
c
c======================== DDM End =======================
c
c-----Sum up fluxes in east and west sides
c
            if (flux1.gt.0.0) then
              fluxtmp(ispc,7,k)  = fluxtmp(ispc,7,k)  + flux1*dtuse
            else
              fluxtmp(ispc,8,k)  = fluxtmp(ispc,8,k)  + flux1*dtuse
            endif
            if (flux2.lt.0.0) then
              fluxtmp(ispc,5,k)  = fluxtmp(ispc,5,k)  - flux2*dtuse
            else
              fluxtmp(ispc,6,k)  = fluxtmp(ispc,6,k)  - flux2*dtuse
            endif
  22      continue
c
c======================== Source Apportion Begin =======================
c
c  --- call routine to update the tracer concentrations
c      based on the calculated fluxes ---
c
          if( ltrace .AND. tectyp .NE. RTRAC ) then
             call xfluxsa(ncol,nrow,nlay,ntotsp,
     &                ptconc(isaptr),i1+1,i2-1,j,k,dy,depth,mapscl,
     &                   fluxo3,cnco3,fluxnox,cncnox,fluxvoc,cncvoc)
          endif
c
c======================== Source Apportion End =======================
c
c  --- next row, non-parallelized loop ---
c
  10    continue
c
c  --- next layer ---
c
  21  continue
  20  continue
c
c  --- end of parallelized loop ---
c
c
      tcpu = dtime(tarray2)
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
c
      if (xyordr.eq.0) goto 300
c
c-----Advection in y-direction
c
 200  continue
c
      write(*,'(a20,$)') 'y advection ......'
      write(iout,'(a20,$)') 'y advection ......'
c
cdbg      write(*,*) 'before checkconc in y advec'
cdbg      call checkconc(conc,ncol,nrow,nlay,nspc,species,3)
cdbg      write(*,*) 'after checkconc in y advec'
c
      do 40 k = 1,nlay
        dtuse = deltat/nadv(k)
        do 41 istep = 1,nadv(k)
c
        do 30 i = 2,ncol-1
          j1 = 1
          j2 = nrow
          if (igrd.eq.1) then
            if (jbeg(i).eq.-999) goto 30
            j1 = jbeg(i) - 1
            j2 = jend(i) + 1
          endif
          l = 0
          do j = j1,j2-1
            l = l + 1
            v1d(l) = 2.*windv(i,j,k)/(mapscl(i,j+1) + mapscl(i,j))
            m1d(l) = mapscl(i,j)*mapscl(i,j)
          enddo
          l = l + 1
          v1d(l) = windv(i,j2,k)
          m1d(l) = mapscl(i,j2)*mapscl(i,j2)
          nn = j2 - j1 + 1
c
c======================== Source Apportion Begin =======================
c
c  --- initialize the tracer fluxes/concs to zero ----
c
         do j=1,MX1D
             fluxo3(j) = 0.
             fluxvoc(j) = 0.
             fluxnox(j) = 0.
             cnco3(j) = 0.
             cncvoc(j) = 0.
             cncnox(j) = 0.
         enddo
c
c========================= Source Apportion End ========================
c
c         do 42 ispc = 1,nspc
c
          do 42 ii = 1,nspc
            ispc = order(ii)
c
            l = 0
            do j = j1,j2
              l = l + 1
c
c For number concentration, unit is #/cm3 * m * m, but eventually,
c it will be converted back to the original unit at the end of 
c this routine.
c
              c1d(l) = conc(i,j,k,ispc)*dx(j)*depth(i,j,k)
            enddo
            if (igrd.eq.1 .and. v1d(1).lt.0.) c1d(1) = c1d(2)
            if (igrd.eq.1 .and. v1d(nn-1).gt.0.) c1d(nn) = c1d(nn-1)
c
c======================== DDM Begin =======================
c
            if( lddm ) then
               do iddm=1,nddmsp
                 l = 0
                 do j = j1,j2
                   l = l + 1
                   c1d0(l) = c1d(l)
                   ioff = iptddm(ispc)+iddm-1
                   sen1d(l,iddm) = sens(i,j,k,ioff)*dx(j)*depth(i,j,k)
                 enddo
                 if (igrd.eq.1 .and. v1d(1).lt.0.) 
     &                                     sen1d(1,iddm) = sen1d(2,iddm)
                 if (igrd.eq.1 .and. v1d(nn-1).gt.0.) 
     &                                 sen1d(nn,iddm) = sen1d(nn-1,iddm)
               enddo
            endif
c
c======================== DDM End =======================
c
cdbg            call mnratio(ii,ngas,naero,c1d,c2d,j1,j2,MX1D,mn,i,k,2) ! by jgj 8/13/06
            call mnratio(ii,c1d,c2d,j1,j2,mn,order,nspc,i,k,2) ! by jgj 8/13/06
c
c           In case of mass species, use the flxarr of the number species.
c           by jgj 8/13/06
c
            if (ii.gt.ngas) then
               spec=mod(ii-ngas,naero)
               if (spec.eq.1) then
                  do jj = 1,nn-1
                     nflxarr(jj) = 0
                  enddo
               else
                  call hadvppm2(ii,ngas,naero,MX1D,mn,nflxarr,
     &             nn,dtuse,dy,c1d,m1d,flxarr,flux1,flux2,
     &             saflux,fc1,fc2)
                  goto 500
               endif
            endif
c
            if (iadvct.eq.2) then
              call hadvbot(nn,dtuse,dy,c1d,v1d,m1d,flxarr,flux1,
     &                                 flux2,saflux,fpc,fmc,fc1,fc2)
            elseif( iadvct .eq. 3) then
              call hadvppm(nn,dtuse,dy,c1d,v1d,m1d,flxarr,flux1,
     &                                         flux2,saflux,fc1,fc2)
            endif
c
c           Save flxarr of number species by jgj 8/14/06
c
            if (ii.gt.ngas) then
               spec=mod(ii-ngas,naero)
               if (spec.eq.1) then
                  do jj = 1,nn-1
                     nflxarr(jj) = flxarr(jj)
                  enddo
               endif
            endif
c
 500        continue
c
c======================== DDM Begin =======================
c
c
            if( lddm .AND. iadvct .eq. 2 ) then
              call bottddm(nn,dtuse,dy,sen1d,nddmsp,c1d0,c1d,
     &                                            fpc,fmc,v1d,m1d)
            endif
c
c======================== DDM End =======================
c
c
c======================== Process Analysis Begin ====================================
c
            if( lipr ) then
c
c-----Change from Y-direction horizontal advection
c 
                l = 1
                do j=j1+1,j2-1
                  l = l+1
c
                  if( ipa_cel(i,j,k) .GT. 0 ) then
                     ipa_idx = ipa_cel(i,j,k)
c
c-----Flux at south boundary
c
                     cipr(IPR_SADV, ipa_idx, ispc) =
     &                   cipr(IPR_SADV, ipa_idx, ispc) +
     &                                         fc1(l)/dx(j)/depth(i,j,k)
c
c-----Flux at north boundary
c
                     cipr(IPR_NADV, ipa_idx, ispc) =
     &                   cipr(IPR_NADV, ipa_idx, ispc) +
     &                                        fc2(l)/dx(j)/depth(i,j,k)
c
                  endif
                enddo
            endif
c
c========================= Process Analysis End =====================================
c
c
c======================== Source Apportion Begin =======================
c
            if( ltrace .AND. tectyp .NE. RTRAC ) then
               l = 0
               do j = j1,j2
                 l = l + 1
                 if( lo3sp(ispc) ) then
                     fluxo3(l) = fluxo3(l) + saflux(l)
                     cnco3(l) = cnco3(l) + conc(i,j,k,ispc)
                 else if( lnoxsp(ispc) ) then
                     fluxnox(l) = fluxnox(l) + saflux(l)
                     cncnox(l) = cncnox(l) + conc(i,j,k,ispc)
                 else if( lvocsp(ispc) ) then
                     fluxvoc(l) = fluxvoc(l) + saflux(l) * crbnum(ispc)
                     cncvoc(l) = cncvoc(l) + conc(i,j,k,ispc) *
     &                                                   crbnum(ispc)
                 endif
              enddo
            endif
c
c======================== Source Apportion End =======================
c
            l = 1
            do j = j1+1,j2-1
              l = l+1
              conc(i,j,k,ispc) = c1d(l)/dx(j)/depth(i,j,k)
c
cjgj For debuging of the negative concentration
c
cdbg              if (conc(i,j,k,ispc).lt.0) then
cdbg                 write(*,*)'Negative concentration occurs during 
cdbg     &		y-advection'
cdbg                 write(*,*)'Species=',spname(ispc)
cdbg                 write(*,*)'Location i, j, k=',i ,j ,k
cdbg                 write(*,*)'Concentration',conc(i,j,k,ispc)
cdbg                 pause
cdbg              endif
            enddo
c
c======================== DDM Begin =======================
c
            if( lddm ) then
               do iddm=1,nddmsp
                 l = 1
                 do j = j1+1,j2-1
                   l = l + 1
                   ioff = iptddm(ispc)+iddm-1
                   sens(i,j,k,ioff) = sen1d(l,iddm)/dx(j)/depth(i,j,k)
                 enddo
               enddo
            endif
c
c======================== DDM End =======================
c
c
c-----Sum up fluxes in north and south sides
c
            if (flux1.gt.0.0) then
              fluxtmp(ispc,3,k)  = fluxtmp(ispc,3,k)  + flux1*dtuse
            else
              fluxtmp(ispc,4,k)  = fluxtmp(ispc,4,k)  + flux1*dtuse
            endif
            if (flux2.lt.0.0) then
              fluxtmp(ispc,1,k)  = fluxtmp(ispc,1,k)  - flux2*dtuse
            else
              fluxtmp(ispc,2,k)  = fluxtmp(ispc,2,k)  - flux2*dtuse
            endif
  42      continue
c
c======================== Source Apportion Begin =======================
c
c  --- call routine to update the tracer concentrations
c      based on the calculated fluxes ---
c
           if( ltrace .AND. tectyp .NE. RTRAC ) then
              call yfluxsa(ncol,nrow,nlay,ntotsp,
     &                ptconc(isaptr),i,j1+1,j2-1,k,dx,depth,mapscl,
     &                     fluxo3,cnco3,fluxnox,cncnox,fluxvoc,cncvoc)
           endif
c
c======================== Source Apportion End =======================
c
  30    continue
c
c  --- next layer, end of parallelized loop ---
c
  41  continue
  40  continue
c
      tcpu = dtime(tarray2)
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
c
      if (xyordr.eq.0) goto 100
c
 300  continue
c
c  ---- put fluxes into global array ----
c
      do i=1,nspc
        do j=1,8
          do k=1,nlay
            fluxes(i,j) = fluxes(i,j) + fluxtmp(i,j,k)
          enddo
        enddo
      enddo
c
      call flush(6)
      call flush(iout)
      return
      end
