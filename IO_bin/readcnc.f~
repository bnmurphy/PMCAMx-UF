      subroutine readcnc
c 
c     
c-----CAMx v4.02 030709
c
c     READCNC operates in two modes:
c        1) when LRSTRT = F, reads and cycles through the AIRQUALITY file
c           to current time/date, initializes coarse grid concentrations,
c           and maps concentrations to any nested grids
c        2) when LRSTRT = T, reads and cycles through the coarse and
c           fine grid INSTANT files (from previous run) to current
c           time/date, and initializes all grid concentrations directly
c     Finally, the routine reads the TOPCON file and initializes CALOFT
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002
c     ENVIRON International Corporation
c           
c     Modifications: 
c        1/7/99    Added conditional to call rdfgcon for only if there are nests
c      10/31/01    Improved assignment of coarse grid IC's based on whether
c                  this is a restart or not
c      04/17/03    Changed the logic so it no longer such a large local array.
c                  It was causing problems with stack size on the latest PGI
c                  compiler.
c  
c     Input arguments: 
c        none
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        INTRPCNC
c        RDFGCON
c             
c     Called by: 
c        CAMx
c
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'filunit.com'
      include 'bndary.com'
      include 'grid.com'
      include 'chmstry.com'
      include 'flags.com'
      include 'section.inc'
      include 'diameters.inc'
      include 'aerpar.inc' ! jgj 12/26/07
c
      character*4 icspec(10)
      dimension cinit(MXCOLA,MXROWA)
      character*10 tpspc
c
      integer ii,jj,kk
      integer ispc,isp
      integer order(nspec)
      integer isund ! check underbar
      real massum
      real nmbr
      real meandp3
      dimension sconc(ncol(1),nrow(1),nlay(1),nspec)
c
c-----Entry point
c
      iunit = iic
      if (lrstrt) iunit = irstc
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
c
c     For setting the order of conversion from con to conc by jgj 12/26/07
c     The order(kk) has a number of gas concentration first.
c     Then it saves mass concentrations at the first size bin.
c     Thereafter it saves number concentration at the size bin.
c     The order(kk) saves the order of species as the same way
c     as increasing the order of size bins.
c
      do kk=1,ngas
         order(kk) = kk
      enddo

      kk=ngas

      do ii=1,nsec
         do jj=1,naero
            kk=kk+1
            if (jj.eq.naero) then
               order(kk) = ngas + (naero-1)*nsec + ii
            else
               order(kk) = ngas + (jj-1)*nsec + ii
            endif
         enddo
      enddo
c
c-----Read through coarse grid concentration records until current time/date
c
      read(iunit,end=900) idat1,tim1,idat2,tim2
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      write(iout,'(a40,2(f7.0,i8.5))') 
     &      'Read initial condition file at ',tim1,idat1,tim2,idat2
      do lread = 1,nicspc
        do k = 1,nz
          read(iunit) idum,(icspec(n),n=1,10), 
     &                ((cinit(i,j),i=1,nx),j=1,ny) 
c
          if ((idat1.lt.date .or. (idat1.eq.date .and. tim1.le.time)) 
     &                 .and. (idat2.gt.date .or. (idat2.eq.date .and. 
     &                                            tim2.gt.time))) then
c
             do 90 lmod = 1,nspec
                lic = licmap(lmod,1)
                if( lic .NE. lread ) goto 90
                do j = 1,ny
                  do i = 1,nx
                    n3d = i + nx*(j - 1) + nx*ny*(k - 1)
                    n4d = n3d + nx*ny*nz*(lmod - 1)
                    if (lrstrt) then
                      conc(n4d) = cinit(i,j)
                    else
                      conc(n4d) = bdnl(lmod)
                      conc(n4d) = amax1(conc(n4d),cinit(i,j))
                    endif
                  enddo
                enddo
  90         continue
          endif
        enddo
      enddo
c
c-----If this is not a restart, interpolate coarse grid concentrations
c     to all fine grids
c
      if (.not.lrstrt) then
        if (ngrid.gt.1) then
          do ip = 1,ngrid
            do ic = 1,nchdrn(ip)
              ig = idchdrn(ic,ip)
              call intrpcnc(nspec,ncol(ip),nrow(ip),nlay(ip),i1(ig),
     &                      j1(ig),nmesh(ig),nmshv(1,ig),ncol(ig),
     &                      nrow(ig),nlay(ig),conc(iptr4d(ip)),
     &                      conc(iptr4d(ig)) )
            enddo
          enddo
        endif
c
c-----Convert from ppm to umol/m3
c
        do igrd = 1,ngrid
          nx = ncol(igrd)
          ny = nrow(igrd)
          nz = nlay(igrd)
          do l = 1,nspec
            do k = 1,nz
              do j = 1,ny
                do i = 1,nx
                  n3d = i + nx*(j - 1) + nx*ny*(k - 1)
                  n4d = n3d + nx*ny*nz*(l - 1)
                  if (l.le.ngas) then
                    convfac = densfac*273./tempk(iptr3d(igrd)-1+n3d)*
     &                        press(iptr3d(igrd)-1+n3d)/1013.
                  else
                    convfac = 1.
                  endif
cbk                  conc(iptr4d(igrd)-1+n4d) = 
cbk     &                AMAX1( bdnl(l),convfac*conc(iptr4d(igrd)-1+n4d) )
                  conc(iptr4d(igrd)-1+n4d) = convfac*
     &                     AMAX1( bdnl(l), conc(iptr4d(igrd)-1+n4d) )
                enddo
              enddo
            enddo
          enddo
        enddo
c     assign_dist will assign number distributions, use in case you do not 
c     have size resolved initial conditions (like in the US). Use numconv if 
c     you want to calculate the number emissions from the mass emissions and 
c     their sizes
        sconc = 0. !always zero here so that subroutine redistribute all of the mass
c        call assign_ndist(conc(1),ncol(1),nrow(1),nlay(1),sconc,nspec)
        call numconv(conc(1),ncol(1),nrow(1),nlay(1),sconc,nspec)
      else
c
c-----Otherwise, read fine grid concentrations from fine grid restart file
c
        if( ngrid .GT. 1 ) call rdfgcon(idat1,tim1)
      endif
c
c-----Read TOPCON file and initialize top concentrations
c
      do l = 1,nspec
        caloft(l) = 0.
      enddo
c
      write(idiag,*)
 200  read(itopc,'(a10,f10.0)',end=901) tpspc,ctin
        if (tpspc.eq.'HNO2      ') tpspc = 'HONO      '
        if (tpspc.eq.'HCHO      ' .and. kHCHO.eq.nspec+1)
     &                                        tpspc = 'FORM      '
      do l = 1,nspec
        if (tpspc.eq.spname(l)) then
          write(idiag,'(a,a10,a)') 
     &          'Topcon species ',tpspc,' read from file'
          caloft(l) = ctin
          goto 200
        endif
      enddo
      write(idiag,'(a,a10)') 'Topcon species not on internal list',
     &                        tpspc
      goto 200
c
c-----Use BDNL value for species not on CALOFT file
c
 901  continue
cjgj      do l = 1,nspec
      do l = 1,ngas
        caloft(l) = amax1(caloft(l), bdnl(l))
      enddo
c
c     In case a caloft is less than bdnl, caloft is set to bdnl if caloft
c     is not number concentration. For number conc, bdnl is not applied
c     It is calculated from sum of mass concentration corresponded size
c     bins. This is to avoid setting very tiny values to caloft so that 
c     we do not have M/N ratio problem in TOMAS
c                                                  12/26/07 jgj
      if (ngas.lt.nspec) then
        ispc=ngas
        do ii=1,nsec
           isp = ispc
           inum = order(ispc+naero)
           massum = 0
           do jj=1,naero
              ispc = ispc + 1
              l = order(ispc)
              if (jj.ne.naero) then ! To avoid number conc.
                 caloft(l) = amax1(caloft(l), bdnl(l))
                 isund=INDEX(spname(l),'_')
                 if (spname(l)(1:isund-1).ne.'PH2O') then ! To avoid water
                    massum = massum + caloft(l)
                 endif
              endif
           enddo
           meandp3 = (dsec_i(ii)*dsec_i(ii+1))**(3./2.)
           nmbr = massum/(pi6*meandp3*rho)*1.0e12 ! # cm-3
cdbg           write(*,*)'nmbr=',nmbr,' massum=',massum,' pi6=',pi6
cdbg           write(*,*)'meandp3=', meandp3,' rho=',rho
cdbg           write(*,*)'dsec_i=',dsec_i
cdbg           pause
           caloft(inum)=nmbr
        enddo

cdbg        do l=ngas+1,nspec
cjgj        caloft(l) = amax1(caloft(l), bdnl(l))
cdbg           caloft(l) = caloft(l)
cdbg        enddo
      endif
     
      goto 999
c
c-----End of IC file reached
c
 900  write(iout,'(//,a)') 'ERROR in READCNC:'
      write(iout,*)'End of IC file'
      write(iout,*)'Make sure initial condition file contains the ',
     &                                   'simulation beginning hour.'
      call camxerr()
c
 999  return
      end
