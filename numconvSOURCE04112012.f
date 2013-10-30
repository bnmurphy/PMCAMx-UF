
c
c     Omninus subroutines
c
c     1. numconv  : generates number concentrations
c     2. saveconc : save concentrations before emissions
c     3. numcheck : check number conc.s in a specific cell and write M/N
c     4. mnratios : check M/N ratios whether they are within 0.8 to 2.2
c
c
      subroutine numconv(conc,ncol,nrow,nlay,sconc,nspc,iflag)
c
c
c     NUMCONV generates number concentrations from mass concentrations.
c     In this version, emissions, BC, IC, and top conditions were 
c     redistributed offline. LNP 03282012
c     
c     Arguments:
c     conc ; concentrations including aerosol and gas
c     ncol ; a number of columns of whole grids
c     nrow ; a number of rows of whole grids
c     nlay ; a number of layers of whole grids
c     tmass ; total mass concentration before doing individual prcesses
c     iflag ; 0 = initial condition
c           ; 1 = boundary condition
c           ; 2 = emissions
c
c     Called by:
c     CAMx
c     emistrns
c
c-----Variable explanation
c
      include 'camx.prm'
      include 'camx.com'
      include 'section.inc'
      include 'chmstry.com'
      include 'aerpar.inc'
      include 'diameters.inc'
c
c-----Variables
c
      integer k, j, i, ii ! counters
      integer isec, iarspc ! section #, aerosol species # 
      integer n3d, n4d ! parameters in conc matrix
      integer iflag ! different depending on calling subroutines
      integer isund
      real mass ! [=] ug/m3, total mass of each size section
      real number ! [=] #/cm3, aerosol number concentration
      real meandp3 ! [=] a cube of mean diameter of particle
      dimension conc(ncol,nrow,nlay,nspc)
      dimension sconc(ncol,nrow,nlay,nspc)
c
c
c-----Generate number conc.
c
      do k = 1,nlay
         do j = 1,nrow
            do i = 1,ncol
               do isec = 1, nsec
                  mass = 0.
                  do iarspc = 1, naero-1
                     lmod = ngas + (iarspc-1)*nsec + isec
cdbg                     write(*,*)'spname(lmod)=',spname(lmod) !debug
                     isund=INDEX(spname(lmod),'_')
                     if (spname(lmod)(1:isund-1).ne.'PH2O') then
                        ! To avoid H2O
                           mass = mass + conc(i,j,k,lmod) 
                     endif
                  enddo
cdbg                  pause
                  lmod = ngas + (naero-1)*nsec + isec ! number conc.
                  meandp3 = (dsec_i(isec)*dsec_i(isec+1))**(3./2.)
                  number = mass/(pi6*meandp3*rho)*1.0d12 ! # cm-3
                  conc(i,j,k,lmod) = number
               enddo
            enddo !ncol
         enddo !nlow
      enddo !nlay

      return
      end
c
c---------------------------------------------------------------------
c
      subroutine saveconc(conc,ncol,nrow,nlay,nspc,sconc)
c
c
c     SAVECONC saves mass concentrations before redistributing.
c
c     Arguments:
c     conc ; concentrations including aerosol and gas
c     ncol ; a number of columes of whole grids
c     nrow ; a number of rows of whole grids
c     nlay ; a number of layers of whole grids
c     nspc ; a number of species
c     iflag ; 0 = Before emiss
c           ; 1 = After emiss
c
c     Called by:
c     CAMx
c     emistrns
c
c-----Variable explanation
c
      include 'camx.prm'
      include 'camx.com'
      include 'section.inc'
      include 'chmstry.com'
      include 'aerpar.inc'
      include 'diameters.inc'
c
c-----Variables
c
      integer k, j, i, ii, jj ! counters
      integer isec, iarspc ! section #, aerosol species # 
      integer n3d, n4d ! parameters in conc matrix
cdbg      integer iflag ! different depending on calling subroutines
      real mass ! [=] ug/m3, total mass of each size section
      real mass2 ! [=] ug/m3, total mass for organic species
      real mass3 ! [=] ug/m3, total mass for number concentratios
      real number ! [=] #/cm3, aerosol number concentration
      real meandp3 ! [=] a cube of mean diameter of particle
      real dist(14)
      dimension conc(ncol,nrow,nlay,nspc)
      dimension sconc(ncol,nrow,nlay,nspc)
c
c-----Generate number conc.
c
      do k = 1,nlay
         do j = 1,nrow
            do i = 1,ncol
               do isec = 1, nsec
                  do iarspc = 1, naero-1 
                         ! Exclude number conc.
                     lmod = ngas + (iarspc-1)*nsec + isec
                     sconc(i,j,k,lmod)=conc(i,j,k,lmod)
                  enddo
               enddo
            enddo
         enddo
       enddo
       return
       end
c
c---------------------------------------------------------------------
c
      subroutine numcheck(conc,ncol,nrow,nlay,nspc)
c
c
c     NUMCHECK checks number concentrations in a specific cell.
c
c     Arguments:
c     conc ; concentrations including aerosol and gas
c     ncol ; a number of columes of whole grids
c     nrow ; a number of rows of whole grids
c     nlay ; a number of layers of whole grids
c     tmass ; total mass concentration before doing individual prcesses
c
c     Called by:
c     CAMx
c     emstrns
c
      include 'camx.prm'
      include 'camx.com'
      include 'section.inc'
      include 'chmstry.com'
      include 'aerpar.inc'
      include 'diameters.inc'
c
c-----Variables
c
      integer k, j, i, isec, iarspc
      integer n3d, n4d
      integer isund
      real mass, number
      real meandp3
      dimension conc(ncol,nrow,nlay,nspc),
     &          xk(nsec) ! boundaries
c
c-----Generate number conc.
c
      !Pittsburgh
c      k=1
c      j=51
c      i=65

      !A problem cell
      ipcell=33
      jpcell=18
      kpcell=1
c
c Setting xk's
c
      xk(1) = 3.75315e-25 ![=]kg for rho of particles = 1.4e+12 ug/m3
      do isec = 2, nsec
         xk(isec) = xk(isec-1) * 2 
      enddo

      write(*,*)'In numcheck'

      do i=ipcell-1,ipcell+1
        do j=jpcell-1,jpcell+1
          do k=kpcell,kpcell
            do isec = 1, nsec-2 !To avoid the bins greater than 10um
              mass = 0.
              do iarspc = 1, naero-1 !To avoid number 
                 lmod = ngas + (iarspc-1)*nsec + isec
                 isund=INDEX(spname(lmod),'_')
                 if (spname(lmod)(1:isund-1).ne.'PH2O') then ! To avoid H2O
                    mass=mass+conc(i,j,k,lmod) 
                 endif
              enddo
              lmod = ngas + (naero-1)*nsec + isec
              number = conc(i,j,k,lmod)
              write(*,*)'Coordinate (i,j,k)',i,j,k
              write(*,*)'Species=',spname(lmod)
              write(*,*)'mass, number=',mass, number
              write(*,*)'M/N ratio',mass/number*1.0d-15/xk(isec)
            enddo
          enddo
        enddo
      enddo
c
cdbg      pause
c
      return
      end
c
c---------------------------------------------------------------------
      subroutine mnratios(conc,ncol,nrow,nlay,nspc,flag)
c
c     Arguments:
c     conc ; concentrations including aerosol and gas
c     ncol ; a number of columes of whole grids
c     nrow ; a number of rows of whole grids
c     nlay ; a number of layers of whole grids
c     tmass ; total mass concentration before doing individual prcesses
c
c     Called by:
c     CAMx
c     emistrns
c
c-----Variable explanation
c
      include 'camx.prm'
      include 'camx.com'
      include 'section.inc'
      include 'chmstry.com'
      include 'aerpar.inc'
      include 'diameters.inc'
c
c-----Variables
c
      integer k, j, i ! counters
      integer isec, iarspc ! section #, aerosol species # 
      integer n3d, n4d ! parameters in conc matrix
      integer iflag ! different depending on calling subroutines
      integer lmod, lmod2 ! the order of species in a given grid cell
      integer flag ! flag to know where it is called
      integer isund ! check underbar
      real mass ! [=] ug/m3, total mass of each size section
      real number ! [=] #/cm3, aerosol number concentration
      real meandp3 ! [=] a cube of mean diameter of particle
      dimension conc(ncol,nrow,nlay,nspc),
     &          rmass(ncol,nrow,nlay,nsec), 
     &          rmnratio(ncol,nrow,nlay,nsec), ! [=] kg/particle
     &          xk(nsec) ! boundaries
c
c Setting xk's
c
      xk(1) = 3.75315e-25 ![=]kg for rho of particles = 1.4e+12 ug/m3
      do isec = 2, nsec
         xk(isec) = xk(isec-1) * 2 
      enddo
c
      write(*,*)'mnratio by',flag
      do k = 1,nlay
         do j = 1, nrow
            do i = 1, ncol
               do isec = 1, nsec-2 ! neglect the size bin above 10 um
                                   ! because of lack of aerosol microphysics
                                   ! and aqueous chemistry
                  lmod2 = ngas + (naero-1)*nsec + isec
                  do iarspc = 1, naero-1 
                     ! Count mass species except number conc.
                     lmod = ngas + (iarspc-1)*nsec + isec
                     isund=INDEX(spname(lmod),'_')
                     if (spname(lmod)(1:isund-1).ne.'PH2O') then
                        ! To avoid H2O
                        rmass(i,j,k,isec) = rmass(i,j,k,isec) +  
     &                                     conc(i,j,k,lmod)
                     endif
                  enddo
                  rmnratio(i,j,k,isec) = rmass(i,j,k,isec) /
     &               conc(i,j,k,lmod2) * 1.0d-15 ! kg/particle
                  rmnratio(i,j,k,isec) = rmnratio(i,j,k,isec) /
     &               xk(isec)
                  if ((rmnratio(i,j,k,isec).lt.0.8).or.
     &              (rmnratio(i,j,k,isec).gt.2.2)) then
                    write(*,*)'Warning a size section =', isec
                    write(*,*)'in a grid cell of i=',i,' j=',j,' k=',k
                    write(*,*)'is out of boundary'
                    write(*,*)'M/N=',rmnratio(i,j,k,isec)
                    write(*,*)'Mass (ug/m3) =',rmass(i,j,k,isec)
                    write(*,*)'Number (#/cm3) =',conc(i,j,k,lmod2)
                    STOP
                  endif
                enddo
             enddo
         enddo
       enddo


       return
       end
