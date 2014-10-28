
      subroutine checkconc(conc,ncol,nrow,nlay,nspc,species,iflag)

      include "camx.prm"
      include "bndary.com"
      include "chmstry.com"
      include "section.inc"
c
      parameter (eps = 1.0d-20)
c
      dimension conc(ncol,nrow,nlay,nspc)
      integer iflag ! A flag shows where this checkconc is called
      integer isec, iarspc ! section #, aerosol species #
      integer lmods ! lmod for printing
      integer iarspcs ! iarspc for printing
      integer isecs ! isec for printing
      character*10 species(nspc) ! by jgj 7/18/06

      do k=1,nlay
         do j=1,nrow
            do i=1,ncol
               do isec = 1, nsec
                  do iarspc =1, naero 
                     lmod = ngas +(iarspc-1)*nsec + isec 
                                                   ! water is not accounted
                     if (conc(i,j,k,lmod).lt.0) then
cdbg                       if (conc(i,j,k,lmod).gt.-eps) then ! trucation error
cdbg                         conc(i,j,k,lmod)=eps
cdbg                       else
                         if (iarspc .eq. naero) then ! Numb. conc.
                           print*,'Negative number conc.'
                           print*,'Flag=',iflag
                           print*,'Coordinate i, j, k=',i, j, k 
                           print*,'Species=',spname(lmod) 
                           print*,'Concentration',conc(i,j,k,lmod)
                           conc(i,j,k,lmod) = eps
                           STOP
                         else
                           write(*,*)'Negative concentration occurs'
                           write(*,*)'Flag=',iflag
                           write(*,*)'Species=',spname(lmod)
                           write(*,*)'Location i,j,k=',i, j, k
                           write(*,*)'Concentration',conc(i,j,k,lmod)
                           do iarspcs =1, naero
                             lmods = ngas +(iarspcs-1)*nsec + isec
                             print*,'Species=',spname(lmods)
                             print*,'Concentration=',conc(i,j,k,lmods)
                           enddo
                           do isecs =1, nsec
                             lmods = ngas +(iarspc-1)*nsec + isecs
                             print*,'Species=',spname(lmods)
                             print*,'Concentration=',conc(i,j,k,lmods)
                           enddo
                           conc(i,j,k,lmod) = eps
                           STOP
cdbg                         endif
                       endif
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
