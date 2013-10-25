
c********************************************************************
c         THIS SUBROUTINE SET EXPERIMENTAL CONDITION
c********************************************************************
c
      subroutine experiment(ygas,icoag,icond,inucl)

      implicit none

      include 'aervaria.inc'
c  
      integer bin
      integer i
      integer icoag ! a coagulation flag
      integer icond ! a condensation flag
      integer inucl ! a nucleation flag
      integer j

      real ygas(ngas)
      real N ! added number concentration of OM.
c
      do i=1, ngas       ! make all gas concentrations equal to zero.
         ygas(i)=0.
      enddo 
c  
      do j=1, icomp      ! make all gas concentrations equal to zero.
         Gc(j)=0.
      enddo
c
c     Diameter of upper boundary of 11th particle equals to 10.159 nm.
c     It is hard wired thing to sellect 12th. Additional works are required. 
c
      do i=12, ibins     
         Nk(i)=0.     
         do j=1, icomp
            Mk(i,j)=0.
         enddo
      enddo
c   
c     Set monodispersed OM material in 0.5 um.
c
      N = 0. ! 10 cm-3
      bin = 28 ! 28th mean diamter is 0.46 um. It is the closest one.
               ! I need improve this.
      Nk(bin) = N * boxvol ! The number concetration goes to 28th bin.
      Mk(bin,srtorg) = 1.4 * xk(bin) * Nk(bin)
c
c     Turn off condensation and nucleation.
c     
      icond=0
      inucl=0
c
      return
      end
