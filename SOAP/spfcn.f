      subroutine spfcn(nsol,sx,f,par,npar)
      implicit none
c
c SPFCN unpacks the variables from array par and calls SPFCN1
c
c Called by HYBRD within the SOAP algorithm
c
c Variable declarations
c
      integer nsol, npar
      real sx(nsol),f(nsol),par(npar)
c
      integer icsat, ictot, imw, icpre, imwpre
c
c Entry point
c
      icsat = 1
      ictot = nsol+1
      imw   = 2*nsol+1
      icpre = 3*nsol+1
      imwpre= 3*nsol+2
c
      call spfcn1(nsol,sx,f,par(icsat),par(ictot),par(imw),
     &            par(icpre),par(imwpre))
c
      return
      end
