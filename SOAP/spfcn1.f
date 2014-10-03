      subroutine spfcn1(nsol,sx,f,scsat,sctot,smw,cpre,mwpre)
      implicit none
c
c SPFCN1 evaluates the functions defining the organic aerosol 
c phase mole fractions
c
c Called by SPFCN within the SOAP algorithm
c
c Variable declarations
c
      integer nsol
      real sx(nsol),f(nsol)
      real scsat(nsol),sctot(nsol),smw(nsol)
      real cpre,mwpre
c
      real sumtot,smprod,sumx,xb,xx,xp
      integer i
c
c Entry point
c
      sumtot=0.0
      smprod=0.0
      sumx=0.0      
      xx = cpre/mwpre
c
c THIS IS SYSTEM OF EQUATIONS TO SOLVE (BY MAKING THEM = 0)
c
      do i=1,nsol
         xb = sctot(i)/smw(i)
         sumtot = sumtot+xb
         smprod = smprod + sx(i)*(scsat(i)/smw(i))
         sumx = sumx + sx(i)
      enddo
c
      do i=1,nsol-1
         f(i)=sx(i)*(sumtot + xx - smprod + (scsat(i)/smw(i)))
     &                                         - sctot(i)/smw(i)
      enddo
c
c     singularity check - bkoo (10/28/03)
      if (xx.eq.0.0) then
        f(nsol) = sumx - 1.0
      else
        xp = sumtot - smprod + xx
        if (xp.gt.-1.e-20 .and. xp.lt.1.e-20) then
          f(nsol) = 1.0
        else
          f(nsol) = sumx + xx/xp - 1.0
        endif
      endif
c
      return
      end
