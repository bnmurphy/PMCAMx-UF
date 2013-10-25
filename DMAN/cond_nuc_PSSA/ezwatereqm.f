
C     **************************************************
C     *  ezwatereqm                                    *
C     **************************************************

C     WRITTEN BY Peter Adams, March 2000

C     This routine uses the current RH to calculate how much water is 
C     in equilibrium with the aerosol.  Aerosol water concentrations 
C     are assumed to be in equilibrium at all times and the array of 
C     concentrations is updated accordingly.

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE ezwatereqm(Mke)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Mke(ibins,icomp)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer k
      double precision so4mass, naclmass
      double precision wrso4, wrnacl
      double precision rhe

      double precision waterso4, waternacl
      external waterso4, waternacl

c      double precision test

C     VARIABLE COMMENTS...

C     This version of the routine works for sulfate and sea salt
C     particles.  They are assumed to be externally mixed and their
C     associated water is added up to get total aerosol water.
C     wr is the ratio of wet mass to dry mass of a particle.  Instead
C     of calling a thermodynamic equilibrium code, this routine uses a
C     simple curve fits to estimate wr based on the current humidity.
C     The curve fit is based on ISORROPIA results for ammonium bisulfate
C     at 273 K and sea salt at 273 K.

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      rhe=rh*100.d0
c      test = 0.8d0*100.d0
c      print*,'rh',rh
c      print*,'rhe',rhe
c      print*,'test',test
      if (rhe .gt. 99.d0) rhe=99.d0
      if (rhe .lt. 1.d0) rhe=1.d0

c      print*,'rhe',rhe

      wrso4=waterso4(rhe)
cjgj      wrnacl=waternacl(rhe) ! No NaCl in DMAN.

c      print*,'wrso4',wrso4
c      print*,'rh',rh

      do k=1,ibins

         so4mass=Mke(k,srtso4)*1.2  !1.2 converts kg so4 to kg nh4hso4
cjgj         naclmass=Mke(k,srtna)      !already as kg nacl - no conv necessary

cjgj         Mke(k,srth2o)=so4mass*(wrso4-1.d0)+naclmass*(wrnacl-1.d0)
         Mke(k,srth2o)=so4mass*(wrso4-1.d0)

      enddo


      RETURN
      END












