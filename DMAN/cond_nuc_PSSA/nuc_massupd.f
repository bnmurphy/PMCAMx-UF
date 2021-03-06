C     **************************************************
C     *  nuc_massupd                                   *
C     **************************************************

C     Written by Jan Julin, November 2014

C     This subroutine updates the number of particles and mass
C     due to nucleation. The calculation was previously done within
C     the subroutine "nucleation"

C-----INPUTS-----------------------------------------------------------

C     Initial values of
C     =================

C     Nkif(ibins) - number of particles per size bin in grid cell
C     Mkif(ibins, icomp) - mass of a given species per size bin/grid cell
C     Gcif(icomp-1) - amount (kg/grid cell) of all species present in the
C                    gas phase except water
C     dt - total model time step to be taken (s)
C     fn - nucleation rate (cm-3 s-1)
C     rnuc - critical radius (nm)
C     mfrac(icomp) - fraction of the nucleated mass to be added to the different 
C                    species

C-----OUTPUTS-----------------------------------------------------------

C     Nkif, Mkif, Gcif - values are updated 

      SUBROUTINE nuc_massupd(Nkif,Mkif,Gcif,nuc_bin,dt,fn,rnuc,mfrac)

      IMPLICIT NONE

      include 'sizecode.COM'

C-----INPUTS AND OUTPUTS------------------------------------------------

      double precision,intent(inout) :: Nkif(ibins), Mkif(ibins, icomp), Gcif(icomp-1)
      integer nuc_bin
      double precision dt
      double precision fn    ! nucleation rate [cm-3 s-1]
      double precision rnuc  ! critical radius [nm]
      double precision mfrac(icomp)

C-----OTHER VARIABLES-------------------------------------------------

      double precision mnuc     !mass of nucleation
      double precision pi

      parameter (pi=3.14159)

C---------------------------------------------------------------------      

      !currently Gc is updated neither here or in nucleation.f but we carry it around anyway
      Gcif(srtso4)=Gcif(srtso4)

      !calculate the mass of the nucleated particles
      mnuc = (4.d0/3.d0*pi*(rnuc*1D-9)**3)*1350.d0
      if (mnuc.lt.xk(1))then !this should never happen anymore since rnuc>=0.4 nm
         print*,'mnuc < xk(1) in nucleation routines',mnuc
         stop
      endif
      !Determine to which bin the mass is added
      nuc_bin = 1
      do while (mnuc .gt. xk(nuc_bin+1))
         nuc_bin = nuc_bin + 1
      enddo
      
      !update mass based on nucleation rate and mfrac
      Mkif(nuc_bin,srtso4) = Mkif(nuc_bin,srtso4)+mfrac(srtso4)*fn*mnuc*boxvol*dt
      Mkif(nuc_bin,srtna) = Mkif(nuc_bin,srtna)+mfrac(srtna)*fn*mnuc*boxvol*dt
      Mkif(nuc_bin,srtnh4) = Mkif(nuc_bin,srtnh4)+mfrac(srtnh4)*fn*mnuc*boxvol*dt
      Mkif(nuc_bin,srth2o) = Mkif(nuc_bin,srth2o)+mfrac(srth2o)*fn*mnuc*boxvol*dt

      !update number distribution
      Nkif(nuc_bin) = Nkif(nuc_bin)+fn*boxvol*dt

      RETURN
      END
