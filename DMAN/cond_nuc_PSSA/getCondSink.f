
C     **************************************************
C     *  getCondSink                                   *
C     **************************************************

C     WRITTEN BY Jeff Pierce, May 2007

C     This subroutine calculates the condensation sink (first order loss
C     rate of condensing gases) from the aerosol size distribution.

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nk(ibins) - number of particles per size bin in grid cell
C     Nnuc - number of particles per size bin in grid cell
C     Mnuc - mass of given species in nucleation pseudo-bin (kg/grid cell)
C     Mk(ibins, icomp) - mass of a given species per size bin/grid cell
C     spec - number of the species we are finding the condensation sink for

C-----OUTPUTS-----------------------------------------------------------

C     CS - condensation sink [s^-1]
C     sinkfrac(ibins) - fraction of condensation sink from a bin

      SUBROUTINE getCondSink(Nko,Mko,spec,CS,sinkfrac)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nko(ibins), Mko(ibins, icomp)
      double precision CS, sinkfrac(ibins)
      integer spec

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      double precision pi, R    ! pi and gas constant (J/mol K)
      double precision mu                  !viscosity of air (kg/m s)
      double precision mfp                 !mean free path of air molecule (m)
      double precision Neps     !tolerance for number
      double precision density  !density [kg m^-3]
      double precision mp       !mass per particle [kg]
      double precision Dpk(ibins) !diameter of particle [m]
      double precision Kn       !Knudson number
      double precision beta(ibins) !non-continuum correction factor
      double precision Mktot    !total mass in bin [kg]
	double precision fudge	  ! fudge factor for Dp of pseudo bin
	                          ! because J is proportional to diameter^2
	                          ! in the kinetic regime
      double precision CSeps      ! Minimum CS for testing
      real Di       !diffusivity of gas in air (m2/s)

C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------
      double precision aerodens_PSSA
      real gasdiff
      external aerodens_PSSA
      external gasdiff

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(Neps=1.0d-10)
	parameter(fudge=1.2d0)
      parameter(CSeps=1.0d-20)

C-----CODE--------------------------------------------------------------

C     get some parameters
      mu=2.5277e-7*temp**0.75302
      mfp=2.0*mu/(pres*sqrt(8.0*0.0289/(pi*R*temp)))  !S&P eqn 8.6
      Di=gasdiff(temp,pres,98.0,Sv(spec))
C     get size dependent values
      do k=1,ibins
         if (Nko(k) .gt. Neps) then
            Mktot=0.d0
            do j=1,icomp
                  Mktot=Mktot+Mko(k,j)
            enddo
cdbg            print*,'k=',k,'Mktot=',Mktot
Ckpc  Density should be changed due to more species involed.             
            density=aerodens_PSSA(Mko(k,srtso4),0.d0,
     &              Mko(k,srtnh4),Mko(k,srtna),Mko(k,srth2o))  !assume bisulfate            
            mp=Mktot/Nko(k)
         else
            !nothing in this bin - set to "typical value"
            density=1500.
            mp=1.4*xk(k)
         endif
         Dpk(k)=((mp/density)*(6./pi))**(0.333)
         Kn=2.0*mfp/Dpk(k)                             !S&P eqn 11.35 (text)
         beta(k)=(1.+Kn)/(1.+2.*Kn*(1.+Kn)/alpha(spec)) !S&P eqn 11.35
      enddo      
         
      

C     get condensation sink
      CS = 0.d0
cdbg      print*,'CS components in getCondSink'
cdbg      print*,'Dpk=',Dpk
cdbg      print*,'Nko=',Nko
cdbg      print*,'beta=',beta
      do k=1,ibins
         CS = CS + Dpk(k)*Nko(k)*beta(k)
      enddo
cdbg      print*,'CS=',CS
      do k=1,ibins
         sinkfrac(k) = Dpk(k)*Nko(k)*beta(k)/CS
      enddo
      CS = 2.d0*pi*Di*CS/(boxvol*1D-6)  
      if (CS.le.CSeps) then
         print*,'CS is less than 1.0d-20 s-1.'
         print*,'CS=',CS
         print*,'Di=',Di
         print*,'Dpk',Dpk
         print*,'beta',beta
         print*,'Nko',Nko
         print*,'Mko',Mko
         print*,'boxvol',boxvol
         print*,'xk',xk
         print*,'A program is stopped.'
         stop
       endif
      return
      end
