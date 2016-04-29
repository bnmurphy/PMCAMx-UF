
C     **************************************************
C     *  ezcond                                        *
C     **************************************************

C     WRITTEN BY Jeff Pierce, May 2007

C     This subroutine takes a given amount of mass and condenses it
C     across the bins accordingly.  

C     ADD MORE HERE!

C-----INPUTS------------------------------------------------------------

C     Initial values of
C     =================

C     Nki(ibins) - number of particles per size bin in grid cell
C     Mki(ibins, icomp) - mass of a given species per size bin/grid cell [kg]
C     mcond - mass of species to condense [kg/grid cell]
C     spec - the number of the species to condense

C-----OUTPUTS-----------------------------------------------------------

C     Nkf, Mkf - same as above, but final values

      SUBROUTINE ezcond(Nki,Mki,mcond,spec,Nkf,Mkf,ichm,jchm,kchm)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nki(ibins), Mki(ibins, icomp)
      double precision Nkf(ibins), Mkf(ibins, icomp)
      double precision mcond
      integer spec
      integer ichm ! i coordinate in PMCAMx
      integer jchm ! j coordinate in PMCAMx
      integer kchm ! k coordinate in PMCAMx

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j,k,c           ! counters
      double precision pi, R    ! pi and gas constant (J/mol K)
      double precision CS       ! condensation sink [s^-1]
      double precision sinkfrac(ibins+1) ! fraction of CS in size bin
      double precision Nk1(ibins), Mk1(ibins, icomp)
      double precision Nk2(ibins), Mk2(ibins, icomp)
      double precision madd     ! mass to add to each bin [kg]
      double precision maddp(ibins)    ! mass to add per particle [kg]
      double precision mconds ! mass to add per step [kg]
      integer nsteps            ! number of condensation steps necessary
      integer floor, ceil       ! the floor and ceiling (temporary)
      double precision eps     ! small number
      double precision tdt      !the value 2/3
      double precision mpo,mpw  !dry and "wet" mass of particle
      double precision WR       !wet ratio
      double precision tau(ibins) !driving force for condensation
      double precision totsinkfrac ! total sink fraction not including nuc bin
      double precision CSeps ! lower limit for condensation sink
      double precision fracch(ibins,icomp)
      double precision totch

C     VARIABLE COMMENTS...

C-----EXTERNAL FUNCTIONS------------------------------------------------


C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(pi=3.141592654, R=8.314) !pi and gas constant (J/mol K)
      parameter(eps=1.d-20)
      parameter(CSeps=1.d-20)

C-----CODE--------------------------------------------------------------

      tdt=2.d0/3.d0

      ! initialize variables
      do k=1,ibins
         Nk1(k)=Nki(k)
         do j=1,icomp
            Mk1(k,j)=Mki(k,j)
         enddo
      enddo

cJJ      call mnfix_PSSA(Nk1,Mk1,ichm,jchm,kchm,1)
      call mnfix_PSSA(Nk1,Mk1,ichm,jchm,kchm)

      ! get the sink fractions
      call getCondSink(Nk1,Mk1,spec,CS,sinkfrac)

      ! make sure that condensation sink isn't too small
      if (CS.lt.CSeps) then ! just make particles in first bin
         Mkf(1,spec) = Mk1(1,spec) + mcond
         Nkf(1) = Nk1(1) + mcond/sqrt(xk(1)*xk(2))
         do j=1,icomp
            if (icomp.ne.spec) then
               Mkf(1,j) = Mk1(1,j)
            endif
         enddo
         do k=2,ibins
            Nkf(k) = Nk1(k)
            do j=1,icomp
               Mkf(k,j) = Mk1(k,j)
            enddo
         enddo
cdbg                print*,'After no condensation in ezcond -3'
cdbg                print*,'Nk=',Nkf
cdbg                print*,'Mk=',Mkf
cdbg                pause
         return
      endif
c      print*,'CS',CS
c      print*,'sinkfrac',sinkfrac
c      print*,'mcond',mcond

      ! determine how much mass to add to each size bin
      ! also determine how many condensation steps we need
      totsinkfrac = 0.d0
      do k=1,ibins
         totsinkfrac = totsinkfrac + sinkfrac(k) 
                     ! get sink frac total not including nuc bin
      enddo
      nsteps = 1
      do k=1,ibins
         if (sinkfrac(k).lt.1.0D-20)then
            madd = 0.d0
         else
            madd = mcond*sinkfrac(k)/totsinkfrac
         endif
         mpo=0.0
         do j=1,icomp-idiag
            mpo=mpo + Mk1(k,j)
         enddo
         floor = int(madd*2.0/mpo)
         ceil = floor + 1
         nsteps = max(nsteps,ceil) 
                     ! don't let the mass increase by more than 10%
      enddo
cdbg      print*,'nsteps',nsteps

      ! mass to condense each step
      mconds = mcond/nsteps

      ! do steps of condensation
      do i=1,nsteps
         if (i.ne.1) then
            call getCondSink(Nk1,Mk1,spec,CS,sinkfrac)
            totsinkfrac = 0.d0
            do k=1,ibins
               totsinkfrac = totsinkfrac + sinkfrac(k) 
                     ! get sink frac total not including nuc bin
            enddo
         endif      

         do k=1,ibins
            mpo=0.0
            mpw=0.0
            !WIN'S CODE MODIFICATION 6/19/06
            !THIS MUST CHANGED WITH THE NEW dmdt_int.f
            if (Nk1(k) .gt. 0.d0) then            
            do j=1,icomp-idiag
               mpo = mpo+Mk1(k,j)  !accumulate dry mass
            enddo
            do j=1,icomp
               mpw = mpw+Mk1(k,j)  ! have wet mass include amso4
            enddo
            WR = mpw/mpo  !WR = wet ratio = total mass/dry mass
cd            if (Nk1(k) .gt. 0.d0) then
               maddp(k) = mconds*sinkfrac(k)/totsinkfrac/Nk1(k)
               mpw=mpw/Nk1(k)
cdbg               print*,'k',k,'mpw',mpw,'maddp(k)',maddp(k),'WR',WR
               tau(k)=1.5d0*((mpw+maddp(k)*WR)**tdt-mpw**tdt)  
                              !added WR to moxid term (win, 5/15/06)
            else
               !nothing in this bin - set tau to zero
               tau(k)=0.d0
               maddp(k) = 0.d0
            endif
         enddo

cJJ         call mnfix_PSSA(Nk1,Mk1,ichm,jchm,kchm,2)
         call mnfix_PSSA(Nk1,Mk1,ichm,jchm,kchm)

         !adjustment by jgj from a comparison with analytic solutions
cdbg         do k=1,ibins
cdbg           tau(k)=1.4*tau(k)
cdbg         enddo

         ! do condensation
         call tmcond(tau,xk,Mk1,Nk1,Mk2,Nk2,spec,maddp)
         totch=0.0
         do k=1,ibins
            do j=1,icomp
               fracch(k,j)=(Mk2(k,j)-Mk1(k,j))
               totch = totch + (Mk2(k,j)-Mk1(k,j))
            enddo
         enddo
c         print*,'fracch',fracch,'totch',totch

         if (i.ne.nsteps)then
            do k=1,ibins
               Nk1(k)=Nk2(k)
               do j=1,icomp
                  Mk1(k,j)=Mk2(k,j)
               enddo
            enddo            
         endif

      enddo

      do k=1,ibins
         Nkf(k)=Nk2(k)
         do j=1,icomp
            Mkf(k,j)=Mk2(k,j)
         enddo
      enddo



      RETURN
      END
