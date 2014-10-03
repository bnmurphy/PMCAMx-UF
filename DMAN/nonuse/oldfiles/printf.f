
C     *************************************************
C     * printf                                        *
C     *************************************************

C     WRITTEN BY JaeGun Jung, November 2007

C     This subroutine generates output files printing variables 
C     in following files.

C     Every time step
C       1. acidity.out - prints acidity of a particle whose diameter is 22
C         nm, ammonia mixing ratio as ppt, sulfuric acid as ppt, and
C         nucleation rate as particles cm-3            

C     Every 15 minutes
C       2. d_mdist.out - prints dM/dlogDp of Organic Matter (OM), and
C                      sulfate, and ammonium
C       3. Mk.out - prints mass concentrations of sulfate, OM, and 
C                  ammonium as ug m-3
C       4. distn.dat - a "MAIN OUTPUT" file, prints number concentration 
C                     as dN/dlogDp
C       5. s_tot.out - prints total surface area um2 cm-3
C       6. orgmass.out - print total mass concentrations of organic mass,
C                       sulfate, and ammonium
C       7. sulfuric_ammonia.out - print sulfuric acid as molec cm-3 
C                                 ammonia as ppt
C       8. Ntot.out - prints total number concentration below cutpoint 
C       9. dndlogdp.out - prints dN/dlogDp for coagulation or condensation
C                        tests

C-----INPUTS------------------------------------------------------------

C     Initial variables of
C     ====================

C     pn(nsect) - number concentrations as particles per cubic cm
C     ygas(ngas) - mixing ratios of gas species for input and output
C     deltat - time step as second
C     tstart - a beginning time of simulation
C     tend - an end time of simulation
C     s_tot - total surface of particles as um2 per cubic cm
C     time - simulation time of running main routine

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE printf(deltat, pn, tend, time, tstart, ygas) 

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'
      include 'IO.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real deltat
      real pn(nsect) ! [=] particles cm-3
      real tend
      real time
      real tstart
      real ygas(ngas) ! ygas [=] ppt

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer cutpoint ! set as the order of size bin of 10 nm paritcle
      integer i, j, z  ! counters for loops

      real cn ! Condensation Nuclei whose diameter is less than 0.1 um
      real cvt  ! a unit converter from kg/cm3 to ug/m3
      real dmdlog(nsect,icomp) ! dM/dlogDp values as ug per cubic m
      real dndlog(nsect) ! dN/dlogDp values as particles per cubic cm
      real dndlog4xk(ibins) 
         ! dN/dlogDp values for TOMAS variables as partilces per cubic cm
      real dummy(20) ! dummy variables
      real h2so4 ! sulfuric acid conc. molec cm-3
      real remainder ! time control
      real meanDp ! mean diameter below cutpoint
      real NiDpi ! a summation of number concentrations,Nk(ibins) times 
                 ! their mean diameters, dpmean4xk(ibins) to calculate 
                 ! a mean diameter of whole paritcles
      real Ntot ! Ntot is total number concentration below cutpoint.
      real R    ! gas constant J/mol-K
      real rnuc ! diameter of nuclei [nm]
      real s_tot ! total surface um2/cm3
      real traf(nsect) ! for traffic correction

      double precision fn ! nucleation rate [cm-3 s-1]
      double precision density
      double precision Mtot(icomp)
      
C-----EXTERNAL FUNCTIONS------------------------------------------------

      double precision aerodens_PSSA
      external aerodens_PSSA

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (R=8.314) !J/mol-K
      parameter (cutpoint = 11) ! corresponds to 10 nm

C-----CODE-------------------------------------------------------------

      !Initialization of variables
      fn = 0.d0

      !Set a conversion factor
      cvt = 1./boxvol*1.0e15 ! conversion from kg/cm3 to ug/m3

C-----PRINT every stepsize.

      ! 1. acidity.out
cacidity      if (ievaluation.eq.1) then ! will be affected
      if (time .eq. tstart) then
        write(iout+10,*)'time     22nm    NH3[ppt]  H2SO4[ppt]   fn'
      endif
      ! Calculate nucleation rate
      h2so4=Gc(srtso4)/boxvol*1000.0/gmw(srtso4)*6.022e+23
      ygas(mgsvi)=1.0e+21*R*temp/(pres*boxvol*gmw(srtso4))
     &              *Gc(srtso4)
      ygas(mgnh3)=1.0e+21*R*temp/(pres*boxvol*gmw(srtnh3))
     &              *Gc(srtnh3)

      if (h2so4.gt.1.0d4) then
        if ((ygas(mgnh3).gt.0.1).and.tern_nuc.eq.1) then
           call napa_nucl(dble(temp),dble(rh),dble(h2so4),
     &            dble(ygas(mgnh3)),fn,rnuc)
        elseif (bin_nuc .eq. 1) then
           call vehk_nucl(dble(temp),dble(rh),dble(h2so4),fn,rnuc)
        endif
      endif
      write(iout+10,100)time,Mk(15,srtnh3)/Mk(15,srtso4)/0.375*2
     &	,ygas(mgnh3),ygas(mgsvi),fn 
 100  format(5(1x,e9.4))
      fn=0.0d0 ! clear after printing out
cacidity      endif 

C-----PRINT distribution every 15 mins.

      ! 2. d_mdist.out
      if (ievaluation.eq.1) then 
        remainder=mod(time,0.25)
        if ((remainder .le. deltat/2.) .or. ((0.25-remainder)
     &    .le. deltat/2.).or.(time .eq. tstart)
     &    .or. ((time+deltat).gt.tend)) then

          if (time .eq. tstart) then
            do i = 1, nsect
              write(iout,*)dpbound(i)
            enddo
            call Nk2pn(pn,cn)
          endif

          write(iout,*)'time=',time

          do i=1, nsect-1
            dndlog(i)=pn(i)/log10(dpbound(i+1)/dpbound(i))
            do j=1,icomp-1
              dmdlog(i,j)=pm(i,j)/
     &          log10(dpbound(i+1)/dpbound(i))
            enddo
            write(iout,*)dmdlog(i,srtorg),dmdlog(i,srtso4)
     & 	      ,dmdlog(i,srtnh3)
          enddo
        endif
      endif

      ! 3. Mk.out
        !Print first boundaries of particles
cacidity      if (ievaluation .eq. 1) then
        if (time .eq. tstart) then
          do i=1, ibins+1
            density=aerodens_PSSA(Mk(i,srtso4),0.0,
     &            Mk(i,srtnh3),0,0,Mk(i,srth2o))
            dp4xk(i)=(xk(i)*1.0e18/pi6/density)**0.3333
                    ! get diameters from TOMAS variables,xk
            write(iout+3,*)dp4xk(i) ! save boundaries in Mk.out
          enddo

          do i=1,ibins+1
            dpmean4xk(i)=sqrt(dp4xk(i)*dp4xk(i+1))
                    ! get mean diameters of xk
          enddo
        endif

        remainder=mod(time,0.25)
        if ((remainder .le. deltat/2.) .or. ((0.25-remainder)
     &    .le. deltat/2.).or.(time .eq. tstart)
     &    .or. ((time+deltat).gt.tend)) then
          write(iout+3,*)'time=',time
          do i=1, ibins
            write(iout+3,*) Mk(i,srtso4)*cvt, Mk(i,srtorg)*cvt, 
     &	      Mk(i,srtnh3)*cvt ![=] ug m-3
          enddo
        endif
cacidity      endif

C-----Print every interval set by kprint

      ! kprint is 60 and deltat, that is time step for main loop,
      ! is 15 sec. So 'dist?.dat' files are updated every 900sec.
      ! Please do not misunderstand timestep as 900sec.

      ! 4. distn.dat
      if (kcount .ge. kprint) then
        write(iout+1,*) 'TIME=', time*3600
cdbg        do i=1, nsect-1
cdbg          traf(i)=0.0
cdbg        enddo
        !Do diesel external mixing.
cdbg        call ext_diesel(pn,time,traf)
        !diploy to dndlog
        call Nk2pn(pn,cn)
        do i=1, nsect-1
          dndlog(i) = pn(i)/log10(dpbound(i+1)/dpbound(i))
cdiesel          dndlog(i) = (pn(i)+traf(i))/log10(dpbound(i+1)/dpbound(i))
        enddo
        !Set dummy variables
        do z=1,20
          dummy(z)=0.0
        enddo
        do i=1,nsect
          write(iout+1,101) dpmean(i),dndlog(i),(dummy(z),z=1,7)
        enddo
      endif
 101  format(22(1x,e9.4))

      if (kcount .ge. kprint) then
        call Nk2pn(pn,cn)
        ! 5. s_tot.out
        do i=1,nsect
          s_tot=s_tot+pn(i)*3.14*dpmean(i)**2
        enddo
        write(iout+15,*)time,s_tot
        s_tot=0.

        ! 6. orgmass.out
        call Mtot_sub()
        write(iout+14,102) time,Mtot(srtorg),Mtot(srtnh3),Mtot(srtso4)
 102    format(8g10.5)

        ! 7. sulfuric_ammonia.out
        h2so4=Gc(srtso4)/boxvol*1000.0/gmw(srtso4)*6.022e+23
        ygas(mgnh3)=1.0e+21*R*temp/(pres*boxvol*gmw(srtnh3))
     &              *Gc(srtnh3)
        write(iout+17,*)time,h2so4,ygas(mgnh3)
        kcount = 0
      endif

      ! 8. Ntot.out
      if ((ievaluation.eq.1).and.(iexptest.eq.1)) then
        if (kcount .ge. kprint) then
           Ntot=0
           do i=1, cutpoint ! It is set now 11, which corresponds
             Ntot=Ntot+Nk(i)  ! to 10 nm. 
           enddo           
           write(iout+12,*)time,Ntot/boxvol
        endif
      endif


C-----Print mean diameter below cutpoint at 12:00

      if ((ievaluation.eq.1).and.(iexptest.eq.1)) then
        if ((time.le.(itcommon+deltat/2.))
     &       .and.(itcommon.eq.12)) then
          NiDpi = 0
          do i=1, cutpoint
            NiDpi = NiDpi + Nk(i) * dpmean4xk(i)
          enddo
          meanDp = NiDpi / Ntot
          write(*,*)'meanDp =',meanDp
        endif   
      endif

C-----Print for analytic solutions

      ! 9. dndlogdp.out
      if (ievaluation .eq. 1) then
        ! First print boundaries of diameters
        if (time .eq. tstart) then
          if ((icoag_test .eq. 1) .or. (icond_test .eq. 1)) then
            do i=1,ibins+1
              write(iout+2,*)dp4xk(i)
            enddo
          endif
        endif

      ! For coagulation tests, print every 3hr.
        if (icoag_test .eq. 1) then
          remainder=0.
          remainder=mod(time,3.)
          if (( remainder .le. deltat/2.) .or. ((3.-remainder) .le.
     &      deltat/2.) .or. (time .eq. tstart) .or. ((time+deltat)
     &      .gt.tend)) then
            write(iout+2,*)'time=',time
            do i=1,ibins
              dndlog4xk(i)=Nk(i)/log10(dp4xk(i+1)/dp4xk(i))/boxvol
              write(iout+2,*)dndlog4xk(i)
            enddo
          endif
        endif

      !For condensation test, print at 0:00, 0:03, 0:10 and 0:40.
        if (icond_test .eq. 1) then
          if ((time .eq. tstart).or.((time.le.3./60.+deltat/2.).and.
     &      (time.gt.3./60.-deltat/2.)).or.((time.le.10./60.+deltat/2.)
     &      .and.(time.gt.10./60.-deltat/2.)).or.((time.le.40./60.+
     &      deltat/2.).and.(time.gt.40./60.-deltat/2.))) then
            write(iout+2,*)'time=',time
            do i=1,ibins
              dndlog4xk(i)=Nk(i)/log10(dp4xk(i+1)
     &          /dp4xk(i))/boxvol
               write(iout+2,*)dndlog4xk(i)
            enddo
          endif
        endif
      endif



      RETURN
      END
