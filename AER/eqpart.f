c=========================================================================
c  03/09/03: bkoo
c            - modified to eliminate organics (eqparto now deals with organics)
c            - moved call to wdiameter
c            - added call to eqparto
c  03/07/03: bkoo
c            - commented out last call to step (redundant)
c            - brought initial call to step out of IF statement
c              so that step can be called when HYBR is selected
c  12/05/02: tmg
c            - eliminated call to subroutine eqneut
c  10/22/02: tmg
c            - added parameter to subroutine step because not all sections
c              are always in use
c  01/25/02: tmg
c            - use dry and wet basis for diameter as appropriate
c  09/20/00: bkoo
c            - combine eqpart & eqparth
c  06/08/00: bkoo
c            - add coagulation & nucleation
c  Apr 2000: bkoo
c            - put a counter to avoid infinity loop
c            - fix negative concentrations
c            - need codes for organics evaporation
c=========================================================================
c
cgy supress error messages to units 39 and 90 (6/12/00)
c
c  EQILIBRIUM PARTITIONING OF CONDENSING OR EVAPORATING SPECIES
C  ACROSS A SIZE AND COMPOSITION RESOLVED AEROSOL DISTRIBUTION
C
C  THIS SUBROUTINE ASSUMES THAT GAS-AEROSOL EQUILIBRIUM IS ACHIEVED 
C  IN FOR THE SUM OF THE SIZE SECTIONAL COMPOSITIONS.
C
C  WRITTEN BY KEVIN P. CAPALDO
c  DEPARTMENT OF CHEMICAL ENGINEERING
C  CARNEGIE MELLON UNIVERSITY
C  JUNE 1998
C
c   SUBROUTINES REQUIRED
C     ISRPIA        EQUILIBRIUM THERMODYNAMIC SUBROUTINE ISORROPIA
C     ORGANICS      GIVES THE EQUILIBRIUM VAPOR PRESSURE OF ORGANIC SPECIES
C     STEP          CALLS ISRPIA IN REVERSE MODE FOR EACH SECTION
C     
C   VARIABLES
c     T             TIME IN SECONDS
C     Q(NTOTALx2)   ARRAY OF AEROSOL SPECIES FOR EACH SECTION (UG/M3) 
C                   THEN GASES (PPM)
C     DQDT(NTOTAL)  THE AEROSOL AND GAS DERIVATIVES (PER SECOND)
C     DQ(NSP)       THE TOTAL MASS OF EACH SPECIES TRANSFERED 
C                   FROM THE GAS TO THE AEROSOL PHASE (uMOLES/m3)
C          THE FOLLOWING TWO DIMENSIONAL ARRAYS ARE NECESARY FOR THE 
C          LINKING OF THE WEIGHTING FACTORS FOR NH4, NO3 AND CL.
C     FRQ(NSEC,nsp)   FRACTION OF DQ TRANSFERED TO EACH AEROSOL SECTION
C
      SUBROUTINE EQPART(t,q)


      include 'dynamic.inc'
      INCLUDE 'equaer.inc'                  ! ISORROPIA declarations

cbk      REAL*8 DQ(NSP),FRQ(nsec,nsp),WI(5),q0(3),q(ntotal)
cbk      real*8 qsav(ntotal), DQsav(nsp), accom(nsp)
c$$$      REAL*8 DQ(nexti),FRQ(nsec,nexti),WI(5),q0(3),q(ntotal) ! bkoo (03/09/03)
c$$$      real*8 qsav(nexti*nsec), DQsav(nexti), accom(nexti)    ! bkoo (03/09/03)
c$$$      real*8 qq, frq0(nsec,nexti) ! bkoo (10/07/03)
      REAL*8 DQ(nsp),FRQ(nsec,nsp),WI(NCOMP),q0(4),q(ntotal) ! JJ (05/08/15)
      real*8 qsav(ntotal), DQsav(nsp), accom(nsp)    ! JJ (05/08/15)
      real*8 qq, frq0(nsec,nsp) ! JJ (05/08/15)
      logical done
      real*8 WT(NCOMP),Gas(4),Aerliq(16),Aersld(22)
      real*8 Other(9),CNTRL(2)
      character SCASI*15

c__cf      if(aerm.eq.'EQUI') then                                     ! cf
cc                                                                     !
cc     STEP 1/3: CALCULATE NUCLEATION RATE FOR THE WHOLE STEP          !
cc                                                                     !
c       call nucl(q)                                                   !
cc                                                                     !
cc     STEP 2/3: CALCULATE COAGULATION RATE FOR THE WHOLE STEP         !
cc                                                                     !
c       call coagul(q)                                                 !
c                                                                      !
ccbk       call step(nsec,q) ! tmg (10/22/02)                          !
c      endif                                                           !
c      call step(nsecx2,q) ! bkoo (03/07/03)                           !
c      call wdiameter(q) ! bkoo (03/09/03)                             !
c      call eqparto(t,q) ! bkoo (03/09/03)                             !
c__cf                                                                  ! cf
C     STEP 1:  DETERMINE BULK EQUILIBRIUM
c
C  INPUT:
c     IPROB = 1=Reverse prob, 0=Foreward prob
C  1. [WI]
C     DOUBLE PRECISION array of length [NCOMP=5]. [NCOMP] is defined in
C     include file 'isrpia.inc'.
C     Total aerosol concentrations, expressed in moles/m3.
C     WI(1) - sodium,   expressed as [Na]
C     WI(2) - sulfate,  expressed as [H2SO4]
C     WI(3) - ammonium, expressed as [NH3]
C     WI(4) - nitrate,  expressed as [HNO3]
C     WI(5) - chloride, expressed as [HCl]
C
C  2. [RHI]
C     DOUBLE PRECISION variable.
C     Ambient relative humidity expressed on a (0,1) scale.
C
C  3. [TEMPI]
C     DOUBLE PRECISION variable. 
C     Ambient temperature expressed in Kelvins. 
C
C *** CONVERT INPUT CONCENTRATIONS TO moles/m3 **************************
C
      nsecx2 = nsec                         ! PMCAMx-UF does not define nsecx2 elsewhere/ JJ 13/07/15
      ! PMCAMx-UF does not calculate qn elsewhere/ JJ 13/07/15
      do i=1,nsec
        qt(i)=0.d0
        do ki=2,nsp-1 ! use dry basis first, tmg (01/23/02)
          qt(i)=qt(i)+q((i-1)*nsp+ki) ! total mass per section(ug/m3)
        enddo  
        qn(i)=qt(i)/dsec(i)**3 ! calculate 0th moment(number of particles)
      enddo

      ng = nsp*nsecx2                       ! gases
      prod=rgas*temp/(1.01325D5*pres)       ! conversion from ppm to umoles/m3
                                            ! pres (bkoo, 06/09/00)

      ! initialize WI and WT /JJ
      WI=0.d0
      WT=0.d0
      !initialize other isorropia arrays
      Gas=0.d0
      Aerliq=0.d0
      Aersld=0.d0

      WI(1) = 0.0D0
      WI(2) = Q(ng+ih2SO4) / PROD*1.0d-6
      WI(3) = q(ng+iNH3)   / PROD*1.0d-6
      WI(4) = q(ng+ihNO3)  / PROD*1.0d-6
c      WI(5) = q(ng+ihCL)   / PROD*1.0d-6
      WI(5) = 0.0D0                                                   ! cf
      WI(9) = q(nq+idma) / PROD*1.0d-6  ! amines
      do k=1,nsecx2                         ! aerosols
         nn=(k-1)*nsp
c         WI(1) = wi(1)+ q(nn+KNa) /emw(KNa) *1.D-6                   ! cf
         WI(2) = wi(2)+ q(nn+KSO4)/emw(KSO4)*1.D-6
         WI(3) = wi(3)+ q(nn+KNH4)/emw(KNH4)*1.D-6
         WI(4) = wi(4)+ q(nn+KNO3)/emw(KNO3)*1.D-6
c         WI(5) = wi(5)+ q(nn+KCL) /emw(KCL) *1.D-6                   ! cf
         WI(9) = wi(9)+ q(nn+kpami)/emw(kpami)*1.D-6 ! amines
      enddo
      RHI=rh
      TEMPI=temp
      CNTRL(1)=0.d0 !0=forward problem
      CNTRL(2)=0.d0 !0=aerosol can have both liquid and solid phases
c      IPROB = 0
C
C *** CALL ISORROPIA+ ***************************************************
C
c      CALL ISRPIA ( WI, RHI, TEMPI, IPROB )
      CALL ISOROPIA( WI, RHI, TEMPI, CNTRL,
     &     WT, GAS, AERLIQ, AERSLD, SCASI, OTHER)

c     initialize DQ & ACCOM arrays - bkoo (05/24/01)
cbk      do i=1,nsp
c      do i=1,nexti ! bkoo (03/09/03)
      do i=1,nsp ! JJ (08/15) 
         dq(i)    = 0.0d0
         accom(i) = 1.0d0
      enddo
C
C   STEP 2:  DETERMINE THE TOTAL TRANSFER BETWEEN GAS AND AEROSOL
C
c      DQ(KH2O)=0.0D0                        ! WATER DONE LATER
c      DQ(KNA) =0.0D0                        ! SODIUM IS IN AEROSOL PHASE ONLY
      DQ(KSO4)=Q(ng+ih2so4)/prod            ! all H2SO4 condenses
      DQ(KNH4)=Q(NG+INH3)  /PROD - dmax1(Gas(1)*1.0d6, 0.d0) ! bkoo (10/07/03)
      DQ(KNO3)=Q(NG+IHNO3) /PROD - dmax1(Gas(2)*1.0d6, 0.d0) ! bkoo (10/07/03)
c_cf      DQ(KCL) =Q(NG+IHCL)  /PROD - dmax1(GHCL *1.0d6, 0.d0) ! bkoo (10/07/03)     ! cf
      DQ(KCL) =0.0D0                                                              !   cf
      DQ(kpami)=Q(ng+idma)/prod - max(Gas(4)*1.0d6, 0.d0)  ! JJ (08/15)
C    DO ORGANICS  **ASSUME PS HAS NO SIZE OR COMPOSITION DEPENDANCE **
cbk  removed - bkoo (03/09/03)
cbk      DO IOG = 1,NORG-1
cbk         DQ(NEXTI+IOG)=(Q(NG+IHCL+IOG)-PS(1,IHCL+IOG)/
cbk     &                 (1.01325d-1*pres))/prod	! pres (bkoo, 06/09/00)
cbk      ENDDO
c
c   STEP 2B: PARTITION DQ's NH4, NO3, AND CL INTO NH4NO3, NH4CL, AND
C            EXCESS NH4 (which is associated with SO4)
C            DQ(KNO3) => now represents the transfer of NH4NO3
C            DQ(KCl) => now represents the transfer of NH4Cl
C            DQ(KNH4) => now represents the transfer of NH4 with SO4
c
      DQ(KNH4)=DQ(KNH4)-DQ(KNO3)-DQ(KCL)
c
c     Save initial aerosol concentrations
      qsav=q  ! qsav and q are now the same size /JJ (08/15) 
cbk      do i=1,ntotalx2
cbk        qsav(i)=q(i)
cJJ      do i=1,nsecx2                          ! bkoo (03/09/03)
cJJ       do ii=1,nexti                         ! bkoo (03/09/03)
cJJ        qsav((i-1)*nexti+ii)=q((i-1)*nsp+ii) ! bkoo (03/09/03)
cJJ       enddo                                 ! bkoo (03/09/03)
cJJ      enddo
c
c     Its possible to not be able to reach equilibrium for NH4Cl and NH4NO3
c     due to diferences between aerosol sectional vs. bulk composition 
      do i=1,4
        q0(i)=0.0d0
      enddo
      do isec=1,nsecx2
        q0(1)=q0(1)+q((isec-1)*nsp+kno3)    ! initial aerosol NO3
        q0(2)=q0(2)+q((isec-1)*nsp+knh4)    ! initial aerosol NH4
c        q0(3)=q0(3)+q((isec-1)*nsp+kcl)     ! initial aerosol Cl     !  cf
        q0(4)=q0(4)+q((isec-1)*nsp+kpami)   !initial aerosol dma
      enddo
C
c   STEP 3: DETERMINE THE RELATIVE RATES OF MASS TRANSFER FOR EACH SECTION
c       if we assume composition changes between size sections are not 
c       important, then the rate of transfer is proportional to the mass
c       mass transfer rate dependance on particle size

! We keep the dma accommodation coefficient on purpose as 1.0 /JJ

cbk   removed - bkoo (03/09/03)
cbk      accom(KH2O)=1.0
cbk      accom(KNa)=1.0
      accom(KSO4)=delta(ih2so4)
      accom(kno3)=delta(ihno3)
      accom(knh4)=delta(inh3)
      accom(kcl) =delta(ihcl)
cbk   removed - bkoo (03/09/03)
cbk      do isp=1,norg
cbk        accom(kcl+isp)=delta(ihcl+isp)
cbk      enddo
cbk      do isp=1,ninert
cbk        accom(kcl+norg+isp)=0.1
cbk      enddo

cbk      call wdiameter(q) ! tmg (01/25/02)

      rlambda=0.065d0
c
c calculate factors
c
      do isp=1,nsp   ! frq goes to nsp now /JJ
cJJ      do isp=1,nexti ! bkoo (03/09/03)
         frqtot = 0.0
         do isec = 1,nsecx2
            frq(isec,isp) = qn(isec)*
     &           dsec(isec)/(1.0+rlambda/(accom(isp)*dsec(isec)))
            frqtot = frqtot + frq(isec,isp)
         enddo
c
c normalize
c
         do isec = 1,nsecx2
            frq(isec,isp) = frq(isec,isp)/frqtot
            frq0(isec,isp) = frq(isec,isp) ! save frq - bkoo (10/07/03)
         enddo
      enddo

      iter=0 ! counter for escaping out of an infinity loop (bkoo: Apr, 2000)
 80   iter=iter+1 ! moved - bkoo (10/07/03)
c     save dq's
      do i=1,nsp
cJJ      do i=1,nexti ! bkoo (03/09/03)
         dqsav(i)=dq(i)
      enddo
c
c   FIRST condense all condensing species
c
      do isp = 1,nsp
cJJ      do isp = 1,nexti ! bkoo (03/09/03)
         if(dq(isp).gt.0.0d0) then
            do isec=1,nsecx2
               INDX=(ISEC-1)*NSP
               q(indx+isp)=q(indx+isp)+frq(isec,isp)*dq(isp)*emw(isp)
c     special cases for remaping of dq
               if(isp.eq.kNO3.or.isp.eq.kCl) then
                  q(indx+kNH4)= q(indx+knh4)+frq(isec,isp)*dq(isp)*emw(knh4)
               endif
            enddo
            dq(isp)=0.0d0
         endif
      enddo
c
c   SECOND evaporate all evaporating species
c   only NH4NO3 and NH4Cl can evaporate
c
C     CHECK FOR COMPLETE EVAPORATION
c$$$ 100  frt=1.0d0
c$$$      frtcl = 0.0 ! bkoo (10/07/03)
c$$$      isp=kCl ! Cl first b/c NH4Cl forms before NH4NO3 when Na exists
c$$$      do m=1,2
c$$$       if(dq(isp).lt.0.0d0) then            ! evaporating
c$$$        do isec=1,nsecx2
c$$$         INDX=(ISEC-1)*NSP
c$$$         dqfx=DQ(ISP)*FRQ(ISEC,isp)
c$$$         IF(-dqfx.gt.tinys) then            ! evaporating significantly
c$$$          if(Q(INDX+ISP).LT.-dqfx*emw(isp)) then ! not enough NO3 or Cl
c$$$           frtq=-q(indx+isp)/dqfx/emw(isp)
c$$$           if(frtq.lt.frt) then
c$$$            frt=frtq
c$$$            ispsav=isp 	                    ! species that evaporates first
c$$$            isecsav=isec                    ! section that evaporates first
c$$$            if(q(indx+isp).lt.tinys) goto 150
c$$$           endif
c$$$          endif
c$$$          qq = q(indx+knh4) ! bkoo (10/07/03)
c$$$          if(isp.eq.KNO3)   ! bkoo (10/07/03)
c$$$     &              qq = qq + frtcl*frq(isec,KCL)*dq(KCL)*emw(knh4)
c$$$cbk          if(q(indx+knh4).lt.-dqfx*emw(kNH4)) then ! not enough NH4+
c$$$cbk           frtq=-q(indx+knh4)/dqfx/emw(kNH4)
c$$$          if(qq.lt.-dqfx*emw(kNH4)) then    ! not enough NH4+
c$$$           frtq=-qq/dqfx/emw(kNH4)
c$$$           if(frtq.lt.frt) then
c$$$            frt=frtq
c$$$            ispsav=isp                      ! species that evaporates first
c$$$            isecsav=isec                    ! section that evaporates first
c$$$            if(q(indx+knh4).lt.tinys) goto 150
c$$$           endif
c$$$          endif
c$$$         endif
c$$$        enddo
c$$$        if(isp.eq.KCL) frtcl = frt ! bkoo (10/07/03)
c$$$       endif
c$$$       isp=kNO3
c$$$      enddo
c$$$c partition species up to evaporation point
c$$$      done=.true.
c$$$      do m=1,2
c$$$       do isec=1,nsecx2
c$$$        INDX=(ISEC-1)*NSP
c$$$        q(indx+isp)=q(indx+isp)+frt*frq(isec,isp)*dq(isp)*emw(isp)
c$$$        q(indx+kNH4)=q(indx+knh4)+frt*frq(isec,isp)*dq(isp)*emw(knh4)
c$$$       enddo
c$$$       dq(isp)=(1-frt)*dq(isp)
c$$$       if(abs(dq(isp)).gt.1e-5) done=.false. ! tolerance for solution
c$$$       isp=kCl
c$$$      enddo
c$$$c check for complete solution
c$$$      if(done) then
c$$$       goto 200
c$$$      endif
c$$$c
c$$$c adjust frq to account for the totally evaporated species
c$$$c
c$$$ 150  frqtot=0.0 ! bkoo (10/07/03) ...
c$$$      do isec=1,nsecx2
c$$$        if(isec.eq.isecsav) frq(isec,ispsav)=0.0d0
c$$$        frqtot = frqtot + frq(isec,ispsav)
c$$$      enddo
c$$$      if(frqtot.gt.tinys) then ! renormalize
c$$$        do isec=1,nsecx2
c$$$          frq(isec,ispsav)=frq(isec,ispsav)/frqtot
c$$$        enddo
c$$$      else
c$$$        if(iter.gt.itmaxeq) then ! moved - bkoo (10/07/03)
c$$$         if(min(dq(KNO3),dq(KCL)).lt.-0.3) then                               ! bkoo_dbg
c$$$         call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)                    ! bkoo_dbg
c$$$         write(*,'(A8,2E15.5,4I4)')'EQUI-F: ',dq(KNO3),dq(KCL),ichm,jchm,kchm ! bkoo_dbg
c$$$         endif                                                                ! bkoo_dbg
c$$$         goto 200
c$$$        endif
c$$$        goto 180 ! can't evaporate given species from any section
c$$$      endif
c$$$
c$$$      goto 100
c$$$c
c$$$c   step 4: when sectional compositions prevent achievement of equilibrium
c$$$c           gas phase concentrations we need to adjust the non limiting
c$$$c           species so that that species does achieve equilibrium.
c$$$c           For example:  if NH4Cl cannot evaporate because a section 
c$$$c           runs out of NH4 and this is the only particle with any Cl
c$$$c           then KNH4Cl will not equal [NH3]final [HCL]final becuase 
c$$$c           this equilibrium could not be achieved.  However the 
c$$$c           equilibrium transfer flux of NH4NO3 was calculated assuming 
c$$$c           that all of the NH4Cl would evaporate.  Since it didn't 
c$$$c           KNH4NO3 will not equal [NH3]final [HNO3]final.  To correct
c$$$c           in this case we will evaporate enough NH4NO3 to establish
c$$$c           equilibrium. 
c$$$c
c$$$c           only evaporating species can encounter this problem
c$$$c           only if NH4 becomes zero in a section can this problem occur
c$$$c
c$$$ 180  if(ispsav.eq.kno3) then
c$$$         isp2=kcl                     ! correcting species
c$$$         g2=ghcl*1.0d6                ! original equilibrium for HCL (umol/m3)
c$$$      else
c$$$         isp2=kno3                    ! correcting species
c$$$         g2=ghno3*1.0d6               ! original equilibrium for HNO3 (umol/m3)
c$$$      endif
c$$$      bb=gnh3*1.0d6+dq(ispsav)+g2
c$$$      dq(isp2)=(-bb+sqrt(bb**2-4*dq(ispsav)*g2))/2.0d0
c$$$
c$$$      dq(isp2) = dmax1(dq(isp2), dq(knh4)+dqsav(ispsav)-dq(ispsav) ! bkoo (10/07/03)
c$$$     &                          +dqsav(isp2)-q(ng+inh3)/prod)
c$$$
c$$$      if(iter.gt.itmaxeq+1) dq(isp2)=0.0d0 ! bkoo (11/14/01)
c$$$c reset aerosol concentrations and dq's for NO3, NH4 and Cl
c$$$      do isec=1,nsecx2
c$$$         indx=(isec-1)*nsp
c$$$cbk         q(indx+kno3)=qsav(indx+kno3) ! initial aerosol no3
c$$$cbk         q(indx+knh4)=qsav(indx+knh4) ! initial aerosol nh4
c$$$cbk         q(indx+kcl)=qsav(indx+kcl)   ! initial aerosol Cl
c$$$         indx2=(isec-1)*nsp          ! bkoo (03/09/03) /JJ (08/15)
c$$$         q(indx+kno3)=qsav(indx2+kno3) ! bkoo (03/09/03)
c$$$         q(indx+knh4)=qsav(indx2+knh4) ! bkoo (03/09/03)
c$$$         q(indx+kcl)=qsav(indx2+kcl)   ! bkoo (03/09/03)
c$$$      enddo
c$$$c restore frq - bkoo (10/07/03)
c$$$c      do isp=1,nexti
c$$$      do isp=1,nsp
c$$$         do isec=1,nsecx2
c$$$            frq(isec,isp) = frq0(isec,isp)
c$$$         enddo
c$$$      enddo
c$$$c adjust DQ
c$$$      do i=1,nsp
c$$$cJJ      do i=1,nexti ! bkoo (03/09/03)
c$$$         if(i.eq.ispsav) then
c$$$            dq(i)=dqsav(i)-dq(i)
c$$$         elseif(i.eq.isp2) then
c$$$            dq(i)=dqsav(i)-dq(i)
c$$$         elseif(i.eq.knh4) then
c$$$            dq(i)=dqsav(i)
c$$$         else
c$$$            dq(i)=0.0d0
c$$$         endif
c$$$      enddo
c$$$      goto 80
c$$$c      endif
C      
C   STEP 4: ASSIGN EQUILIBRIUM VALUES TO GASES
c   Its possible to not be able to reach equilibrium for NH4Cl and NH4NO3
c   do to diferences between aerosol sectional vs. bulk composition 
C   So we treat NH3, HNO3, and HCl conservatively
C
 200  Q(NG+IH2SO4)=0.0D0              ! all sulfuric acid condenses
      Q(NG+IHNO3)= Q(NG+IHNO3)+q0(1)*PROD/GMW(IHNO3)
      Q(NG+INH3) = Q(NG+INH3) +q0(2)*PROD/GMW(INH3)
      Q(NG+IHCL) = Q(NG+IHCL) +q0(3)*PROD/GMW(IHCL)
      Q(ng+idma) = q(ng+idma) +q0(4)*prod/gmw(idma) ! dma /JJ (08/15)
      do isec=1,nsecx2
       Q(NG+IHNO3)= Q(NG+IHNO3)-Q((ISEC-1)*NSP+KNO3)*PROD/GMW(IHNO3)
       Q(NG+INH3) = Q(NG+INH3) -Q((ISEC-1)*NSP+KNH4)*PROD/GMW(INH3)
       Q(NG+IHCL) = Q(NG+IHCL) -Q((ISEC-1)*NSP+KCL) *PROD/GMW(IHCL)
       q(ng+idma) = q(ng+idma) -q((isec-1)*nsp+kpami)*prod/gmw(idma) ! dma /JJ (08/15)
      enddo
c     correct negative NH3 - bkoo (10/07/03)
      dq(KNH4) = q(ng+inh3) / prod ! umol/m3
      if(dq(KNH4).lt.-tinys) then
        iter = 0
 300    frt = 1.d0
        do isec = 1,nsecx2
          indx=(isec-1)*nsp       
          if(frq0(isec,KNH4).gt.0.d0) frt = dmax1( dmin1(q(indx+KNH4)
     &             /(-dq(KNH4)*frq0(isec,KNH4)*emw(KNH4)),frt), 0.d0)
        enddo
        frqtot = 0.d0
        do isec = 1,nsecx2
          indx=(isec-1)*nsp
          q(indx+KNH4) = dmax1(q(indx+KNH4)
     &                  +frt*dq(KNH4)*frq0(isec,KNH4)*emw(KNH4),0.d0)
          if(q(indx+KNH4).lt.tinys) frq0(isec,KNH4) = 0.d0
          frqtot = frqtot + frq0(isec,KNH4)
        enddo
        q(ng+inh3) = q(ng+inh3) - frt*dq(KNH4)*prod
        ! check if we should evaporate more
        dq(KNH4) = (1.d0 - frt) * dq(KNH4)
        if(dq(KNH4).lt.-tinys) then
          if(frqtot.gt.tinys) then ! we have sections which are not empty    
            if(iter.le.itmaxeq) then ! check infinite loop
              iter = iter + 1
              do isec = 1,nsecx2
                frq0(isec,KNH4) = frq0(isec,KNH4) / frqtot
              enddo
              goto 300
            endif
          endif
          ! we need to evaporate more to achieve equilibrium
          ! but we completely evaporate the species in all sections
          ! or exceed itermax
          write(*,*)'EQUI-F2: dq(KNH4)=',dq(KNH4),' iter=',iter ! bkoo_dbg
        endif
      endif
C   DO ORGANICS  **ASSUME PS HAS NO COMPOSITION DEPENDANCE **
cbk removed - bkoo (03/09/03)
cbk      DO IOG = 1,NORG-1
cbk       Q(NG+IHCL+IOG) =PS(1,IHCL+NORG)/(1.01325d-1*pres) ! pres (bkoo, 06/09/00)
cbk      ENDDO
ccf
ccf      call negchk(t,q,nsecx2)
ccf
ccf      call eqneut(0, q) ! tmg (12/05/02)
c
cc_cf      if(aerm.eq.'EQUI') then ! (tmg,01/31/02)                       ! cf
C                                                                         !
C   STEP 5: CALL EQUILIBRIUM CODE FOR EACH SECTION TO DETERMINE WATER     !
C                                                                         !
cbk        CALL STEP(nsec,Q)  ! tmg (10/22/02) bkoo (03/07/03)            !
C                                                                         !
C   STEP 6: Determine new diameter                                        !
C                                                                         !
cc_cf        call ddiameter(q)                                            ! cf
c      endif
C
      RETURN
      END     

c ----------------------------------------------------------------------------
c neutralize the acidic or basicness of the aerosol in each section
c
c is   : 0 (false) for equilibrium section
c        1 (true)  for dynamic     section
c q    : aerosol [ug/m3] and gas [ppm] species array
c        dimension - ntotal  for equilibrium section of EQUI &
c                                dynamic     section of MADM
c                    ntotale for equilibrium section of HYBR
c                    ntotald for dynamic     section of HYBR
c
c CHTOT: Total charge balance [umole/m3]
c
cc_cf      subroutine eqneut(is, q)
c                                                                     ! cf
c      include 'dynamic.inc'
c
c      real*8      q(ntotal),chtot,chmove,chno3,chcl,chnh4,prod
c      integer     nsx,ng,na,chspc(nexti)
c      integer     jani,jcat ! bkoo (10/07/03)
c      integer     kerr ! bkoo (01/06/04)
c
cc charges of H2O, Na, SO4, NO3, NH4, Cl
c      data chspc / 0, 1, -2, -1, 1, -1 /
c
c      nsx  = nsecx * is + nsecx2 * (1 - is)
c      ng   = nsp * nsx
c      prod = rgas*temp/(1.01325d5*pres) ! PPM / prod = UMOLE/M3
c
cc      if(is) call negchk(0.d0,q,nsecx)
c
c      do i = 1, nsx
c        na = (i-1)*nsp
c        chtot = 0.0d0
c        jani = 0 ! bkoo (10/07/03)
c        jcat = 0 ! bkoo (10/07/03)
c        do j = 2, nexti ! the first element is H2O
c          chtot = chtot + DBLE(chspc(j)) * q(na+j) / emw(j)
c          if(chspc(j).lt.0 .and. q(na+j).gt.tinys) jani = 1 ! bkoo (10/07/03)
c          if(chspc(j).gt.0 .and. q(na+j).gt.tinys) jcat = 1 ! bkoo (10/07/03)
c        enddo
cc BASIC:
c        if(chtot .gt. 0.0d0 .and. jani .eq. 0) then ! bkoo (10/07/03)
c          chno3 = (q(ng+IHNO3)/prod)
c          chcl  = (q(ng+IHCL) /prod)
c          chnh4 = (q(na+KNH4) /emw(KNH4))
cc first move HCl in
c          chmove     = max(min(chtot, chcl),0.0d0)
c          q(na+KCL)  = q(na+KCL)  + chmove * emw(KCL)
c          q(ng+IHCL) = max(q(ng+IHCL) - chmove * prod,0.0d0)
c          chtot      = chtot - chmove
cc then move HNO3 in
c          chmove     = max(min(chtot, chno3),0.0d0)
c          q(na+KNO3) = q(na+KNO3) + chmove * emw(KNO3)
c          q(ng+IHNO3)= max(q(ng+IHNO3)- chmove * prod,0.0d0)
c          chtot      = chtot - chmove
cc then move NH4 out
c          chmove     = max(min(chtot, chnh4),0.0d0)
c          q(na+KNH4) = max(q(na+KNH4) - chmove * emw(KNH4),0.0d0)
c          q(ng+INH3) = q(ng+INH3) + chmove * prod
c          chtot      = chtot - chmove
c        endif                                       ! bkoo (10/07/03)
cc ACIDIC:
c        if(chtot .lt. 0.0d0 .and. jcat .eq. 0) then ! bkoo (10/07/03)
c          chtot = -1.0d0 * chtot
c          chno3 = (q(na+KNO3)/emw(KNO3))
c          chcl  = (q(na+KCL) /emw(KCL))
c          chnh4 = (q(ng+INH3)/prod)
cc first move NH3 in
c          chmove     = max(min(chtot, chnh4),0.0d0)
c          q(na+KNH4) = q(na+KNH4) + chmove * emw(KNH4)
c          q(ng+INH3) = max(q(ng+INH3) - chmove * prod,0.0d0)
c          chtot      = chtot - chmove
cc then move HCl out
c          chmove     = max(min(chtot, chcl),0.0d0)
c          q(na+KCL)  = max(q(na+KCL)  - chmove * emw(KCL),0.0d0)
c          q(ng+IHCL) = q(ng+IHCL) + chmove * prod
c          chtot      = chtot - chmove
cc then move HNO3 out
c          chmove     = max(min(chtot, chno3),0.0d0)
c          q(na+KNO3) = max(q(na+KNO3) - chmove * emw(KNO3),0.0d0)
c          q(ng+IHNO3)= q(ng+IHNO3)+ chmove * prod
c          chtot      = chtot - chmove
c        endif
c      enddo
cc
cc     check the charge balance in each section - bkoo (01/06/04)
cc
c      do i = 1, nsx
c        na = (i-1)*nsp
c        kerr = 0
c        if (q(na+KNA).gt.tinys .or. q(na+KNH4).gt.tinys) then
c          if (q(na+KSO4).le.tinys .and. q(na+KNO3).le.tinys .and.
c     &        q(na+KCL).le.tinys) kerr = 1
ccbk          if ( (2.*q(na+KSO4)/emw(KSO4)+q(na+KNO3)/emw(KNO3)+q(na+KCL)
ccbk     &            /emw(KCL))/(q(na+KNA)/emw(KNA)+q(na+KNH4)/emw(KNH4))
ccbk     &            .gt.0.2e10) kerr = 1
c        else
c          if (q(na+KNO3).gt.tinys .or. q(na+KCL).gt.tinys) kerr = 1
c        endif
c        if (kerr.ne.0) then
c          call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)
cctmg          write(iout,'(//,A)')'CHECK in EQNEUT'
cctmg          write(iout,*)' q(',i,'): ',(q(na+j),j=2,nexti)
cctmg          write(iout,*)' igrd,i,j,k: ',igrdchm,ichm,jchm,kchm
c        endif
c      enddo
c
c      return
c      end                                                            ! cf
c
