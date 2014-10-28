      subroutine emistrns(igrd)
c
c-----CAMx v4.02 030709
c
c     EMISTRNS performs the following tasks for one grid:
c        1.  Determines mass-conserving vertical velocity parameters
c        2.  Updates concentrations due to emissions
c        3.  Grows PiG puffs and adds dumped mass to grid
c        4.  Performs 3-D transport
c        5.  Performs wet scavenging 
c        6.  Updates met fields to current time step
c        7.  Performs 3-D diffusion
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003 
c     ENVIRON International Corporation
c            
c     Modifications:  
c        01/30/02    --gwilson--  Added code for RTRAC probing tool
c
c     Input arguments:
c        igrd                grid index
c
c     Output arguments:
c        none
c
c     Routines called:
c        ZRATES
c        EMISS
c        PIGDRIVE
c        XYADVEC
c        ZADVEC
c        WETDEP
c        UPDTMET
c        DIFFUS
c
c     Called by:
c        CAMx
c        NESTING
c  
      
      integer, intent(in):: igrd
      
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'grid.com'
      include "ptemiss.com"
      include "bndary.com"
      include 'chmstry.com'
      include 'flags.com'
      include 'filunit.com'
      include 'section.inc'
      include 'ahomap.com'
c
c======================== Process Analysis Begin ====================================
c
      include "procan.com"
c
c========================= Process Analysis End =====================================
c
c
c======================== Source Apportion Begin =======================
c
      include 'tracer.com'
      include 'rtracchm.com'
c
      real*8 ardum(MXTRSP), ptdum(MXTRSP), fluxdum(MXTRSP*11)
      real depdum(MXCOLA,MXROWA,3*MXSPEC)
      data ardum /MXTRSP*0.0/
      data ptdum /MXTRSP*0.0/
      data fluxdum /MXTRSP*0.0,MXTRSP*0.0,MXTRSP*0.0,MXTRSP*0.0,
     &              MXTRSP*0.0,MXTRSP*0.0,MXTRSP*0.0,MXTRSP*0.0,
     &              MXTRSP*0.0,MXTRSP*0.0,MXTRSP*0.0/
      dimension sconc(ncol(1),nrow(1),nlay(1),nspec)
                  ! For saving concentrations before calling emiss
      
c
c========================= Source Apportion End ========================
c
c-----Entry point
c
      write(*,'(a,i3,a,f7.1)') 'Processing grid:', igrd,
     &                         ' Timestep: ',deltat(igrd)
      write(iout,'(a,i3,a,f7.1)') 'Processing grid:', igrd,
     &                            ' Timestep: ',deltat(igrd)
c
c-----Determine mass-conserving vertical velocity parameters
c

      write(*,'(a20,$)') 'zrates ......'
      write(iout,'(a20,$)') 'zrates ......'
cdbg      call checkconc(conc(iptr4d(igrd)),ncol,nrow,nlay,nspec,spname,10)
      call zrates(igrd,xyordr,ncol(igrd),nrow(igrd),nlay(igrd),
     &            nadv(1,igrd),deltat(igrd),deltax(1,igrd),
     &            deltay(igrd),depth(iptr3d(igrd)),
     &            phpt(iptr3d(igrd)),pppt(iptr3d(igrd)),
     &            ptpt(iptr3d(igrd)),windu(iptr3d(igrd)),
     &            windv(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &            press(iptr3d(igrd)),mapscl(iptr2d(igrd)),dilut,entrn)
      tcpu = dtime(tarray2)
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
cdbg      call checkconc(conc(iptr4d(igrd)),ncol,nrow,nlay,nspec,spname,11)
      call flush(6)
      call flush(iout)
c
c-----Inject emissions
c
      write(*,'(a20,$)') 'emiss ......'
      write(iout,'(a20,$)') 'emiss ......'
      tcpu = dtime(tarray2)
cdbg      do k=1,14
cdbg         n3d = 2+97*(2-1)+97*90*(k-1)
cdbg         n4d = n3d+97*90*14*(230-1)
cdbg         write(*,*)conc(n4d)
cdbg      enddo
cdbg      pause
c
cdbg      write(*,*)'before saveconc'
      call saveconc(conc(iptr4d(igrd)),ncol,nrow,nlay,nspec,sconc)
cdbg      write(*,*)'after saveconc'
c
cdbg      call checkconc(conc(iptr4d(igrd)),ncol,nrow,nlay,nspec,spname,12)
      call emiss(igrd,kno,kno2,nspec,narspc,nptspc,larmap(1,igrd),
     &           lptmap,nosrc(igrd),idsrc(1,igrd),isrc(1,igrd),
     &           jsrc(1,igrd),ncol(igrd),nrow(igrd),nlay(igrd),
     &           deltat(igrd),deltax(1,igrd),deltay(igrd),
     &           mapscl(iptr2d(igrd)),
     &           height(iptr3d(igrd)),depth(iptr3d(igrd)),
     &           windu(iptr3d(igrd)),windv(iptr3d(igrd)),
     &           tempk(iptr3d(igrd)),press(iptr3d(igrd)),
     &           aremis(iptrem(igrd)),
     &           ptemis,
     &           armass(1,igrd),ptmass(1,igrd),
     &           conc(iptr4d(igrd)),ipacl_3d(iptr3d(igrd)) )
c
cdbg      call checkconc(conc(iptr4d(igrd)),ncol,nrow,nlay,nspec,spname,13)
cdbg      write(*,*)'before numconv'
      call numconv(conc(1),ncol(1),nrow(1),nlay(1),sconc,nspec,3)
cdbg      write(*,*)'after numconv after emiss'
cdbg      call checkconc(conc(iptr4d(igrd)),ncol,nrow,nlay,nspec,spname,14)
                           ! by jgj 2/17/06
c======================== Source Apportion Begin =======================
c
c  --- call routine with tracer species arrays, send dummy
c      arguemnts for the total mass arrays ---
c
c
      if( ltrace .OR. lddm ) then
         call emiss(igrd,kno,kno2,ntotsp,ntotsp,ntotsp,lsamap,lsamap,
     &              nosrc(igrd),idsrc(1,igrd),isrc(1,igrd),jsrc(1,igrd),
     &              ncol(igrd),nrow(igrd),nlay(igrd),deltat(igrd),
     &              deltax(1,igrd),deltay(igrd),mapscl(iptr2d(igrd)),
     &              height(iptr3d(igrd)),
     &              depth(iptr3d(igrd)),windu(iptr3d(igrd)),
     &              windv(iptr3d(igrd)),
     &              tempk(iptr3d(igrd)),press(iptr3d(igrd)),
     &              saemis(ipsa2d(igrd)),
     &              sapnts,ardum,ptdum,
     &              ptconc(ipsa3d(igrd)),ipacl_3d(iptr3d(igrd)) )
      endif
c
c========================= Source Apportion End ========================
c
      tcpu = dtime(tarray2)
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
      call flush(6)
      call flush(iout)
c
c-----Perform 2-D transport
c
c
cdbg      call mnratios(conc(1),ncol(1),nrow(1),nlay(1),nspec,11)
      call xyadvec(igrd,xyordr,ncol(igrd),nrow(igrd),nlay(igrd),nspec,
     &             MAX(1,ntotsp),nadv(1,igrd),
     &             deltat(igrd),deltax(1,igrd),deltay(igrd),
     &             windu(iptr3d(igrd)),
     &             windv(iptr3d(igrd)),depth(iptr3d(igrd)),
     &             mapscl(iptr2d(igrd)),conc(iptr4d(igrd)),
     &             fluxes(1,igrd),ptconc(MAX(1,ipsa3d(igrd))),
cjgj     &             tarray2,ipsa3d(igrd),ipacl_3d(iptr3d(igrd)) )
     &             tarray2,ipsa3d(igrd),ipacl_3d(iptr3d(igrd)),spname)
cdbg      call mnratios(conc(1),ncol(1),nrow(1),nlay(1),nspec,12)
      call flush(6)
      call flush(iout)
c
c-----Perform vertical transport
c
cdbg      do k=1,14
cdbg         n3d = 2+97*(2-1)+97*90*(k-1)
cdbg         n4d = n3d+97*90*14*(230-1)
cdbg         write(*,*)conc(n4d)
cdbg      enddo
cdbg      pause
c
      call zadvec(.FALSE.,igrd,ncol(igrd),nrow(igrd),nlay(igrd),nspec,
     &             MAX(1,ntotsp),deltat(igrd),deltax(1,igrd),
     &             deltay(igrd),
     &             densfac,idfin(iptr2d(igrd)),caloft,
     &             depth(iptr3d(igrd)),entrn,dilut,
     &             tempk(iptr3d(igrd)),press(iptr3d(igrd)),
     &             spname,conc(iptr4d(igrd)),
     &             fluxes(1,igrd),ptloft,ptconc(MAX(ipsa3d(igrd),1)),
     &             tarray2,ipacl_2d(iptr2d(igrd)),
     &                                        ipacl_3d(iptr3d(igrd)) )
cdbg      call mnratios(conc(1),ncol(1),nrow(1),nlay(1),nspec,13)
      call flush(6)
      call flush(iout)
c
c
c======================== Source Apportion Begin =======================
c
      if( ltrace ) then
c
         if( tectyp .EQ. RTRAC ) then
           call xyadvec(igrd,xyordr,ncol(igrd),nrow(igrd),nlay(igrd),
     &             ntotsp,1,nadv(1,igrd),deltat(igrd),deltax(1,igrd),
     &             deltay(igrd),windu(iptr3d(igrd)),windv(iptr3d(igrd)),
     &             depth(iptr3d(igrd)),mapscl(iptr2d(igrd)),
     &             ptconc(ipsa3d(igrd)),fluxdum,ptconc(1),
     &             tarray2,ipsa3d(igrd),ipacl_3d(iptr3d(igrd)) )
         endif
c
c----- call routine to advect the timing tracers horizontally ---
c
         if( npttim .GT. 0 ) call timadv(igrd,xyordr,ncol(igrd),
     &             nrow(igrd),nlay(igrd),nspec,
     &             deltat(igrd),deltax(1,igrd),deltay(igrd),
     &             windu(iptr3d(igrd)),
     &             windv(iptr3d(igrd)),depth(iptr3d(igrd)),
     &             mapscl(iptr2d(igrd)),ptconc(ipsa3d(igrd)))
c
c  --- call routine with tracer species arrays, send dummy
c      arguments for the total mass arrays ---
c
         call zadvec(.TRUE.,igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &                ntotsp,1,deltat(igrd),deltax(1,igrd),deltay(igrd),
     &                densfac,idfin(iptr2d(igrd)),ptloft,
     &                depth(iptr3d(igrd)),entrn,dilut,
     &                tempk(iptr3d(igrd)),press(iptr3d(igrd)),
     &                ptname, ptconc(ipsa3d(igrd)),
     &                fluxdum,ptloft,ptconc(1),
     &                tarray2,ipacl_2d(iptr2d(igrd)),
     &                ipacl_3d(iptr3d(igrd)) )
      endif
c
c========================= Source Apportion End ========================
cdbg      write(*,*)'after source apportion'
cdbg      call numcheck(conc,ncol(igrd),nrow(igrd),nlay(igrd),nspec)
c
c-----Perform wet scavenging
c
      if (lwet) then
        write(*,'(a20,$)') 'wetdep  ......'
        write(iout,'(a20,$)') 'wetdep  ......'
        call wetdep(igrd,ncol(igrd),nrow(igrd),nlay(igrd),nspec,
     &              deltat(igrd),deltax(1,igrd),deltay(igrd),
     &              depth(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &              press(iptr3d(igrd)),cwc(iptr3d(igrd)),
     &              pwc(iptr3d(igrd)),densfac,idfin(iptr2d(igrd)),
     &              conc(iptr4d(igrd)),fluxes(1,igrd),
     &              depfld(iptrdp(igrd)),dtout,ipacl_3d(iptr3d(igrd)) )
cdbg        call mnratios(conc(1),ncol(1),nrow(1),nlay(1),nspec,14)
        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
        call flush(6)
        call flush(iout)
c
c======================== Source Apportion Begin =======================
c
         if( ltrace .AND. tectyp .EQ. RTRAC ) then
           call wetdeprt(igrd,ncol(igrd),nrow(igrd),nlay(igrd),nrtrac,
     &              deltat(igrd),deltax(1,igrd),deltay(igrd),
     &              depth(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &              press(iptr3d(igrd)),cwc(iptr3d(igrd)),
     &              pwc(iptr3d(igrd)),densfac,idfin(iptr2d(igrd)),
     &              ptconc(ipsa3d(igrd)))
         endif
c
c======================== Source Apportion End =======================
c
cdbg      write(*,*)'after wet deposition'
cdbg      call numcheck(conc,ncol(igrd),nrow(igrd),nlay(igrd),nspec)
      endif
c
c-----Update vertical grid and environmental fields for this timestep
c
      write(*,'(a20,$)') 'updtmet ......'
      write(iout,'(a20,$)') 'updtmet ......'
      call updtmet(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &             ngas,densfac,deltat(igrd),phpt(iptr3d(igrd)),
     &             height(iptr3d(igrd)),depth(iptr3d(igrd)),
     &             pppt(iptr3d(igrd)),press(iptr3d(igrd)),
     &             pupt(iptr3d(igrd)),windu(iptr3d(igrd)),
     &             pvpt(iptr3d(igrd)),
     &             windv(iptr3d(igrd)),pspt(iptr2d(igrd)),
     &             tsurf(iptr2d(igrd)),
     &             ptpt(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &             conc(iptr4d(igrd)) )
      tcpu = dtime(tarray2)
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
      call flush(6)
      call flush(iout)
c
c-----Perform 3-D Diffusion
c
      call diffus(igrd,ncol(igrd),nrow(igrd),nlay(igrd),nspec,
     &          MAX(1,ntotsp),deltat(igrd),deltax(1,igrd),deltay(igrd),
     &          idfin(iptr2d(igrd)),vdep(iptrem(igrd)),
     &          rkx(iptr3d(igrd)),rky(iptr3d(igrd)),
     &          rkv(iptr3d(igrd)),depth(iptr3d(igrd)),
     &          tempk(iptr3d(igrd)),press(iptr3d(igrd)),
     &          mapscl(iptr2d(igrd)),conc(iptr4d(igrd)),
     &          fluxes(1,igrd),depfld(iptrdp(igrd)),
     &          ptconc(MAX(1,ipsa3d(igrd))),tarray2,
     &          '  z diffusion ......','     x/y diff ......',
     &          ipacl_2d(iptr2d(igrd)),ipacl_3d(iptr3d(igrd)) )
cdbg      call mnratios(conc(1),ncol(1),nrow(1),nlay(1),nspec,15)
      call flush(6)
      call flush(iout)
c
c======================== Source Apportion Begin =======================
c
      if( ltrace ) then
c
c   --- call routine to calculate the depostion velocities 
c       for the tracer species ---
c
         call filvdsa(ncol(igrd),nrow(igrd),nlay(igrd),nspec,ntotsp,
     &                     conc(iptr4d(igrd)),vdep(iptrem(igrd)),ptvdep,
     &                                             vdeprt(ipsa2d(igrd)))
c
c  --- call routine with tracer species arrays, send dummy
c      arguemnts for the total mass arrays ---
c
         call diffus(igrd,ncol(igrd),nrow(igrd),nlay(igrd),ntotsp,1,
     &               deltat(igrd),deltax(1,igrd),deltay(igrd),
     &               idfin(iptr2d(igrd)),ptvdep,rkx(iptr3d(igrd)),
     &               rky(iptr3d(igrd)),rkv(iptr3d(igrd)),
     &               depth(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &               press(iptr3d(igrd)),mapscl(iptr2d(igrd)),
     &               ptconc(ipsa3d(igrd)),fluxdum,depdum,ptconc(1),
     &               tarray2,
     &               '  SA z diffus ......','  SA x/y diff ......',
     &               ipacl_2d(iptr2d(igrd)),ipacl_3d(iptr3d(igrd)) )
      endif
c
c========================= Source Apportion End ========================
c
c
c-----Perform PiG evolution
c
      if( ipigflg .EQ. GRESPIG ) then
        write(*,'(a20,$)') 'pigdrive (GREASD) ......'
        write(iout,'(a20,$)') 'pigdrive (GREASD) ......'
        call flush(6)
        call gresdriv(igrd,iptr2d(igrd),ncol(igrd),nrow(igrd),
     &                nlay(igrd),itzon,deltat(igrd),
     &                deltax(1,igrd),deltay(igrd),mapscl(iptr2d(igrd)),
     &                height(iptr3d(igrd)),
     &                rkv(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &                press(iptr3d(igrd)),water(iptr3d(igrd)),
     &                windu(iptr3d(igrd)),windv(iptr3d(igrd)),
     &                cldtrns(iptr3d(igrd)),fcloud(iptr3d(igrd)),
     &                cellat(iptr2d(igrd)),cellon(iptr2d(igrd)),
     &                conc(iptr4d(igrd)),pigdump(1,igrd),ipsa3d(igrd),
     &                ipacl_3d(iptr3d(igrd)) )
        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
        call flush(6)
        call flush(iout)
      else if( ipigflg .EQ. IRONPIG ) then
        write(*,'(a20,$)') 'pigdrive (IRON) ......'
        write(iout,'(a20,$)') 'pigdrive (IRON) ......'
        call irondriv(igrd,iptr2d(igrd),ncol(igrd),nrow(igrd),
     &                nlay(igrd),lradm,itzon,deltat(igrd),
     &                deltax(1,igrd),deltay(igrd),mapscl(iptr2d(igrd)),
     &                height(iptr3d(igrd)),
     &                rkv(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &                press(iptr3d(igrd)),water(iptr3d(igrd)),
     &                windu(iptr3d(igrd)),windv(iptr3d(igrd)),
     &                cldtrns(iptr3d(igrd)),fcloud(iptr3d(igrd)),
     &                cellat(iptr2d(igrd)),cellon(iptr2d(igrd)),
     &                conc(iptr4d(igrd)),cncrad(iptrad(igrd)),
     &                pigdump(1,igrd))
        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
        call flush(6)
        call flush(iout)
      endif
cadded from ver 6 by Elham
c  --- call routine to do haze adjustment ---
c
c     added by LA
c      do ispc = 1,600
c         indx = 2 + 150 + (ispc-1)*150*162
c         write(*,*)'bef drtuv ispc =',ispc
c         write(*,*)'conc(i=2,j=2,k=1) =',conc(indx)
c      enddo   
c     end LA
         if( lchem ) then
c          print*,'ELHAM,lchem:', lchem
c           print *,'Elham,iproc_id:',iproc_id
           if( iproc_id .LE. 1 ) then
              write(*,'(a20,$)') 'TUV ......'
           endif
           write(iout,'(a20,$)') 'TUV ......'
           call flush(6)
           call drvtuv(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &                 nspec,time,date,itzon,
     &                 idfin(iptr2d(igrd)),
     &                 cellat(iptr2d(igrd)),cellon(iptr2d(igrd)),
     &                 height(iptr3d(igrd)),press(iptr3d(igrd)), 
     &                 tsurf(iptr2d(igrd)),tempk(iptr3d(igrd)),
     &                 water(iptr3d(igrd)),conc(iptr4d(igrd)),
     &                 cod(iptr3d(igrd)),icdalb(iptr2d(igrd)),
     &                 cldtrns(iptr3d(igrd)))
      
c      print *, 'igrd:',igrd,'mzp',mzp
c      print *,'iptr2d:',iptr2d     
c        if(cldtrns(iptr3d(igrd)).gt.1)
c      print *,'Elham.emistrns.cldtrns:',cldtrns(56+(129-1)*150) ! Hyytila grid cell
c      print *, 'Elham.emistrns.time:',time
c        endif
c      print *, 'Elham:date:',date 

    
c           if( iproc_id .LE. 1 ) then
c              write(*,'(a)') '   Done'
c           endif
           write(iout,'(a)') '   Done'
           call flush(6)
           call flush(iout)
c     added by LA
c           write(*,*)'Gone through drvtuv'
c           do ispc = 1,600
c              indx = 2 + 150 + (ispc-1)*150*162
c              write(*,*)'aft drtuv ispc =',ispc
c              write(*,*)'conc(i=2,j=2,k=1) =',conc(indx)
c           enddo   
c     end LA    
        endif
c
cdbg      write(*,*)'at the end of emistrns'
cdbg      call numcheck(conc,ncol(igrd),nrow(igrd),nlay(igrd),nspec)
      return
      end
