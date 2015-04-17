      subroutine wrtcon(iflag,tim2,idat2,iunit,incf,c_ncf_avrg,nox,noy,noz,
     &                  nsptmp, cncfld, iunitNuc, Jnucfld,
     &                  cellon, cellat, height)

      USE io_ezcdf
c
c-----CAMx v4.03 031205
c 
c     WRTCON writes average and instantaneous concentration fields; for
c     average files, optionally writes only layer 1; for instantaneous files,
c     writes all layers.
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        1/20/99   Grid cell size on file should be meters for all cartesian
c                  projections (UTM, LCP, PSP)
c        10/24/01  Removed BSWAP and converted integer strings to character*4
c 
c     Input arguments:
c        iflag               output type flag (0=average/1=instantaneous)
c        tim2                output time (HHMM)
c        idat2               output date (YYJJJ)
c        iunit               output unit
c        nox                 number of cells in x-direction
c        noy                 number of cells in y-direction
c        noz                 number of layers
c        nsptmp              number of species
c        cncfld              concentration field to output (ppm or umol/m3)
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        CAMx
c 
      include 'camx.prm'
      include 'camx.com'
      include 'grid.com'
      include 'chmstry.com'
      include 'flags.com'
c
      character*20 Snuc(2)
      real Jnucfld(nox,noy,noz,2)
      integer iunitNuc, incf

      integer  iflag, idat2, iunit
      character*4 ispec(10,MXSPEC), ifile(10), note(60)
      character*200 c_ncf_avrg, cfil
      integer nox, noy, noz, nvar
      real cncfld(nox,noy,noz,nsptmp), tim2
      real cellon(nox,noy), cellat(nox,noy), height(nox, noy, noz)
      character*10 cnfil
c
c-----Data statements
c
      
      INTEGER :: ji, jj, jk

      INTEGER :: id_fil, id_var, nsptmp, lz, lt, lct, lx, ly 
      REAL,DIMENSION(nox,noy) :: vlon,  vlat
      REAL,DIMENSION(noz) :: vdpth
      REAL, DIMENSION(24) :: vtime
      REAL, DIMENSION(nox,noy,noz) :: x3d
      REAL :: vflag
      CHARACTER(len=20) :: cvarlon, cvarlat, cvardpth
      CHARACTER(len=20), DIMENSION(nsptmp+2) :: cvar, cunit, cln
      CHARACTER(len=20) :: cvartime, cvar1, cun_z, cun_t
 

c
c-----Data statements
c
      data cnfil /'INSTANT   '/
      data nseg,izero /1,0/
      data zero /0./

      data Snuc /'Ternary','Binary'/
c
c-----Entry point
c
c-----Determine time/date range
c
      idat1 = idat2 
      etim = ANINT(tim2)/100.
      if( dtout .GE. 60.0 ) then
          btim = ANINT( 1000*(etim - ANINT(dtout)/60.) )/1000.
      else
          btim = ANINT( 1000*(etim - ANINT(dtout)/100.) )/1000.
      endif
      if (btim.lt.0.) then 
        btim = btim + 24. 
        idat1 = idat1 - 1 
      endif 
      idat3 = idat2
      etim3 = etim + 0.1
      if (etim3.gt.24.) then 
        etim3 = 24. - etim3 
        idat3 = idat3 + 1 
        if( MOD(idat3,1000) .GT. 365 ) then
            if( MOD(INT(idat3/1000),4) .EQ. 0 ) then
               if( MOD(idat3,1000) .EQ. 367 )
     &                     idat3 = (INT(idat3/1000)+1)*1000 + 1
            else
               idat3 = (INT(idat3/1000)+1)*1000 + 1
            endif
         endif
      endif 
c
c-----For instantaneous files, rewind file and write the header
c
      nlayer = noz
      if (iflag.eq.1) then
        read(cnfil,'(10a1)') (ifile(n),n=1,10)
        read(runmsg(1:60),'(60a1)') (note(n),n=1,60)
        if (.NOT.llatlon) then 
          orgx = 1000.*xorg 
          orgy = 1000.*yorg 
          dx = 1000.*delx 
          dy = 1000.*dely 
          izone = 0
          if (lutm) izone = iuzon 
        else 
          orgx = xorg 
          orgy = yorg 
          dx = delx 
          dy = dely 
          izone = 0 
        endif
        do l=1,nsptmp
          read(spname(l),'(10a1)') (ispec(n,l),n=1,10) 
        enddo
c
        rewind(iunit)
        write(iunit) ifile,note,nseg,nsptmp,idat2,etim,idat3,etim3
        write(iunit) zero,zero,izone,orgx,orgy,dx,dy,
     &               nox,noy,noz,izero,izero,zero,zero,zero
        write(iunit) izero,izero,nox,noy
        write(iunit) ((ispec(n,l),n=1,10),l=1,nsptmp)
        write(iunit) idat2,etim,idat3,etim3

      else

        !Prepare to Write Average Concentrations to Binary File
        do l = 1,nsptmp
          read(spname(lavmap(l)),'(10a1)') (ispec(n,l),n=1,10)
        enddo
        if (.not.l3davg) nlayer = 1
        write(iunit) idat1,btim,idat2,etim

        !------------------------NETCDF--------------------------------!
        !Prepare to Write Data to NETCDF File
        lx  = nox
        ly  = noy
        lz  = noz
        if (.not.l3davg) lz    = 1

        lt  = 24
        lct = int(tim2/100.)  !Current Time step (hours)
	if (lct.eq.0) lct = 24

        vlon   = cellon  !2D array of longitude
        vlat   = cellat  !2D array of latitude
        vdpth  = sum(sum(height,1),1)/(nox*noy)  !1D array of heights

        do it = 1,lt
          vtime(it) = real(it,8)
        enddo

        cfil     = c_ncf_avrg
        cvarlon  = 'lon'
        cvarlat  = 'lat'
        cvardpth = 'depth'
        cvartime = 'time'
        vglag  = 999.0
        cun_z  = 'm'
        cun_t  = 'hr'

        !First Initialize The NETCDF File
        if (lct.eq.1) then
          linit = 1   !initialization

          !Attributes for Variables
          nvar = nsptmp+2
          do l = 1,nvar   !Add 2 to make room for Nucleation Rates

            if (l.gt.nsptmp) then    !Variable is a Nucleation Rate
              cvar(l)  = Snuc(l-nsptmp)      !Species shortname
              cunit(l) = 'N cm-3 s-1'
              cln(l)   = cvar(l) 

            elseif (lavmap(l).le.ngas) then  !Species is a gas
              cvar(l)  = spname(lavmap(l))  !Species shortname
              cunit(l) = 'ppbv'   !Units
              cln(l)   = cvar(l)  !Long Name

            else    !Species is an Aerosol
              cvar(l)  = spname(lavmap(l))  !Species shortname
              cunit(l) = 'ug m-3' !Units
              cln(l)   = cvar(l)  !Long Name
            endif 
          enddo

          CALL P3D_T_irr_init(incf, id_var, linit, lx, ly, lz, lt, lct, 
     &       vlon, vlat, vdpth(1:lz), 
     &       vtime, cfil, cvarlon, cvarlat, cvardpth, cvartime,  
     &       nvar, cvar, cunit, cln, vflag, cun_z, cun_t) 
        endif  
      
        !Add species concentrations
        do l = 1,nsptmp+2 !Loop over Variables
          id_var = l      !Temporary Species identifier 
                          !  It will be overwritten in P3D_T_irr
          if (l.le.nsptmp) then  !Variable is a Species
            x3d(:,:,:) = cncfld(:,:,:,l)  !Conc. Array
            cvar1 = spname(lavmap(l))      !Species shortname

          else                   !Variable is a Nucleation Rate
            x3d(:,:,:) = Jnucfld(:,:,:,l-nsptmp)  !Nuc. Array
            cvar1 = Snuc(l-nsptmp)
	  endif

          !Write Array
          CALL P3D_T_irr_var(incf, id_var, lx, ly, lz, lt, lct, 
     &       x3d, cfil, cvar1)
        enddo

	!Sync NetCDF File on Disk with Memory Buffer
        call disp_err(NF_SYNC(incf), 'WRTCON','AVERAGE FILE','SYNCING')     
        !------------------------NETCDF--------------------------------!
      endif


!-----Write gridded concentration field for instantaneous concentrations
!     Or Binary Average File Depending on Value of iflag
      if (iflag.eq.1) then
         do l = 1,nsptmp
            do k = 1,nlayer
               write(iunit) nseg,(ispec(n,l),n=1,10),
     &               ((cncfld(i,j,k,l),i=1,nox),j=1,noy)
            enddo
         enddo
      endif

!-----Write Gridded Nucleation Rates
c$$$      if (iflag.eq.0) then
c$$$        write (iunitNuc) idat1,btim,idat2,etim
c$$$	do l = 1,2
c$$$	  do k = 1,nlayer
c$$$	    write(iunitNuc) nseg,Snuc(l),
c$$$     &              ((Jnucfld(i,j,k,l),i=1,nox),j=1,noy)
c$$$          enddo 
c$$$        enddo
c$$$      endif
      

      return
      end
