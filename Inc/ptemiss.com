c-----CAMx v4.02 030709
c  
c     PTEMISS.COM contains all information regarding elevated point sources
c                            
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c            
c     Modifications:  
c        none
c 
c-----------------------------------------------------------------------
c     Variables for point source parameters:
c
c     nptsrc  -- number of point source stacks
c     xstk    -- stack x-location (km or degrees)
c     ystk    -- stack y-location (km or degrees)
c     lpiglet -- flag indicating a PiG point source
c     idstk   -- stack ID
c     hstk    -- stack height (meters)
c     dstk    -- stack diameter (meters)
c     tstk    -- stack exit temperature (K)
c     vstk    -- stack exit velocity (meters/s)
c     effph   -- effective plume height (meters)
c                normally ignored; if negative, PLUMERIS is overridden
c     ptemis  -- point source emission rates (mol/s)
c     ididstk -- stack id pointer for NetCDF I/O
c     idtstk  -- stack temperature pointer for NetCDF I/O
c     idvstk  -- stack velocity pointer for NetCDF I/O
c     locpnt  -- stack location pointer for NetCDF I/O
c-----------------------------------------------------------------------
c
      character*20 idstk(MXPTSRC)
      integer      nptsrc
      real         xstk(MXPTSRC,MXGRID)
      real         ystk(MXPTSRC,MXGRID)
      real         hstk(MXPTSRC)
      real         dstk(MXPTSRC)
      real         tstk(MXPTSRC)
      real         vstk(MXPTSRC)
      real         effph(MXPTSRC)
      real         ptemis(MXPTSRC,MXSPEC)
      logical      lpiglet(MXPTSRC)
      integer      ididstk
      integer      idtstk
      integer      idvstk
      integer      locpnt(MXHRS+1)
c
      common /ptechr/ idstk
      common /ptemiss/ nptsrc, xstk, ystk, lpiglet, hstk,
     &                 dstk, tstk, vstk, effph, ptemis, ididstk, idtstk,
     &                 idvstk, locpnt

c     for scaling the point emissions
c     nh3_scale  -- scaling factors for NH3
c     so2_scale  -- scaling factors for SO2
c     nox_scale  -- scaling factors for NOx
c     voc_scale  -- scaling factors for VOCs
c     pm25_scale -- scaling factors for PM 2.5
c     pm10_scale -- scaling factors for PM between 2.5 and 10 um

      real nh3_scale(150,162)
      real so2_scale(150,162)
      real nox_scale(150,162)
      real voc_scale(150,162)
      real pm25_scale(150,162)
      real pm10_scale(150,162)

      common /emission_scale/ nh3_scale,so2_scale,nox_scale,voc_scale,
     &                        pm25_scale,pm10_scale