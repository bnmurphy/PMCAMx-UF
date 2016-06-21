      subroutine readpt
c
c-----CAMx v4.02 030709
c
c     READPT reads the time-variant records of the point source file and
c     cycles through to current time/date to load point emission rates
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c     05/01/03  Time span of emissions must now match emiss update interval
c
c     Input arguments:
c        none
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
      include 'filunit.com'
      include 'ptemiss.com'
      include 'chmstry.com'
      include 'flags.com'
      include 'grid.com'
c
      character*4 ptspec(10)
      real flowrat(MXPTSRC)
c
      !indexes for bin 18 of various aerosol species
      integer kptno3_18
      integer kptpoc_18
      integer kptpec_18
      integer kptcrst_18
      integer kptso4_18
      integer lpt
c
      data pi /3.1415927/
c
c-----Entry point
c
      !set aero species pointers
      do lpt=1,nptspc
         if (spname(lptmap(lpt)).eq.'PNO3_18   ') kptno3_18=lpt
         if (spname(lptmap(lpt)).eq.'POC_18    ') kptpoc_18=lpt
         if (spname(lptmap(lpt)).eq.'PEC_18    ') kptpec_18=lpt
         if (spname(lptmap(lpt)).eq.'CRST_18   ') kptcrst_18=lpt
         if (spname(lptmap(lpt)).eq.'PSO4_18   ') kptso4_18=lpt
      end do
c
      kount = 1
 100  continue
      read(iptem,end=900) idat1,tim1,idat2,tim2 
      ichktm1 = NINT( 1000*(tim1) )
      ichktm2 = NINT( 1000*(tim2) )
      if( NINT(tim2) .EQ. 0 ) ichktm2 = 24000
      ichkems = NINT( 1000*(dtems/60.) )
      if( (ichktm2 - ichktm1) .NE. ichkems ) then
          write(iout,'(//,a)')'ERROR in READPT:'
          write(iout,*) 'Time interval in point emissions file does'
          write(iout,*)  ' not match emissions update time interval.'
          write(iout,*) '   Beginning Date/Time (HHMM): ',idat1,tim1
          write(iout,*) '   Ending Date/Time    (HHMM): ',idat2,tim2
          write(iout,*) '   Emiss Input interval (min): ',dtems
          call camxerr()
      endif
      if(NINT(tim2) .EQ. 0) then
        tim2 = 24.
        idat2 = idat2 - 1
      endif
      tim1 = 100.*tim1 
      tim2 = 100.*tim2 
      read(iptem) idum,npts 
      if (npts.ne.nptsrc) then 
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,*)'Number of points read from point source file',
     &               ' differs from header value'
        write(iout,*)'Record: ',npts,' Header: ',nptsrc 
        call camxerr()
      endif 
      read (iptem) (idum,idum,idum,flowrat(n),effph(n),n=1,npts)
      do ll = 1,nptspc 
        read(iptem) idum,(ptspec(i),i=1,10),(ptemis(n,ll),n=1,npts) 
      enddo
      write(iout,'(a40,2(f7.0,i8.5))')
     &       'Read point source file at ',tim1,idat1,tim2,idat2
c
      !scaling emissions
      !aerosol species
      do n=1,npts
         if (ptemis(n,kptno3_18).gt.0.) cycle !fire source
         do lpt=kptpec_18,kptpec_18+17 !bins 18-35
            ptemis(n,lpt)=pm25_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
         end do
         do lpt=kptpoc_18,kptpoc_18+17 !bins 18-35
            ptemis(n,lpt)=pm25_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
         end do
         do lpt=kptcrst_18,kptcrst_18+13 !bins 18-31 (highest anthro crustal matter bin)
            ptemis(n,lpt)=pm25_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
         end do
         do lpt=kptso4_18,kptso4_18+17 !bins 18-35
            ptemis(n,lpt)=pm25_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
         end do
         do lpt=kptpec_18+18,kptpec_18+23 !bins 36-41
            ptemis(n,lpt)=pm10_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
         end do
         do lpt=kptpoc_18+18,kptpoc_18+23 !bins 36-41
            ptemis(n,lpt)=pm10_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
         end do
      end do
      !gas species
      do lpt=1,nptspc
         do n=1,npts
            if (ptemis(n,kptno3_18).gt.0.) cycle !fire source
            if (spname(lptmap(lpt)).eq.'NH3       ') then
               ptemis(n,lpt)=nh3_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
            else if (spname(lptmap(lpt)).eq.'SO2       ') then
               ptemis(n,lpt)=so2_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
            else if (spname(lptmap(lpt)).eq.'NO        ') then
               ptemis(n,lpt)=nox_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
            else if (spname(lptmap(lpt)).eq.'NO2       ') then
               ptemis(n,lpt)=nox_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
            else if (spname(lptmap(lpt)).eq.'ALK1      '.or.
     &              spname(lptmap(lpt)).eq.'ALK2      '.or.
     &              spname(lptmap(lpt)).eq.'ALK3      '.or.
     &              spname(lptmap(lpt)).eq.'ALK4      '.or.
     &              spname(lptmap(lpt)).eq.'ALK5      '.or.
     &              spname(lptmap(lpt)).eq.'OLE1      '.or.
     &              spname(lptmap(lpt)).eq.'OLE2      '.or.
     &              spname(lptmap(lpt)).eq.'ARO1      '.or.
     &              spname(lptmap(lpt)).eq.'ARO2      '.or.
     &              spname(lptmap(lpt)).eq.'ETHE      '.or.
     &              spname(lptmap(lpt)).eq.'MEK       '.or.
     &              spname(lptmap(lpt)).eq.'HCHO      '.or.
     &              spname(lptmap(lpt)).eq.'CCHO      '.or.
     &              spname(lptmap(lpt)).eq.'RCHO      '.or.
     &              spname(lptmap(lpt)).eq.'MEOH      '.or.
     &              spname(lptmap(lpt)).eq.'PROD      ') then
               ptemis(n,lpt)=voc_scale(isrc(n,1),jsrc(n,1))*ptemis(n,lpt)
            end if
         end do
      end do
c
c-----Check times only if LE1DAY = T, otherwise check both time and date
c
      if (le1day) then
        if (tim1.le.time .and. tim2.gt.time) goto 200
        if (tim1.gt.time) goto 900
      else
        if ((idat1.lt.date .or. (idat1.eq.date .and. tim1.le.time))
     &    .and. (idat2.gt.date .or. (idat2.eq.date .and. tim2.gt.time)))
     &    goto 200
      endif
      goto 100
c 
c-----Convert emission rates from moles/(dtems-hours) to moles/s for gas 
c     or g/(dtems-hours) to g/s for aero species 
c 
 200  do 10 l = 1,nptspc
        do n = 1,npts 
          ptemis(n,l) = ptemis(n,l)/(60.*dtems) 
        enddo 
 10   continue 
c 
c-----Convert flow rate to new exit velocity in m/s  
c 
      do n = 1,npts 
        if (flowrat(n).gt.0. .AND. dstk(n) .NE. 0. )
     &    vstk(n) = flowrat(n)/(3600.*pi*(abs(dstk(n))/2)**2) 
      enddo 
      goto 999
c
c-----End of file reached; if 1-day emissions requested, rewind and read 
c     through header once more.  Otherwise, report error and stop
c
 900  continue
      if (le1day) then
        if (kount.ge.2) then
          write(iout,'(//,a)')'ERROR in READPT:'
          write(iout,*)'Cannot match model time with point source time.'
          call camxerr()
        endif
        rewind(iptem)
        read(iptem) idum 
        read(iptem) dum  
        read(iptem) idum  
        read(iptem) idum 
        read(iptem) idum 
        read(iptem) dum
        kount = kount + 1
        goto 100
      else
        write(iout,'(//,a)')'ERROR in READPT:'
        write(iout,*)'End of point source file reached.  Make sure the '
        write(iout,*)
     &            'file is for the correct day and contains all hours.'
        call camxerr()
      endif
c
 999  return
      end
