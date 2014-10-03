c***********************************************************************
c
c      STARTUP - open IO files
c
c***********************************************************************
c
       subroutine startup(date)
       implicit none
c
       include 'aervaria.inc'
       include 'IO.inc'
c
       integer isblk
c
       character*6 date
       character*99 fname1, fname2, fname3, fname4
c
c   OPEN I/O UNITS 
c     
c            Order:
c            - INPUT, OUTPUT and ERROR TRACKING FILES come sequetianlly.
c            - File I/O unit number follows an alphabetic order, but 
c              sometimes old files come first. The similar files may
c              come together.
c
c-----------------------------------------------------------------------
c
c     INPUT FILES...
c
c-----------------------------------------------------------------------
c     aerdist.dat, qi, qi2 and qi21;
c      read initial number distributions [dndlog]. dat, qi and qi2*
c      correspond the date of 072701, 091202 and 090802. qi21 starts
c      at 6:30 AM. Others starts at midnight.
c
       fname1='aerdist.'//date
       isblk=INDEX(fname1,' ')
       open(unit=idat,     file=dir(1)//fname1(1:isblk-1), status='old')
cdate      open(unit=idat,     file=dir(1)//'aerdist.dat', status='old')
cdate      open(unit=idat,     file=dir(1)//'aerdist.qi', status='old')
cdate      open(unit=idat,     file=dir(1)//'aerdist.qi2', status='old')
c      open(unit=idat,     file=dir(1)//'aerdist.qi21', status='old')
c
c     met.dat, qi and qi2;
c      read temperature [degree C], relative humidity, UV [W/m2] and 
c      windspeed [m/s]. dat, qi and qi2 correspond the dates of 072701
c      , 091202 and 090802. Time interval is set as every 15 minutes.
c
       fname2='met.'//date
       isblk=INDEX(fname2,' ')
       open(unit=idat+1,     file=dir(1)//fname2(1:isblk-1), 
     & status='old')
cdate      open(unit=idat+1,     file=dir(1)//'met.dat', status='old')
cdate      open(unit=idat+1,     file=dir(1)//'met.qi', status='old')
cdate      open(unit=idat+1,	file=dir(1)//'met.qi2', status='old')
c
c     nh3.dat;
c      reads ammonia concentrations [ppt] every hour. Ammonia
c      concentration is reproduced from measured total (PM2.5 + gas) 
c      sulfate, ammonium and nitrate.
c
       if((date.eq.'091202').or.(date.eq.'090802'))goto10
         ! If date is Zhang et al.(2002), skip open NH3.
       fname3='nh3.'//date
       isblk=INDEX(fname3,' ')
       open(unit=idat+2,     file=dir(1)//fname3(1:isblk-1), 
     & status='old')
 10    continue
c      open(unit=idat+2,     file=dir(1)//'nh3.dat', status='old')
c
c     so2.dat, qi and qi2;
c      read SO2 mixing ratio [ppb]. dat, qi and qi2 correspond the
c      date of 072701, 091202 and 090802. Time interval is 10 mins.
c
       fname4='so2.'//date
       isblk=INDEX(fname4,' ')
       open(unit=idat+3,     file=dir(1)//fname4(1:isblk-1), 
     & status='old')
cdate      open(unit=idat+3,     file=dir(1)//'so2.dat', status='old')
cdate      open(unit=idat+3,     file=dir(1)//'so2.qi', status='old')
cdate      open(unit=idat+3      file=dir(1)//'so2.qi2', status='old')
c
c     nh3ratio.dat;
c      reads nh3 change ratio from begining. These ratios is provided
c      from Rob Pinder (rwp@andrew.cmu.edu). The final result was not
c      good enough to reporduce Sept. 2002 simulation. So constant
c      ammonia emission rates were used.
c
c      open(unit=idat+6,    file=dir(1)//'nh3ratio.dat', status='old')
c
c     aerdist1.cgl and cnd;
c      initial Nk for coagulation and condensation tests. These files
c      are used in initconv.f.
c
      open(unit=idat+7,     file=dir(1)//'aerdist1.cgl',status='old')
      open(unit=idat+8,     file=dir(1)//'aerdist1.cnd',status='old')
c-----------------------------------------------------------------------
c
c     OUTPUT FILES...
c
c-----------------------------------------------------------------------
c
c     1. Size resolved outputs
c
c     d_mdist.out;
c      prints mass distribution of OM, sulfate and ammonium every 15 
c      mins (dmdlog).
c
      open(unit=iout,     file=dir(2)//'d_mdist.out', status='unknown')
c
c     distn.dat; 
c      prints number distribution (dndlog). A print time interval is 
c      adjusted by kprint. The print time interval is a product of 
c      timestep (deltat) and kprint.
c
      open(unit=iout+1,     file=dir(2)//'distn.dat', status='unknown')
c
c     dndlogdp.out;
c      prints dndlogdp's which are seperated by TOMAS' Nk values. It is
c      designed to print only when icoag_test or icond_test is on. 
c      Time interval is every 15 mins.
c
      if ((icoag_test.eq.1).or.(icond_test.eq.1)) then
        open(unit=iout+2,  file=dir(2)//'dndlogdp.out', status=
     & 'unknown') 
      endif
c
c     Mk.out;
c      prints mass distribution sulfate, OM and ammonium every 15 
c      mins [ug/m3]. It is same with d_mdist.out but unit and 
c      expression are different.
c
      open(unit=iout+3,	 file=dir(2)//'Mk.out', status='unknown')
c
c
c     2. Total properties
c
c     acidity.out;
c      prints acidity of a particle whose diameter is 22 nm, ammonia
c      concentration [ppt], sulfuric acid concentratin [ppt] and
c      nucleation rate as particles cm-3. It writes every
c      timestep unlike nucrate in order to get high resolution.
c
      open(unit=iout+10,  file=dir(2)//'acidity.out',status='unknown')
c
c     ccn.out;
c      prints cn, cn20 and ccn every timestep (deltat).
c
cdel      open(unit=iout+11,  file=dir(2)//'ccn.out',status='unknown')
c
c     Ntot.out;
c      prints total number concentration [particles/cm3] below cutpoint
c      every timestep.
      open(unit=iout+12,  file=dir(2)//'Ntot.out', status='unknown') 
c
c     nucrate.out;
c      prints sulfuric acid [ppt], temperature, relative 
c      humidity ammonia mixing ratio [ppt] every 15 minutes.
c
cdel      open(unit=iout+13,  file=dir(2)//'nucrate.out', status='unknown')
c     
c     orgmass.out;
c      prints total mass concentrations of organic, sulfate, and ammonium.
c      Time interval is every 15 mins.
c
      open(unit=iout+14,   file=dir(2)//'orgmass.out', status='unknown')
c
c     s_tot.out;
c      prints total surface area [um2/cm3] every 15 mins.
c
      open(unit=iout+15,     file=dir(2)//'s_tot.out', status='unknown')
c       
c     PM25.out;
c      prints PM2.5 every 15 mins.
c
cdel      open(unit=iout+16,     file=dir(2)//'PM25.out',  status='unknown')
c
cc     sulfuric.out;
c     sulfuric_ammonia.out
c      print sulfuric acid as ppt and molecules cm-3
c
c      open(unit=iout+17,     file=dir(2)//'sulfuric.out',
      open(unit=iout+17,     file=dir(2)//'sulfuric_ammonia.out',
     &  status='unknown')
c
c     3. Only once writing
c
c     distavg.out;
c      prints Dpmean(um), dNavg/dlogDp and Navg(cm-3) at the end of 
c      simulation
c
cdel      open(unit=iout+20,    file=dir(2)//'distavg.out',status='unknown')
c-----------------------------------------------------------------------
c
c     ERROR TRACKING FILES...
c
c-----------------------------------------------------------------------
c     coag.trk;
c      tracks mass conc. change between before and after coagulation
c      process. The first and second columns are sulfuric acid and
c      ammonia conc. change. It is written every 15 mins.
c
      open(unit=itrk,     file=dir(3)//'coag.trk', status='unknown')
c
c     errmnfix.trk;
c      tracks error inside of mnfix and write in mnfix subroutine not
c      in main routine.      
c
      open(unit=itrk+1,   file=dir(3)//'errmnfix.trk', status='unknown'
     &   ,access='append')
c     
c      tracks errors every timestep in main routine. The errors are
c      mnfixdiagm, mnfixdiagn, conddiagn and conddiagm. The detail 
c      meanings are described in aervaria.inc. It is designed to check
c      whether the changes by mnfix is less enough compared to changes
c      by condensation or coagulation.
c
       open(unit=itrk+2,  file=dir(3)//'mnfixerr.trk',status='unknown') 
c
c     so4conditr.trk;
c      tracks iteration times inside of condensation subroutine every
c      timestep. It also prints reason of iteration with size bin 
c      order and specie number and increase amount by condensation as
c      loss of gas concentration [ppt]. It is designed to find optimum
c      timestep.
      open(unit=itrk+3,	file=dir(3)//'so4conditr.trk', status='unknown') 
c-----------------------------------------------------------------------
c
      return
      end
