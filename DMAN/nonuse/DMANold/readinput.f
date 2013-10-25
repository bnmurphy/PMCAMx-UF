
C     *************************************************
C     * readinput - part 1. initinput                 *
C     *************************************************

C     WRITTEN BY JaeGun Jung, November 2007

C     This subroutine consists of two parts. The first one is initinput,
C     which reads input data from input files from measurement before
C     a main loop. The second part is readinput, which read the input data
C     in the main loop. Variables are updated as new input data appear.

C-----INPUTS------------------------------------------------------------


C-----OUTPUTS-----------------------------------------------------------

C     mettime - meteorology measured time, hr
C     tempdat - measured temperature, K
C     rhdat   - measured relative humidity
C     raddat  - measured UV radiation, W m-2
C     winddat - measured wind speed, m s-1

C     so2time - SO2 measured time,hr
C     so2dat  - measured SO2, ppb

C     nh3time - NH3 time, hr
C     nh3dat  - NH3 mixing ratio, ppt

C     metterms- a number of lines in meteorology input files
C     so2terms- a number of lines in SO2 input files
C     nh3terms- a number of lines in NH3 input files

C     date    - input date
C     tend    - the end time of simulation

 
      SUBROUTINE initinput(mettime,tempdat,rhdat,raddat,winddat,so2time,
     & so2dat,nh3time,nh3dat,metterms,so2terms,nh3terms,date,tend)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'
      include 'IO.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer metterms
      integer so2terms
      integer nh3terms

      real mettime(inmax)
      real tempdat(inmax)
      real rhdat(inmax)
      real raddat(inmax)
      real winddat(inmax)

      real so2time(inmax)
      real so2dat(inmax)

      real nh3time(intime)
      real nh3dat(intime)

      real tend

      character*6 date ! simulation date

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i

C-----EXTERNAL FUNCTIONS-----------------------------------------------


C-----ADJUSTABLE PARAMETERS--------------------------------------------


C-----CODE--------------------------------------------------------------

      ! Skip read measurement data if any test flag is on.
      if ((icoag_test .eq. 1) .or. (icond_test .eq. 1)) then
         return
      endif

      ! Read meteorology data
      i = 1
      read(idat+1,*)mettime(i),tempdat(i),rhdat(i),raddat(i),winddat(i)
      do while (mettime(i) .lt. tend)
        i = i+1
        read(idat+1,*)mettime(i),tempdat(i),rhdat(i),raddat(i),
     &   winddat(i)
      enddo
      metterms = i

      ! Read SO2 mixing ratio, ppb
      i = 1
      read(idat+3,*)so2time(i),so2dat(i)
      do while (so2time(i) .lt. tend)
        i = i + 1
        read(idat+3,*)so2time(i),so2dat(i)
      enddo
      so2terms = i
    
      ! If ammonia mixing ratio is calculated from emission or out-setting,
      ! then skip reading ammonia.
      if ((date .eq. '091202') .or. (date .eq. '090802')) goto 10

      ! Read ammonia mixing ratio, ppt
      i = 1
      read(idat+2,*)nh3time(i),nh3dat(i)
      do while (nh3time(i) .lt. tend)
         i = i + 1
         read(idat+2,*)nh3time(i),nh3dat(i)
      enddo
      nh3terms = i

 10   continue


      RETURN
      END

C=======================================================================
C
C=======================================================================

C     *************************************************
C     * readinput - part 2. readinput                 *
C     *************************************************

C-----INPUTS------------------------------------------------------------

C     time    - a main loop time

C     mettime - meteorology measured time, hr
C     tempdat - measured temperature, K
C     rhdat   - measured relative humidity
C     raddat  - measured UV radiation, W m-2
C     winddat - measured wind speed, m s-1

C     so2time - SO2 measured time,hr
C     so2dat  - measured SO2, ppb

C     nh3time - NH3 time, hr
C     nh3dat  - NH3 mixing ratio, ppt

C     metterms- a number of lines in meteorology input files
C     so2terms- a number of lines in SO2 input files
C     nh3terms- a number of lines in NH3 input files

C     date    - input date

C-----OUTPUTS-----------------------------------------------------------
 
C     ygas    - mixing ratios of gas species

      SUBROUTINE readinput(time,mettime,tempdat,rhdat,raddat,winddat,
     & so2time,so2dat,nh3time,nh3dat,metterms,so2terms,nh3terms,metcount
     & ,so2count,nh3count,ygas,date)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'
      include 'IO.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer metterms
      integer so2terms
      integer nh3terms
      integer metcount
      integer so2count
      integer nh3count

      real time

      real mettime(inmax)
      real tempdat(inmax)
      real rhdat(inmax)
      real raddat(inmax)
      real winddat(inmax)

      real so2time(inmax)
      real so2dat(inmax)

      real nh3time(intime)
      real nh3dat(intime)
   
      real ygas(ngas)

      character*6 date ! simulation date

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i

C-----EXTERNAL FUNCTIONS-----------------------------------------------


C-----ADJUSTABLE PARAMETERS--------------------------------------------


C-----CODE--------------------------------------------------------------

      ! Set constant values if any test flag is on.
      if ((icoag_test .eq. 1) .or. (icond_test .eq. 1)) then
         temp = 298.15
         rh = 0.5
         rad = 0.0
         goto 10
      endif

      ! Read meteorology data
      do while ((time .ge. mettime(metcount+1)) .and. (metcount .le. 
     & metterms))
        metcount = metcount + 1
        temp = tempdat(metcount) + 273.15
        rh = rhdat(metcount) * 1.0e-2
        rad = raddat(metcount)
      enddo

      ! Read SO2 mixing ratio, ppb
      do while ((time .ge. so2time(so2count+1)) .and. (so2count .le. 
     & so2terms))
        so2count = so2count + 1
        ygas(mgso2) = so2dat(so2count) * 1.0e+3
      enddo
    
      ! If ammonia mixing ratio is calculated from emission or out-setting,
      ! then skip reading ammonia.
      if ((date .eq. '091202') .or. (date .eq. '090802')) goto 10

      ! Read ammonia mixing ratio, ppt
      do while ((time . ge. nh3time(nh3count + 1)) .and. (nh3count .lt.
     & nh3terms))
        nh3count = nh3count + 1
        ygas(mgnh3) = nh3dat(nh3count)
      enddo

 10   continue


      RETURN
      END
