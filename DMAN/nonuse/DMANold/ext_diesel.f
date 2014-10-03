
C     *************************************************
C     * ext_diesel                                    *
C     *************************************************

C     RE-WRITTEN BY JaeGun Jung, November 2007

C     This subroutine is designed for external mixing of 
C     diese. First calculate total amount of number to
C     be distributed log-normally, then distribute by calling
C     logdist. The amount of number of particle in each size 
C     bin is decided in the subroutine. traf variables are
C     the values which will be added particle distribution.

C-----INPUTS------------------------------------------------------------

C     Initial variables
C     =================

C     pn(nsect) - number concentrations as particles per cubic cm
C     t - time of a main routine
C     traf - the amount of particle added in each section

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE ext_diesel(pn, t, traf)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      real pn(nsect) ! [=] particles per cm3
      real t ! time of a main routine
      real traf(nsect) ! [=] particles per cm3

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i

      real traffic
      real tnumb
      real dptr
      real strfc

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

       if(nsect.eq.1) then
         write(6,*) 'monodisperse particles'
         return
       endif
c
        if (iwknd.eq.0) then
c          if (t.le.6) then
c            tnumb=4000
          if(t.le.9) then
            tnumb=2.0d3+9.0d3*(t-3.0d0)/6.0d0
          elseif(t.le.15) then
            tnumb=5.5d3+(15.-t)/9.*5.5d3
c          elseif(t.le.22) then
c            tnumb=6000
          else
            tnumb=2.0d3+(24.0-t)/9.0d0*3.5d3
          endif
        else
          if (t.le.6) then
            tnumb=3000
          elseif(t.le.10) then
            tnumb=5000
          elseif(t.le.16) then
            tnumb=8000
          elseif(t.le.22) then
            tnumb=5000
          else
            tnumb=3000
          endif
        endif

C       Add traffic signal as external particles if iemis2=1
C       Can also add external traffic signal in Matlab script
        traffic = tnumb*FLOAT(iemis2)

        dptr = 0.01035                     ! mean diameter in um (was 0.7)
        strfc = 0.6                      ! std
        call logdist(traffic,dptr,strfc,traf)
       do i=1,nsect



       RETURN
       END
