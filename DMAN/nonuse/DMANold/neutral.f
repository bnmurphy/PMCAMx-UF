
C     **************************************************
C     * neutral                                        *
C     **************************************************
     
C     WRITTEN BY JaeGun Jung, November 2007

C     This subroutine is called by main routine not to make basic 
C     particles, because currently condensation scheme does not
C     allow basic particles. Most cases belong to small number (less
C     than frac*) of moderate (65% - 3.0%) basic particles

C                         *frac is currently set as 1.0E-15 11/25/07

C-----INPUTS------------------------------------------------------------

C     iflag - a flag shows where this subroutine is called after
C           - 1st condensation = 1, 1st coagulation = 2
C           -       nucleation = 3, 2nd coagulation = 4
C           - 2nd condensation = 5

C-----OUTPUTS-----------------------------------------------------------

C     Mk(ibins,icomp) - mass concentration

      SUBROUTINE neutral(iflag)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer iflag

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,j
      real ratio
      real Neps
      real eps
      real frac
      double precision Nktot ! total number concentration

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (Neps = 1.0e-20)
cdbg      parameter (eps  = 3.0e-2)  ! 1.5% of neutral ratio, 2
      parameter (eps  = 6.0e-2)  ! 3.0% of neutral ratio, 2
      parameter (frac = 1.0e-7) ! allowable number fraction to pass 

C-----CODE--------------------------------------------------------------

      !Calculate total number concentration
      do i=1, ibins
        Nktot=Nktot+Nk(i)
      enddo

C-----Check acidity of particles
      do i=1,ibins
        if (Nk(i).gt.Neps)then
          ratio=Mk(i,srtnh3)/Mk(i,srtso4)*(2.0/0.375) ! ratio=2 is neutral

          !Molar ratio of ammonium to sulfate is bigger than 3.3
          if (ratio .gt. 3.3) then 
            if (Nk(i)/Nktot.gt.frac)then ! Nk is bigger than frac of Ntot
              write(*,*)'Ammonium takes more than 65 % of its neutral
& value and its number on the size is bigger than frac of total number'
              write(*,*)'Over neutralization,iflag',iflag
              write(*,*)'Section=',i,'Molar ratio=',ratio
              write(*,*)'Number frac',Nk(i)/Nktot,' of the section',i
              STOP
            else
              write(*,*)'Ammonium takes more than 65 %, but correspond-
&ing Nk consists less than or equal to frac of total number conc.'
              write(*,*)'Over neutralization,iflag',iflag
              write(*,*)'Section=',i,'Molar ratio=',ratio
              write(*,*)'Number frac',Nk(i)/Nktot,' of the section',i
              Mk(i,srtnh3)=0.375*Mk(i,srtso4)
            endif

          !Molar ratio of ammonium to sulfate is less than or equal to 3.3
          else
            !Number is bigger than frac
            if (Nk(i)/Nktot.gt.frac)then ! Nk is bigger than frac of Ntot
              if (ratio .gt. 2+eps)then ! ratio is bigger than frac
                write(*,*)'Ammonium takes more than 3.0 % of its neutral
& value and less than 65%, and corresponding Nk consists more than frac 
& of total number conc'
                write(*,*)'Over neutralization,iflag',iflag
                write(*,*)'Section=',i,'Molar ratio=',ratio
                write(*,*)'Number frac',Nk(i)/Nktot,' of the section',i
cdbg                pause !debug
                Mk(i,srtnh3)=0.375*Mk(i,srtso4)
              else
                if (ratio .gt. 2.0) then
                  Mk(i,srtnh3)=0.375*Mk(i,srtso4) ! Truncation error
                endif
              endif

            !Number is less than or equal to frac
            else
              if (ratio .gt. 2+eps)then
                Mk(i,srtnh3)=0.375*Mk(i,srtso4) ! Probably truncation error 
                                   ! Small number of moderate basic particles 
              else
                if (ratio .gt. 2.0) then
                  Mk(i,srtnh3)=0.375*Mk(i,srtso4) ! Truncation error
                endif
              endif
            endif
          endif
        endif
      enddo


      
      RETURN
      END
