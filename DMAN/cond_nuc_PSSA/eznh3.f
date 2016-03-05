
C     **************************************************
C     *  eznh3                                         *
C     **************************************************

C     WRITTEN BY JaeGun Jung, March 2008

C     This is a subroutine calling eznh3eqm 

C-----INPUTS------------------------------------------------------------

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE eznh3(Gci,Mki,Nki,Gcf,Mkf,Nkf,ichm,jchm,kchm)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Gci(icomp)
      double precision Mki(ibins,icomp)
      double precision Nki(ibins)
      double precision Gcf(icomp)
      double precision Mkf(ibins,icomp)
      double precision Nkf(ibins)
      integer ichm, jchm, kchm

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer k,j
      double precision tot_nh3  !total kmoles of ammonia
      double precision tot_so4  !total kmoles of so4
      double precision sfrac    !fraction of sulfate that is in that bin

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      call eznh3eqm(Gci,Mki)
      call mnfix_PSSA(Nki,Mki,ichm,jchm,kchm,3)
          ! adjust average mass in each size bin within boundaries

      do k=1, ibins
         do j=1, icomp
            Mkf(k,j)=Mki(k,j)
         enddo
         Nkf(k)=Nki(k)
      enddo
      Gcf(srtso4)=Gci(srtso4)
      Gcf(srtnh4)=Gci(srtnh4)


      RETURN
      END












