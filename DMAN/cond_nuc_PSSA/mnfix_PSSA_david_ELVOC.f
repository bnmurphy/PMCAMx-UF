
C     **************************************************
C     *  mnfix_PSSA                                    *
C     **************************************************

C     WRITTEN BY Peter Adams, September 2000

C     This subroutine examines the mass and number distributions and
C     determines if any bins have an average mass outside their normal
C     range.  I have seen this happen because the GCM advection seems
C     to treat the mass and number tracers inconsistently.  If any bins
C     are out of range, I shift some mass and number to a new bin in
C     a way that conserves both.

C     re-written by Jeff Pierce, August 2007

C-----INPUTS------------------------------------------------------------

C     Nkx and Mkx are the number and mass distributions

C-----OUTPUTS-----------------------------------------------------------

      SUBROUTINE mnfix_PSSA(Nkx,Mkx,ichm,jchm,kchm,pq)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

      include 'sizecode.COM'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      double precision Nkx(ibins), Mkx(ibins,icomp)
      integer ichm ! i coordinate in PMCAMx
      integer jchm ! j coordinate in PMCAMx
      integer kchm ! k coordinate in PMCAMx
      integer pq !flaqs error 

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer k,j,L,jj        !counters
      integer is, js          !for printing Nkx, Mkx by jgj 01/27/08
      double precision tot_mass ! total dry mass
      double precision xk_hi, xk_lo ! geometric mean mass of bins that mass is moving to
      double precision xk_hi1, xk_hi2 ! used as tests to see if mass is way to low or high
      double precision Neps,Meps
      double precision tmpvar ! temporary variable
      double precision frac_lo_n, frac_lo_m ! fraction of the number and mass going to the lower bin
      double precision Mktot !To calculate fraction of negative mass

C     VARIABLE COMMENTS...

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter(Neps=1.0d-10,Meps=1.0d-32)

C-----CODE--------------------------------------------------------------

      !Check variables
cdbg      print*,'at the beginning of mnfix_PSSA'
cdbg      print*,'xk'
cdbg      do k=1,ibins
cdbg         print*,xk(k)
cdbg      enddo

cdbg      print*,'Nkx'
cdbg      do k=1,ibins
cdbg         print*,Nkx(k)
cdbg      enddo

cdbg      print*,'Mkx'
cdbg      do j=1,icomp
cdbg         print*,'j=',j
cdbg         do k=1,ibins
cdbg            print*,Mkx(k,j)
cdbg         enddo
cdbg      enddo
cdbg      pause

cdbg      write(*,*)'mnfix - chkpt 1'

      ! first check for negative tracers
      do k=1,ibins
         if (Nkx(k).lt.Neps)then
            if (Nkx(k).gt.-1.0d10)then
               Nkx(k)=Neps
               do j=1,icomp
                  if (j.eq.srtso4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5*0.727273
                  elseif (j.eq.srtinrt)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.   
                  elseif (j.eq.srtsoa1)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa2)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa3)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa5)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.   
                  elseif (j.eq.srtnh4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5*(1.0-0.727273)
                  else
                     Mkx(k,j)=0.d0      
                  endif
               enddo
            else
               print*,'Negative tracer in mnfix',pq
               print*,'Coord=',ichm,jchm,kchm
               print*,'Nkx='
               do is=1, ibins
                  write(*,*)Nkx(is)
               enddo
               print*,'Mkx='
               do js=1, icomp
                  write(*,*)'species=',js
                  do is=1, ibins
                     write(*,*)Mkx(is,js)
                  enddo
               enddo
               STOP
            endif
         endif
         tot_mass=0.d0
         do j=1,icomp-idiag
            tot_mass=tot_mass+Mkx(k,j)
         enddo
         if (tot_mass.lt.Meps) then
            Nkx(k)=Neps
            do jj=1,icomp
                  if (jj.eq.srtso4)then
                     Mkx(k,jj)=Neps*1.4*xk(k)*0.5*0.727273
                  elseif (jj.eq.srtinrt)then
                     Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.   
                  elseif (jj.eq.srtsoa1)then
                     Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                  elseif (jj.eq.srtsoa2)then
                     Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                  elseif (jj.eq.srtsoa3)then
                     Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                  elseif (jj.eq.srtsoa4)then
                     Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                  elseif (jj.eq.srtsoa5)then
                     Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.   
                  elseif (jj.eq.srtnh4)then
                     Mkx(k,jj)=Neps*1.4*xk(k)*0.5*(1.0-0.727273)
                  else
                      Mkx(k,jj)=0.d0      
                  endif
            enddo
         endif            
cdbg         write(*,*)'mnfix - chkpt 2'
         do j=1,icomp-1 !-idiag
            if (Mkx(k,j).lt.0.d0)then
               if (Mkx(k,j).gt.-1.0d0)then
                  Nkx(k)=Neps
                  do jj=1,icomp
                   if (jj.eq.srtso4)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5*0.727273
                   elseif (jj.eq.srtinrt)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.   
                   elseif (jj.eq.srtsoa1)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                   elseif (jj.eq.srtsoa2)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                   elseif (jj.eq.srtsoa3)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                   elseif (jj.eq.srtsoa4)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                   elseif (jj.eq.srtsoa5)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.   
                   elseif (jj.eq.srtnh4)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5*(1.0-0.727273)
                   else
                      Mkx(k,jj)=0.d0      
                   endif
                  enddo
               else
                  write(*,*)'Negative tracer in mnfix',pq
                  write(*,*)'Coord=',ichm,jchm,kchm
                  write(*,*)'Nkx='
                  do is=1, ibins
                     write(*,*)Nkx(is)
                  enddo
                  write(*,*)'Mkx='
                  do js=1, icomp
                     write(*,*)'species=',js
                     do is=1, ibins
                        write(*,*)Mkx(is,js)
                     enddo
                  enddo
c
cdbg                  write(*,*)'mnfix - chkpt 3'
                  Mktot=0.d0
                  do jj=1,icomp-idiag
                     Mktot=Mktot+Mkx(k,jj)
                  enddo
cdbg                  write(*,*)'mnfix - chkpt 4'
c
                  !If the negative is less than 1.0% of total mass without
                  !H2O, take an eps treatment.
                  if (abs((Mkx(k,j)/Mktot)).lt.1.0d-2) then
                     Nkx(k)=Neps
                     do jj=1,icomp
                     if (jj.eq.srtso4)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5*0.727273
                     elseif (jj.eq.srtinrt)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.   
                     elseif (jj.eq.srtsoa1)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                     elseif (jj.eq.srtsoa2)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                     elseif (jj.eq.srtsoa3)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                     elseif (jj.eq.srtsoa4)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6.
                     elseif (jj.eq.srtsoa5)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5/6. 
                     elseif (jj.eq.srtnh4)then
                      Mkx(k,jj)=Neps*1.4*xk(k)*0.5*(1.0-0.727273)
                     else
                      Mkx(k,jj)=0.d0        
                     endif
                     enddo
                  else
                     write(*,*)'bin=',k,'species=',j
                     STOP 'Negative tracer takes more than 1.0%'
                  endif
               endif
            endif
         enddo
      enddo

      do k=1,ibins
         tot_mass=0.0d0
         do j=1,icomp-idiag
            tot_mass=tot_mass+Mkx(k,j)
         enddo
cJJ         if (tot_mass/Nkx(k).gt.xk(k+1).or.
cJJ     &        tot_mass/Nkx(k).lt.xk(k)) then
         if (tot_mass.gt.xk(k+1)*Nkx(k).or.
     &        tot_mass.lt.xk(k)*Nkx(k)) then
c            print*,'out of bounts in mnfix, fixing'
c            print*,k,'AVG',tot_mass/Nkx(k),'lo',xk(k),'hi',xk(k+1)
            ! figure out which bins to redistribute to
            xk_hi1 = sqrt(xk(2)*xk(1))
            xk_hi2 = sqrt(xk(ibins+1)*xk(ibins))
            if (xk_hi1.gt.tot_mass/Nkx(k)) then
               !mass per particle very low
               !conserve mass at expense of number
               tmpvar = Nkx(k)
               Nkx(k)=Neps
               Nkx(1)=Nkx(1)+tot_mass/sqrt(xk(2)*xk(1))
               do j=1,icomp
                  tmpvar=Mkx(k,j)
                  if (j.eq.srtso4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5*0.727273
                  elseif (j.eq.srtinrt)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.   
                  elseif (j.eq.srtsoa1)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa2)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa3)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa5)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.   
                  elseif (j.eq.srtnh4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5*(1.0-0.727273)
                  else
                     Mkx(k,j)=0.d0      
                  endif

                  Mkx(1,j)=Mkx(1,j)+tmpvar
               enddo
            elseif (xk_hi2.lt.tot_mass/Nkx(k)) then
               !mass per particle very high
               !conserve mass at expernse of number
               Nkx(k)=Neps
               Nkx(ibins)=Nkx(ibins)+tot_mass/sqrt(xk(ibins+1)*
     &                    xk(ibins))
               do j=1,icomp
                  tmpvar=Mkx(k,j)

                  if (j.eq.srtso4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5*0.727273
                  elseif (j.eq.srtinrt)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.   
                  elseif (j.eq.srtsoa1)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa2)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa3)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa5)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6. 
                  elseif (j.eq.srtnh4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5*(1.0-0.727273)
                  else
                     Mkx(k,j)=0.d0      
                  endif

                  Mkx(ibins,j)=Mkx(ibins,j)+tmpvar
               enddo               
            else ! mass of particle is somewhere within the bins
               L = 2
               xk_hi = sqrt(xk(L+1)*xk(L))
cJJ              do while (xk_hi .lt. tot_mass/Nkx(k))
               do while (xk_hi*Nkx(k) .lt. tot_mass)
                  L=L+1
                  xk_hi = sqrt(xk(L+1)*xk(L))
               enddo
               xk_lo = sqrt(xk(L)*xk(L-1))
                         ! figure out how much of the number to distribute to the lower bin
               frac_lo_n = (tot_mass - Nkx(k)*xk_hi)/
     &              (Nkx(k)*(xk_lo-xk_hi))
               frac_lo_m = frac_lo_n*Nkx(k)*xk_lo/tot_mass
c               print*,'frac_lo_n',frac_lo_n
c               print*,'frac_lo_m',frac_lo_m
               tmpvar = Nkx(k)
               Nkx(k) = Neps
               Nkx(L-1) = Nkx(L-1) + frac_lo_n*tmpvar
               Nkx(L) = Nkx(L) + (1-frac_lo_n)*tmpvar
               do j=1,icomp
                  tmpvar = Mkx(k,j)

                  if (j.eq.srtso4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5*0.727273
                  elseif (j.eq.srtinrt)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.   
                  elseif (j.eq.srtsoa1)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa2)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa3)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtsoa5)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5/6.
                  elseif (j.eq.srtnh4)then
                     Mkx(k,j)=Neps*1.4*xk(k)*0.5*(1.0-0.727273)
                  else
                     Mkx(k,j)=0.d0      
                  endif

                  Mkx(L-1,j) = Mkx(L-1,j) + frac_lo_m*tmpvar
                  Mkx(L,j) = Mkx(L,j) + (1-frac_lo_m)*tmpvar
               enddo
               tot_mass=0.0d0
               do j=1,icomp-idiag
                  tot_mass=tot_mass+Mkx(k,j)
               enddo
            endif 
         endif
      enddo

      do k=1,ibins
         tot_mass=0.0d0
         do j=1,icomp-idiag
            tot_mass=tot_mass+Mkx(k,j)
         enddo
cJJ         if (tot_mass/Nkx(k).gt.xk(k+1).or.
cJJ     &        tot_mass/Nkx(k).lt.xk(k)) then
         if (tot_mass.gt.xk(k+1)*Nkx(k).or.
     &        tot_mass.lt.xk(k)*Nkx(k)) then
            print*,'ERROR in mnfix'
            print*,'bin',k,'lo',xk(k),'hi',xk(k+1),'avg',tot_mass/Nkx(k)
            print*,'tot_mass',tot_mass
            print*,'Nkx(k)',Nkx(k)
            print*,'Mkx'
            do j=1,icomp
               print*,Mkx(k,j)
            enddo
            STOP
         endif
      enddo

      !Check variables
cdbg      print*,'at the end of mnfix_PSSA'
cdbg      print*,'xk'
cdbg      do k=1,ibins
cdbg         print*,xk(k)
cdbg      enddo

cdbg      print*,'Nkx'
cdbg      do k=1,ibins
cdbg         print*,Nkx(k)
cdbg      enddo

cdbg      print*,'Mkx'
cdbg      do j=1,icomp
cdbg         print*,'j=',j
cdbg         do k=1,ibins
cdbg            print*,Mkx(k,j)
cdbg         enddo
cdbg      enddo
cdbg      pause



      return
      end
