
C     *************************************************
C     * lubksb                                        *
C     *************************************************

C     WRITTEN BY JaeGun Jung, November 2007

C     This subroutine solves the set of n linear equations A * X = B.
C     Here a is input, not as the matrix A but rather as its LU decom-
C     position, determined by the routine ludcmp. indx is input as the
C     permutation vector returned by ludcmp. b(n) is input as the right
C     hand side vector B, and returns with the solution vector X. a,
C     n, np, and indx are not modified by this routine and can be left
C     in place for successive calls with different right-hand sides b.
C     This routine takes into account the possibility that b will begin
C     with many zero elements, so it is efficient for use in matrix 
C     inversion.

C-----INPUTS------------------------------------------------------------     

C     (variable) - (description)
C     (variable) - (description)

C-----OUTPUTS-----------------------------------------------------------

C     (variable) - (description)
C     (variable) - (description)

      SUBROUTINE lubksb(a,n,np,indx,b)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

cdbg      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer n,np,indx(n)
      real a(np,np), b(n)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer i,ii,j,ll
      real sum

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

C-----CODE--------------------------------------------------------------

      !Initialization
      ii=0             !When ii is set to a positive value, it will become
                       !the index of the first nonvanishing element of b.
      !Solve equations !We now do the forward substitution, equation
      do i=1, n        !(2.3.6). The only new wrinkle is to unscramble
        ll=indx(i)     !the permutation as we go.
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0) then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        elseif (sum.ne.0) then
          ii=i         !A non zero element was encountered, so from now on
        endif          !we will have to do the sums in the loop above.
        b(i)=sum
      enddo
      do i=n, 1, -1    !Now we do the backsubstitution, equation (2.3.7).
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,j)!Store a component of the solution vector X.
      enddo
  




      RETURN           !All done~
      END
