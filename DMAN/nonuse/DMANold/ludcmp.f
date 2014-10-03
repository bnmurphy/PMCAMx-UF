
C     *************************************************
C     * ludcmp                                        *
C     *************************************************

C     WRITTEN BY JaeGun Jung, November 2007

C     This subroutine replaces a given matrix a(np,np) by the LU decomp-
C     osition of a rowwise permutation of itself with physical dimension 
C     np by np. a and n are input. a is output, arranged as in equation
C     (2.3.14) in ref.; indx(n) is an output vector that records the row 
C     permutation effected by the partial pivoting; d is output as (+/-) 1 
C     depending on whether the number of row interchanges was even or odd,
C     respectively. This routine is used in combination with lubksb to 
C     solve linear equations or invert a matrix.

C     Referrence: Chapter 2.3 in Numerical Recipes in Fortran 77

C-----INPUTS------------------------------------------------------------     

C     (variable) - (description)
C     (variable) - (description)

C-----OUTPUTS-----------------------------------------------------------

C     (variable) - (description)
C     (variable) - (description)

      SUBROUTINE ludcmp(a,n,np,indx,d)

      IMPLICIT NONE

C-----INCLUDE FILES-----------------------------------------------------

cdbg      include 'aervaria.inc'

C-----ARGUMENT DECLARATIONS---------------------------------------------

      integer n, np, indx(n)
      real d, a(np,np)

C-----VARIABLE DECLARATIONS---------------------------------------------

      integer NMAX
      integer i,imax,j,k
      real TINY
      real aamax, dum, sum

      parameter (NMAX=500)     !Largest expected n

      real vv(NMAX)            !vv stores the implicit scaling of each row.

C-----EXTERNAL FUNCTIONS------------------------------------------------

C-----ADJUSTABLE PARAMETERS---------------------------------------------

      parameter (TINY=1.0e-20) !a small number

C-----CODE--------------------------------------------------------------

      !Initialization
      d=1                              !No row interchanges yet
      do i=1, n                        !Loop over rows to get the implicit
        aamax=0.                       !scaling information.
        do j=1, n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
                                       !No nonzero largest element.
        vv(i)=1./aamax                 !Save the scaling.
      enddo

      !Do LU Decomposition
      do j=1, n           !This is the loop over columns of Crout's method.
        do i=1, j-1       !This is equation (2.3.12) except for i=j
          sum=a(i,j)
          do k=1, i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0.          !Initialize for the search for largest pivot element.
        do i=j, n         !This is i=j of equation (2.3.12) and i=j+1...N
          sum=a(i,j)      !of equation (2.3.13),
          do k=1, j-1
            sum=sum-a(i,j)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*abs(sum)           !Figure of merit for the pivot.
          if (dum.ge.aamax) then       !Is it better than the best so far?
            imax=i
            aamax=dum
          endif
        enddo
        if (j.ne.imax)then             !Do we need to interchange rows?
          do k=1, n                    !Yes, do so...
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(i,k)=dum
          enddo
          d=-d                         !...and change the parity of d.
          vv(imax)=vv(j)               !Also interchange the scale factor.
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
          !If the pivot element is zero the matrix is singular (at least
          !to the precision of the algorithm). For some applications on
          !singular matrices, it is desirable to substitute TINY to zero.
        if (j.ne.n) then   !Now, finally, divide by the pivot element.
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)+dum
          enddo
        endif
      enddo                !Go back for the next column in the reduction.





      RETURN
      END
