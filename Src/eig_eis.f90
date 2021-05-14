!
!*******************************************************************************
!
subroutine tri_eig(nm,n,d,e,z,ierror)
!
!*******************************************************************************
!
!  Discussion:
!
!---  This subroutine is a rewrite and update of the eispack routine tql2 
!---  Revised 31 August 2020 to remove obsolete FORTRAN and to make
!---  compatible with FORTRAN 90 
!
!----------    Original EISPACK coments   --------------------------------
!
!     this subroutine is a translation of the algol procedure tql2,
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
!     wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierror-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierror is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag_90 for  dsqrt(a*a + b*b) .
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated august 1983.
!--------------     End of original EISPACK comments   ------------------
!
!  Date:
!
!    31 August 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   implicit none
   integer(kind=4), intent(in) :: nm, n
   real(kind=8), intent(out) :: d(n), e(n), z(nm,n)
   integer(kind=4), intent(out) :: ierror
!-----------------------------------------------------------------------------
   integer(kind=4) :: i, j, k, l, m, ii, l1, l2, mml
   real(kind=8) :: c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1, tst2
!
!-------------    External functions
!
   real(kind=8) :: pythag_90
!
!*******************************************************************************
!
   ierror = 0
   if(n == 1)return   !   Nothing to be done if n = 1
!
   do i = 2, n
      e(i-1) = e(i)
   end do
!
   f = 0.0d0
   tst1 = 0.0d0
   e(n) = 0.0d0
!
   do l = 1, n
      j = 0
      h = abs(d(l)) + abs(e(l))
      if(tst1 < h) tst1 = h
!---   Look for small sub-diagonal element
      do m = l, n
         tst2 = tst1 + abs(e(m))
!         if (tst2 == tst1)exit
         if(abs(tst2 - tst1) < 1.0d-8)exit
      end do
!---   An eigenvalue is found
      if(m == l)then
         d(l) = d(l) + f
         cycle
      end if
!---   Do up to 30 iterations
      do j = 1, 30
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g)/(2.0d0*e(l))
         r = pythag_90(p,1.0d0)
         d(l) = e(l)/(p + sign(r,p))
         d(l1) = e(l)*(p + sign(r,p))
         dl1 = d(l1)
         h = g - d(l)
!
         if(l2 <= n)then
            do i = l2, n
               d(i) = d(i) - h
            end do
        end if
!
        f = f + h
!
!---   QL transformation
!
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
         do ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c*e(i)
            h = c*p
            r = pythag_90(p,e(i))
            e(i+1) = s*r
            s = e(i) / r
            c = p / r
            p = c*d(i) - s*g
            d(i+1) = h + s*(c*g + s*d(i))
!
!---   Form eigenvectors
!
            do k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s*z(k,i) + c*h
               z(k,i) = c*z(k,i) - s*h
            end do
         end do
!
         p = -s*s2*c3*el1*e(l)/dl1
         e(l) = s*p
         d(l) = c*p
         tst2 = tst1 + abs(e(l))
         if(tst2 <= tst1)exit
      end do
!
!----   Too many iterations and failed to find eigenvalues
!----   return with error flag set
!
      if(j == 30)then
         ierror = 1
         return
      end if
      d(l) = d(l) + f
   end do
!
!---   Order eigenvalues and eigenvectors
!
   do ii = 2, n
      i = ii - 1
      k = i
      p = d(i)
!
      do j = ii, n
         if (d(j) < p)then
            k = j
            p = d(j)
         end if
      end do
!
      if(k == i) cycle
      d(k) = d(i)
      d(i) = p
      do j = 1, n
         p = z(j,i)
         z(j,i) = z(j,k)
         z(j,k) = p
      end do
   end do
   return
end subroutine tri_eig
!
!*******************************************************************************
!
real(kind=8) function pythag_90(a,b)
   implicit none
   real(kind=8), intent(in) :: a, b
!
!*******************************************************************************
!
!  Discussion:
!
!----   Conversion of EISPACK routine pythag to FORTRAN 90
!----   Removed outdated FORTRAN and updated to be compatible 
!----   with FORTRAN 90
!
!
!--------------     Original EISPACK comments   ------------------
!
!     finds dsqrt(a**2+b**2) without overflow or destructive underflow
!
!--------------     End of original EISPACK comments   ------------------
!
!
!  Date:
!
!    31 August 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   real(kind=8) :: p, r, s, t, u
   pythag_90 = 0.0d0
   p = max(abs(a),abs(b))
   if(abs(p - 0.0d0) < 1.0d-30) return
   r = (min(abs(a),abs(b))/p)**2
   t = 4.0d0 + r
   do while(abs(t - 4.0d0) >= 1.0d-8)
      s = r/t
      u = 1.0d0 + 2.0d0*s
      p = u*p
      r = (s/u)**2*r
      t = 4.0d0 + r
   end do
   pythag_90 = p 
   return
end function pythag_90
