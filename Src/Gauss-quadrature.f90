
subroutine gauss_quad(nt,kind,alpha,beta,x,w)
!
!*******************************************************************************
!
!  Discussion:
!
!  FORTRAN 90 subroutine to calculate abscissas and weights for Gauss-quadruature
!  for a wide range of orthogonal polynomials ortogonalized over a given range
!  [a,b] and weight function w(x).
!
!  Inputs:
!         nt    = Number of Gauss-quadruature points, [integer(kind = 4)]
!         kind  = Type of Quadrature rule,  [integer(kind = 4)]
!
!               = 1,   Legendre,             [-1,1],     w(x) = 1.0
!               = 2,   Associated Laguerre,  [0,inf],    w(x) = exp(-x)*x^alpha
!               = 3,   Hermite,              [-inf,inf], w(x) = exp(-x^2)
!               = 4,   Jacobi,               [-1,1],     w(x) = (1-x)^alpha*(1+x)^beta
!               = 5,   Chebyshev (Type I),   [-1,1],     w(x) = (1-x^2)^(-0.5)
!               = 6,   Chebyshev (type II),  [-1,1],     w(x) = (1-x^2)^(0.5)
!               = 7,   Gegenbauer,           [-1,1],     w(x) = (1-x^2)^(alpha-0.5)
!
!         alpha = Value of alpha, if needed, [real(kind = 8)]
!         beta  = Value of beta, if needed, [real(kind = 8)]
!
!   Outputs:
!         x(1:nt) = abscissas for Gauss-quadrature rule, [real(kind=8)]
!         w(1:nt) = weights for Gauss-quadrature rule, [real(kind=8)]
!
!   Note that Chebyshev Type II is a special case of Gegenbauer with alpha = 1
!
!   The algorithm is based on the three-term recurrence relation for each
!   polynomial. The terms are found for the normalized polynonmials, i.e.,
!
!   F(n,alpha,beta,x) = p(n,alpha,beta,x)/sqrt(N(n,alpha,beta)),
!
!   where n is the order of the polynomial, and p(n,alpha,beta) denotes 
!   the standard textbook definition with p(0,alpha,beta) = 1.0, and 
!   N(n,alpha,beta) is the norm integral defined as
!
!   N(n,alpha,beta) = int[a,b] dx w(x)*[p(n,alpha,beta,x)]^2.
!
!   The three-term recurrence relations can be written as (with fixed 
!   alpha and beta)
!
!   a(n-1)*F(n-1,x) + b(n)*F(n,x) + a(n)*F(n+1,x) = x*F(n,x)
!
!   Using the notation
!   <g,p> = int[a,b] dx w(x)*g(x)*F(n,x).
!
!   b(n) = <x*F(n,alpha,beta),F(n,alpha,beta)>
!
!   The coefficients a(n) = k(n)/k(n+1)
!
!   where k(n) is the leading coefficient for F(n), that is
!
!   F(n) = k(n)*x^n + k(n-1)*x^(n-1) + ...
!
!   The abscicca x(k) is the kth eigenvalue of the tridigonal matrix
!
!   b(0)   a(0)      0      0      0                  0       0
!   a(0)   b(1)   a(1)      0      0                  0       0
!      0   a(1)   b(2)   a(2)      0                  0       0
!      0      0   a(2)   b(3)   a(3)                  0       0
!   
!                                                b(n-1)  a(n-1)
!      0      0      0      0      0             a(n-1)    b(n)
!
!   The weights are obtained from the 1st component of the kth eigenvector, 
!   vec(1,k), via:
!
!   w(k) = N(0,alpha,beta)*vec(1,k)^2/(sum_k vec(1,k)^2)
!
!   The eigenvalues and eigenvectors are obtained from a modified version of 
!   the EISPACK routine tql2, named tri_eig
!
!   Dependencies:
!     jacobi_norm
!     tri_eig
!
!  Date:
!
!    10 March 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   implicit none
   integer(kind=4), intent(in) :: nt, kind
   real(kind=8), intent(in) :: alpha, beta
   real(kind=8), intent(out) :: x(nt), w(nt)

   integer(kind=4) :: i, n
   real(kind=8) :: pi

   real(kind=8) :: xn, two_xn

   real(kind=8), allocatable :: d(:), e(:), z(:,:)
   real(kind=8) :: norm, int_norm
   integer(kind=4) :: ierror
!
!-----    External functions
!
   real(kind=8) :: poly_a, poly_b, poly_norm
!
!*******************************************************************************
!
   pi = 2.0d0*asin(1.0d0)
!
!----   Allocate and initalize arrays for the tri-diagonal matrix
!
   if(.not.allocated(d))allocate(d(nt))
   if(.not.allocated(e))allocate(e(nt))
   if(.not.allocated(z))allocate(z(nt,nt))
   d(1:nt) = 0.0d0
   e(1:nt) = 0.0d0
   z(1:nt,1:nt) = 0.0d0
!
!----   Set up tridiagonal matrix to be diagonalized for each integration type
!
   do i = 1, nt
      n = i - 1                           !   integer order of the polynomial, p(n,alpha,beta)
      xn = real(n,kind=8)
      two_xn = 2.0d0*xn
      d(i) = poly_b(n, kind, alpha, beta)
      if(i < nt)e(i+1)= poly_a(n, kind, alpha, beta)
      z(i,i) = 1.0d0
   end do
!
!   Find eigenvalues and eigenvectors of the tri-diagonal matrix
!
   call tri_eig(nt, nt, d, e, z, ierror)
!
!----    Find norm to set the weights
!
   norm = 0.0d0
   do i = 1, nt
      norm = norm + z(1,i)**2
   end do
!
!----    Set norm for the weights - integral with w(x)*[p(0,alpha,beta)]^2*dx
!
   int_norm = poly_norm(0,kind,alpha,beta)
!
!----   Put abscissas and weights into output arrays
! 
   do i = 1, nt
      x(i) = d(i)
      w(i) = int_norm*z(1,i)**2/norm
   end do
!
!----   Delete temporary data 
!
   if(allocated(d))deallocate(d)
   if(allocated(e))deallocate(e)
   if(allocated(z))deallocate(z)
!
!----   All done
!
   return
!   
end subroutine gauss_quad

!
!
subroutine tri_eig(nm,n,d,e,z,ierror)
!
   implicit none
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
         if (tst2 == tst1)exit
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
   if (p == 0.0d0) return
   r = (min(abs(a),abs(b))/p)**2
   t = 4.0d0 + r
   do while(t /= 4.0d0)
      s = r/t
      u = 1.0d0 + 2.0d0*s
      p = u*p
      r = (s/u)**2*r
      t = 4.0d0 + r
   end do
   pythag_90 = p 
   return
end function pythag_90
!
real(kind=8) function poly_b(n, kind, alpha, beta)
!*******************************************************************************
!
!  Discussion:
!
!  FORTRAN 90 function to calculate recurssion coefficient b(n), where 
!   the three-term recurrence relations are written as (with fixed 
!   alpha and beta) for the normalized orthogonal polynomials as:
!
!   a(n-1)*F(n-1,x) + b(n)*F(n,x) + a(n)*F(n+1,x) = x*F(n,x)
!
!  Inputs:
!         n     = order of the polynomial, [integer(kind = 4)]
!         kind  = Type of polynomial,  [integer(kind = 4)]
!
!               = 1,   Legendre,             [-1,1],     w(x) = 1.0
!               = 2,   Associated Laguerre,  [0,inf],    w(x) = exp(-x)*x^alpha
!               = 3,   Hermite,              [-inf,inf], w(x) = exp(-x^2)
!               = 4,   Jacobi,               [-1,1],     w(x) = (1-x)^alpha*(1+x)^beta
!               = 5,   Chebyshev (Type I),   [-1,1],     w(x) = (1-x^2)^(-0.5)
!               = 6,   Chebyshev (type II),  [-1,1],     w(x) = (1-x^2)^(0.5)
!               = 7,   Gegenbauer,           [-1,1],     w(x) = (1-x^2)^(alpha-0.5)
!
!         alpha = Value of alpha, if needed, [real(kind = 8)]
!         beta  = Value of beta, if needed, [real(kind = 8)]
!
!  Date:
!
!    10 March 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
    implicit none
    integer(kind=4), intent(in) :: n, kind
    real(kind=8), intent(in) :: alpha, beta
!
!--------------------------------------------------------------------
!
    real(kind=8) :: xn, two_xn, xLn
!
!---------   External functions
!
!*************************************************************************
!
    poly_b = 0.0d0
    xn = real(n,kind=8)
    two_xn = 2.0d0*xn
    if(kind == 1)then                       !   Legendre
       poly_b = 0.0d0
    elseif(kind == 2)then                   !   Associated Laguerre
       poly_b = two_xn + alpha + 1.0d0
    elseif(kind == 3)then                   !   Hermite
       poly_b = 0.0d0
    elseif(kind == 4)then                   !   Jacobi
       xLn = two_xn + alpha + beta
       poly_b = (beta**2 - alpha**2)/(xLn*(xLn+2.0d0))
    elseif(kind == 5)then                   !   Chebyshev (Type I)
       poly_b = 0.0d0
    elseif(kind == 6)then                   !   Chebyshev (Type II)
       poly_b = 0.0d0
    elseif(kind == 7)then                   !   Gegenbauer
       if(alpha < 1.0d0)stop 'Error in Gengenbauer polynomilas, alpha < 1'
       poly_b = 0.0d0
    end if
    return
end function poly_b

real(kind=8) function poly_a(n, kind, alpha, beta)
!*******************************************************************************
!
!  Discussion:
!
!  FORTRAN 90 function to calculate recurssion coefficient a(n), where 
!   the three-term recurrence relations are written as (with fixed 
!   alpha and beta) for the normalized orthogonal polynomials as:
!
!   a(n-1)*F(n-1,x) + b(n)*F(n,x) + a(n)*F(n+1,x) = x*F(n,x)
!
!  Inputs:
!         n     = order of the polynomial, [integer(kind = 4)]
!         kind  = Type of polynomial,  [integer(kind = 4)]
!
!               = 1,   Legendre,             [-1,1],     w(x) = 1.0
!               = 2,   Associated Laguerre,  [0,inf],    w(x) = exp(-x)*x^alpha
!               = 3,   Hermite,              [-inf,inf], w(x) = exp(-x^2)
!               = 4,   Jacobi,               [-1,1],     w(x) = (1-x)^alpha*(1+x)^beta
!               = 5,   Chebyshev (Type I),   [-1,1],     w(x) = (1-x^2)^(-0.5)
!               = 6,   Chebyshev (type II),  [-1,1],     w(x) = (1-x^2)^(0.5)
!               = 7,   Gegenbauer,           [-1,1],     w(x) = (1-x^2)^(alpha-0.5)
!
!         alpha = Value of alpha, if needed, [real(kind = 8)]
!         beta  = Value of beta, if needed, [real(kind = 8)]
!
!  Date:
!
!    10 March 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
    implicit none
    integer(kind=4), intent(in) :: n, kind
    real(kind=8), intent(in) :: alpha, beta
!
!--------------------------------------------------------------------
!
    real(kind=8) :: xn, two_xn, XLn, temp
!
!---------   External functions
    real(kind=8) :: poly_norm
!
!*************************************************************************
!
    poly_a = 0.0d0
    xn = real(n,kind=8)
    two_xn = 2.0d0*xn
    if(kind == 1)then                       !   Legendre
       poly_a = (xn + 1.0d0)/sqrt((two_xn + 3.0d0)*(two_xn + 1.0d0))
    elseif(kind == 2)then                   !   Associated Laguerre
       poly_a = -sqrt((xn+1.0d0)*(xn+alpha+1.0d0))
    elseif(kind == 3)then                   !   Hermite
       poly_a = sqrt((xn + 1.0d0)/2.0d0)
    elseif(kind == 4)then                   !   Jacobi
       xLn = two_xn + alpha + beta
       temp = log_gamma(alpha + beta + two_xn + 1.0d0) + log_gamma(alpha + beta + xn + 2.0d0) -   &
              log_gamma(alpha + beta + xn + 1.0d0) - log_gamma(alpha + beta + two_xn + 3.0d0)
       temp = exp(temp)*2.0d0*(xn + 1.0d0)*sqrt(poly_norm(n+1,kind,alpha,beta)/poly_norm(n,kind,alpha,beta))
       poly_a = temp
    elseif(kind == 5)then                   !   Chebyshev (Type I)
       if(n == 0)then
          poly_a = 1.0d0/sqrt(2.0d0)
       else
          poly_a = 0.5d0
       end if
    elseif(kind == 6)then                   !   Chebyshev (Type II)
       poly_a = 0.5d0
    elseif(kind == 7)then                   !   Gegenbauer
       if(alpha < 1.0d0)stop 'Error in Gengenbauer polynomilas, alpha < 1'
       poly_a = 0.5d0*sqrt((xn+1.0d0)*(xn+2.0d0*alpha)/((xn+alpha)*(xn+alpha+1.0d0)))
    end if
    return
end function poly_a
!
!
real(kind=8) function poly_norm(n, kind, alpha, beta)
!
!*******************************************************************************
!
!  Discussion:
!
!  FORTRAN 90 function to calculate the normalization of orthogonal Polynomials
!
!   int[-1,1] w(x) [P(n,alpha,beta)]^2 dx
!  Inputs:
!         n     = order of the polynomial, [integer(kind = 4)]
!         kind  = Type of polynomial,  [integer(kind = 4)]
!
!               = 1,   Legendre,             [-1,1],     w(x) = 1.0
!               = 2,   Associated Laguerre,  [0,inf],    w(x) = exp(-x)*x^alpha
!               = 3,   Hermite,              [-inf,inf], w(x) = exp(-x^2)
!               = 4,   Jacobi,               [-1,1],     w(x) = (1-x)^alpha*(1+x)^beta
!               = 5,   Chebyshev (Type I),   [-1,1],     w(x) = (1-x^2)^(-0.5)
!               = 6,   Chebyshev (type II),  [-1,1],     w(x) = (1-x^2)^(0.5)
!               = 7,   Gegenbauer,           [-1,1],     w(x) = (1-x^2)^(alpha-0.5)
!
!         alpha = Value of alpha, if needed, [real(kind = 8)]
!         beta  = Value of beta, if needed, [real(kind = 8)]
!
!  Date:
!
!    10 March 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
    implicit none
    integer(kind=4), intent(in) :: n, kind
    real(kind=8), intent(in) :: alpha, beta
!
!--------------------------------------------------------------------
!
    real(kind=8) :: factor, xn, two_xn, pi
!
!*************************************************************************
!
    pi = 2.0d0*asin(1.0d0)
    poly_norm = 0.0d0
    xn = real(n,kind=8)
    two_xn = 2.0d0*xn
    if(kind == 1)then                       !   Legendre
       poly_norm = 2.0d0/(two_xn + 1.0d0)
    elseif(kind == 2)then                   !   Associated Laguerre
       poly_norm = exp(log_gamma(xn + alpha + 1.0d0) - log_gamma(xn + 1.0d0))
    elseif(kind == 3)then                   !   Hermite
       poly_norm = 2.0d0**xn*sqrt(pi)*exp(log_gamma(xn + 1.0d0))
    elseif(kind == 4)then                   !   Jacobi
       factor = log_gamma(xn + alpha + 1.0d0) + log_gamma(xn + beta + 1.0d0) -  &
                log_gamma(xn + 1.0d0) - log_gamma(xn + alpha + beta + 1.0d0)
       factor = factor + (alpha + beta + 1.0d0)*log(2.0d0) -                    &
                log(2.0d0*xn + alpha + beta + 1.0d0)
       poly_norm = exp(factor)
    elseif(kind == 5)then                   !   Chebyshev (Type I)
       if(n == 0)then
          poly_norm = pi
       else
          poly_norm = 0.5d0*pi
       end if
    elseif(kind == 6)then                   !   Chebyshev (Type II)
       poly_norm = 0.5d0*pi
    elseif(kind == 7)then                   !   Gegenbauer
       factor = log_gamma(xn + 2.0d0*alpha) - log_gamma(xn + 1.0d0) - 2.0d0*log_gamma(alpha)
       poly_norm = pi*2.0d0**(1.0d0 - 2.0d0*alpha)*exp(factor)/(xn + alpha)
    end if
    return
end function poly_norm
!
!
!
real(kind=8) function poly(n, kind, alpha, beta, x)
!*******************************************************************************
!
!  Discussion:
!
!  FORTRAN 90 function to calculate values at a given value of x for the 
!   orthogonal polynomial p(n,alpha,beta,x). 
!   Uses the three-term recurrence relation (with fixed 
!   alpha and beta) for the normalized orthogonal polynomials F(n,alpha,beta,x) as:
!
!   a(n-1)*F(n-1,x) + b(n)*F(n,x) + a(n)*F(n+1,x) = x*F(n,x)
!
!   or
!
!   F(n) = ((x - b(n-1)*F(n-1)) - a(n-2)*F(n-2))/a(n-1) 
!
!   with
!
!   F(n,alpha,beta,x) = p(n,alpha,beta,x)/sqrt(N(n,alpha,beta)),
!
!  Inputs:
!         n     = order of the polynomial, [integer(kind = 4)]
!         kind  = Type of polynomial,  [integer(kind = 4)]
!
!               = 1,   Legendre,             [-1,1],     w(x) = 1.0
!               = 2,   Associated Laguerre,  [0,inf],    w(x) = exp(-x)*x^alpha
!               = 3,   Hermite,              [-inf,inf], w(x) = exp(-x^2)
!               = 4,   Jacobi,               [-1,1],     w(x) = (1-x)^alpha*(1+x)^beta
!               = 5,   Chebyshev (Type I),   [-1,1],     w(x) = (1-x^2)^(-0.5)
!               = 6,   Chebyshev (type II),  [-1,1],     w(x) = (1-x^2)^(0.5)
!               = 7,   Gegenbauer,           [-1,1],     w(x) = (1-x^2)^(alpha-0.5)
!
!         alpha = Value of alpha, if needed, [real(kind = 8)]
!         beta  = Value of beta, if needed, [real(kind = 8)]
!
!  Date:
!
!    10 March 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   implicit none
   integer(kind=4), intent(in) :: n, kind
   real(kind=8), intent(in) :: alpha, beta
   real(kind=8), intent(in) :: x
!----------------------------------------------------------------
   integer(kind=4) :: nn
   real(kind=8) :: p0, p1, p2
!
!-----    External functions
!
   real(kind=8) :: poly_a, poly_b, poly_norm
!
!*******************************************************************************
!
   p0 = 1.0d0
   p1 = 0.0d0
!---   If n=1, trivial, and exit
   poly = p0
   if(n == 0)return
   if(kind == 1)then
      p1 = x
   elseif(kind == 2)then
      p1 = (-x + alpha + 1.0d0)
   elseif(kind == 3)then
      p1 = 2.0d0*x
   elseif(kind == 4)then
      p1 = (alpha + 1.0d0) + (alpha + beta + 2.0d0)*0.5d0*(x-1.0d0)
   elseif(kind == 5)then
      p1 = x
   elseif(kind == 6)then
      p1 = 2.0d0*x
   elseif(kind == 7)then
      p1 = 2.0d0*alpha*x
   end if
!---   If n=1, trivial, and exit
   poly = p1
   if(n == 1)return
!---   First and second normalized polynomials
   p0 = p0/sqrt(poly_norm(0,kind,alpha,beta))
   p1 = p1/sqrt(poly_norm(1,kind,alpha,beta))
   p2 = 0.0d0
!---   Build up to n via recurrence relation
   do nn = 2, n
      p2 = (x - poly_b(nn-1,kind,alpha,beta))*p1 -            &
           poly_a(nn-2,kind,alpha,beta)*p0
      p2 = p2/poly_a(nn-1,kind,alpha,beta)
      p0 = p1
      p1 = p2
   end do
!---   Undo normalization and return
   poly = p2*sqrt(poly_norm(n,kind,alpha,beta))
   return
end function poly

