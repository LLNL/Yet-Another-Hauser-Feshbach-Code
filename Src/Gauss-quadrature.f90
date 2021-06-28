!
!*******************************************************************************
!
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
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        jacobi_norm
!        tri_eig
!
!     External functions:
!
!       real(kind=8) :: poly_a
!       real(kind=8) :: poly_b
!       real(kind=8) :: poly_norm
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
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
   real(kind=8) :: poly_a
   real(kind=8) :: poly_b 
   real(kind=8) :: poly_norm
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
!*******************************************************************************
!
real(kind=8) function poly_b(n, kind, alpha, beta)
!
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
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
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
!--------------------------------------------------------------------
    real(kind=8) :: xn, two_xn, xLn
!*************************************************************************
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
       poly_b = 0.0d0
       if(n == 0 .and. abs(beta**2 - alpha**2) < 1.0d-10)return
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
!
!*******************************************************************************
!
real(kind=8) function poly_a(n, kind, alpha, beta)
!
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
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        jacobi_norm
!        tri_eig
!
!     External functions:
!
!        real(kind=8) :: poly_norm
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
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
!--------------------------------------------------------------------
    real(kind=8) :: xn, two_xn, XLn, temp
!---------   External functions
    real(kind=8) :: poly_norm
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
!*******************************************************************************
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
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
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
!--------------------------------------------------------------------
    real(kind=8) :: factor, xn, two_xn, pi
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
!*******************************************************************************
!
real(kind=8) function poly(n, kind, alpha, beta, x)
!
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
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: poly_a
!        real(kind=8) :: poly_b
!        real(kind=8) :: poly_norm
!
!     MPI routines:
!
!        None
!
!  Licensing:
!
!    SPDX-License-Identifier: MIT 
!
!  Date:
!
!    11 May 2021
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
   real(kind=8) :: poly_a
   real(kind=8) :: poly_b
   real(kind=8) :: poly_norm
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

