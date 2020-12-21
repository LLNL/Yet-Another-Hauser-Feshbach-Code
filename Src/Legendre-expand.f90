subroutine Legendre_expand(jx_max,xval,Ang_Dist,L_max,Ang_L)
!
!*******************************************************************************
!
!  Discussion:
!
!  FORTRAN 90 subroutine to find the coefficients for an expaansion into
!  Legendre polynomials
!
!  Ang_Dist(x) = Sum Ang_L(L)*poly(L,1,alpha, beta,x)
!
!  Ang_L(L) = int(-1,1) Ang_Dist(x)*poly(L,1,alpha,beta,x)*dx*(2L+1)/2.
!
!  Note, input array is interpolated to be put onto the grid for the 
!  Gauss-Legendre integration points
!
!  Date:
!
!    20 March 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   implicit none 
   integer(kind=4), intent(in) :: jx_max
   real(kind=8), intent(in) :: xval(jx_max),Ang_Dist(jx_max)
   integer(kind=4), intent(in) :: L_max
   real(kind=8), intent(out) :: Ang_L(0:L_max)
!------------------------------------------------------
   integer(kind=4) :: L, LL_max
   real(kind=8) :: sum, value
   real(kind=8) :: alpha, beta
   real(kind=8) :: xnorm
!   integer(kind=4) :: jx
   integer(kind=4) :: n
   integer(kind=4) :: n_gleg
   real(kind=8), allocatable :: x_gleg(:), w_gleg(:)
!-------------   External function calls
   real(kind=8) :: Poly
   real(kind=8) :: interp
!------------------------------------------------------
!
!----   Set up data for Legendre integration
!----   Can't have more polynomials than number of data points in Ang_Dist
!----   The array Ang_Dist is assumed to be on a fixed grid
!----   with values xval(jx) and fixed spacing dx with dx = xval(2) - xval(1)
!----   First number of computed coefficients is LL_max
!----   set LL_max = L_max, but if L_max > jx_max, adjust
!----   Then the same for number of Gauss-Legendre points. Must be <= jx_max
!----   In general, preferably less
! 
   LL_max = L_max
   if(L_max > jx_max)LL_max = jx_max - 1
   n_gleg = jx_max 

   Ang_L(0) = 0.5d0
   if(LL_max == 0)return    !   All done, just one value and must be 0.5
   
   if(.not.allocated(x_gleg))allocate(x_gleg(n_gleg))
   if(.not.allocated(w_gleg))allocate(w_gleg(n_gleg))

   alpha = 0.0d0
   beta = 0.0d0
   call gauss_quad(n_gleg, 1, alpha, beta, x_gleg, w_gleg)

   xnorm = 0.0d0
   do n = 1, n_gleg
      value = interp(x_gleg(n),jx_max,xval,Ang_Dist)
      xnorm = xnorm + value*w_gleg(n)
   end do

   Ang_L(0:LL_max) = 0.0d0

   do L = 0, L_max
      sum = 0.0d0
      do n = 1, n_gleg
         value = interp(x_gleg(n),jx_max,xval,Ang_Dist)
         sum = sum + value*w_gleg(n)*poly(L,1,alpha,beta,x_gleg(n))
      end do
!      Ang_L(L) = sum*0.5d0*(2.0d0*real(L,kind=8)+1.0d0)
      Ang_L(L) = sum*0.5d0*(2.0d0*real(L,kind=8)+1.0d0)/xnorm
   end do
   if(allocated(x_gleg))deallocate(x_gleg)
   if(allocated(w_gleg))deallocate(w_gleg)
   return
end subroutine Legendre_expand
