!
!*******************************************************************************
!
real(kind=8) function interp(x, num, xx, yy)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function interpolates the array yy, which is evalauted on 
!    the equally spaced array xx at the value x.
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
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
   use variable_kinds
   implicit none
   real(kind=8), intent(in) :: x
   integer(kind=4), intent(in) :: num
   real(kind=8), intent(in) :: xx(num), yy(num)
!--------------------------------------------------------
   integer(kind=4) :: idex
   real(kind=8) :: dxx
   real(kind=8) :: p
!--------------------------------------------------------
   dxx = xx(2) - xx(1)
   idex = nint((x - xx(1))/dxx) + 1
   if(idex == 1) then
      p = (x - xx(1))/dxx
      interp = (1.0d0 - p)*yy(1) + p*yy(2)
      return
   elseif (idex == num)then
      p = (x - xx(num-1))/dxx
      interp = (1.0d0 - p)*yy(num-1) + p*yy(num)
      return
   else
      p = (x - xx(idex))/dxx
      interp = p*(p-1.0d0)/2.0d0*yy(idex-1) + (1.0d0-p**2)*yy(idex) +  &
               p*(p+1.0d0)/2.0d0*yy(idex+1)
      return
   end if 
   return
end function interp

