!
!*******************************************************************************
!
!  Discussion:
!
!    This function compares two real(kind=8) variables x and y for equality
!    within a tolerance of 1.0d-7.
!    It returns .false. if abs(x-y) > 1.0d-7, and .true. otherwise
!
!  Reference:
!
!  S. Hilaire, Ch. lagrange, and A. J. Koning, Ann. Phys. 306, 209 (2003)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL version 2 license. 
!
!  Date:
!
!    25 September 2019
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
logical function real8_equal(x,y)
   implicit none
   real(kind=8), intent(in) :: x, y
   real8_equal = .true.
   if(abs(x-y) > 1.0d-7)real8_equal = .false.
   return
end function real8_equal
