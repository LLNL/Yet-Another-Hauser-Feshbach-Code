!
!*******************************************************************************
!
!  Discussion:
!
!    This function compares two real(kind=8) variables x and y for equality
!    within a tolerance of 1.0d-7.
!    It returns .false. if abs(x-y) > 1.0d-7, and .true. otherwise
!    replacement for x == y, unsafe compare of real variables
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
logical function real8_equal(x,y)
   implicit none
   real(kind=8), intent(in) :: x, y
!-------------------------------------------------------------------------------
   real8_equal = .true.
   if(abs(x-y) > 1.0d-7)real8_equal = .false.
   return
end function real8_equal
