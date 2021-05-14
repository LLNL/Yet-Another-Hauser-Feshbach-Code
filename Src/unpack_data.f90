 
subroutine unpack_data(Ix_f, ip_f, n_f, idb, l, iss, itemp)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine unpacks data in the integer(kind=4) itemp
!    that was packed into it using the subroutine pack_data
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
   integer(kind=4), intent(out) :: Ix_f, ip_f, n_f, idb, l, iss
   integer(kind=4), intent(in) :: itemp


   Ix_f = iand(itemp,2**6-1)
   ip_f = iand(ishft(itemp,-6),1)
   n_f = iand(ishft(itemp,-7),2**14-1)
   idb = iand(ishft(itemp,-21),1)
   l = iand(ishft(itemp,-22),2**6-1)
   iss = iand(ishft(itemp,-28),2**5-1)

   return
end subroutine unpack_data 
 
