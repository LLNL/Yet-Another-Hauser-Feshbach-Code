!
!*******************************************************************************
!
subroutine pack_data(Ix_f, ip_f, n_f, idb, l, iss, itemp)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine packs integer data defining aspects of a decay
!    into a single integer(kind=4) itemp
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
   integer(kind=4), intent(in) :: Ix_f, ip_f, n_f, idb, l, iss
   integer(kind=4), intent(out) :: itemp

   itemp = Ix_f
   itemp = ior(itemp,ishft(ip_f,6))
   itemp = ior(itemp,ishft(n_f,7))
   itemp = ior(itemp,ishft(idb,21))
   itemp = ior(itemp,ishft(l,22))
   itemp = ior(itemp,ishft(iss,28))
   return
end subroutine pack_data 
