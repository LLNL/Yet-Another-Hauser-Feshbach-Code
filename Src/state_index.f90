!
!*******************************************************************************
!
integer function state_index(inuc, spin, parity, energy)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function searches through the discrete states to find the index
!    of the state that best matches the requested spin, parity, and energy
!
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        logical :: real8_equal
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
   use nodeinfo
   use nuclei
   implicit none
   integer(kind=4), intent(in) :: inuc
   real(kind=8), intent(in) :: spin
   real(kind=8), intent(in) :: parity
   real(kind=8), intent(in) :: energy
!------------------------------------------------------------
   integer(kind=4) :: i
   real(kind=8) :: e_diff
!------------   External functions   ------------------------
   logical :: real8_equal
!------------------------------------------------------------
!----   Start with -1. Returning -1 means the state was not found
   state_index = -1
   do i = 1, nucleus(inuc)%num_discrete       !   these states need to be below E_cut
      e_diff = abs(energy - nucleus(inuc)%state(i)%energy)
      if(real8_equal(spin,nucleus(inuc)%state(i)%spin) .and.               &
         real8_equal(parity,nucleus(inuc)%state(i)%parity) .and.           &
         e_diff <= 1.0d-3)then
         state_index = i
         exit
      end if
   end do
   return
end function state_index
