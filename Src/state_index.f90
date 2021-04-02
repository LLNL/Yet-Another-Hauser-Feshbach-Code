integer function state_index(inuc, spin, parity, energy)
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
!   write(6,*)i,spin,spin,nucleus(inuc)%state(i)%spin
!   write(6,*)parity,nucleus(inuc)%state(i)%parity
!   write(6,*)e_diff
      if(real8_equal(spin,nucleus(inuc)%state(i)%spin) .and.               &
         real8_equal(parity,nucleus(inuc)%state(i)%parity) .and.           &
         e_diff <= 1.0d-3)then
         state_index = i
         exit
      end if
   end do
   return
end function state_index
