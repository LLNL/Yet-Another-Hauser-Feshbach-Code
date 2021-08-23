!
!*****************************************************************************80
!
subroutine remove_states(inuc, num, state_map)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine examines discrete states and removes them based on the marker
!    marker state_map(i)=-1.  
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
!*****************************************************************************80
!
   use nodeinfo
   use nuclei
   implicit none
   integer(kind=4), intent(in) :: inuc
   integer(kind=4), intent(in) :: num
   integer(kind=4), intent(inout) :: state_map(num)
!---------------------------------------------------------------------------
   integer(kind=4) :: m, n
   integer(kind=4) :: ibranch
   integer(kind=4) :: nnd
!---------------------------------------------------------------------------
!---   Before removing states, check on ecut. If levels below ecut are being 
!---   removed, we may have to adjust ecut
!---   Loop over levels up to ecut (ncut), pick maximum 
   nnd = 0   
   do n = 1, nucleus(inuc)%num_discrete
      if(state_map(n) > 0)then
         nnd = nnd + 1
         if(state_map(n) /= n)then     !   put state n into state nnd
            nucleus(inuc)%state(nnd)%energy = nucleus(inuc)%state(n)%energy
            nucleus(inuc)%state(nnd)%spin = nucleus(inuc)%state(n)%spin
            nucleus(inuc)%state(nnd)%parity = nucleus(inuc)%state(n)%parity
            nucleus(inuc)%state(nnd)%t12 = nucleus(inuc)%state(n)%t12
            nucleus(inuc)%state(nnd)%isomer = nucleus(inuc)%state(n)%isomer
            nucleus(inuc)%state(nnd)%shift = nucleus(inuc)%state(n)%shift
            nucleus(inuc)%state(nnd)%level_float = nucleus(inuc)%state(n)%level_float
            nucleus(inuc)%state(nnd)%iband = nucleus(inuc)%state(n)%iband
            nucleus(inuc)%state(nnd)%nbranch = nucleus(inuc)%state(n)%nbranch
            nucleus(inuc)%state(nnd)%exit_lab = nucleus(inuc)%state(n)%exit_lab
            nucleus(inuc)%state(nnd)%level_float = nucleus(inuc)%state(n)%level_float
            nucleus(inuc)%state(nnd)%state_modified = nucleus(inuc)%state(n)%state_modified
            if(allocated(nucleus(inuc)%state(nnd)%ibranch))                            &
               deallocate(nucleus(inuc)%state(nnd)%ibranch)
            allocate(nucleus(inuc)%state(nnd)%ibranch(nucleus(inuc)%state(nnd)%nbranch))
            if(allocated(nucleus(inuc)%state(nnd)%branch))                             &
               deallocate(nucleus(inuc)%state(nnd)%branch)
            allocate(nucleus(inuc)%state(nnd)%branch(nucleus(inuc)%state(nnd)%nbranch))
            if(allocated(nucleus(inuc)%state(nnd)%egamma))                             &
               deallocate(nucleus(inuc)%state(nnd)%egamma)
            allocate(nucleus(inuc)%state(nnd)%egamma(nucleus(inuc)%state(nnd)%nbranch))
            if(allocated(nucleus(inuc)%state(nnd)%p_gamma))                            &
               deallocate(nucleus(inuc)%state(nnd)%p_gamma)
            allocate(nucleus(inuc)%state(nnd)%p_gamma(nucleus(inuc)%state(nnd)%nbranch))
            if(allocated(nucleus(inuc)%state(nnd)%p_ic))                               &
               deallocate(nucleus(inuc)%state(nnd)%p_ic)
            allocate(nucleus(inuc)%state(nnd)%p_ic(nucleus(inuc)%state(nnd)%nbranch))
            if(allocated(nucleus(inuc)%state(nnd)%cs))                                 &
               deallocate(nucleus(inuc)%state(nnd)%cs)
            allocate(nucleus(inuc)%state(nnd)%cs(nucleus(inuc)%state(nnd)%nbranch))
            if(allocated(nucleus(inuc)%state(nnd)%branch_modified))                                 &
               deallocate(nucleus(inuc)%state(nnd)%branch_modified)
            allocate(nucleus(inuc)%state(nnd)%branch_modified(nucleus(inuc)%state(nnd)%nbranch))
            do m = 1, nucleus(inuc)%state(n)%nbranch
               ibranch = nucleus(inuc)%state(n)%ibranch(m)
               nucleus(inuc)%state(nnd)%ibranch(m) = state_map(ibranch)
               nucleus(inuc)%state(nnd)%egamma(m) = nucleus(inuc)%state(n)%egamma(m)
               nucleus(inuc)%state(nnd)%branch(m) = nucleus(inuc)%state(n)%branch(m)
               nucleus(inuc)%state(nnd)%p_gamma(m) = nucleus(inuc)%state(n)%p_gamma(m)
               nucleus(inuc)%state(nnd)%p_ic(m) = nucleus(inuc)%state(n)%p_ic(m)
               nucleus(inuc)%state(nnd)%branch_modified(m) = nucleus(inuc)%state(n)%branch_modified(m)
               nucleus(inuc)%state(nnd)%cs(m) = 0.0d0
            end do
            nucleus(inuc)%state(nnd)%num_decay = nucleus(inuc)%state(n)%num_decay
            if(allocated(nucleus(inuc)%state(nnd)%decay_prob))                         &
               deallocate(nucleus(inuc)%state(nnd)%decay_prob)
            allocate(nucleus(inuc)%state(nnd)%decay_prob(nucleus(inuc)%state(nnd)%num_decay))
            if(allocated(nucleus(inuc)%state(nnd)%decay_type))                         &
               deallocate(nucleus(inuc)%state(nnd)%decay_type)
            allocate(nucleus(inuc)%state(nnd)%decay_type(nucleus(inuc)%state(nnd)%num_decay))
            if(allocated(nucleus(inuc)%state(nnd)%decay_to))                           &
               deallocate(nucleus(inuc)%state(nnd)%decay_to)
            allocate(nucleus(inuc)%state(nnd)%decay_to(nucleus(inuc)%state(nnd)%num_decay))
            do m = 1, nucleus(inuc)%state(n)%num_decay
               nucleus(inuc)%state(nnd)%decay_prob(m) = nucleus(inuc)%state(n)%decay_prob(m)
               nucleus(inuc)%state(nnd)%decay_type(m) = nucleus(inuc)%state(n)%decay_type(m)
               nucleus(inuc)%state(nnd)%decay_to(m) = nucleus(inuc)%state(n)%decay_to(m)
            end do
         end if
      end if
   end do
   do n = nnd + 1, nucleus(inuc)%num_discrete
      if(allocated(nucleus(inuc)%state(n)%ibranch))deallocate(nucleus(inuc)%state(n)%ibranch)
      if(allocated(nucleus(inuc)%state(n)%branch))deallocate(nucleus(inuc)%state(n)%branch)
      if(allocated(nucleus(inuc)%state(n)%egamma))deallocate(nucleus(inuc)%state(n)%egamma)
      if(allocated(nucleus(inuc)%state(n)%p_gamma))deallocate(nucleus(inuc)%state(n)%p_gamma)
      if(allocated(nucleus(inuc)%state(n)%p_ic))deallocate(nucleus(inuc)%state(n)%p_ic)
      if(allocated(nucleus(inuc)%state(n)%cs))deallocate(nucleus(inuc)%state(n)%cs)
      if(allocated(nucleus(inuc)%state(n)%branch_modified))deallocate(nucleus(inuc)%state(n)%branch_modified)
   end do
!---------  Once states have been removed, reset the number of discretes
   nucleus(inuc)%num_discrete = nnd

   return
end subroutine remove_states

