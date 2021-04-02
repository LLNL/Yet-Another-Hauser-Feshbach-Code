!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine examines the incident energy grid and fixes issues,
!    such as photon energies being too small, mapping onto the nergy grid,
!    and removing redundant incident energies 
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
!*****************************************************************************80
!
subroutine fix_incident_energies(iproj, rel_factor)
   use variable_kinds
   use options
   use nodeinfo
   use useful_data
   use nuclei
   use particles_def
   implicit none
   integer(kind=4), intent(in) :: iproj
   real(kind=8), intent(in) :: rel_factor
!-------------------------------------------------------------------------
   integer(kind=4) :: num_energies
   integer(kind=4) :: j, in, in2, in3
   integer(kind=4) :: ishift
   logical :: e_in_problem
   real(kind=8) :: e_in, e_in2, e_rel, e_x
!------------    External functions
   integer(kind=4) :: find_ibin
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------                                                                   +
!--------   Now that energy grids are set up, remap incident projectile   +
!--------   energies to fit on excitation energy grid                     +
!------                                                                   +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   num_energies = projectile%num_e


!---    If incident photons, remove low-energy below E_cut
   e_in_problem = .false.
   if(iproj == 0)then
      in = 1
      do while(in <= num_energies)
         e_in = projectile%energy(in)
         if(e_in >= nucleus(1)%e_grid(1))then   !  Energy is good
            in = in + 1
            cycle
         else                                   !  It is not, remove it, push rest down
            do in2 = in + 1, num_energies
               projectile%energy(in2-1) = projectile%energy(in2)
            end do
            num_energies = num_energies - 1
            e_in_problem = .true.
         end if
      end do
      if(e_in_problem .and. iproc == 0)then
         write(6,*)'**************************************************************'
         write(6,*)'* Warning, some incident photon energies were too low, and   *'
         write(6,*)'* will be removed from this calculation.                     *'
         write(6,*)'**************************************************************'
      end if
   end if

   e_in_problem = .false.
   do in = 1, num_energies
      e_in = projectile%energy(in)
      if(.not.pop_calc)then
!   write(6,*)in, e_in
         e_rel = e_in*rel_factor
         e_x = e_rel + nucleus(1)%sep_e(projectile%particle_type)
         j = find_ibin(E_x, 1)
         if(e_x < nucleus(1)%sep_e(projectile%particle_type) +     &
               0.5d0*nucleus(1)%delta_e(1))cycle
!
!----   Find bin associated with this energy
!
         j = find_ibin(E_x, 1)
!   write(6,*)in, e_in, j, e_x, nucleus(1)%sep_e(projectile%particle_type)
!   write(6,*)nucleus(1)%e_grid(j+1)
!   write(6,*)nucleus(1)%e_grid(j)
!   write(6,*)nucleus(1)%e_grid(j-1)
         if(j > 0)then
            e_x = nucleus(1)%e_grid(j)
            e_rel = e_x - nucleus(1)%sep_e(projectile%particle_type)
            e_in = e_rel/rel_factor
            projectile%energy(in) = e_in
         else
            e_x = nucleus(1)%e_grid(1)
            e_rel = e_x - nucleus(1)%sep_e(projectile%particle_type)
            e_in = e_rel/rel_factor
            projectile%energy(in) = e_in
            e_in_problem = .true.
         end if
      else
         e_x = e_in
         j = find_ibin(e_x,1)
         e_x = nucleus(1)%e_grid(j)
         projectile%energy(in) = e_in
      end if
   end do
   if(iproj > 0 .and. e_in_problem .and. iproc == 0)then
      write(6,*)'**************************************************************'
      write(6,*)'* Warning, compound excitation extends to low energy with    *'
      write(6,*)'* with low level density. Calculations will only be          *'
      write(6,*)'* for incident energies that populate the continuous energy  *'
      write(6,*)'* bins.                                                      *'
      write(6,*)'* Continue at your own risk!                                 *'
      write(6,*)'**************************************************************'
   end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Now look to see if we have to remove any energies from the list  +
!------   due to potential overlap - no goto's                             +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   in = 1
   do while(in < num_energies)
      e_in = projectile%energy(in)
      if(iproj > 0 .and. e_in < 0.5d0*de)then         !   keep it if energy is less than de/2
         in = in + 1
         cycle
      end if
      in2 = in + 1
      do while(in2 <= num_energies)
         e_in2 = projectile%energy(in2)
         if(abs(e_in - e_in2) < 1.0d-4)then          !  next energy is the same
            do in3 = in2 + 1, num_energies           !  push energies in list down
               ishift = in3 - in2 - 1
               projectile%energy(in2+ishift) = projectile%energy(in3)
            end do
            num_energies = num_energies - 1     !  reduce number of energies in list
         else
            in2 = in2 + 1
         end if
      end do
      in = in + 1
   end do

!  do in = 1, num_energies
!      write(6,*)in,projectile%energy(in)
!   end do

! stop

   projectile%num_e = num_energies

!   if(allocated(test_e))deallocate(test_e)
!     do in = 1, num_energies
!        write(6,*)in,projectile%energy(in)
!     end do

!     stop
   return

end subroutine fix_incident_energies

!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine examines the incident energy grid and fixes issues,
!    such as photon energies being too small, mapping onto the nergy grid,
!    and removing redundant incident energies 
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
!*****************************************************************************80
!
subroutine fix_pop_energies
   use variable_kinds
   use options
   use nodeinfo
   use useful_data
   use nuclei
   use particles_def
   implicit none
!-------------------------------------------------------------------------
   integer(kind=4) :: num_energies
   integer(kind=4) :: i, j, k, n
   integer(kind=4) :: in, in2, in3
   integer(kind=4) :: ishift
   logical :: e_in_problem
   real(kind=8) :: e_in, e_in2
!------------    External functions
   integer(kind=4) :: find_ibin
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------                                                                   +
!--------   Now that energy grids are set up, remap incident projectile   +
!--------   energies to fit on excitation energy grid                     +
!------                                                                   +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   num_energies = projectile%num_e

!---    If incident photons, remove low-energy below E_cut

   e_in_problem = .false.

   do in = 1, num_energies
      e_in = projectile%energy(in)
      j = find_ibin(e_in,1)
      if(j > 0)then
         e_in = nucleus(1)%e_grid(j)
         projectile%energy(in) = e_in
         Pop_data(in)%Ex_pop = e_in
      else
         do in2 = in + 1, num_energies
            projectile%energy(in2-1) = projectile%energy(in2)
!---   We also have to do all the pop data
            k = in2-1
            Pop_data(k)%Ex_pop = Pop_data(in2)%Ex_pop
            Pop_data(k)%dEx_pop = Pop_data(in2)%dEx_pop
            Pop_data(k)%num_pop = Pop_data(in2)%num_pop
            n = Pop_data(in2)%num_pop
            if(allocated(Pop_data(k)%j_pop))deallocate(Pop_data(k)%j_pop)
            if(allocated(Pop_data(k)%par_pop))deallocate(Pop_data(k)%par_pop)
            if(allocated(Pop_data(k)%bin_pop))deallocate(Pop_data(k)%bin_pop)
            if(.not.allocated(Pop_data(k)%j_pop))allocate(Pop_data(k)%j_pop(n))
            if(.not.allocated(Pop_data(k)%par_pop))allocate(Pop_data(k)%par_pop(n))
            if(.not.allocated(Pop_data(k)%bin_pop))allocate(Pop_data(k)%bin_pop(n))
            do i = 1, Pop_data(in2)%num_pop
               Pop_data(k)%j_pop(i) = Pop_data(in2)%j_pop(i)
               Pop_data(k)%par_pop(i) = Pop_data(in2)%par_pop(i)
               Pop_data(k)%bin_pop(i) = Pop_data(in2)%bin_pop(i)
           end do
         end do
         num_energies = num_energies - 1
         e_in_problem = .true.
      end if
   end do
   if(e_in_problem .and. iproc == 0)then
      write(6,'(''**************************************************************'')')
      write(6,'(''* Warning, some population excitation energies extend to     *'')')
      write(6,'(''* very low energy with with a low level density.             *'')')
      write(6,'(''* Excitation energies must be above E_cut = '',f10.5''       *'')')nucleus(1)%level_ecut
      write(6,'(''* Energies requested below this value will be removed        *'')')
      write(6,'(''* --------      Continue at your own risk!      -------------*'')')
      write(6,'(''**************************************************************'')')
   end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Now look to see if we have to remove any energies from the list  +
!------   due to potential overlap - no goto's                             +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   in = 1
   do while(in < num_energies)
      e_in = projectile%energy(in)
      in2 = in + 1
      do while(in2 <= num_energies)
         e_in2 = projectile%energy(in2)
         if(abs(e_in - e_in2) < 1.0d-4)then          !  next energy is the same
            do in3 = in2 + 1, num_energies           !  push energies in list down
               ishift = in3 - in2 - 1
               k = in2 + ishift
               projectile%energy(k) = projectile%energy(in3)
!---   We also have to do all the pop data
               Pop_data(k)%Ex_pop = Pop_data(in3)%Ex_pop
               Pop_data(k)%dEx_pop = Pop_data(in3)%dEx_pop
               Pop_data(k)%num_pop = Pop_data(in3)%num_pop
               n = Pop_data(in3)%num_pop
               if(allocated(Pop_data(k)%j_pop))deallocate(Pop_data(k)%j_pop)
               if(allocated(Pop_data(k)%par_pop))deallocate(Pop_data(k)%par_pop)
               if(allocated(Pop_data(k)%bin_pop))deallocate(Pop_data(k)%bin_pop)
               if(.not.allocated(Pop_data(k)%j_pop))allocate(Pop_data(k)%j_pop(n))
               if(.not.allocated(Pop_data(k)%par_pop))allocate(Pop_data(k)%par_pop(n))
               if(.not.allocated(Pop_data(k)%bin_pop))allocate(Pop_data(k)%bin_pop(n))
               do i = 1, Pop_data(in3)%num_pop
                  Pop_data(k)%j_pop(i) = Pop_data(in3)%j_pop(i)
                  Pop_data(k)%par_pop(i) = Pop_data(in3)%par_pop(i)
                  Pop_data(k)%bin_pop(i) = Pop_data(in3)%bin_pop(i)
              end do
            end do
            num_energies = num_energies - 1     !  reduce number of energies in list
         else
            in2 = in2 + 1
         end if
      end do
      in = in + 1
   end do

   projectile%num_e = num_energies
   num_pop_e = num_energies

   if(num_energies == 0)then
     if(iproc == 0)write(6,*)'num_energies = 0! Calculation will be terminated'
     call MPI_Abort(icomm,101,ierr)
   end if

   return

end subroutine fix_pop_energies


