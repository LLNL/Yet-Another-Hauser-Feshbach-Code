!
!*****************************************************************************80
!
subroutine Setup_Energy_bins(de)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine sets up the energy bins for the continuous part of the 
!    energy spectrum
!
!   Dependencies:
!
!     Modules:
!
!        options
!        nuclei
!        nodeinfo
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: delta_e_value
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
   use options
   use nuclei
   use nodeinfo
   implicit none
   real(kind=8), intent(in) :: de
!---------------------------------------------------------------------------
   integer(kind=4) :: i, j, iproj
   integer(kind=4) :: nbin
   real(kind=8) :: e_grid, delta_e, delta_e_old
   real(kind=8) :: ex_min, ex_base
   real(kind=8) :: Ex_max, Ex_calc
   real(kind=8) :: S_part
   real(kind=8) :: dep
   real(kind=8) :: E_cut
   logical :: finished
!-------------     External functions              -------------------------
   real(kind=8) :: delta_e_value
!---------------------------------------------------------------------------
   iproj = projectile%particle_type
   S_part = 0.0d0
   if(iproj >= 0)S_part = nucleus(1)%sep_e(iproj)
   if(.not. use_unequal_bins)then
      do i = 1, num_comp                  !    loop over nuclei
         Ex_max = nucleus(i)%Ex_max       !   max energy
         S_part = nucleus(i)%sep_e(1)
         E_cut = nucleus(i)%level_ecut
         if(E_cut + de < S_part)then
            Ex_calc = min(Ex_max,S_part) - E_cut - de   !   start at E_cut
            nbin = max(int(Ex_calc/de),1)               !   number of bins
            dep = Ex_calc/real(nbin,kind=8)             !   new energy step
            nbin = nint((Ex_max - E_cut)/dep) + 1       !   actual number of bins
            nucleus(i)%nbin = nbin
         else
            Ex_calc = (Ex_max - E_cut)
            if(Ex_calc > 0.0d0)then
               nbin = max(int(Ex_calc/de),1)               !   number of bins
               dep = Ex_calc/real(nbin,kind=8)             !   new energy step
               nbin = nint((Ex_max - E_cut)/dep) + 1       !   actual number of bins
               nucleus(i)%nbin = nbin
            else
               dep = de
               nbin = 1
               nucleus(i)%nbin = nbin
            end if
         end if

         if(nbin > 0)then
            allocate(nucleus(i)%e_grid(nbin))
            allocate(nucleus(i)%delta_e(nbin))
         end if
!-------------------   Set up energy grids for the nucleus
         nbin = nucleus(i)%nbin
         if(nbin == 0)cycle
         do j = 1, nbin
            nucleus(i)%e_grid(j) = E_cut + dep*real(j-1) + 0.5d0*dep
            nucleus(i)%delta_e(j) = dep
         end do
        nucleus(i)%Ex_max = nucleus(i)%e_grid(nbin)
      end do
   else
      if(iproc == 0)then
         write(6,*)'***************************************************************'
         write(6,*)'*        Using UNEQUAL BINS                                   *'
         write(6,*)'***************************************************************'
      end if
      do i = 1, num_comp                  !    loop over nuclei
         Ex_max = nucleus(i)%Ex_max       !   max energy
         S_part = nucleus(i)%sep_e(1)
         E_cut = nucleus(i)%level_ecut
         Ex_calc = min(Ex_max,S_part) - E_cut - de   !   start at E_cut
         nbin = max(int(Ex_calc/de),1)               !   number of bins
         dep = Ex_calc/real(nbin,kind=8)             !   new energy step
!-----    Now that we have dep, count how many bins will be needed
         Ex_min = E_cut + 0.5d0*dep
         Ex_base = nucleus(i)%sep_e(iproj) + 2.0d0
         nbin = 0
         delta_e = delta_e_value(Ex_min, Ex_base, dep)
         delta_e_old = delta_e
         finished = .false.
         e_grid = E_cut - 0.5d0*delta_e
!----   Cycle until e_grid > Ex_max
         do while(.not. finished)
            e_grid = e_grid + 0.5d0*delta_e_old + 0.5d0*delta_e
            nbin = nbin + 1
            delta_e_old = delta_e
            delta_e = delta_e_value(e_grid, ex_base, dep)
            if(e_grid > Ex_max)finished = .true.
         end do
         nucleus(i)%nbin = nbin
         if(nbin > 0)then
            allocate(nucleus(i)%e_grid(nbin))
            allocate(nucleus(i)%delta_e(nbin))
         end if
         delta_e = delta_e_value(Ex_min, Ex_base, dep)
         delta_e_old = delta_e
         e_grid = E_cut - 0.5d0*delta_e
!-------------------   Set up energy grids 
         do j = 1, nucleus(i)%nbin
            e_grid = e_grid + 0.5d0*delta_e_old + 0.5d0*delta_e
            nucleus(i)%e_grid(j) = e_grid
            nucleus(i)%delta_e(j) = delta_e
            delta_e_old = delta_e
            delta_e = delta_e_value(e_grid, ex_base, dep)
         end do
         nucleus(i)%Ex_max = nucleus(i)%e_grid(nbin)
      end do
   end if

   return
end subroutine Setup_Energy_bins
