!
!*****************************************************************************80
!
subroutine memory_used
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine calculates and prints out the amount of memory used
!    in the calculation, nucleus arrays, etc.
!
!   Dependencies:
!
!        variable_kinds
!        options
!        constants
!        nodeinfo
!        nuclei
!        Channel_info
!        Scatter_info
!        pre_equilibrium_no_1
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
   use variable_kinds
   use options
   use constants
   use nodeinfo
   use nuclei
   use Channel_info
   use Scatter_info
   use pre_equilibrium_no_1
   implicit none
!------------------------------------------------------------------------
   integer(kind=4) :: icomp
   integer(kind=4) :: i, j, ip, n, nn, nbin, k, lx
   integer(kind=4) :: inuc
   integer(kind=4) :: num_s
   real(kind=8) :: mem
   real(kind=8) :: mem_tot
   real(kind=8) :: mem_icomp
   real(kind=8) :: mem_bins
   real(kind=8) :: mem_preeq
   real(kind=8) :: mem_channels
   real(kind=8) :: mem_HF
   real(kind=8) :: mem_scatter
!
!-------------------------------------------------------------------------+
!---------                                                                +
!---------      Cacluate memory being used for each compound nucleus      +
!---------                                                                +
!-------------------------------------------------------------------------+
      mem_tot=0.0d0
      mem_icomp = 0.0d0
      mem_preeq = 0.0d0
      mem_channels = 0.0d0
      mem_scatter = 0.0d0
      mem_HF = 0.0d0
      if(print_me)write(6,*)
      if(print_me)write(6,'(''*************************************************************'')')
      do icomp = 1, num_comp
         mem_icomp = 0.0d0
         mem_bins = 0.0d0
         if(allocated(nucleus(icomp)%PREEQ_pop))then
            mem = real(size(nucleus(icomp)%PREEQ_pop)*8,kind=8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%PREEQ_cs))then
            mem = real(size(nucleus(icomp)%PREEQ_cs)*8,kind=8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%PREEQ_part_cs))then
            mem = real(size(nucleus(icomp)%PREEQ_part_cs)*8,kind=8)
            mem_icomp = mem_icomp + mem
         end if
!---   Electric strength function data 
         if(allocated(nucleus(icomp)%EL_mode))then
            do lx = 1, e_l_max
               mem = mem + 8.0d0
               do k = 1, max_num_gsf
                  mem = mem + 32.0d0 + 8.0d0
                  mem_icomp = mem_icomp + real(nucleus(icomp)%EL_mode(lx)%gsf(k)%num_data,kind=8)*16.0d0
               end do
            end do
         end if
!---   Magnetic strength function data 
         if(allocated(nucleus(icomp)%ML_mode))then
            do lx = 1, m_l_max
               mem = mem + 8.0d0
               do k = 1, max_num_gsf
                  mem = mem + 32.0d0 + 8.0d0
                  mem_icomp = mem_icomp + real(nucleus(icomp)%ML_mode(lx)%gsf(k)%num_data,kind=8)*16.0d0
               end do
            end do
         end if

!-------   Size of derived type bins

         mem_bins = 0.0d0
         if(allocated(nucleus(icomp)%bins))then
            nbin = nucleus(icomp)%nbin
            do n = 1, nbin
               mem_bins = mem_bins + 8.0d0
               mem_bins = mem_bins + 8.0d0
               do ip = 0, 1
                  do j = 0, nucleus(icomp)%j_max
                     mem_icomp = mem_icomp + 8.0d0                                    !   rho
                     mem_icomp = mem_icomp + 8.0d0                                    !   pop
                     mem_icomp = mem_icomp + 8.0d0                                    !   HF_den
                     mem_icomp = mem_icomp + 4.0d0                                    !   num_decay
                     mem_bins = mem_bins + real(size(nucleus(icomp)%bins(j,ip,n)%decay_to)*4,kind=8)
                     mem_bins = mem_bins + real(size(nucleus(icomp)%bins(j,ip,n)%decay_particle)*4,kind=8)
                     mem_bins = mem_bins + real(size(nucleus(icomp)%bins(j,ip,n)%HF_prob)*8,kind=8)
                     mem_bins = mem_bins + real(size(nucleus(icomp)%bins(j,ip,n)%HF_trans)*8,kind=8)
!-----     Data attached to Nuke_decay
                     do nn = 1, nucleus(icomp)%num_decay + 1         !  nuke_decay element
                        mem_bins = mem_bins + 4.0d0                 !  num_decay
                        if(allocated(nucleus(icomp)%bins(j,ip,n)%nuke_decay(nn)%decay_prob))               &
                           mem_bins = mem_bins + real(size(nucleus(icomp)%bins(j,ip,n)%nuke_decay(nn)%decay_prob)*8,kind=8)
                        if(allocated(nucleus(icomp)%bins(j,ip,n)%nuke_decay(nn)%decay_list))               &
                           mem_bins = mem_bins + real(size(nucleus(icomp)%bins(j,ip,n)%nuke_decay(nn)%decay_list)*4,kind=8)
                     end do
                  end do
               end do
            end do
         end if
         mem_icomp = mem_icomp

         mem_tot = mem_tot + mem_icomp + mem_bins
         mem_HF = mem_HF + mem_bins
         if(print_me)write(6,'(''Memory used in data for compound nucleus#'',i3,''      = '',f14.2,'' kbytes'')')  &
            icomp,mem_icomp/1.0d3
         if(print_me)write(6,'(''Memory used in HF decays for compound nucleus#'',i3,'' = '',f14.2,'' kbytes'')')  &
            icomp,mem_bins/1.0d3
         if(print_me)write(6,'(''Total memory used for compound nucleus#'',i3,''        = '',f14.2,'' kbytes'')')  &
            icomp,(mem_icomp+mem_bins)/1.0d3
         if(print_me)write(6,'(''--------------------------------------------------------------------'')')
      end do
      if(print_me)write(6,'(''Memory used for all HF decay arrays     = '',f14.3,'' kbytes'')')mem_HF/1.0d3
      if(print_me)write(6,'(''Total memory used for all nuclei arrays = '',f14.3,'' kbytes'')')mem_tot/1.0d3
!-------------------------------------------------------------------------+
!---  Pre-equilibrium arrays
!-------------------------------------------------------------------------+
      mem_preeq = 0.0d0
      if(.not. pop_calc .and. PREEQ_Model == 1)then
         mem_preeq = mem_preeq + real(size(dWk)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Wk)*8,kind=8)
         mem_preeq = mem_preeq + real(size(W)*8,kind=8)
         mem_preeq = mem_preeq + real(size(tau)*8,kind=8)
         mem_preeq = mem_preeq + real(size(taup)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Pre_prob)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Gam_p)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Gam_n)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Gam_pp)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Gam_pn)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Gam_0_pn)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Gam_0_np)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Lamb_p_p)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Lamb_p_n)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Lamb_0_pn)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Lamb_0_np)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Vwell)*8,kind=8)
         mem_preeq = mem_preeq + real(size(Vwell_g)*8,kind=8)
         if(print_me)write(6,'(''Memory used for pre-equilibrium model   = '',f14.2,'' kbytes'')')mem_preeq/1.0d3
      end if

      mem_tot = mem_tot + mem_preeq
!-------------------------------------------------------------------------+
!-----   Memory in Exit_Channel derived type
!-------------------------------------------------------------------------+
      mem_channels = 0.0d0
      if(allocated(Exit_Channel))then
         do i = 1, num_channels
            inuc = Exit_Channel(i)%Final_nucleus
            num_s = Exit_Channel(i)%num_cs
            if(allocated(Exit_Channel(i)%Channel_cs))                                            &
               mem_channels = mem_channels + real(size(Exit_Channel(i)%Channel_cs)*8,kind=8)
            if(allocated(Exit_Channel(i)%StateLabel))                                            &
               mem_channels = mem_channels + real(size(Exit_Channel(i)%StateLabel)*4,kind=8)
            if(allocated(Exit_Channel(i)%decay_particles))                                       &
               mem_channels = mem_channels + real(size(Exit_Channel(i)%decay_particles)*4,kind=8)
            if(track_gammas)then
               do n = 1, nucleus(inuc)%num_discrete
                  if(allocated(Exit_Channel(i)%state(n)%cs))                                     &
                     mem_channels = mem_channels + real(size(Exit_channel(i)%state(n)%cs)*8,kind=8)
               end do
            end if
            do j = -1, num_s
               do k = 0, 6
                  if(allocated(Exit_Channel(i)%Spect(k,j)%E_spec))                               &
                     mem_channels = mem_channels + real(size(Exit_Channel(i)%Spect(k,j)%E_spec)*8,kind=8)
                  if(allocated(Exit_Channel(i)%Spect(k,j)%E_count))                              &
                     mem_channels = mem_channels + real(size(Exit_Channel(i)%Spect(k,j)%E_count)*8,kind=8)
                  if(allocated(Exit_Channel(i)%Spect(k,j)%Ang_Dist))                             &
                     mem_channels = mem_channels + real(size(Exit_Channel(i)%Spect(k,j)%Ang_Dist)*8,kind=8)
                  if(allocated(Exit_Channel(i)%Spect(k,j)%Ang_L))                                &
                     mem_channels = mem_channels + real(size(Exit_Channel(i)%Spect(k,j)%Ang_L)*8,kind=8)
                  if(allocated(Exit_Channel(i)%Spect(k,j)%E_Ang_Dist))                           &
                     mem_channels = mem_channels + real(size(Exit_Channel(i)%Spect(k,j)%E_Ang_Dist)*8,kind=8)
                  if(allocated(Exit_Channel(i)%Spect(k,j)%E_Ang_L))                              &
                     mem_channels = mem_channels + real(size(Exit_Channel(i)%Spect(k,j)%E_Ang_L)*8,kind=8)
                  if(allocated(Exit_Channel(i)%Spect(k,j)%E_Ang_L_max))                          &
                     mem_channels = mem_channels + real(size(Exit_Channel(i)%Spect(k,j)%E_Ang_L_max)*4,kind=8)
               end do
            end do
         end do    
      end if
      if(print_me)write(6,'(''Memory used for Exit_Channel arrays     = '',f14.2,'' kbytes'')')mem_channels/1.0d3
!-------------------------------------------------------------------------+
!------   Memory in Inelastic arrays
!-------------------------------------------------------------------------+
      mem_scatter = 0.0d0
      if(allocated(Inelastic_cs))mem_scatter = mem_scatter + real(size(Inelastic_cs)*8,kind=8)
      if(allocated(Inelastic_count))mem_scatter = mem_scatter + real(size(Inelastic_count)*4,kind=8)
      if(allocated(Inelastic_L_max))mem_scatter = mem_scatter + real(size(Inelastic_L_max)*4,kind=8)
      if(allocated(Inelastic_Ang_L))mem_scatter = mem_scatter + real(size(Inelastic_Ang_L)*8,kind=8)
      if(allocated(Inelastic_Ang_dist))mem_scatter = mem_scatter + real(size(Inelastic_Ang_dist)*8,kind=8)
      if(allocated(Inelastic_Total))mem_scatter = mem_scatter + real(size(Inelastic_Total)*8,kind=8)
      if(print_me)write(6,'(''Memory used for Inelastic arrays        = '',f14.2,'' kbytes'')')mem_scatter/1.0d3
      if(allocated(direct_cs))mem_scatter = mem_scatter + real(size(direct_cs)*8,kind=8)
      if(allocated(direct_Ang))mem_scatter = mem_scatter + real(size(direct_Ang)*8,kind=8)
      if(allocated(direct_tot))mem_scatter = mem_scatter + real(size(direct_tot)*8,kind=8)
      if(allocated(direct_cc))mem_scatter = mem_scatter + real(size(direct_cc)*8,kind=8)
      if(allocated(direct_dwba))mem_scatter = mem_scatter + real(size(direct_dwba)*8,kind=8)
      if(allocated(SE_prob))mem_scatter = mem_scatter + real(size(SE_prob)*8,kind=8)
      if(allocated(Elastic_cs))mem_scatter = mem_scatter + real(size(Elastic_cs)*8,kind=8)
      if(allocated(Elastic_Ang))mem_scatter = mem_scatter + real(size(Elastic_Ang)*8,kind=8)
      if(allocated(SE_cs))mem_scatter = mem_scatter + real(size(SE_cs)*8,kind=8)
      if(allocated(SE_Ang))mem_scatter = mem_scatter + real(size(SE_Ang)*8,kind=8)
      if(allocated(CE_cs))mem_scatter = mem_scatter + real(size(CE_cs)*8,kind=8)
      if(allocated(CE_Ang))mem_scatter = mem_scatter + real(size(CE_Ang)*8,kind=8)
      if(allocated(QE_cs))mem_scatter = mem_scatter + real(size(QE_cs)*8,kind=8)
      if(allocated(QE_Spectrum))mem_scatter = mem_scatter + real(size(QE_Spectrum)*8,kind=8)
      if(allocated(QE_Ang))mem_scatter = mem_scatter + real(size(QE_Ang)*8,kind=8)
      if(print_me)write(6,'(''Memory used for all Scattering arrays   = '',f14.2,'' kbytes'')')mem_scatter/1.0d3
     



      mem_tot = mem_tot + mem_channels + mem_scatter
      if(print_me)write(6,'(''--------------------------------------------------------------------'')')

      if(print_me)write(6,'(''Memory for this calculation is at least = '',f14.3,'' kbytes'')')mem_tot/1.0d3
      if(print_me)write(6,'(''*************************************************************'')')
      if(print_me)write(6,*)

   return
end subroutine memory_used
