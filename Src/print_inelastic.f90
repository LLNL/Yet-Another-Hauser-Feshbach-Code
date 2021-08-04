!
!*****************************************************************************80
!
subroutine print_inelastic(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,     &
                           num_energies, Ang_L_max, max_jx_50, delta_jx_50,                &
                           cs_threshold, nstates, absorption_cs, Inelastic_cs,             &
                           Inelastic_Ang_L, Inelastic_L_max, Inelastic_cont_cs,            &
                           reaction_cs, write_error)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write all the inelastic cross sections to the library 
!    directory
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        constants
!        nodeinfo
!        nuclei
!        Channel_info
!        particles_def
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        poly
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
   use particles_def
   implicit none  
   integer(kind=4), intent(in):: itarget
   integer(kind=4), intent(in):: istate
   integer(kind=4), intent(in):: ilab
   character(len=100), intent(in) :: file_lab
   integer(kind=4), intent(in):: ilib_dir
   character(len=132), intent(in) :: lib_dir
   character(len=1), intent(in) :: ch_par(0:1)
   integer(kind=4), intent(in) :: num_energies
   integer(kind=4), intent(in) :: Ang_L_max
   integer(kind=4), intent(in) :: max_jx_50
   real(kind=8), intent(in) :: delta_Jx_50
   real(kind=8), intent(in) :: cs_threshold
   integer(kind=4), intent(in) :: nstates
   real(kind=8), intent(in) :: absorption_cs(num_energies)
   real(kind=8), intent(in) :: Inelastic_cs(0:nstates,num_energies)
   real(kind=8), intent(in) :: Inelastic_Ang_L(0:Ang_L_max,0:nstates,num_energies)
   integer(kind=4), intent(in) :: Inelastic_L_max(0:nstates,num_energies)
   real(kind=8), intent(in) :: reaction_cs(num_energies)
   real(kind=8), intent(in) :: Inelastic_cont_cs(0:nstates,num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   real(kind=8) :: Sum_Inelastic
   integer(kind=4) :: iproj
   integer(kind=4) :: inuc
   integer(kind=4) :: i, j, n
   integer(kind=4) :: ipi, ipf, in, jx, L, Max_L
   real(kind=8) :: x
   real(kind=8) :: e_in
   real(kind=8) :: xnorm
   real(kind=8) :: Temp
   real(kind=8) :: alf, bet
   real(kind=8) :: sum
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   character(len=132) :: directory
   character(len=132) :: outfile
   character(len=1) :: quote
   real(kind=8) :: thresh
   logical :: print_cs
   real(kind=8), allocatable :: Ang_Dist(:)
!----------------------------------------------------------------------
   real(kind=8) :: poly
!----------------------------------------------------------------------

   write_error = .false.
   quote = "'"
   iproj = projectile%particle_type

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   idir = idir + 1
   directory(idir:idir) = particle(iproj)%label
   idir = idir + 1
   directory(idir:idir) = '/'

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Inelastic to discrete states    ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   open(unit=100, file = directory(1:idir)//'check_inelastic.dat',status = 'unknown')
!   do in = 1, num_energies
!      sum = 0.0d0
!      do j = 1, nucleus(itarget)%num_discrete
!         sum = sum + Inelastic_cont_cs(j,in)
!      end do
!      e_in = projectile%energy(in)
!      write(100,'(1x,4(3x,1pe16.7))')e_in, sum*cs_scale, Inelastic_cs(0,in)*cs_scale
!   end do
!   close(unit=100)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Now inelastic transitions to continuous bins with first     ---+
!---------   gamma decay to discrete state, cross section on entry into  ---+
!---------   the knwon discrete spectrum                                 ---+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do j = 1, nucleus(itarget)%num_discrete

      ifile = 26
      outfile(1:8) = 'Channel_'
      outfile(9:9) = particle(iproj)%label
      outfile(10:26) = '_cs_entry_to_L000'
      l = j - 1
      if(l < 10)then
         write(outfile(ifile:ifile),'(i1)')l
      elseif(l >= 10 .and. l < 100)then
         write(outfile(ifile-1:ifile),'(i2)')l
      elseif(l >= 100 .and. l < 1000)then
         write(outfile(ifile-2:ifile),'(i3)')l
      else
         write(6,*)'too many discrete states inelastic_cs ',l
         write_error = .true.
         return
      end if
!**********************************************************************************
!-----   Check if this state has appreciable cross section. If not, skip     -----*
!-----   Threshold is based on the number of samples, and the smallest       -----*
!-----   cross section observable with one sample, i.e., abs_cs/num_mc_samp  -----*
!-----   We use 10 times this as a cut off. Cross ection below this will be  -----*
!-----   unreliable due to limitations in Monte Carlo sampling               -----*
!**********************************************************************************
      print_cs = .false.
      do in = 1, num_energies
         thresh = 20.0d0*absorption_cs(in)/real(num_mc_samp,kind=8)
         if(Inelastic_cont_cs(j,in) > thresh)print_cs = .true.
      end do
      if(.not. print_cs)cycle
!**********************************************************************************

      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
      file_lab2(1:20) = ' '
      ilab2 = ilab
      file_lab2(1:ilab2) = file_lab(1:ilab)
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = particle(iproj)%label
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = quote
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = ')'
      write(100,'(''# '',a20)')file_lab2
      ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
      write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
            nucleus(itarget)%state(istate)%energy
      write(100,'(''# Inelastic cross section.'')')
      write(100,'(''# Contains cross section of inelastic transitions to continuous enrgy bins '')')
      write(100,'(''# on entry into the discrete spectrum.'')')
      write(100,'(''# Note subsequent gamma decays from the final state are not tracked.'')')
      write(100,'(''#'')')
      write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
      write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
      inuc = itarget
      write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
      ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
      write(100,'(''# Final state  = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
      write(100,'(''#'')')
      write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
      write(100,'(''#'',2(''   ----------------''))')
      do in = 1, num_energies
         e_in = projectile%energy(in)
         write(100,'(1x,4(3x,1pe16.7))')e_in, Inelastic_cont_cs(j,in)*cs_scale
      end do
      close(unit=100)
   end do


   if(.not.allocated(Ang_Dist))allocate(Ang_Dist(0:max_jx_50))

   do j = 1, nucleus(itarget)%num_discrete
      if(j == istate)cycle
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    cross section v. E_in           ------------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Start with primary decays of projectile directly to a       ---+
!---------   discrete state                                              ---+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ifile = 17
      outfile(1:ifile) = 'Inelastic_cs_L000'
      l = j - 1
      if(l < 10)then
         write(outfile(ifile:ifile),'(i1)')l
      elseif(l >= 10 .and. l < 100)then
         write(outfile(ifile-1:ifile),'(i2)')l
      elseif(l >= 100 .and. l < 1000)then
         write(outfile(ifile-2:ifile),'(i3)')l
      else
         write(6,*)'too many discrete states inelastic_cs ',l
         write_error = .true.
         return
      end if
!**********************************************************************************
!-----   Check if this state has appreciable cross section. If not, skip     -----*
!-----   Threshold is based on the number of samples, and the smallest       -----*
!-----   cross section observable with one sample, i.e., abs_cs/num_mc_samp  -----*
!-----   We use 10 times this as a cut off. Cross ection below this will be  -----*
!-----   unreliable due to limitations in Monte Carlo sampling               -----*
!**********************************************************************************
      print_cs = .false.
      do in = 1, num_energies
         thresh = 20.0d0*absorption_cs(in)/real(num_mc_samp,kind=8)
         if(Inelastic_cs(j,in) > thresh)print_cs = .true.
      end do
      if(.not. print_cs)cycle
!**********************************************************************************

      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
      file_lab2(1:20) = ' '
      ilab2 = ilab
      file_lab2(1:ilab2) = file_lab(1:ilab)
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = particle(iproj)%label
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = quote
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = ')'
      write(100,'(''# '',a20)')file_lab2
      ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
      write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
            nucleus(itarget)%state(istate)%energy
      write(100,'(''# Inelastic cross section.'')')
      write(100,'(''# Contains only primary decays of projectile particle to discrete final state.'')')
      write(100,'(''# Note subsequent gamma decays from the final state are not tracked.'')')
      write(100,'(''#'')')
      write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
      write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
      inuc = itarget
      write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
      ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
      write(100,'(''# Final state  = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
      write(100,'(''#'')')
      write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
      write(100,'(''#'',2(''   ----------------''))')
      do in = 1, num_energies
         e_in = projectile%energy(in)
         write(100,'(1x,4(3x,1pe16.7))')e_in, Inelastic_cs(j,in)*cs_scale
      end do
      close(unit=100)

      if(xs_only)cycle
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Angular distributions           ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ifile = 18
      outfile(1:ifile) = 'Inelastic_Ang_L000'
      l = j - 1
      if(l < 10)then
         write(outfile(ifile:ifile),'(i1)')l
      elseif(l >= 10 .and. l < 100)then
         write(outfile(ifile-1:ifile),'(i2)')l
      elseif(l >= 100 .and. l < 1000)then
         write(outfile(ifile-2:ifile),'(i3)')l
      else
         write(6,*)'too many discrete states inelastic_ANG ',l
         write_error = .true.
         return
      end if
      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
      ilab2 = ilab
      file_lab2(1:ilab2) = file_lab(1:ilab)
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = particle(iproj)%label
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = quote
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = ')'
      write(100,'(''# '',a20)')file_lab2
      ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
      write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
            nucleus(itarget)%state(istate)%energy
      write(100,'(''# Inelastic Angular Distribution data '')')
      write(100,'(''# Contains only primary decays of the projectile particle to discrete final state.'')')
      write(100,'(''# Note subsequent gamma decays from the final state are not tracked.'')')
      write(100,'(''# The probability density function P(x|E_in) is printed, x = cos(theta)'')')
      write(100,'(''# This distribution is normalized using the trapazoidial rule'')')
      write(100,'(''#'')')
      write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
      write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
      inuc = itarget
      write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
      ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
      write(100,'(''# Final state  = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
      write(100,'(''# Frame = COM (center-of-momentum)'')')
      write(100,'(''#'')')
      do in = 1, num_energies
         e_in = projectile%energy(in)
         write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                 e_in,Inelastic_cs(j,in),cs_units
         if(Inelastic_cs(j,in) < cs_threshold)cycle
         write(100,'(''#         E_in           x=cos(theta)         P(x|E_in)'')')
         write(100,'(''#'',3(''   ----------------''))')
         xnorm = 0.0d0
         do jx = 0, max_jx_50
            x = real(jx,kind=8)*delta_jx_50 - 1.0d0
            if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.99999d0
            alf = 0.0d0
            bet = 0.0d0
            sum = 0.0d0
            do L = 0, Inelastic_L_max(j,in)
               sum = sum + Inelastic_Ang_L(L,j,in)*poly(L,1,alf,bet,x)
            end do
            Ang_Dist(jx) = sum
            if(jx >= 1)xnorm = xnorm + (Ang_Dist(jx-1) + Ang_Dist(jx))*delta_jx_50*0.5d0
         end do
         do jx = 0, max_jx_50
            x = real(jx,kind=8)*delta_jx_50 - 1.0d0
            if(xnorm > 1.0d-8)Ang_Dist(jx) = Ang_Dist(jx)/xnorm
            write(100,'(1x,3(3x,1pe16.7))')e_in, x, Ang_Dist(jx)
         end do
      end do
      close(unit=100)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Angular distributions Legendre Coefficients        -----+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ifile = 22
      outfile(1:ifile) = 'Inelastic_Ang_Leg_L000'
      l = j - 1
      if(l < 10)then
         write(outfile(ifile:ifile),'(i1)')l
      elseif(l >= 10 .and. l < 100)then
         write(outfile(ifile-1:ifile),'(i2)')l
      elseif(l >= 100 .and. l < 1000)then
         write(outfile(ifile-2:ifile),'(i3)')l
      else
         write(6,*)'too many discrete states Inelastic_ANG_Leg ',l
         write_error = .true.
         return
      end if
      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
      ilab2 = ilab
      file_lab2(1:ilab2) = file_lab(1:ilab)
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = particle(iproj)%label
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = quote
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = ')'
      write(100,'(''# '',a20)')file_lab2
      ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
      write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
            nucleus(itarget)%state(istate)%energy
      write(100,'(''# Inelastic Angular Distribution data '')')
      write(100,'(''# Contains only primary decays to discrete final state.'')')
      write(100,'(''# Note subsequent gamma decays from the final state are not tracked.'')')
      write(100,'(''# Coefficients of Legendre expansion'')')
      write(100,'(''# Ang_dist(x) = Sum_L a(L)*P(L,x)'')')
      write(100,'(''#'')')
      write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
      write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
      inuc = itarget
      write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
      ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
      write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
            j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
      write(100,'(''# Frame = COM (center-of-momentum)'')')
      write(100,'(''#'')')

      do in = 1, num_energies
         e_in = projectile%energy(in)
         write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                 e_in,Inelastic_cs(j,in),cs_units
         if(Inelastic_cs(j,in) < cs_threshold)cycle
         write(100,'(''#         E_in            L            a(L)'')')
         write(100,'(''#   ----------------     ---     ----------------'')')
!         write(100,'(''#'')')
         xnorm = 2.0d0*Inelastic_Ang_L(0,j,in)
         Max_L = 0
         do L = 0, Inelastic_L_max(j,in)
            if(abs(Inelastic_Ang_L(L,j,in)) > 1.0d-6 .and. L > Max_L)Max_L = L
         end do
         do L = 0, Max_L
            if(xnorm > 1.0d-8)Temp = Inelastic_Ang_L(L,j,in)/xnorm
            write(100,'(1x,3x,1pe16.7,3x,i5,3x,1pe16.7)')e_in, L, Temp
         end do 
      end do
      close(unit=100)
!------------------------------------------------------------------------------------
   end do
   if(allocated(Ang_Dist))deallocate(Ang_Dist)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    cross section v. E_in           ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ifile = 18
      outfile(1:ifile) = 'Inelastic_cs_Total'
!**********************************************************************************
!-----   Check if this state has appreciable cross section. If not, skip     -----*
!-----   Threshold is based on the number of samples, and the smallest       -----*
!-----   cross section observable with one sample, i.e., abs_cs/num_mc_samp  -----*
!-----   We use 10 times this as a cut off. Cross ection below this will be  -----*
!-----   unreliable due to limitations in Monte Carlo sampling               -----*
!**********************************************************************************
!**********************************************************************************

      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
      file_lab2(1:20) = ' '
      ilab2 = ilab
      file_lab2(1:ilab2) = file_lab(1:ilab)
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = particle(iproj)%label
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = quote
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = ')'
      write(100,'(''# '',a20)')file_lab2
      ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
      write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
            nucleus(itarget)%state(istate)%energy
      write(100,'(''# Inelastic cross section summed over all inelastic channels'')')
      write(100,'(''#'')')
      write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
      write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
      inuc = itarget
      write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
      ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
      write(100,'(''#'')')
      write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
      write(100,'(''#'',2(''   ----------------''))')
      do in = 1, num_energies
         e_in = projectile%energy(in)

         Sum_Inelastic = 0.0d0
         do j = 1, nucleus(itarget)%num_discrete
            if(j == istate)cycle
            Sum_Inelastic = Sum_Inelastic + Inelastic_cs(j,in)   !   already has units of b
!   write(6,*)j,in,Inelastic_cs(J,in), Sum_Inelastic
         end do
         do i = 1, num_channels
            inuc = Exit_Channel(i)%Final_nucleus
            if(inuc == itarget)then
               do n = 1, Exit_Channel(i)%num_cs
                  Sum_Inelastic = Sum_Inelastic +                         &
                     Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)  ! does not yet have units
!   write(6,*)Sum_Inelastic, Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)
               end do
            end if
         end do

         write(100,'(1x,4(3x,1pe16.7))')e_in, Sum_Inelastic*cs_scale
      end do
      close(unit=100)

   return
end subroutine print_inelastic
