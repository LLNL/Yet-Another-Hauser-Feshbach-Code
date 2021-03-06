!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine checks that the directory structure for the output
!    libraries and events exists, and if not creates them
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
subroutine check_directories(ntar, target_label, ilib_dir, lib_dir)
   use variable_kinds
   use options
   use constants
   use nodeinfo
   use nuclei
   use Channel_info
   use particles_def
   implicit none
   integer(kind=4), intent(in) :: ntar
   character(len=5), intent(in) :: target_label
   integer(kind=4), intent(out) :: ilib_dir
   character(len=132), intent(out) :: lib_dir   
!--------------------------------------------------------------------------
   integer(kind=4) :: i
   integer(kind=4) :: if_check, ilast
   integer(kind=4) ::icmd
   logical :: f_exist
   character(len=132) :: file_check
   character(len=132) :: unix_cmd
!--------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------------------   Set up directories for output libraries

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   write(6,*)
   write(6,*)'Checking on directory structure for Data Libraries'
   f_exist = .false.
   lib_dir(1:ntar) = target_label(1:ntar)
   ilib_dir = ntar
   inquire(file = lib_dir(1:ilib_dir), exist = f_exist)
   if(.not. f_exist)then
      icmd = 6
      unix_cmd(1:icmd) = 'mkdir '
      icmd = icmd + 1
      unix_cmd(icmd:icmd + ilib_dir) = lib_dir(1:ilib_dir)
      icmd = icmd + ilib_dir
      write(6,*)unix_cmd(1:icmd)
      call system(unix_cmd(1:icmd))
   end if
   ilib_dir = ilib_dir + 1
   lib_dir(ilib_dir:ilib_dir) = '/'
   if(.not. pop_calc)then
      ilib_dir = ilib_dir + 1
      lib_dir(ilib_dir:ilib_dir) = particle(projectile%particle_type)%label
   else
      ilib_dir = ilib_dir + 1
      lib_dir(ilib_dir:ilib_dir+2) = 'Pop'
      ilib_dir = ilib_dir + 2
   end if
   inquire(file = lib_dir(1:ilib_dir), exist = f_exist)
   if(.not. f_exist)then
      icmd = 6
      unix_cmd(1:icmd) = 'mkdir '
      icmd = icmd + 1
      unix_cmd(icmd:icmd + ilib_dir) = lib_dir(1:ilib_dir)
      icmd = icmd + ilib_dir
      write(6,*)unix_cmd(1:icmd)
      call system(unix_cmd(1:icmd))
   end if
   if(.not. event_generator)then
      if(fission)then
         if_check = ilib_dir
         file_check(1:ilib_dir) = lib_dir
         if_check = if_check + 1
         file_check(if_check:if_check) = '/'
         if_check = if_check + 1
         file_check(if_check:if_check) = 'f'
         inquire(file = file_check(1:if_check), exist = f_exist)
         if(.not. f_exist)then
            icmd = 6
            unix_cmd(1:icmd) = 'mkdir '
            icmd = icmd + 1
            unix_cmd(icmd:icmd + if_check) = file_check(1:if_check)
            icmd = icmd + if_check
            write(6,*)unix_cmd(1:icmd)
            call system(unix_cmd(1:icmd))
         end if
      end if
      do i = 1, num_channels
         if_check = ilib_dir
         file_check(1:ilib_dir) = lib_dir
         if_check = if_check + 1
         file_check(if_check:if_check) = '/'
         if_check = if_check + 1
         ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
         file_check(if_check:if_check+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         if_check = if_check + ilast
         inquire(file = file_check(1:if_check), exist = f_exist)
         if(.not. f_exist)then
            icmd = 6
            unix_cmd(1:icmd) = 'mkdir '
            icmd = icmd + 1
            unix_cmd(icmd:icmd + if_check) = file_check(1:if_check)
            icmd = icmd + if_check
            write(6,*)unix_cmd(1:icmd)
             call system(unix_cmd(1:icmd))
         end if
      end do
   end if
   if_check = ilib_dir
   file_check(1:ilib_dir) = lib_dir
   if_check = if_check + 1
   file_check(if_check:if_check) = '/'
   if_check = if_check + 1
   ilast = 6
   file_check(if_check:if_check+ilast) = 'Spectra'
   if_check = if_check + ilast
   inquire(file = file_check(1:if_check), exist = f_exist)
   if(.not. f_exist)then
      icmd = 6
      unix_cmd(1:icmd) = 'mkdir '
      icmd = icmd + 1
      unix_cmd(icmd:icmd + if_check) = file_check(1:if_check)
      icmd = icmd + if_check
      write(6,*)unix_cmd(1:icmd)
      call system(unix_cmd(1:icmd))
   end if
!-----   Direct printing of event to their own directory


   if(dump_events)then
      if_check = ilib_dir
      file_check(1:ilib_dir) = lib_dir
      if_check = if_check + 1
      file_check(if_check:if_check) = '/'
      if_check = if_check + 1
      ilast = 10
      file_check(if_check:if_check+ilast) = 'Event-files'
      if_check = if_check + ilast
      inquire(file = file_check(1:if_check), exist = f_exist)
      if(.not. f_exist)then
         icmd = 6
         unix_cmd(1:icmd) = 'mkdir '
         icmd = icmd + 1
         unix_cmd(icmd:icmd + if_check) = file_check(1:if_check)
         icmd = icmd + if_check
         write(6,*)unix_cmd(1:icmd)
         call system(unix_cmd(1:icmd))
      end if
   end if
   return
end subroutine check_directories
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write data on nuclear states: properties and how they 
!    decay to the library directory
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
subroutine print_nuke_data(num_comp, ilib_dir, lib_dir)
   use variable_kinds
   use options
   use constants
   use nodeinfo
   use nuclei
   use Channel_info
   use particles_def
   implicit none
   integer(kind=4), intent(in) :: num_comp
   integer(kind=4), intent(in):: ilib_dir
   character(len=132), intent(in) :: lib_dir
!----------------------------------------------------------------------
   integer(kind=4) :: i, iA, k, m
   integer(kind=4) :: idir, ifile
   character(len=132) :: directory
   character(len=132) :: outfile
   character(len=14) iso_label

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Core Directory     -------------------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Decay pattern for each for each nucleus   --------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   do i = 1,num_comp
      iA = nucleus(i)%A
      ifile = 0
      if(iA < 10)then
         ifile = 1
         write(outfile(1:ifile),'(i1)')iA
      elseif(iA < 99)then
         ifile = 2
         write(outfile(1:ifile),'(i2)')iA
      elseif(iA < 1000)then
         ifile = 3
         write(outfile(1:ifile),'(i3)')iA
      end if
      if(nucleus(i)%atomic_symbol(1:1) == ' ')then
         outfile(ifile+1:ifile+1) = nucleus(i)%atomic_symbol(2:2)
         ifile = ifile + 1
      else
         outfile(ifile+1:ifile+2) = nucleus(i)%atomic_symbol(1:2)
         ifile = ifile + 2
      end if
      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'-decay-properties.dat',status = 'unknown')            
      write(100,'(''#Decay properties of discrete states for'',a5)')outfile(1:ifile)
      write(100,'(''#Mass = '',1pe23.16,'' MeV'')')nucleus(i)%mass
      write(100,'(''#Mass = '',1pe23.16,'' amu'')')nucleus(i)%mass/mass_u
      write(100,'(''#AMU = '',1pe23.16,'' MeV'')')mass_u
      write(100,'(''#'')')
      do k = 1, nucleus(i)%num_discrete, 1
         write(100,'(''    #       Ex     J   par #decay     Isomer'')')
         write(100,'(''  ---     ------  --- ----  -----  --------------'')')
         iso_label = 'Isomer = False'
         if(nucleus(i)%state(k)%isomer)iso_label = 'Isomer = True '
         write(100,'(1x,i4,1x,f10.4,2(1x,f4.1),1x,i5,3x,a14)')k-1,nucleus(i)%state(k)%energy,    &
                     nucleus(i)%state(k)%spin,nucleus(i)%state(k)%parity,                        &
                     nucleus(i)%state(k)%nbranch, iso_label
         write(100,'(''#      i --->   f            branch      prob_gamma         prob_ic'')')
         write(100,'(''#    ---      ---     -------------   -------------   -------------'')')
         do m = 1, nucleus(i)%state(k)%nbranch
            write(100,'(4x,i4,'' --->'',i4,4(3x,1pe15.7))')k-1,nucleus(i)%state(k)%ibranch(m)-1, &
                                                           nucleus(i)%state(k)%branch(m),        &
                                                           nucleus(i)%state(k)%p_gamma(m),       &
                                                           nucleus(i)%state(k)%p_ic(m)
         end do
         write(100,'(''#'')')
      end do
   end do
  return
end subroutine print_nuke_data
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write reaction, absorption,and total cross sections 
!    to the library directory
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
subroutine print_reaction_cs(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,   &
                             num_energies, reaction_cs, absorption_cs, SE_cs, write_error)
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
   real(kind=8), intent(in) :: reaction_cs(num_energies)
   real(kind=8), intent(in) :: absorption_cs(num_energies)
   real(kind=8), intent(in) :: SE_cs(num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: ipi, in
   real(kind=8) :: e_in, cs
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   character(len=132) :: directory
   character(len=132) :: outfile
!----------------------------------------------------------------------

   write_error = .false.

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Total cross section          ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   iproj = projectile%particle_type
   if(iproj == 1)then
      ifile = 8
      outfile(1:ifile) = 'Total_cs'
      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
      file_lab2(1:20) = ' '
      ilab2 = ilab
      file_lab2(1:ilab2) = file_lab(1:ilab)
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = 'X'
      ilab2 = ilab2 + 1
      file_lab2(ilab2:ilab2) = ')'
      write(100,'(''# '',a20)')file_lab2
      ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
      write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
            istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                        &
            nucleus(itarget)%state(istate)%energy
      write(100,'(''# Total Cross Section data '')')
      write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
      write(100,'(''#'',2(''   ----------------''))')
      do in = 1, num_energies
         e_in=projectile%energy(in)
         cs = (reaction_cs(in) + SE_cs(in))*cs_scale
         write(100,'(1x,4(3x,1pe16.7))')e_in, cs
      end do
      close(unit=100)
   end if
         
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Reaction cross section          ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ifile = 11
   outfile(1:ifile) = 'Reaction_cs'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = 'X'
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
         nucleus(itarget)%state(istate)%energy
   write(100,'(''# Reaction Cross Section data = Total - Elastic'')')
   write(100,'(''# Reaction Cross Section data = Absorption + direct'')')
   write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
   write(100,'(''#'',2(''   ----------------''))')
   do in = 1, num_energies
      e_in=projectile%energy(in)
      cs = reaction_cs(in)*cs_scale
      write(100,'(1x,4(3x,1pe16.7))')e_in, cs
   end do
   close(unit=100)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Absorption cross section          ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ifile = 13
   outfile(1:ifile) = 'Absorption_cs'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = 'X'
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
         nucleus(itarget)%state(istate)%energy
   write(100,'(''# Absorption Cross Section data'')')
   write(100,'(''# Absorption Cross Section data = Reaction - Direct - Compound Elastic'')')
   write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
   write(100,'(''#'',2(''   ----------------''))')
   do in = 1, num_energies
      e_in=projectile%energy(in)
      cs = absorption_cs(in)*cs_scale
         write(100,'(1x,4(3x,1pe16.7))')e_in, cs
   end do
   close(unit=100)
   return
end subroutine print_reaction_cs
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write pre-equilibrium cross section to the library 
!    directory
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
subroutine print_preeq_cs(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,   &
                          num_energies, reaction_cs, preeq_css, write_error)
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
   real(kind=8), intent(in) :: reaction_cs(num_energies)
   real(kind=8), intent(in) :: preeq_css(0:6,1:num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: k, ipi, in
   real(kind=8) :: e_in, cs
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   character(len=132) :: directory
   character(len=132) :: outfile
!----------------------------------------------------------------------
   write_error = .false.
   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   ifile = 17
   outfile(1:ifile) = 'Preequilibrium_cs'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = 'X'
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
         nucleus(itarget)%state(istate)%energy
   write(100,'(''# Preequilibrium Cross section. Units = '',a2)')
   write(100,'(''#'',3x,6x,''E_in'',6x,3x,4x,''Tot('',a2,'')'',7(3x,6x,a1,''('',a2,'')'',5x))')  &
      cs_units,(particle(k)%label,cs_units,k = 0,6)
   write(100,'(''#'',9(''   ----------------''))')
   do in = 1, num_energies
      e_in = projectile%energy(in)
      cs = 0.0d0
      do k = 0, 6
         cs = cs + preeq_css(k,in)
      end do
      cs = cs*reaction_cs(in)*cs_scale
      write(100,'(1x,9(3x,1pe16.7))')e_in, cs, (preeq_css(k,in)*reaction_cs(in)*cs_scale, k = 0, 6)
   end do
   close(unit=100)
   return
end subroutine print_preeq_cs
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write direct, coupled-channles cross section to the library 
!    directory
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
subroutine print_direct_cs(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,   &
                           num_energies, direct_cc, direct_dwba, direct_tot, write_error)
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
   real(kind=8), intent(in) :: direct_cc(1:num_energies)
   real(kind=8), intent(in) :: direct_dwba(1:num_energies)
   real(kind=8), intent(in) :: direct_tot(1:num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: ipi, in
   real(kind=8) :: e_in
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   character(len=132) :: directory
   character(len=132) :: outfile
!----------------------------------------------------------------------
   write_error = .false.
   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   ifile = 9
   outfile(1:ifile) = 'Direct_cs'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = 'X'
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                           &
         nucleus(itarget)%state(istate)%energy
   write(100,'(''# Direct Cross Section data'')')
   write(100,'(''# Reaction Cross Section data'')')
   write(100,1100)cs_units,cs_units,cs_units
 1100    format('#         E_in      ',10x,'cc(',a2,')',3x,8x,'dwba(',a2,')',3x,              &
          9x,'tot(',a2,')')
   write(100,'(''#'',4(''   ----------------''))')
   do in = 1, num_energies
      e_in = projectile%energy(in)
!      cs = absorption_cs(in)*cs_scale
      write(100,'(1x,4(3x,1pe16.7))')e_in, direct_cc(in)*cs_scale, direct_dwba(in)*cs_scale,  &
              direct_tot(in)*cs_scale
   end do
   close(unit=100)
   return
end subroutine print_direct_cs 
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write fission cross section to the library 
!    directory
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
subroutine print_fission_cs(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,   &
                            num_energies, reaction_cs, fission_cs,                        &
                            num_comp, Fiss_J_avg, Fiss_J_var, Fiss_tally, write_error)
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
   real(kind=8), intent(in) :: reaction_cs(num_energies)
   real(kind=8), intent(in) :: fission_cs(num_energies)
   integer(kind=4), intent(in) :: num_comp
   real(kind=8), intent(inout) :: Fiss_J_avg(num_comp,num_energies)
   real(kind=8), intent(inout) :: Fiss_J_var(num_comp,num_energies)
   real(kind=8), intent(inout) :: Fiss_tally(num_comp,num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: ipi, in, icomp
   real(kind=8) :: e_in
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   character(len=132) :: directory
   character(len=132) :: outfile
!----------------------------------------------------------------------
   write_error = .false.

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   idir = idir + 1
   directory(idir:idir) = 'f'
   idir = idir + 1
   directory(idir:idir) = '/'
   ifile = 10
   outfile(1:ifile) = 'fission_cs'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = 'f'
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                           &
         nucleus(itarget)%state(istate)%energy
   write(100,'(''# Fission cross section data '')')
   write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
   write(100,'(''#'',2(''   ----------------''))')
   do in = 1, num_energies
      e_in=projectile%energy(in)
      write(100,'(1x,4(3x,1pe16.7))')e_in, fission_cs(in)*reaction_cs(in)*cs_scale
   end do
   close(unit=100)
!--------------   Output data on fission <J> <J**2>
   open(unit=24,file=directory(1:idir)//'Fission-Stats.dat',status='unknown')
   do in = 1, num_energies
      E_in = projectile%energy(in)
      do icomp = 1, num_comp
         if(Fiss_tally(icomp,in) > 0.0d0)then
            Fiss_J_avg(icomp,in) = Fiss_J_avg(icomp,in)/Fiss_tally(icomp,in)
            Fiss_J_var(icomp,in) = Fiss_J_var(icomp,in)/Fiss_tally(icomp,in)
            Fiss_J_var(icomp,in) = sqrt(abs(Fiss_J_var(icomp,in) - Fiss_J_avg(icomp,in)**2))
         end if
      end do
      write(24,'(25(1x,f10.5))')E_in,(Fiss_J_avg(icomp,in),Fiss_J_var(icomp,in),icomp = 1, num_comp)
   end do
   close(unit=24)
   return
end subroutine print_fission_cs
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write all the elastic cross sections to the library 
!    directory
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
subroutine print_elastic(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,     &
                         num_energies, Ang_L_max, max_jx_100, delta_Jx_100,              &
                         cs_threshold, SE_cs, SE_Ang, Elastic_cs, Elastic_Ang,           &
                         nstates, Inelastic_cs, Inelastic_Ang_L, Inelastic_L_max, write_error)
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
   integer(kind=4), intent(in) :: max_jx_100
   real(kind=8), intent(in) :: delta_Jx_100
   real(kind=8), intent(in) :: cs_threshold
   real(kind=8), intent(in) :: SE_cs(num_energies)
   real(kind=8), intent(in) :: SE_Ang(0:Ang_L_max,num_energies)
   real(kind=8), intent(inout) :: Elastic_cs(num_energies)
   real(kind=8), intent(inout) :: Elastic_Ang(0:Ang_L_max,num_energies)
   integer(kind=4), intent(in) :: nstates
   real(kind=8), intent(in) :: Inelastic_cs(0:nstates,num_energies)
   real(kind=8), intent(in) :: Inelastic_Ang_L(0:Ang_L_max,0:nstates,num_energies)
   integer(kind=4), intent(in) :: Inelastic_L_max(0:nstates,num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: ipi, ipf, in, j, jx, L
   real(kind=8) :: xA
   real(kind=8) :: x
   real(kind=8) :: e_in, cs
   real(kind=8) :: xnorm
   real(kind=8) :: Temp
   real(kind=8) :: alf, bet
   real(kind=8) :: sum
   real(kind=8) :: comp
   real(kind=8) :: P_L
   real(kind=8) :: shape, Coul, Sig_C
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   character(len=132) :: directory
   character(len=132) :: outfile
   real(kind=8), allocatable :: xvalue(:)
   real(kind=8), allocatable :: Ang_Dist(:)
!----------------------------------------------------------------------
   real(kind=8) :: poly
!----------------------------------------------------------------------
   write_error = .false.
   iproj = projectile%particle_type

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   idir = idir + 1
   directory(idir:idir) = particle(iproj)%label
   idir = idir + 1
   directory(idir:idir) = '/'
                                      ! don't have with population calculation
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Shape Elastic scattering cross section  ----------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ifile = 16
   outfile(1:ifile) = 'Shape_Elastic_cs'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = particle(iproj)%label
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   j = istate
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                           &
         nucleus(itarget)%state(istate)%energy
   if(iproj == 1)write(100,'(''# Shape Elastic cross section data '')')
   if(iproj > 1)write(100,'(''# <Shape Elastic/Rutherford> '')')
   write(100,'(''#'')')
   write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
   write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   if(iproj > 0 .and. iproj < 7)then
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
   else
      write(100,'(''# Population decay, no projectile'')')
   end if
   write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
         j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   if(iproj == 1)then
      write(100,'(''#          E_in              xs'',''('',a2,'')'')')cs_units
   elseif(iproj > 1)then
      write(100,'(''#          E_in         <OM/Rutherford>'')')
   end if
   write(100,'(''#'',2(''   ----------------''))')
   do in = 1, num_energies
      cs = SE_cs(in)
      e_in = projectile%energy(in)
      if(iproj == 1)then
         write(100,'(1x,4(3x,1pe16.7))')e_in, cs*cs_scale
      elseif(iproj > 1)then
         write(100,'(1x,4(3x,1pe16.7))')e_in, cs*0.5d0
      end if
   end do
   close(unit=100)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------   Compound Elastic crosss section     -------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ifile = 19
   outfile(1:ifile) = 'Compound_Elastic_cs'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = particle(iproj)%label
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   j = istate
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
         nucleus(itarget)%state(istate)%energy
   write(100,'(''# Compound Elastic cross section data '')')
   write(100,'(''#'')')
   write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
   write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   if(iproj > 0 .and. iproj < 7)then
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
   else
      write(100,'(''# Population decay, no projectile'')')
   end if
   write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
      j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''#            E_in              xs'',''('',a2,'')'')')cs_units
   write(100,'(''#'',2(''   ----------------''))')
   do in = 1, num_energies
      cs = Inelastic_cs(istate,in)
      e_in=projectile%energy(in)
      write(100,'(1x,4(3x,1pe16.7))')e_in, cs*cs_scale
   end do
   close(unit=100)

   if(xs_only)return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------     Angular distributions     ---------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(.not.allocated(Ang_Dist))allocate(Ang_Dist(0:max_jx_100))
   ifile = 17
   outfile(1:ifile) = 'Shape_Elastic_Ang'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = particle(iproj)%label
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   j = istate
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
         nucleus(itarget)%state(istate)%energy
   if(iproj == 1)write(100,'(''# Shape Elastic Angular Distribution data '')')
   if(iproj == 1)write(100,'(''# Angular Distribution normalized to unity '')')
   if(iproj > 1)write(100,'(''# Elastic Scattering Angular Distribution data '')')
   if(iproj > 1)write(100,'(''# Calculated as ratio to Rutherford: OM/Rutherford '')')
   write(100,'(''#'')')
   write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
   write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   if(iproj > 0 .and. iproj < 7)then
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
   else
      write(100,'(''# Population decay, no projectile'')')
   end if
   write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
            j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''# Frame = COM'')')
   alf = 0.0d0
   bet = 0.0d0
   do in = 1, num_energies
      e_in = projectile%energy(in)
      if(iproj == 1)then
         write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                    e_in,SE_cs(in),cs_units
      elseif(iproj > 1)then
         write(100,'(''# E_in = '',1pe16.7,3x,''OM/Rutherford = '',1pe16.7)')                    &
                    e_in,SE_cs(in)
      end if
      if(iproj == 1)then
         write(100,'(''#         E_in            cos(theta)            Prob'')')
      elseif(iproj > 1)then
         write(100,'(''#         E_in            cos(theta)        OM/Rutherford'')')
      end if
      write(100,'(''#'',3(''   ----------------''))')
      do jx = 0, max_jx_100
         x = real(jx,kind=8)*delta_jx_100 - 1.0d0
         if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.999845
         alf = 0.0d0
         bet = 0.0d0
         sum = 0.0d0
         do L = 0, Ang_L_max
            sum = sum + SE_Ang(L,in)*poly(L,1,alf,bet,x)
         end do
         Ang_Dist(jx) = sum
      end do
      xnorm = 0.0d0
      do jx = 1, max_jx_100
         xnorm = xnorm + (Ang_Dist(jx-1) + Ang_Dist(jx))*delta_jx_100*0.5d0
      end do
      do jx = 0, max_jx_100
         x = real(jx,kind=8)*delta_jx_100 - 1.0d0
         Temp = Ang_Dist(jx)
         if(iproj == 1 .and. xnorm > 1.0d-20)Temp = Temp/xnorm
         write(100,'(1x,3(3x,1pe16.7))')e_in, x, Temp
      end do
   end do
   close(unit=100)
   if(allocated(Ang_Dist))deallocate(Ang_Dist)
!
!---   Also print out Legendre Coefficients
!
   ifile = 21
   outfile(1:ifile) = 'Shape_Elastic_Ang_Leg'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = particle(iproj)%label
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   j = istate
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
         nucleus(itarget)%state(istate)%energy
   if(iproj == 1)write(100,'(''# Shape Elastic Angular Distribution data '')')
   if(iproj == 1)write(100,'(''# Angular Distribution normalized to unity '')')
   if(iproj == 1)write(100,'(''# Angular Distribution normalized to unity '')')
   if(iproj > 1)write(100,'(''# Elastic Scattering Angular Distribution data '')')
   if(iproj > 1)write(100,'(''# Calculated as ratio to Rutherford: OM/Rutherford '')')
   write(100,'(''# Coefficients of Legendre expansion'')')
   write(100,'(''# Ang_dist(x) = Sum_L a(L)*P(L,x)'')')
   write(100,'(''#'')')
   write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
   write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   if(iproj > 0 .and. iproj < 7)then
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
   else
      write(100,'(''# Population decay, no projectile'')')
   end if
   write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
            j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''# Frame = COM'')')

   do in = 1, num_energies
      e_in = projectile%energy(in)
      if(iproj == 1)then
         write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                    e_in,SE_cs(in),cs_units
      elseif(iproj > 1)then
         write(100,'(''# E_in = '',1pe16.7,3x,''OM/Rutherford = '',1pe16.7)')                   &
                          e_in,SE_cs(in)
      end if
      write(100,'(''#         E_in            L            a(L)'')')
      write(100,'(''#   ----------------     ---     ----------------'')')
      xnorm = 2.0d0*SE_Ang(0,in)
      do L = 0, Ang_L_max
         Temp = SE_Ang(L,in)
         if(iproj == 1)Temp = SE_Ang(L,in)/xnorm
         write(100,'(1x,3x,1pe16.7,3x,i5,3x,1pe16.7)')e_in, L, Temp
      end do 
   end do
   close(unit=100)

   if(.not.allocated(Ang_Dist))allocate(Ang_Dist(0:max_jx_100))
   ifile = 20
   outfile(1:ifile) = 'Compound_Elastic_Ang'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = particle(iproj)%label
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   j = istate
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
         nucleus(itarget)%state(istate)%energy
   write(100,'(''# Compound Elastic Angular Distribution data '')')
   write(100,'(''# Distribution normalized to unity'')')
   write(100,'(''#'')')
   write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
   write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   if(iproj > 0 .and. iproj < 7)then
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
   else
      write(100,'(''# Population decay, no projectile'')')
   end if
   write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
      j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''# Frame = COM'')')
   alf = 0.0d0
   bet = 0.0d0
   do in = 1, num_energies
      e_in = projectile%energy(in)
      write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')    &
                 e_in,Inelastic_cs(istate,in),cs_units
      if(Inelastic_cs(istate,in) < cs_threshold)cycle
      write(100,'(''#         E_in            cos(theta)            Prob'')')
      write(100,'(''#'',3(''   ----------------''))')
      xnorm = 0.0d0
      do jx = 0, max_jx_100
         x = real(jx,kind=8)*delta_jx_100 - 1.0d0
         if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.999845
         alf = 0.0d0
         bet = 0.0d0
         sum = 0.0d0
         do L = 0, Inelastic_L_max(istate,in)
            sum = sum + Inelastic_Ang_L(L,istate,in)*poly(L,1,alf,bet,x)
         end do
         Ang_Dist(jx) = sum
         if(jx >= 1)xnorm = (Ang_Dist(jx-1) + Ang_Dist(jx))*delta_jx_100*0.5d0
      end do
      do jx = 0, max_jx_100
         x = real(jx,kind=8)*delta_jx_100 - 1.0d0
         Temp = Ang_Dist(jx)
         if(xnorm > 1.0d-7)Temp = Temp/xnorm
         write(100,'(1x,3(3x,1pe16.7))')e_in, x, Temp
      end do
   end do
   close(unit=100)
   if(allocated(Ang_Dist))deallocate(Ang_Dist)
!
!---   Also print out Legendre Coefficients
!
   ifile = 24
   outfile(1:ifile) = 'Compound_Elastic_Ang_Leg'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = particle(iproj)%label
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   j = istate
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
         nucleus(itarget)%state(istate)%energy
   if(iproj == 1)write(100,'(''# Compound Elastic Angular Distribution data '')')
   if(iproj == 1)write(100,'(''# Angular Distribution normalized to unity '')')
   write(100,'(''# Coefficients of Legendre expansion'')')
   write(100,'(''# Ang_dist(x) = Sum_L a(L)*P(L,x)'')')
   write(100,'(''#'')')
   write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
   write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   if(iproj > 0 .and. iproj < 7)then
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
   else
      write(100,'(''# Population decay, no projectile'')')
   end if
   write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
            j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''# Frame = COM'')')

   do in = 1, num_energies
      e_in = projectile%energy(in)
      write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                 e_in,Inelastic_cs(istate,in),cs_units
      if(Inelastic_cs(istate,in) < cs_threshold)cycle
      if(iproj == 1)then
         write(100,'(''#         E_in            L            a(L)'')')
      elseif(iproj > 1)then
         write(100,'(''#         E_in            L         a(L)('',a2,''/sr)'')')cs_units
      end if
      write(100,'(''#   ----------------     ---     ----------------'')')
      xnorm = 2.0d0*Inelastic_Ang_L(0,istate,in)
      do L = 0, Inelastic_L_max(istate,in)
         if(xnorm > 1.0d-8)Temp = Inelastic_Ang_L(L,istate,in)/xnorm
         if(iproj > 1)Temp = Temp*Inelastic_cs(istate,in)
         write(100,'(1x,3x,1pe16.7,3x,i5,3x,1pe16.7)')e_in, L, Temp
      end do 
   end do
   close(unit=100)

!
!----   Angular Distribution for Elastic
!
   if(.not.allocated(Ang_Dist))allocate(Ang_Dist(0:max_jx_100))
   ifile = 11
   outfile(1:ifile) = 'Elastic_Ang'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = particle(iproj)%label
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   j = istate
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate - 1, nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
         nucleus(itarget)%state(istate)%energy
   write(100,'(''# Elastic Angular Distribution data '')')
   if(iproj ==1)write(100,'(''# Distribution normalized to unity'')')
   if(iproj > 1)write(100,'(''# Distribution given as Elastic/Rutherford'')')
   write(100,'(''#'')')
   write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
   write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   if(iproj > 0 .and. iproj < 7)then
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
   else
      write(100,'(''# Population decay, no projectile'')')
   end if
   write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
          j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''# Frame = COM'')')

   xA = (nucleus(itarget)%mass+nucleus(itarget)%state(istate)%energy)/particle(iproj)%mass
   Coul = (particle(iproj)%Z*nucleus(itarget)%Z*                                                &
          (1.0d0+xA)/xA*fine_structure*hbar_c)**2*0.25d0
   if(.not.allocated(xvalue))allocate(xvalue(0:max_jx_100))
   do jx = 0, max_jx_100
      x = real(jx,kind=8)*delta_jx_100 - 1.0d0
      if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.99999d0
      xvalue(jx) = x
   end do

   Elastic_cs(1:num_energies) = 0.0d0
   sig_C = 0.0d0  
   do in = 1, num_energies
      e_in = projectile%energy(in)
!-------   Shape + Compound Elastic
      xnorm = 0.0d0
      do jx = 0, max_jx_100
         x = real(jx,kind=8)*delta_jx_100 - 1.0d0
         if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.99999d0
         if(iproj > 1)sig_C = Coul/e_in**2/(1.0d0 - x)**2*fmsq_eq_barn                ! Rutherford in b/sr
         alf = 0.0d0
         bet = 0.0d0
         shape = 0.0d0
         comp = 0.0d0
         do L = 0, Ang_L_max
            P_L = poly(L,1,alf,bet,x)
            shape = shape + SE_Ang(L,in)*P_L
            comp = comp + Inelastic_Ang_L(L,istate,in)*P_L
         end do
         shape = shape
         comp = comp*Inelastic_cs(istate,in)
         if(iproj > 1)comp = comp/Sig_C
         Ang_Dist(jx) = shape + comp
         if(jx > 0)xnorm = xnorm + (Ang_Dist(jx-1) + Ang_Dist(jx))*delta_jx_100*0.5d0
      end do
      if(iproj == 1)then                   !  neutrons - cross section
         Elastic_Ang(0:Ang_L_max,in) = 0.0d0
         do L = 0, Ang_L_max
            Temp = SE_Ang(L,in) + Inelastic_Ang_L(L,istate,in)*Inelastic_cs(istate,in)
            Elastic_Ang(L,in) = Temp/(SE_cs(in) + Inelastic_cs(istate,in))
         end do
         Elastic_cs(in) = (SE_cs(in) + Inelastic_cs(istate,in))
      elseif(iproj> 1)then                 !  Charged particles - ratio to Rutherford
         call Legendre_expand(max_jx_100+1, xvalue(0), Ang_Dist(0),                             &
                              Ang_L_max, Elastic_Ang(0,in))
         do L = 0, Ang_L_max
            Elastic_Ang(L,in) = Elastic_Ang(L,in)*xnorm
         end do
         Elastic_cs(in) = xnorm*0.5d0
      end if
      if(iproj == 1)then
         write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                    e_in,(SE_cs(in) + Inelastic_cs(istate,in)),cs_units
      elseif(iproj > 1)then
         write(100,'(''# E_in = '',1pe16.7,3x,'' <EL/Rutherford> = '',1pe16.7)')                &
                    e_in,xnorm
      end if
      if(iproj == 1)then
         write(100,'(''#         E_in            cos(theta)            Prob'')')
      elseif(iproj > 1)then
         write(100,'(''#         E_in            cos(theta)        EL/Rutherford'')')
      end if
      write(100,'(''#'',3(''   ----------------''))')
      do jx = 0, max_jx_100
         x = real(jx,kind=8)*delta_jx_100 - 1.0d0
         if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.99999d0
         Temp = Ang_Dist(jx)
         if(iproj == 1 .and. xnorm > 1.0d-8)Temp = Temp/xnorm
         if(iproj > 1 )Temp = Temp
         write(100,'(1x,3(3x,1pe16.7))')e_in, x, Temp
      end do
   end do
   close(unit=100)
   if(allocated(Ang_Dist))deallocate(Ang_Dist)
   if(allocated(xvalue))deallocate(xvalue)
!----
!----   Moved this below angular distribution to collect ratio for charged particles
!----
   ifile = 10
   outfile(1:ifile) = 'Elastic_cs'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = particle(iproj)%label
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   j = istate
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
         nucleus(itarget)%state(istate)%energy
   write(100,'(''# Total Elastic cross section data '')')
   if(iproj == 1)write(100,'(''# Total Elastic cross section data '')')
   if(iproj > 1)write(100,'(''# <Total Elastic/Rutherford> '')')
   if(iproj > 1)write(100,'(''# Integrated over d(cos(theta)) '')')
   write(100,'(''#'')')
   write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
   write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   if(iproj > 0 .and. iproj < 7)then
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
   else
      write(100,'(''# Population decay, no projectile'')')
   end if
   write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
         j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   if(iproj == 1)then
      write(100,'(''#          E_in              xs'',''('',a2,'')'')')cs_units
   elseif(iproj > 1)then
      write(100,'(''#          E_in         <EL/Rutherford>'')')
   end if
   write(100,'(''#'',2(''   ----------------''))')
   do in = 1, num_energies
      cs = Elastic_cs(in)
      if(iproj == 1) cs = cs*cs_scale
      e_in = projectile%energy(in)
      write(100,'(1x,4(3x,1pe16.7))')e_in, cs
   end do
   close(unit=100)

!
!---   Also print out Legendre Coefficients
!
   ifile = 15
   outfile(1:ifile) = 'Elastic_Ang_Leg'
   open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
   file_lab2(1:20) = ' '
   ilab2 = ilab
   file_lab2(1:ilab2) = file_lab(1:ilab)
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = particle(iproj)%label
   ilab2 = ilab2 + 1
   file_lab2(ilab2:ilab2) = ')'
   j = istate
   write(100,'(''# '',a20)')file_lab2
   ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         j-1, nucleus(itarget)%state(j)%spin, ch_par(ipi),                                       &
         nucleus(itarget)%state(j)%energy
   if(iproj == 1)write(100,'(''# Elastic Angular Distribution data '')')
   if(iproj == 1)write(100,'(''# Angular Distribution normalized to unity '')')
   if(iproj > 1)write(100,'(''# Elastic Scattering Angular Distribution data '')')
   if(iproj > 1)write(100,'(''# Calculated as ratio to Rutherford: OM/Rutherford '')')
   write(100,'(''# Coefficients of Legendre expansion'')')
   write(100,'(''# Ang_dist(x) = Sum_L a(L)*P(L,x)'')')
   write(100,'(''#'')')
   write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
   write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   if(iproj > 0 .and. iproj < 7)then
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
   else
      write(100,'(''# Population decay, no projectile'')')
   end if
   write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
   ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
      j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''# Frame = COM'')')

   if(iproj == 1)then
      write(100,'(''#         E_in            L            a(L)'')')
   elseif(iproj > 1)then
      write(100,'(''#         E_in            L            a(L)'')')
   end if
   write(100,'(''#   ----------------     ---     ----------------'')')
   do in = 1, num_energies
      e_in = projectile%energy(in)
      write(100,'(''#'')')
      do L = 0, Ang_L_max
         write(100,'(1x,3x,1pe16.7,3x,i5,3x,1pe16.7)')e_in, L, Elastic_Ang(L,in)
      end do 
   end do
   close(unit=100)
   return
end subroutine print_elastic
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write all the inelastic cross sections to the library 
!    directory
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
subroutine print_inelastic(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,     &
                           num_energies, Ang_L_max, max_jx_50, delta_jx_50,                &
                           cs_threshold, nstates, absorption_cs, Inelastic_cs,             &
                           Inelastic_Ang_L, Inelastic_L_max, write_error)
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
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: inuc
   integer(kind=4) :: ipi, ipf, in, j, jx, L
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

   if(.not.allocated(Ang_Dist))allocate(Ang_Dist(0:max_jx_50))
   do j = 1, nucleus(itarget)%num_discrete
      if(j == istate)cycle
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    cross section v. E_in           ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
      write(100,'(''# Inelastic cross section'')')
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
      write(100,'(''#'')')
      write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
      write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
      write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
      inuc = itarget
      write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
      ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
      write(100,'(''# Final state  = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
      write(100,'(''# Frame = COM'')')
      write(100,'(''#'')')
      do in = 1, num_energies
         e_in = projectile%energy(in)
         write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                 e_in,Inelastic_cs(j,in),cs_units
         if(Inelastic_cs(j,in) < cs_threshold)cycle
         write(100,'(''#         E_in            cos(theta)            Prob'')')
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
      write(100,'(''# Frame = COM'')')
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
         do L = 0, Inelastic_L_max(j,in)
            if(xnorm > 1.0d-8)Temp = Inelastic_Ang_L(L,j,in)/xnorm
            write(100,'(1x,3x,1pe16.7,3x,i5,3x,1pe16.7)')e_in, L, Temp
         end do 
      end do
      close(unit=100)
!------------------------------------------------------------------------------------
   end do
   if(allocated(Ang_Dist))deallocate(Ang_Dist)
   return
end subroutine print_inelastic
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write exit channel reaction data the library 
!    directory
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
subroutine print_channels(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,     &
                          num_energies, num_e, max_jx_20, delta_jx_20,                    &
                          de_spec, cs_threshold, reaction_cs, write_error)
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
   integer(kind=4), intent(in) :: num_e
   integer(kind=4), intent(in) :: max_jx_20
   real(kind=8), intent(in) :: delta_Jx_20
   real(kind=8), intent(in) :: de_spec
   real(kind=8), intent(in) :: cs_threshold
   real(kind=8), intent(in) :: reaction_cs(num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: inuc
   integer(kind=4) :: ipi, ipf, in, i
   integer(kind=4) :: j, jx
   integer(kind=4) :: k
   integer(kind=4) :: m
   integer(kind=4) :: n, nf, nn
   integer(kind=4) :: L, L_max
   integer(kind=4) :: icc, icc_max
   real(kind=8) :: cs
   real(kind=8) :: ratio
   real(kind=8) :: Q
   real(kind=8) :: x
   real(kind=8) :: e_in, e_out
   real(kind=8) :: xnorm
   real(kind=8) :: Temp
   real(kind=8) :: alf, bet
   real(kind=8) :: sum
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   integer(kind=4) :: ilast
   character(len=132) :: directory
   character(len=132) :: outfile
   character(len=5) :: resid_label
   integer(kind=4) :: nres

   real(kind=8), allocatable :: Temp_Ang_Dist(:,:,:)
   real(kind=8), allocatable :: Ang_Dis(:,:)
!----------------------------------------------------------------------
   real(kind=8) :: poly
!----------------------------------------------------------------------

   write_error = .false.

   iproj = projectile%particle_type

   if(.not.allocated(Temp_Ang_Dist))allocate(Temp_Ang_Dist(0:max_jx_20,0:num_e,0:6))
   do i = 1, num_channels
      inuc = Exit_Channel(i)%Final_nucleus
      directory(1:ilib_dir) = lib_dir(1:ilib_dir)
      idir = ilib_dir + 1
      directory(idir:idir) = '/'
      ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
      directory(idir+1:idir+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
      idir = idir + ilast
      idir = idir + 1
      directory(idir:idir) = '/'

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------    Channel reaction data   -------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do n = 1, Exit_Channel(i)%num_cs
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------    Channel cross section data   --------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         j = Exit_Channel(i)%StateLabel(n)
         Q = Exit_Channel(i)%Q_value - nucleus(inuc)%state(j)%energy +                             &
             nucleus(itarget)%state(istate)%energy

         ifile = 8
         outfile(1:ifile) = 'Channel_'
         outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ifile = ifile + ilast
         outfile(ifile+1:ifile+8) = '_cs_L000'
         ifile = ifile + 8
         l = j - 1
         if(l < 10)then
            write(outfile(ifile:ifile),'(i1)')l
         elseif(l >= 10 .and. l < 100)then
            write(outfile(ifile-1:ifile),'(i2)')l
         elseif(l >= 100 .and. l < 1000)then
            write(outfile(ifile-2:ifile),'(i3)')l
         else
            write(6,*)'too many discrete states Channel_cs ',l
            write_error = .true.
            return
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         ilast = index(Exit_Channel(i)%Channel_Label,' ') - 1
         file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ilab2 = ilab2 + ilast
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''# To level #'',i3)')l
!-----------------------------------
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         inuc = Exit_Channel(i)%Final_nucleus
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
         ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Xs collected in Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
         write(100,'(''#'')')
         if(Exit_channel(i)%num_particles > 0)then
            if(explicit_channels)then
               write(100,'(''# Mass of emitted particles '')')
               do m = 1, Exit_channel(i)%num_particles
                  k = Exit_channel(i)%decay_particles(m)
                  write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                  write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')               &
                     m, particle(k)%mass/mass_u
               end do
            else
               write(100,'(''# Mass of emitted particles '')')
               m = 0
               do k = 1, 6
                  do nn = 1, Exit_channel(i)%num_part(k)
                     m = m + 1
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                        m,particle(k)%mass/mass_u
                  end do
               end do
            end if
         end if
         write(100,'(''#'')')
         write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
         write(100,'(''#'')')
         write(100,'(''# Cross section data '')')
         write(100,1902)cs_units,(particle(k)%label, k = 0, 6)
 1902    format('#         E_in              xs(',a2,')',5x,7(8x,'Mult(',a1,')',4x))
         write(100,'(''#'',9(''   ----------------''))')
         do in = 1, num_energies
            e_in=projectile%energy(in)
            cs = Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)*cs_scale
            write(100,'(1x,9(3x,1pe16.7))')e_in, cs, (Exit_Channel(i)%part_mult(k,n,in),k=0,6)
         end do
         close(unit=100)

         if(xs_only)cycle

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Angular distribution and spectrum  ---------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------   Output P(Eout,x|Ein)   -----------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ifile = 8
         outfile(1:ifile) = 'Channel_'
         outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ifile = ifile + ilast
         outfile(ifile+1:ifile+20) = '_Eout_Ang_G_Ein_L000'
         ifile = ifile + 20
         l = j - 1
         if(l < 10)then
            write(outfile(ifile:ifile),'(i1)')l
         elseif(l >= 10 .and. l < 100)then
            write(outfile(ifile-1:ifile),'(i2)')l
         elseif(l >= 100 .and. l < 1000)then
            write(outfile(ifile-2:ifile),'(i3)')l
         else
            write(6,*)'too many discrete states Channel_Eout_ANG_G_Ein ',l
            stop
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
         file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ilab2 = ilab2 + ilast
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         inuc = Exit_Channel(i)%Final_nucleus
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
         ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Xs collected in Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
         write(100,'(''#'')')

         if(Exit_channel(i)%num_particles > 0)then
            if(explicit_channels)then
               write(100,'(''# Mass of emitted particles '')')
               do m = 1, Exit_channel(i)%num_particles
                  k = Exit_channel(i)%decay_particles(m)
                  write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                  write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')               &
                     m, particle(k)%mass/mass_u
               end do
            else
               write(100,'(''# Mass of emitted particles '')')
               m = 0
               do k = 1, 6
                  do nn = 1, Exit_channel(i)%num_part(k)
                     m = m + 1
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                        m,particle(k)%mass/mass_u
                  end do
                end do
            end if
         end if
         write(100,'(''#'')')
         write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
         write(100,'(''#'')')
         write(100,'(''# Angular Distribution data  P(x,Eout|Ein) '')')
         write(100,'(''# Normalized to unity over integral Int(-1,1) int(0,E_max)  P(x,Eout|Ein)dx,dEout '')')
         write(100,'(''# dx = '',1pe16.7,'' dEout = '',1pe16.7,'' MeV'')')delta_jx_20, de_spec
         write(100,'(''# Frame = COM'')')
         do in = 1, num_energies
            e_in = projectile%energy(in)
            if(e_in < -Q)cycle
            cs = Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)
            write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                    e_in,cs,cs_units
            if(cs < cs_threshold)cycle
            write(100,'(''#'')')
            write(100,1910)(particle(k)%label, k = 0,6)
 1910    format('#         E_in               E_out           cos(theta)   ',                      &
                7(3x,(5x,'Prob(',a1,')',4x)))
            write(100,'(''#'',10(''   ----------------''))')
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------    Find max value of icc to print out
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Normalize the particle spectra and find max icc to print out
            Temp_Ang_Dist(0:max_jx_20,0:num_e,0:6) = 0.0d0
            do k = 0, 6
               xnorm = 0.0d0
               do icc = 0, num_e
                  xnorm = xnorm + Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)*de_spec
               end do
               if(xnorm > 1.0d-8)then
                  do icc = 0, num_e
                     Exit_Channel(i)%Spect(k,n,in)%E_spec(icc) =                                   &
                       Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)/xnorm
                  end do
               else
                  do icc = 0, num_e
                     Exit_Channel(i)%Spect(k,n,in)%E_spec(icc) = 0.0d0
                  end do
               end if
            end do
            icc_max = num_e
            do icc = num_e - 1, 0, -1
               sum = 0.0d0
               do k = 0, 6
                  sum = sum + Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)
               end do
               if(sum >= 1.0d-8)then
                  icc_max = icc + 1
                  exit
               end if
            end do
            if(icc_max < 1)icc_max = 1
!------   Create angular distribution array - normalize over angle and multiply by prob 
            alf = 0.0d0
            bet = 0.0d0
            do k = 0, 6
               do icc = 0, icc_max
                  do jx = 0, max_jx_20
                     x = real(jx,kind=8)*delta_jx_20 - 1.0d0
                     sum = 0.0d0
                     do L = 0, Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_Max(icc)
                        sum = sum + Exit_Channel(i)%Spect(k,n,in)%E_Ang_L(L,icc)*                  &
                                    poly(L,1,alf,bet,x)
                     end do
                     Temp_Ang_Dist(jx,icc,k) = sum
                  end do
                  xnorm = 0.0d0
                  do jx = 0, max_jx_20 - 1
                     xnorm = xnorm + (Temp_Ang_Dist(jx,icc,k) + Temp_Ang_Dist(jx+1,icc,k))*        &
                             delta_jx_20*0.5d0
                  end do
                  if(xnorm > 0.0d0)then
                     do jx = 0, max_jx_20
                        Temp_Ang_Dist(jx,icc,k) = Temp_Ang_Dist(jx,icc,k)/xnorm*                   &
                                                  Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)
                     end do
                  end if
               end do
               xnorm = 0.0d0
               temp = 0.0d0
               do icc = 0, icc_max
                  sum = 0.0d0
                  do jx = 0, max_jx_20 - 1
                     sum = sum + (Temp_Ang_Dist(jx,icc,k) + Temp_Ang_Dist(jx+1,icc,k))*            &
                                  0.5d0*delta_jx_20
                  end do
                  if(Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)> 1.0d-8)                            &
                  xnorm = xnorm + sum*de_spec
                  temp = temp + Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)*de_spec
               end do
            end do

!----------   Print out data 
           do icc = 0, icc_max
               e_out = real(icc,kind=8)*de_spec
               if(icc >0)e_out = e_out - de_spec/2.0d0
               do jx = 0, max_jx_20
                  x = real(jx,kind=8)*delta_jx_20 - 1.0d0
                  write(100,'(1x,10(3x,1pe16.7))')e_in, e_out, x,                                  &
                     (Temp_Ang_Dist(jx,icc,k),k = 0,6)
               end do
            end do
          end do
         close(unit=100)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Angular distributions Legendre Coefficients        -----+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ifile = 8
         outfile(1:ifile) = 'Channel_'
         outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ifile = ifile + ilast
         outfile(ifile+1:ifile+24) = '_Eout_Ang_Leg_G_Ein_L000'
         ifile = ifile + 24
         l = j - 1
         if(l < 10)then
            write(outfile(ifile:ifile),'(i1)')l
         elseif(l >= 10 .and. l < 100)then
            write(outfile(ifile-1:ifile),'(i2)')l
         elseif(l >= 100 .and. l < 1000)then
            write(outfile(ifile-2:ifile),'(i3)')l
         else
            write(6,*)'too many discrete states Channel_Eout_ANG_Leg_G_Ein ',l
            stop
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
         file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ilab2 = ilab2 + ilast
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         inuc = Exit_Channel(i)%Final_nucleus
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
         ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Xs collected in Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
         write(100,'(''#'')')
         if(Exit_channel(i)%num_particles > 0)then
            if(explicit_channels)then
               write(100,'(''# Mass of emitted particles '')')
               do m = 1, Exit_channel(i)%num_particles
                  k = Exit_channel(i)%decay_particles(m)
                  write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                  write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')               &
                     m, particle(k)%mass/mass_u
               end do
            else
               write(100,'(''# Mass of emitted particles '')')
               m = 0
               do k = 1, 6
                  do nn = 1, Exit_channel(i)%num_part(k)
                     m = m + 1
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                        m,particle(k)%mass/mass_u
                  end do
                end do
            end if
         end if
         write(100,'(''#'')')
         write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
         write(100,'(''#'')')
         write(100,'(''# Angular Distribution data  P(x,Eout|Ein) '')')
         write(100,'(''# Normalized to unity over integral Int(-1,1) int(0,E_max)  P(x,Eout|Ein)dx,dEout '')')
         write(100,'(''# dx = '',1pe16.7,'' dEout = '',1pe16.7,'' MeV'')')delta_jx_20, de_spec
         write(100,'(''# Frame = COM'')')

         do in = 1, num_energies
            e_in = projectile%energy(in)
            icc_max = num_e
            cs = Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)
            write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                    e_in,cs,cs_units
            if(cs < cs_threshold)cycle
            do icc = num_e - 1, 0, -1
               sum = 0.0d0
               do k = 0, 6
                  sum = sum + Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)
               end do
               if(sum >= 1.0d-8)then
                  icc_max = icc + 1
                  exit
               end if
            end do
            if(icc_max < 1)icc_max = 1
            do icc = 0, icc_max
               E_out = real(icc,kind=8)*de_spec
               if(icc > 0)E_out = E_out - 0.5d0*de_spec
               write(100,'(''#'')')
               write(100,'(''#'',36x,''Prob(Eout)'',7(3x,1pe16.7))')                               &
                  (Exit_Channel(i)%Spect(k,n,in)%E_spec(icc),k=0,6)
               write(100,'(''#'')')
         write(100,2110)(particle(k)%label, k = 0,6)
 2110    format('#         E_in               E_out           L',7(3x,(5x,'a(L)(',a1,')',4x)))
         write(100,'(''#   ----------------   ----------------     ---'',7(3x,''----------------''))')
               L_max = 0
               do k = 0, 6
                  if(Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_Max(icc) > L_max)                       &
                     L_max = Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_Max(icc)
               end do
               do L = 0, L_max
                 write(100,'(1x,3x,1pe16.7,3x,1pe16.7,3x,i5,7(3x,1pe16.7))')e_in, e_out, L,        &
                    (Exit_Channel(i)%Spect(k,n,in)%E_Ang_L(L,icc),k=0,6)
               end do
            end do 
         end do

         close(unit=100)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------   Output P(x|Ein)   ------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(.not.allocated(Ang_dis))allocate(Ang_dis(0:max_jx_20,0:6))
         ifile = 8
         outfile(1:ifile) = 'Channel_'
         outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ifile = ifile + ilast
         outfile(ifile+1:ifile+15) = '_Ang_G_Ein_L000'
         ifile = ifile + 15
         l = j - 1
         if(l < 10)then
            write(outfile(ifile:ifile),'(i1)')l
         elseif(l >= 10 .and. l < 100)then
            write(outfile(ifile-1:ifile),'(i2)')l
         elseif(l >= 100 .and. l < 1000)then
            write(outfile(ifile-2:ifile),'(i3)')l
         else
            write(6,*)'too many discrete states Channel_ANG_G_Ein ',l
            stop
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
         file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ilab2 = ilab2 + ilast
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         inuc = Exit_Channel(i)%Final_nucleus
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
         ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Xs collected in Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
         write(100,'(''#'')')
         if(Exit_channel(i)%num_particles > 0)then
            if(explicit_channels)then
               write(100,'(''# Mass of emitted particles '')')
               do m = 1, Exit_channel(i)%num_particles
                  k = Exit_channel(i)%decay_particles(m)
                  write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                  write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')               &
                     m, particle(k)%mass/mass_u
               end do
            else
               write(100,'(''# Mass of emitted particles '')')
               m = 0
               do k = 1, 6
                  do nn = 1, Exit_channel(i)%num_part(k)
                     m = m + 1
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                        m,particle(k)%mass/mass_u
                  end do
                end do
            end if
         end if
         write(100,'(''#'')')
         write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
         write(100,'(''#'')')
         write(100,'(''# Angular Distribution data  P(x|Ein) '')')
         write(100,'(''# Normalized to unity over integral Int(-1,1) P(x|Ein)dx '')')
         write(100,'(''# dx = '',1pe16.7)')delta_jx_20
         write(100,'(''# Frame = COM'')')
         Ang_dis(0:max_jx_20,0:6) = 0.0d0
         do in = 1, num_energies
            e_in=projectile%energy(in)
            if(e_in < -Q)cycle
            cs = Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)
            write(100,'(''#'')')
            write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                    e_in,cs,cs_units
            if(cs < cs_threshold)cycle
             write(100,1903)(particle(k)%label, k = 0,6)
 1903       format('#         E_in             cos(theta)   ',3x,7(3x,(5x,'Prob(',a1,')',4x)))
            write(100,'(''#'',9(''   ----------------''))')
!----    Calculate Angular distribution again and then integrate over energy
            alf = 0.0d0
            bet = 0.0d0
            do k = 0, 6
               do jx = 0, max_jx_20
                  Ang_Dis(jx,k) = 0.0d0
                  sum = 0.0d0
                  do icc = 0, num_e
                     x = real(jx,kind=8)*delta_jx_20 - 1.0d0
                     do L = 0, Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_Max(icc)
                        sum = sum + Exit_Channel(i)%Spect(k,n,in)%E_Ang_L(L,icc)*                  &
                                    poly(L,1,alf,bet,x)*                                           &
                                    Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)*de_spec
                     end do
                  end do
                  Ang_dis(jx,k) = sum
               end do
               xnorm = 0.0d0
               do jx = 1, max_jx_20
                  xnorm = xnorm + (Ang_dis(jx-1,k)+Ang_dis(jx,k))*delta_jx_20*0.5d0
               end do
               if(xnorm > 0.0d0)then
                  do jx = 0, max_jx_20
                     Ang_dis(jx,k) = Ang_dis(jx,k)/xnorm
                  end do
               end if
            end do
            do jx = 0, max_jx_20
               x = real(jx,kind=8)*delta_jx_20 - 1.0d0
               write(100,'(1x,9(3x,1pe16.7))')e_in, x, (Ang_dis(jx,k),k = 0,6)
            end do
         end do
         close(unit=100)
         if(allocated(Ang_dis))deallocate(Ang_dis)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------   Output P(Eout|Ein)   ------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ifile = 8
         outfile(1:ifile) = 'Channel_'
         outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ifile = ifile + ilast
         outfile(ifile+1:ifile+16) = '_Eout_G_Ein_L000'
         ifile = ifile + 16
         l = j - 1
         if(l < 10)then
            write(outfile(ifile:ifile),'(i1)')l
         elseif(l >= 10 .and. l < 100)then
            write(outfile(ifile-1:ifile),'(i2)')l
         elseif(l >= 100 .and. l < 1000)then
            write(outfile(ifile-2:ifile),'(i3)')l
         else
            write(6,*)'too many discrete states Channel_G_Ein ',l
            stop
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
         file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ilab2 = ilab2 + ilast
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         inuc = Exit_Channel(i)%Final_nucleus
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
         ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Xs collected in Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
         write(100,'(''#'')')
         if(Exit_channel(i)%num_particles > 0)then
            if(explicit_channels)then
               write(100,'(''# Mass of emitted particles '')')
               do m = 1, Exit_channel(i)%num_particles
                  k = Exit_channel(i)%decay_particles(m)
                  write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                  write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')               &
                    m, particle(k)%mass/mass_u
               end do
            else
               write(100,'(''# Mass of emitted particles '')')
               m = 0
               do k = 1, 6
                  do nn = 1, Exit_channel(i)%num_part(k)
                     m = m + 1
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                         m,particle(k)%mass/mass_u
                  end do
               end do
            end if
         end if
         write(100,'(''#'')')
         write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
         write(100,'(''#'')')
         write(100,'(''# Emission Spectrum data  P(Eout|Ein) '')')
         write(100,'(''# Normalized to unity over integral Int(0,Emax) P(Eout|Ein)dEout '')')
         write(100,'(''# dEout = '',1pe16.7)')de_spec
         write(100,'(''# Frame = COM'')')
         do in = 1, num_energies
            e_in=projectile%energy(in)
            if(e_in < -Q)cycle
            write(100,'(''#'')')
            cs = Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)
            write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                    e_in,cs,cs_units
            if(cs < cs_threshold)cycle
         write(100,1904)(particle(k)%label, k = 0,6)
  1904   format('#         E_in               E_out      ',7(3x,(5x,'Prob(',a1,')',4x)))
         write(100,'(''#'',9(''   ----------------''))')
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------    Find max value of icc to print out
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----
            icc_max = num_e
            do icc = num_e - 1, 0, -1
               sum = 0.0d0
               do k = 0, 6
                  sum = sum + Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)
               end do
               if(sum >= 1.0d-8)then
                  icc_max = icc + 1
                  exit
               end if
            end do
            if(icc_max < 1)icc_max = 1
            do icc = 0, icc_max
               e_out = real(icc,kind=8)*de_spec
               if(icc > 0)e_out = e_out - de_spec/2
               write(100,'(1x,9(3x,1pe16.7))')e_in, e_out,                                         &
                  (Exit_Channel(i)%Spect(k,n,in)%E_spec(icc),k = 0,6)
            end do
         end do
         close(unit=100)
      end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------    Channel fission cross section data   ------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(fission)then
         ifile = 8
         outfile(1:ifile) = 'Channel_'
         outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ifile = ifile + ilast
         outfile(ifile+1:ifile+11) = '_fission_cs'
         ifile = ifile + 11
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
         file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ilab2 = ilab2 + ilast
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''# To level = N/A'')')
         write(100,'(''# Cross section data '')')
         write(100,'(''#         E_in              xs('',a2,'')'')')cs_units
         write(100,'(''#'',2(''   ----------------''))')
         do in = 1, num_energies
            e_in = projectile%energy(in)
            cs = Exit_Channel(i)%Channel_cs(0,in)*reaction_cs(in)*cs_scale
            write(100,'(1x,9(3x,1pe16.7))')e_in,cs
         end do
         close(unit=100)
      end if
   end do
   if(allocated(Temp_Ang_Dist))deallocate(Temp_Ang_Dist)

!
!-------------------------------------------------------------------------------------------
!-----   Print out gammas for each channel if track_gammas == .true. for each energy
!-------------------------------------------------------------------------------------------
!
   if(track_gammas)then
      do i = 1, num_channels
         inuc = Exit_Channel(i)%Final_nucleus
         call nucleus_label(inuc,nres,resid_label)
         directory(1:ilib_dir) = lib_dir(1:ilib_dir)
         idir = ilib_dir + 1
         directory(idir:idir) = '/'
         ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
         directory(idir+1:idir+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         idir = idir + ilast
         idir = idir + 1
         directory(idir:idir) = '/'
         Q = Exit_Channel(i)%Q_value
         Q = -Q
         ifile = 8
         outfile(1:ifile) = 'Channel_'
         outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ifile = ifile + ilast
         outfile(ifile+1:ifile+13) = '_Gammas_G_Ein'
         ifile = ifile + 13
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
         file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ilab2 = ilab2 + ilast
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''# Mass amu = '',1pe15.7,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
         write(100,'(''#'')')
         if(Exit_channel(i)%num_particles > 0)then
            if(explicit_channels)then
               write(100,'(''# Mass of emitted particles '')')
               do m = 1, Exit_channel(i)%num_particles
                  k = Exit_channel(i)%decay_particles(m)
                  write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                  write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')               &
                     m, particle(k)%mass/mass_u
               end do
            else
               write(100,'(''# Mass of emitted particles '')')
               m = 0
               do k = 1, 6
                  do nn = 1, Exit_channel(i)%num_part(k)
                     m = m + 1
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                        m,particle(k)%mass/mass_u
                  end do
                end do
         end if
         end if
         write(100,'(''#'')')
         write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
         write(100,'(''#'')')
         if(.not. pop_calc)then
            write(100,'(''# Gamma-ray probabilities for each incident energy '')')
         else
            write(100,'(''# Gamma-ray probabilities for each initial excitation energy '')')
         end if


         if(.not. pop_calc)then
            write(100,'(''# Cross section for discrete gamma transitions in '',a5)')resid_label(1:nres)
            write(100,'(''#'')')
         else
            write(100,'(''# Fraction of decay for discrete gamma transitions in '',a5)')           &
                 resid_label(1:nres)
            write(100,'(''#'')')
         end if

         do k = 1, nucleus(inuc)%num_discrete
            do m = 1, nucleus(inuc)%state(k)%nbranch
               nf  = nucleus(inuc)%state(k)%ibranch(m)
               ipi = nint((nucleus(inuc)%state(k)%parity+1.0d0)/2.0d0)
               ipf = nint((nucleus(inuc)%state(nf)%parity+1.0d0)/2.0d0)

               write(100,'(''# Transitions involving gamma ray #'',i4,'' from state #'',i4)')m,k

               write(100,'(''#   i    Jp     f    Jp     e-gamma (MeV)'')')
               write(100,'(''#  --- ----    --- ----     ------------- '')')
               write(100,1911)k-1, nucleus(inuc)%state(k)%spin,ch_par(ipi),                        &
                              nf-1, nucleus(inuc)%state(nf)%spin, ch_par(ipf),                     &
                              nucleus(inuc)%state(k)%egamma(m)
1911        format('#',i4,1x,f4.1,a1,2x,i4,1x,f4.1,a1,8x,f10.5) 

               if(.not. pop_calc)then
                  write(100,1906)cs_units
               else
                  write(100,1920)
               end if
1906        format('#        E_in    i    Jp     f    Jp     e-gamma (MeV)        sigma (',a2,    &
                   ')     sigma(g)/channel')
1920        format('#        E_in    i    Jp     f    Jp     e-gamma (MeV)        Fraction ',     &
                   '      Fraction/channel')
               write(100,1907)
1907        format('#  ----------   --- ----    --- ----     -------------     -------------     -------------')

               do in = 1, num_energies
                  e_in = projectile%energy(in)
                  cs = Exit_channel(i)%state(k)%cs(m,in)*reaction_cs(in)*cs_scale
                  ratio = 0.0d0
                  sum = 0.0d0
                  do nn = 1, Exit_Channel(i)%num_cs
                     sum = sum + Exit_channel(i)%Channel_cs(nn,in)
                  end do
                  if(sum > 1.0d-10)ratio = cs/sum
1908   format(1x,f12.6,1x,i4,1x,f4.1,a1,'  ',i4,1x,f4.1,a1,8x,f10.5,3x,e15.7,3x,e15.7,3x,e15.7) 
                  write(100,1908)e_in, k-1, nucleus(inuc)%state(k)%spin, ch_par(ipi),              &
                                 nf-1, nucleus(inuc)%state(nf)%spin, ch_par(ipf),                  &
                                 nucleus(inuc)%state(k)%egamma(m),                                 &
                                 cs, ratio
               end do
            end do
         end do
         close(unit = 100)
      end do                               
   end if
   close(unit = 13)
   return
end subroutine print_channels
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write exit channel gamma ray data in the library 
!    directory
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
subroutine print_channel_gammas(itarget, ilab, file_lab, ilib_dir, lib_dir,            &
                                ch_par, in, e_in, reaction_cs, write_error)
   use variable_kinds
   use options
   use constants
   use nodeinfo
   use nuclei
   use Channel_info
   use particles_def
   implicit none  
   integer(kind=4), intent(in):: itarget
   integer(kind=4), intent(in):: ilib_dir
   character(len=132), intent(in) :: lib_dir
   character(len=1), intent(in) :: ch_par(0:1)
   integer(kind=4) :: in
   real(kind=8), intent(in) :: e_in
   real(kind=8), intent(in) :: reaction_cs
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: inuc
   integer(kind=4) :: ipi, ipf, i
   integer(kind=4) :: k
   integer(kind=4) :: m
   integer(kind=4) :: nf, nn
   real(kind=8) :: cs
   real(kind=8) :: ratio
   real(kind=8) :: Q
   real(kind=8) :: sum, check_sum
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   integer(kind=4) :: ilast
   character(len=132) :: directory
   character(len=132) :: outfile
   integer(kind=4) :: ilab
   character(len=100) :: file_lab
   character(len=5) :: target_label
   integer(kind=4) :: ntar
   character(len=5) :: resid_label
   integer(kind=4) :: nres
!----------------------------------------------------------------------
   write_error = .false.
   iproj = projectile%particle_type

   call nucleus_label(itarget,ntar,target_label)
   ilab = ntar
   file_lab(1:ntar) = target_label(1:ntar)
   ilab = ilab + 1
   file_lab(ilab:ilab) = '('
   ilab = ilab + 1
   file_lab(ilab:ilab) = particle(iproj)%label
   ilab = ilab + 1
   file_lab(ilab:ilab) = ','

   do i = 1, num_channels
      inuc = Exit_Channel(i)%Final_nucleus
      call nucleus_label(inuc,nres,resid_label)
      directory(1:ilib_dir) = lib_dir(1:ilib_dir)
      idir = ilib_dir + 1
      directory(idir:idir) = '/'
      ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
      directory(idir+1:idir+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
      idir = idir + ilast
      idir = idir + 1
      directory(idir:idir) = '/'

      Q = Exit_Channel(i)%Q_value
      Q = -Q

      ifile = 8
      outfile(1:ifile) = 'Channel_'
      outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
      ifile = ifile + ilast
      outfile(ifile+1:ifile+5) = '_Ein_'
      ifile = ifile + 5
      if(e_in < 10.0)then
         outfile(ifile+1:ifile+1) ='0'
         ifile = ifile + 1
         write(outfile(ifile+1:ifile+6),'(f6.4)')e_in
         ifile = ifile + 6
      else
         write(outfile(ifile+1:ifile+7),'(f7.4)')e_in
         ifile = ifile + 7
      end if
      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'_cs_gammas.dat',status = 'unknown')
      file_lab2(1:20) = ' '
      ilab2 = ilab
      file_lab2(1:ilab2) = file_lab(1:ilab)
      ilab2 = ilab2 + 1
      ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
      file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
      ilab2 = ilab2 + ilast
      file_lab2(ilab2:ilab2) = ')'
      write(100,'(''# '',a20)')file_lab2
!-----------------------------------
      write(100,'(''#'')')
      write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
      write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
      if(iproj > 0 .and. iproj < 7)then
         write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
      else
         write(100,'(''# Population decay, no projectile'')')
      end if
      write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
      write(100,'(''#'')')
      if(Exit_channel(i)%num_particles > 0)then
         if(explicit_channels)then
            write(100,'(''# Mass of emitted particles '')')
            do m = 1, Exit_channel(i)%num_particles
               k = Exit_channel(i)%decay_particles(m)
               write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
               write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')                  &
                     m, particle(k)%mass/mass_u
            end do
         else
            write(100,'(''# Mass of emitted particles '')')
            m = 0
            do k = 1, 6
               do nn = 1, Exit_channel(i)%num_part(k)
                  m = m + 1
                  write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                  write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')                &
                        m, particle(k)%mass/mass_u
               end do
            end do
         end if
      end if

      write(100,'(''#'')')
      write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
      write(100,'(''#'')')
      write(100,'(''# Listing discrete gammas observed in this particle decay channel'')')
      write(100,'(''# Note these are not exclusive to a specific discrete state'')')
      write(100,'(''# Incident Energy = '',f10.5)')E_in
      sum = 0.0d0

      do nn = 1, Exit_Channel(i)%num_cs
         sum = sum + Exit_channel(i)%Channel_cs(nn,in)
      end do

      if(.not. pop_calc)then
         sum = sum*reaction_cs*cs_scale
         write(100,'(''#'')')
         write(100,'(''# Total cross section to this exit decay channel = '',1pe15.7,1x,a2)')sum,cs_units
         write(100,'(''#'')')
      else
         sum = sum*reaction_cs
         write(100,'(''#'')')
         write(100,'(''# Fraction of decay strength to this exit decay channel = '',1pe15.7,1x,a2)')sum
         write(100,'(''#'')')
      end if


      write(100,'(''# Cross section for discrete gamma transitions in '',a5)')resid_label(1:nres)
      write(100,'(''#'')')
      write(100,1900)cs_units
1900           format('#   i    Jp     f    Jp     e-gamma (MeV)        sigma (',a2,')     sigma(g)/channel')
      write(100,'(''#  --- ----    --- ----     -------------     -------------     -------------'')')

      check_sum = 0.0d0
      do k = 1, nucleus(inuc)%num_discrete
         do m = 1, nucleus(inuc)%state(k)%nbranch
            nf  = nucleus(inuc)%state(k)%ibranch(m)
            ipi = nint((nucleus(inuc)%state(k)%parity+1.0d0)/2.0d0)
            ipf = nint((nucleus(inuc)%state(nf)%parity+1.0d0)/2.0d0)
            cs = Exit_channel(i)%state(k)%cs(m,in)*reaction_cs*cs_scale
            ratio = 0.0d0
            if(sum > 1.0d-10)ratio = Exit_channel(i)%state(k)%cs(m,in)/sum
            check_sum = check_sum + ratio 
            write(100,'(1x,i4,1x,f4.1,a1,''  '',i4,1x,f4.1,a1,8x,f10.5,3x,e15.7,3x,e15.7,3x,e15.7)') &
                                                  k, nucleus(inuc)%state(k)%spin, ch_par(ipi),       &
                                                  nf, nucleus(inuc)%state(nf)%spin, ch_par(ipf),     &
                                                  nucleus(inuc)%state(k)%egamma(m),                  &
                                                  cs, ratio
         end do
      end do
      close(unit = 100)
   end do                               
   return
end subroutine print_channel_gammas

