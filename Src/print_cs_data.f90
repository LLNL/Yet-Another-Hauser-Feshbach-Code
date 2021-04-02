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
subroutine print_total_cs(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,   &
                          num_energies, reaction_cs, SE_cs, write_error)
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
   
   return

end subroutine print_total_cs
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
                             num_energies, reaction_cs, write_error)
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
!---------    Reaction cross section          ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   iproj = projectile%particle_type

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

   return

end subroutine print_reaction_cs
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
subroutine print_absorption_cs(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,   &
                               num_energies, absorption_cs, write_error)
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
   real(kind=8), intent(in) :: absorption_cs(num_energies)
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

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   iproj = projectile%particle_type
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
end subroutine print_absorption_cs
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
                            num_nuc, Fiss_J_avg, Fiss_J_var, Fiss_tally, write_error)
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
   integer(kind=4), intent(in) :: num_nuc
   real(kind=8), intent(inout) :: Fiss_J_avg(num_nuc,num_energies)
   real(kind=8), intent(inout) :: Fiss_J_var(num_nuc,num_energies)
   real(kind=8), intent(inout) :: Fiss_tally(num_nuc,num_energies)
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
