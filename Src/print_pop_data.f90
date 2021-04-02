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
subroutine print_pop_data(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,   &
                          num_energies, decay_mult, decay_var, write_error)
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
   real(kind=8), intent(in) :: decay_mult(0:6,num_energies)
   real(kind=8), intent(in) :: decay_var(0:6,num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: k
   integer(kind=4) :: ipi, in
   real(kind=8) :: e_in
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

   ifile = 14
   outfile(1:ifile) = 'Multiplicities'
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
   write(100,'(''# Multiplicity data '')')
   write(100,'(''#''9x,''E_in'',6x,7(3x,6x,''Avg('',a1,'')'',4x,3x,6x,''Var('',a1,'')'',4x))')         &
       (particle(k)%label,particle(k)%label, k= 0,6)
   write(100,'(''#   ----------------'',14(''   ----------------   ----------------''))')
   do in = 1, num_energies
      e_in = projectile%energy(in)
      write(100,'(1x,15(3x,1pe16.7))')e_in, (decay_mult(k,in),decay_var(k,in),k=0,6)
   end do
   close(unit=100)
   return
end subroutine print_pop_data
