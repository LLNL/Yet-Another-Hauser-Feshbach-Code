!
!*****************************************************************************80
!
subroutine check_directories(ntar, target_label, ilib_dir, lib_dir)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine checks that the directory structure for the output
!    libraries and events exists, and if not creates them
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
!        system
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
   if(refresh_library_directories)then
      lib_dir(1:ntar) = target_label(1:ntar)
      ilib_dir = ntar
      icmd = 6
      unix_cmd(1:icmd) = 'rm -r '
      icmd = icmd + 1
      unix_cmd(icmd:icmd + ilib_dir) = lib_dir(1:ilib_dir)
      icmd = icmd + ilib_dir
      write(6,*)unix_cmd(1:icmd)
      call system(unix_cmd(1:icmd))
   end if
   
   icmd = 9
   unix_cmd(1:icmd) = 'mkdir -p '
   icmd = icmd + 1
   unix_cmd(icmd:icmd + ilib_dir) = lib_dir(1:ilib_dir)
   icmd = icmd + ilib_dir
   write(6,*)unix_cmd(1:icmd)
   call system(unix_cmd(1:icmd))

   ilib_dir = ilib_dir + 1
   lib_dir(ilib_dir:ilib_dir) = '/'
   if(.not. pop_calc)then
      lib_dir(ilib_dir+1:ilib_dir+1) = particle(projectile%particle_type)%label
      ilib_dir = ilib_dir + 1
   else
      lib_dir(ilib_dir+1:ilib_dir+3) = 'Pop'
      ilib_dir = ilib_dir + 3
   end if

   icmd = 9
   unix_cmd(1:icmd) = 'mkdir -p '
   icmd = icmd + 1
   unix_cmd(icmd:icmd + ilib_dir) = lib_dir(1:ilib_dir)
   icmd = icmd + ilib_dir
   write(6,*)unix_cmd(1:icmd)
   call system(unix_cmd(1:icmd))

   if(.not. event_generator)then
      if(primary_decay)then
         if_check = ilib_dir
         file_check(1:ilib_dir) = lib_dir
         if_check = if_check + 1
         file_check(if_check:if_check) = '/'
         file_check(if_check+1:if_check+13) = 'Primary-decay'
         if_check = if_check + 13
         icmd = 9
         unix_cmd(1:icmd) = 'mkdir -p '
         unix_cmd(icmd+1:icmd + if_check) = file_check(1:if_check)
         icmd = icmd + if_check
         write(6,*)unix_cmd(1:icmd)
         call system(unix_cmd(1:icmd))
      end if

      if(fission)then
         if_check = ilib_dir
         file_check(1:ilib_dir) = lib_dir
         if_check = if_check + 1
         file_check(if_check:if_check) = '/'
         file_check(if_check+1:if_check+1) = 'f'
         if_check = if_check + 1
         icmd = 9
         unix_cmd(1:icmd) = 'mkdir -p '
         unix_cmd(icmd+1:icmd + if_check) = file_check(1:if_check)
         icmd = icmd + if_check
         write(6,*)unix_cmd(1:icmd)
         call system(unix_cmd(1:icmd))
      end if

      do i = 1, num_channels
         if_check = ilib_dir
         file_check(1:ilib_dir) = lib_dir
         if_check = if_check + 1
         file_check(if_check:if_check) = '/'
         ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
         file_check(if_check+1:if_check+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         if_check = if_check + ilast
         icmd = 9
         unix_cmd(1:icmd) = 'mkdir -p '
         icmd = icmd
         unix_cmd(icmd+1:icmd + if_check) = file_check(1:if_check)
         icmd = icmd + if_check
         write(6,*)unix_cmd(1:icmd)
         call system(unix_cmd(1:icmd))

         if(track_gammas)then
            unix_cmd(icmd+1:icmd+7) = '/gammas'
            icmd = icmd + 7
            write(6,*)unix_cmd(1:icmd)
            call system(unix_cmd(1:icmd))
         end if
      end do
   end if
   if_check = ilib_dir
   file_check(1:ilib_dir) = lib_dir
   if_check = if_check + 1
   file_check(if_check:if_check) = '/'
   file_check(if_check+1:if_check+7) = 'Spectra'
   if_check = if_check + 7
   icmd = 9
   unix_cmd(1:icmd) = 'mkdir -p '
   unix_cmd(icmd+1:icmd + if_check) = file_check(1:if_check)
   icmd = icmd + if_check
   write(6,*)unix_cmd(1:icmd)
   call system(unix_cmd(1:icmd))
!-----   Direct printing of event to their own directory


   if(dump_events)then
      if_check = ilib_dir
      file_check(1:ilib_dir) = lib_dir
      if_check = if_check + 1
      file_check(if_check:if_check) = '/'
      file_check(if_check+1:if_check+11) = 'Event-files'
      if_check = if_check + 11
      icmd = 9
      unix_cmd(1:icmd) = 'mkdir -p '
      unix_cmd(icmd+1:icmd + if_check) = file_check(1:if_check)
      icmd = icmd + if_check
      write(6,*)unix_cmd(1:icmd)
      call system(unix_cmd(1:icmd))
   end if
   return
end subroutine check_directories
