!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine prints populations following the first, primary decay step
!    to the appropriate file in the library directory 
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
subroutine print_primary_decay(ilib_dir, lib_dir, ifile, file_name,        &
                               e_in, reaction_cs)
   use variable_kinds
   use options
   use constants
   use nodeinfo
   use nuclei
   use Channel_info
   use particles_def
   implicit none  
   integer(kind=4), intent(in):: ilib_dir
   character(len=132), intent(in) :: lib_dir
   integer(kind=4), intent(in):: ifile
   character(len=40), intent(in) :: file_name
   real(kind=8), intent(in) :: e_in
   real(kind=8), intent(in) :: reaction_cs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer(kind=4) :: jj, k, ip
   integer(kind=4) :: num_states
   integer(kind=4) :: ilast
   integer(kind=4) :: if1, inuc
   integer(kind=4) :: idir
   real(kind=8) :: energy
!-check-norm   real(kind=8) :: sum
   character(len=132) :: directory
   character(len=132) :: temp_string, fstring
   character(len=132) :: file_char

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   idir = idir + 1
   directory(idir:idir+12) = 'Primary-decay'
   idir = idir + 13
   directory(idir:idir) = '/'


!-check-norm   sum = 0.0d0
   do if1 = 1, nucleus(1)%num_decay
      inuc = nucleus(1)%decay_to(if1)
      ilast = index(nucleus(inuc)%Label,' ')-1
      if(ilast < 0)ilast = 5

      file_char(1:132) = ' '
!      file_char = directory(1:idir)//file_name(1:ifile)//   &
!                  '-'//nucleus(inuc)%Label(1:ilast)//'decay.dat'
      file_char = directory(1:idir)//file_name(1:ifile)//   &
                  "-"//trim(adjustl(nucleus(inuc)%Label))//".dat"
      ilast = index(file_char,' ') - 1
      if(ilast < 1)ilast = 132

!  write(6,*)ilast
!  write(6,*)file_char
!  write(6,*)file_char(1:ilast)

      open(unit = 24, file = trim(adjustl(file_char)), status = 'unknown')
!      open(unit=24,file=                                                                     &
!          directory(1:idir)//file_name(1:ifile)//'_'//nucleus(inuc)%Label(1:ilast)//'.dat',  &
!          status='unknown')
      write(24,'(''#Cross section from the primary decay '')')
      write(24,'(''#Population for ead discrete state and continuous energy bin '')')
      write(24,'(''#Incident energy = '',f10.4)')e_in
      num_states = nucleus(inuc)%ncut
      if(all_discrete_states)num_states = nucleus(inuc)%num_discrete
      write(24,'(''#'')')
      write(24,'(''#Discrete states '')')
      do k = 1, num_states
         write(24,'(''    #       Ex     J   par     cs('',a2,'')'')')cs_units
         write(24,'(''  ---     ------  --- ----  ---------------'')')
         write(24,'(1x,i4,1x,f10.4,2(1x,f4.1),1x,1pe15.7)')                                  &
                     k-1,nucleus(inuc)%state(k)%energy,                                      &
                     nucleus(inuc)%state(k)%spin, nucleus(inuc)%state(k)%parity,             &
                     nucleus(inuc)%state(k)%pop*reaction_cs*cs_scale
!-check-norm         sum = sum + nucleus(inuc)%state(k)%pop
      end do
      ip = 0
      write(24,'(''#'')')
      write(24,'(''#Positive-Parity Continuous Energy Bins '')')
      write(temp_string,*)min(nucleus(inuc)%j_max,60)+1
      fstring = "('#  Energy  ',"//trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"
      write(24,fstring)(real(jj) + nucleus(inuc)%jshift, jj = 0, min(nucleus(inuc)%j_max,60))
      fstring = "('#----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
      write(24,fstring)
      fstring = "(1x,f10.4,"//trim(adjustl(temp_string))//"(1x,1pe15.7))"
      do k = 1, nucleus(inuc)%nbin
         energy = nucleus(inuc)%e_grid(k)
         write(24,fstring)                                                                   &
            energy,(nucleus(inuc)%bins(jj,ip,k)%pop*reaction_cs*cs_scale,jj = 0, min(nucleus(inuc)%j_max,60))
!-check-norm         do jj = 0, min(nucleus(inuc)%j_max,60)
!-check-norm            sum = sum + nucleus(inuc)%bins(jj,ip,k)%pop
!-check-norm         end do
      end do
      ip = 1
      write(24,'(''#'')')
      write(24,'(''#Negative-Parity Continuous Energy Bins '')')
      fstring = "('#  Energy  ',"//trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"
      write(24,fstring)(real(jj) + nucleus(inuc)%jshift, jj = 0, min(nucleus(inuc)%j_max,60))
      fstring= "('#----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
      write(24,fstring)
      fstring= "(1x,f10.4,"//trim(adjustl(temp_string))//"(1x,1pe15.7))"
      do k = 1, nucleus(inuc)%nbin
         energy = nucleus(inuc)%e_grid(k)
         write(24,fstring)                                                                   &
            energy,(nucleus(inuc)%bins(jj,ip,k)%pop*reaction_cs*cs_scale,jj = 0, min(nucleus(inuc)%j_max,60))
!-check-norm         do jj = 0, min(nucleus(inuc)%j_max,60)
!-check-norm            sum = sum + nucleus(inuc)%bins(jj,ip,k)%pop
!-check-norm         end do
      end do
      close(unit=24)
   end do
!-check-norm   write(6,*)'sum = ',sum, reaction_cs

   return
end subroutine print_primary_decay
