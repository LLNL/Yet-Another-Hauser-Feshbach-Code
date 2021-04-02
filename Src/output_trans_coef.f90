!
!*******************************************************************************
!
subroutine output_trans_coef
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine prints out transmission coefficients on the calculation
!    grid
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
!*******************************************************************************
!
   use variable_kinds
   use nodeinfo
   use options
   use print_control
   use useful_data
   use nuclei
   use Channel_info
   use particles_def
   use directory_structure
   use constants
   implicit none
!-------------------------------------------------------------
   integer(kind=4) :: j, k, l
   integer(kind=4) :: i, isp
   real(kind=8) :: energy
   real(kind=8) :: xj, xi
   character(len=132) :: temp_string
   character(len=132) :: fstring
!------------------------------------------------------------------
   write(13,'(''Transmission coefficients for '','//               &
              '''each emitted particle (COM energy)'')')
   do k = 1, 6
      if(.not.particle(k)%in_decay)cycle

      write(temp_string,*)particle(k)%lmax+1

      isp = nint(2.0d0*particle(k)%spin)
      xj = -particle(k)%spin
      do i = 0, isp
         xi = real(i,kind=8)
         write(13,*)
         write(13,'('' Transmission coefficients for '',a8)')        &
               particle(k)%name
         if(xj + xi < 0)then
            write(13,'(''For j = l -'',f4.1)')abs(xj+xi)
         else
            write(13,'(''For j = l +'',f4.1)')xj+xi
         end if
         fstring = "('       T         l=',i3,5x,"                   &
                   //trim(adjustl(temp_string))//"(8x,i3,5x))"
         write(13,fstring)(l,l = 0, particle(k)%lmax)
         fstring = "(1x,'----------',"//trim(adjustl(temp_string))// &
                   "(1x,'---------------'))"
         write(13,fstring)
         fstring = "(1x,f10.5,"//trim(adjustl(temp_string))//"(1x,e15.7))"
         do j = 1, particle(k)%nbin
            energy = real(j,kind=8)*de
            write(13,fstring)energy,                                 &
                (particle(k)%trans(i,l,j),l = 0, particle(k)%lmax)
      end do
      end do
   end do

   flush(13)

   return
end subroutine output_trans_coef
