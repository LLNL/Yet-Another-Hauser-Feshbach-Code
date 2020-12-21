!
!*****************************************************************************80
!
subroutine range_char(line,iend)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This routine finds the last set character in a character(len=132)
!    string line
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
   use variable_kinds
   implicit none
   character(len=132), intent(in) :: line
   integer(kind=4), intent(out) :: iend
!-----------------------------------------------
   integer i
   iend = 0
   do i = 132, 1, -1
      if(line(i:i) /= ' ')then
         iend = i
         exit
      end if
   end do
   return
end subroutine range_char
