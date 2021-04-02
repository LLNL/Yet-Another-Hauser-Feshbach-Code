!
!
!*******************************************************************************
!
real(kind=8) function tco_interpolate(e,nume,e_grid,tco)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine finds the transmission coeeficient at a 
!    specific energy e by interpolating from a list of nume values 
!    in the array tco on a energy grid defined by e_grid 
!    Interpolate using a log-log linear approximation
!    
!    Improved 6 Jan 2021 making use of the constant grid in log(e) to
!    find mid points. Sped up by 6x.
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
   use nodeinfo
   use variable_kinds
   implicit none
   real(kind=8), intent(in) :: e
   integer(kind=4), intent(in) :: nume
   real(kind=8), intent(in) :: e_grid(nume),tco(nume)
!-----------------------------------------------------------------------------
   integer(kind=4) :: i1, i2 
!   real(kind=8) :: tco_1,tco_2
   real(kind=8) :: x1, x2, y1, y2, x, y, a, b 

   real(kind=8) :: delta
!-----------------------------------------------------------------------------
   tco_interpolate = 1.0d-9
   if(e <= 0.0d0)return
   tco_interpolate = tco(1)
   if(e <= e_grid(1))return

   x = log(e)
   x1 = log(e_grid(1))

   delta = log(e_grid(2)) - x1
   i1 = int((x - x1)/delta) + 1
   i2 = i1 + 1

!   if(i1 >= nume)then
!      if(iproc == 0)write(6,*)'Error in log_tco_interpolate e > egrid(nume)'
!      call MPI_Abort(icomm,101,ierr)
!   end if

   x1 = log(e_grid(i1))
   x2 = log(e_grid(i2))
!------   Put in lower threshold
!   tco_1 = max(tco(i1),1.0d-9)
!   tco_2 = max(tco(i2),1.0d-9)
!   y1 = log(tco_1)  
!   y2 = log(tco_2)
   y1 = log(max(tco(i1),1.0d-9))
   y2 = log(max(tco(i2),1.0d-9))
   a = (y2-y1)/(x2-x1)
   b = y1 - a*x1
   y = a*x + b
   tco_interpolate = exp(y)

!  write(6,*)'LOG'
!  write(6,*)i1,i2
!  write(6,*)x,x1,x2,y1,y2
!  write(6,*)a,b,y

  return
end function tco_interpolate

