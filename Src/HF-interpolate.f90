!
!*******************************************************************************
!
real(kind=8) function EM_interpolate(de,energy,l,ndim,str)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function interpolates the E(L) transmission coefficient along a grid
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
   implicit none
   real(kind=8), intent(in) :: de, energy
   integer(kind=4), intent(in) :: l
   integer(kind=4), intent(in) :: ndim
   real(kind=8), intent(in) :: str(0:ndim)
!---------------------------------------------------------------------------
   integer(kind=4) :: i
   real(kind=8) :: ym1,y0,yp1
   real(kind=8) :: det,a,b,c
   real(kind=8) :: xm1,x0,xp1,x
   real(kind=8) :: yinterpolate
   integer(kind=4) :: i1,i2,i3
   i=nint(energy/de)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   8/16/2011 added ndim == 1 to if statement to      +
!----   fix crash when there is only one element.         +
!----   Otherwise we get i1 = -1 and an array out of      +
!----   bounds.                                           +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   i1 = -1
   i2 = -1
   i3 = -1
   if(i == 0 .or. ndim == 1)then
      x0=de
      a=str(1)/x0**(2*l+2)
      EM_interpolate=a*energy**(2*l+2)
      return
   else if(i < ndim)then
      i1=i-1
      i2=i
      i3=i+1
   else if(i == ndim)then
      i1=i-2
      i2=i-1
      i3=i
   end if
   x0=dfloat(i2)*de
   xm1=x0-de
   xp1=x0+de
   ym1 = str(i1)
   y0  = str(i2)
   yp1 = str(i3)
   det=(x0*xp1**2-xp1*x0**2)-xm1*(xp1**2-x0**2)+xm1**2*(xp1-x0)
   a=(ym1*(x0*xp1**2-xp1*x0**2)-xm1*(y0*xp1**2-yp1*x0**2)+xm1**2*(y0*xp1-yp1*x0))/det
   b=((y0*xp1**2-yp1*x0**2)-ym1*(xp1**2-x0**2)+xm1**2*(yp1-y0))/det
   c=((x0*yp1-xp1*y0)-xm1*(yp1-y0)+ym1*(xp1-x0))/det
   x=energy
   yinterpolate=a+b*x+c*x**2
   EM_interpolate=yinterpolate
      if(EM_interpolate < 0.0)then
         write(6,*)'Warning EM_interpolation returned negative strength function'
         write(6,*)i,energy
         write(6,*)str(i-1),str(i),str(i+1)
         write(6,*)a,b,c
         write(6,*)xm1,x0,xp1,x
         write(6,*)str(i-1),a+b*xm1+c*xm1**2
         write(6,*)str(i),a+b*x0+c*x0**2
         write(6,*)str(i+1),a+b*xp1+c*xp1**2
         write(6,*)x,EM_interpolate
         stop
      end if
   return
end function EM_interpolate
!
!*******************************************************************************
!
real(kind=8) function spec_interpolate(de,a,b,c,energy,ndim,str)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function interpolates within the array str
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
   implicit none
   real(kind=8), intent(in) :: de, a, b, c, energy
   integer(kind=4), intent(in) :: ndim
   real(kind=8), intent(in) :: str(0:ndim)
!---------------------------------------------------------------------------
   real(kind=8) :: x,p
   integer(kind=4) :: i
!---------------------------------------------------------------------------
   if(energy < de)then
!------   Do a cubic spline back to the origin
      x=energy
      spec_interpolate=a*x+b*x**2+c*x**3
      return
   end if
   if(energy > dfloat(ndim)*de)then
      spec_interpolate=0.0d0
      return
   end if
   i = int(energy/de)
   if(i == 1)i=2
   p=(energy-dfloat(i)*de)/de
   spec_interpolate=p*(p-1.0d0)*str(i-1)/2.0d0+(1.0d0-p**2)*str(i)+ &
                   p*(p+1.0d0)*str(i+1)/2.0d0
   return
end function spec_interpolate
