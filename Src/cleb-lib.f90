subroutine factorial(lfactc,faclog)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine sets up a stored array of the log of n factorial
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
   implicit none
   integer(kind=4), intent(in) :: lfactc
   real(kind=8), intent(out) :: faclog(lfactc)
   integer(kind=4) :: i
   real(kind=8) :: fn
   faclog(1) = 0.0d0
   faclog(2) = 0.0d0
   fn = 1.0d0
   do i = 3, lfactc
      fn = fn + 1.0d0
      faclog(i) = faclog(i-1) + log(fn)
   end do
   return
end subroutine factorial
!
!*******************************************************************************
!
real(kind=8) function wignerj(j1,j2,j3,m1,m2,m3)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the Wigner 3J symbol
!    Modified from a FORTRAN77 code of unknown origin
!    Converted and updtated to FORTRAN 90 to remove outdated
!    FORTRAN 
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
!      Erich Ormand, LLNL updated to FORTRAN 90 removing computed IF and GOTO
!      statements
!
!*******************************************************************************
!
   implicit none
   real(kind=8), intent(in) :: j1, m1, j2, m2, j3, m3
!-------------------------------------------------------------------
   integer(kind=4), parameter :: num_fac = 200
!
   real(kind=8) :: fact(num_fac)
   real(kind=8) :: x, y, z

   integer(kind=4) :: ia, ib, ic, id, ie
   integer(kind=4) :: i
   integer(kind=4) :: ih, ii, ij, ik, il, im, in
   integer(kind=4) :: ix, iy, iz
   integer(kind=4) :: minn, nin
   real(kind=8) :: s, t, ta, tb, xddd

   real(kind=8) :: phase 
   logical :: frac, triangle
!----  Phase (-1)**i
   phase(i) = (-1.0d0)**i
!----  Check x is an integer
   frac(x) = abs(x-int(x)) > 1.0d-10
!----   Triangle relation to check if angular momenta can couple
   triangle(x,y,z) = frac(x+y+z) .or. x > y+z .or. x < abs(y-z)
!--------------------------------------------------------------
!
   call factorial(num_fac,fact)

   wignerj = 0.0d0
   if((m1 + m2 + m3) /= 0.0d0)return    !   check on z-component 
   if(triangle(j3,j1,j2))return            !   angular momenta can't couple
   if(j1 - abs(m1) < 0)return           !   m1 > j1 - not allowed
   if(j2 - abs(m2) < 0)return           !   m2 > j1 - not allowed
   if(j3 - abs(m3) < 0)return           !   m3 > j - not allowed

   ia = int(j3 - j2 + m1)

   ib = int(j3 - j1 - m2)

   if(ia < 0)then
      minn = -ia
      if(minn + ib < 0)minn  = -ib
   else
      if(ib < 0)then
         minn = -ia
         if(minn + ib < 0)minn  = -ib
      else
         minn = 0
      end if
   end if

   ic = int(j1 - m1)
   id = int(j2 + m2)
   ie = int(j1 + j2 - j3)
   nin = minn
   t = real(phase(minn),kind=8)
   s = t

   ix = 1
   do while(ix > 0)
      minn = minn + 1
      i = int(j2 - j1 + m3)
      iz = ic - minn + 1
      if(iz > 0)then
         iy = id - minn + 1
         if(iy > 0)then
            ix = ie - minn + 1
            if(ix > 0)then
               ta = real(ix*iy*iz,kind=8)
               tb = real(minn*(ia+minn)*(ib+minn),kind=8)
               t = -t*ta/tb
               s = s + t
            end if
         else
            ix = -1
         end if
      else
         ix = -1
      end if
   end do   


   if(abs(s) <= 1.0d-10)then
      wignerj = 0.0d0
      return
   end if

   ih = int(j1 + m1)
   ii = int(j2 - m2)
   ij = int(j3 + m3)
   ik = int(j3 - m3)
   il = int(j1 + j3 - j2)
   im = int(j2 + j3 - j1)
   in = int(j1 + j2 + j3) + 1
   xddd = 0.5d0*(fact(ih+1) + fact(ic+1) + fact(id+1) +                  &
                 fact(ii+1) + fact(ij+1) + fact(ik+1) +                  &
                 fact(ie+1) + fact(im+1) + fact(il+1) - fact(in+1)) -    &
                 (fact(ic-nin+1) + fact(ia+nin+1) + fact(id-nin+1) +     &
                 fact(ib+nin+1) + fact(nin+1) + fact(ie-nin+1))
   wignerj = phase(i)*exp(xddd)*s

   return
end function wignerj
!
!*******************************************************************************
!
real(kind=8) function clebr(j1,m1,j2,m2,j3,m3)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes Clebsch-Gordan coefficient
!    Modified from a FORTRAN77 code of unknown origin
!    Converted and updtated to FORTRAN 90 to remove outdated
!    FORTRAN 
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
!      Erich Ormand, LLNL updated to FORTRAN 90 removing computed IF and GOTO
!      statements
!
!*******************************************************************************
!
   implicit none
   real(kind=8), intent(in) :: j1, j2, j3
   real(kind=8), intent(in) :: m1, m2, m3
!---------------------------------------------------------------------------
   integer(kind=4) :: i
   integer(kind=4) :: iphase
   real(kind=8) :: wignerj
   real(kind=8) :: phase
!----  Phase (-1)**i
   phase(i) = (-1.0d0)**i
!--------------------------------------------------------------

   iphase = nint(-j1 + j2 - m3)

   clebr = phase(iphase)*wignerj(j1,j2,j3,m1,m2,-m3)*sqrt(2.0d0*j3+1.0d0)

   return
end function clebr
!
!*******************************************************************************
!
real(kind=8) function wig6j(j1,j2,j3,j4,j5,j6)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes Wigner 6J symbols defined in Brink and Satcher
!    Modified from a FORTRAN77 code of unknown origin
!    Converted and updtated to FORTRAN 90 to remove outdated
!    FORTRAN 
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
   implicit none
   real(kind=8), intent(in) :: j1, j2, j3, j4, j5, j6
!---------------------------------------------------------------------------
   integer(kind=4) :: i
   integer(kind=4) :: ia, ib, ic, id, ie, ig, ih, m, n, mup
   integer(kind=4) :: it, iu, iv, iw
   real(kind=8) :: xd
   real(kind=8) :: ta, tb, t, s

   real(kind=8) :: x, y, z

   integer(kind=4) :: num_fac
   parameter (num_fac = 200)
   real(kind=8) :: fact(num_fac)
!--------------------------------------------------------------
   real(kind=8) :: phase 
   logical :: frac, triangle
!----  Phase (-1)**i
   phase(i) = (-1.0d0)**i
!----  Check x is an integer
   frac(x) = abs(x-int(x)) > 1.0d-10
!----   Triangle relation to check if angular momenta can couple
   triangle(x,y,z) = frac(x+y+z) .or. x > y+z .or. x < abs(y-z)
!--------------------------------------------------------------
!
   call factorial(num_fac,fact)
!
   wig6j = 0.0d0

   if(triangle(j1,j2,j3))return
   if(triangle(j1,j5,j6))return
   if(triangle(j2,j4,j6))return
   if(triangle(j3,j4,j5))return

   ic = int(j1+j2-j3)
   id = int(j5+j4-j3)
   ie = int(j1+j5-j6)
   ig = int(j2+j4-j6)
   ia = int(j3+j6-j1-j4)
   ib = int(j3+j6-j2-j5)
   ih = int(j1+j2+j5+j4) + 1

   m = min(ih,ic,id,ie,ig)

   if(m < 0)return

   mup = min(ia,ib,0)
   t = phase(m)
   n = m
   s = t

   m = m - 1
   do while((m+mup) >= 0)
      ta = (ia+m+1)*(ib+m+1)*(ih-m)*(m+1)
      tb = (ic-m)*(id-m)*(ie-m)*(ig-m)
      t = -t*ta/tb
      s = s + t
      m = m - 1
   end do

   it = int(j1+j2+j3) + 1
   iu = int(j1+j5+j6) + 1
   iv = int(j2+j4+j6) + 1
   iw = int(j3+j4+j5) + 1
   xd = 0.5d0*(fact(1+ic) + fact(1+ie+ib) + fact(1+ia+ig) +        &
               fact(1+ie) + fact(1+ib+ic) + fact(1+ia+id) +        &
               fact(1+ig) + fact(1+ic+ia) + fact(1+id+ib) +        &
               fact(1+id) + fact(1+ia+ie) + fact(1+ib+ig) -        &
               fact(1+it) - fact(1+iu) - fact(1+iv)-fact(1+iw)) +  &
        fact(1+ih-n) - fact(1+n) - fact(1+ia+n) - fact(1+ib+n) -   &
        fact(1+ic-n) - fact(1+id-n)-fact(1+ie-n)-fact(1+ig-n)
   wig6j = phase(ih-1)*s*exp(xd)

   return
end function wig6j
!
!*******************************************************************************
!
real(kind=8) function wig9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes Wigner 9J symbols defined in Brink and Satcher
!    Modified from a FORTRAN77 code of unknown origin
!    Converted and updtated to FORTRAN 90 to remove outdated
!    FORTRAN 
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
   implicit none
   real(kind=8), intent(in) :: j1, j2, j3, j4, j5, j6, j7, j8, j9
!---------------------------------------------------------------------------
   real(kind=8) :: s
   real(kind=8) :: xa, xb, xc, x
   integer(kind=4) :: k
   integer(kind=4) :: i
   real(kind=8) :: xx, yy, zz
   real(kind=8) :: wig6j
!------------------------------------------------------------------------
   real(kind=8) :: phase 
   logical :: frac, triangle
!----  Phase (-1)**i
   phase(i) = (-1.0d0)**i
!----  Check x is an integer
   frac(x) = abs(xx-int(xx)) > 1.0d-10
!----   Triangle relation to check if angular momenta can couple
   triangle(xx,yy,zz) = frac(xx+yy+zz) .or. xx > yy+zz .or. xx < abs(yy-zz)
!--------------------------------------------------------------
!
   wig9j = 0.0d0

   if(triangle(j1,j2,j3))return
   if(triangle(j4,j5,j6))return
   if(triangle(j3,j6,j9))return
   if(triangle(j1,j4,j7))return
   if(triangle(j2,j5,j8))return
   if(triangle(j7,j8,j9))return
   s = 0.0d0

   xa = abs(j1-j9)
   xb = abs(j4-j8)
   xc = abs(j2-j6)
   x = xa

   if((x - xb) < 0)x = xb
   if((x - xc) < 0)x = xc

   do while((x-j1-j9) <= 0.0d0 .and. (x-j4-j8) <= 0.0d0 .and.       &
            (x-j2-j6) <= 0.0d0 )
      s = s + (2.0d0*x+1.0d0)*wig6j(j1,j9,X,j8,j4,j7)*              &
          wig6j(j2,j6,X,j4,j8,j5)*wig6j(j1,j9,X,j6,j2,j3)
      x = x + 1.0d0
   end do

   k = int(2.0d0*(j1+j2+j4+j6+j8+j9))
   wig9j = phase(k)*s
   return
end function wig9j
!
!*******************************************************************************
!
real(kind=8) function racahr(j1,j2,j3,j4,j5,j6)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes Racah 9J symbols defined in Brink and Satcher
!    Modified from a FORTRAN77 code of unknown origin
!    Converted and updtated to FORTRAN 90 to remove outdated
!    FORTRAN 
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
   implicit none
   real(kind=8), intent(in) :: j1, j2, j3, j4, j5, j6
!---------------------------------------------------------------------------
   integer(kind=4) :: i
   real(kind=8) :: z
   real(kind=8) :: phase
   real(kind=8) :: wig6j
!----   Compute phase (-1)**i
   phase(i) = (-1.0d0)**i
!--------------------------------------------------------------
!
   z = abs(j1+j2+j3+j4)
   i = int(z + 0.5d0)
   racahr = phase(i)*wig6j(j1,j2,j5,j4,j3,j6)
   return
end function racahr


