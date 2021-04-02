
!
!*******************************************************************************
!
!  Discussion:
!
!  FORTRAN 90 function to return 32-bit pseudo-random numbers with a large 
!  period
!  Based on a modification of the linear congruential method
!
!      r(i+1) = mod(A*ran(i),M)
!
!  A = 2946716567, M = 2^63-1
!
!  where the mod operation is computing using Schrage's algorithm based
!  on an approximate factorization
!
!  M = A*Q + R
!
!  Q = [M/A} = 3130050626 and R = mod(M,A) = 1671854155
!
!  Then 
!
!  mod(A*Y,M) = A*mod(Y,Q) - R*[Y/Q]      if >= 0
!               A*mod(Y,Q) - R*[Y/Q] + M  if < 0
!
!  The generator is augmented with Marsaglia shift operartors for a 
!  64-bit random number generator
!
!  Date:
!
!    20 March 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
real(kind=8) function random_32(iseed_32)
   implicit none
   integer, parameter :: int32 = selected_int_kind(9)
   integer(kind=int32), intent(inout) :: iseed_32
   integer(kind=int32), save :: A, M, Q, R
   integer(kind=int32), save :: kk
   real(kind=8), save :: AM
   integer(kind=int32), save :: i1 = -1_int32, i2 = -1_int32
!----------------------------------------------------------
   if(iseed_32 <= 0_int32 .or. i2 < 0_int32)then              !   Initialize the seed
      M = 2147483647_int32
      A = 16807_int32
      Q = 127773_int32
      R = 2836_int32
      AM = nearest(1.0d0,-1.0d0)/M
      i2 = ior(ieor(888889999_int32,abs(iseed_32)),1_int32)
      i1 = ieor(777755555_int32, abs(iseed_32))
      iseed_32 = abs(iseed_32) + 1
   end if

   i1 = ieor(i1,ishft(i1,13))                  !  Apply Marasglia shift sequence
   i1 = ieor(i1,ishft(i1,-17))
   i1 = ieor(i1,ishft(i1,5))
   kk = i2/Q                                   !  Compute mod using Schrage's algorithm
   i2 = A*(i2 - kk*Q) - R*kk
   if(i2 < 0)i2 = i2 + M

   random_32 = AM*ior(iand(M,ieor(i1,i2)),1_int32)      !  Combine generators
   return
end function random_32
