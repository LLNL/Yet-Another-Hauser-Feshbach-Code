!
!*******************************************************************************
!
!  Discussion:
!
!  FORTRAN 90 function to return 64-bit pseudo-random numbers with a large 
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
real(kind=8) function random_64(iseed_64)
   implicit none
   integer, parameter :: int64 = selected_int_kind(18)
   integer(kind=int64), intent(inout) :: iseed_64
   integer(kind=int64), save :: A, M, Q, R
   integer(kind=int64), save :: kk
   real(kind=8), save :: AM
   integer(kind=int64), save :: i1 = -1_int64, i2 = -1_int64
!----------------------------------------------------------
   if(iseed_64 <= 0_int64 .or. i2 < 0_int64)then              !   Initialize the seed
      M = 9223372036854775097_int64
      A = 2946716567_int64
      Q = 3130050626_int64
      R = 1671854155_int64
      AM = nearest(1.0d0,-1.0d0)/M
      i2 = ior(ieor(88888888899999999_int64,abs(iseed_64)),1_int64)
      i1 = ieor(77777777555555555_int64, abs(iseed_64))
      iseed_64 = abs(iseed_64) + 1
   end if

   i1 = ieor(i1,ishft(i1,13))                  !  Apply Marasglia shift sequence
   i1 = ieor(i1,ishft(i1,-7))
   i1 = ieor(i1,ishft(i1,17))
   kk = i2/Q                                   !  Compute mod using Schrage's algorithm
   i2 = A*(i2 - kk*Q) - R*kk
   if(i2 < 0)i2 = i2 + M

   random_64 = AM*ior(iand(M,ieor(i1,i2)),1_int64)      !  Combine generators
   return
end function random_64
