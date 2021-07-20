!
!*****************************************************************************80
!
subroutine EM_str_param(inuc)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to set up E&M strength function parameters for
!    nucleus inuc
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
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: EL_f
!        real(kind=8) :: ML_f
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
   implicit none
!-------------------------------------------------------------------------
   integer(kind=4), intent(in) :: inuc
!-------------------------------------------------------------------------
   integer(kind=4) :: l_radiation
   real(kind=8) :: xZ, xA
   real(kind=8) :: f, f1, ratio
!--------------------External functions--------------------------------
   real(kind=8) :: EL_f
   real(kind=8) :: ML_f
!-------------------------------------------------------------------------
   nucleus(inuc)%lmax_E = e_l_max                     !   maximum Electro-magnetic multipole 
!----------------------------------------------------------------------
   xA = real(nucleus(inuc)%A,kind=8)
   xZ = real(nucleus(inuc)%Z,kind=8)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Magnetic dipole
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   l_radiation = 1
   if(m_l_max >= 1 .and. nucleus(inuc)%ML_mode(l_radiation)%default)then

      nucleus(inuc)%lmax_M = m_l_max

      nucleus(inuc)%ML_mode(l_radiation)%num_gsf = 1
      nucleus(inuc)%ML_mode(l_radiation)%gsf(1)%er = 41.0d0/xa**(1.0d0/3.0d0)
      nucleus(inuc)%ML_mode(l_radiation)%gsf(1)%gr = 4.0d0
      f1 = EL_f(inuc, 1, 7.0d0, nucleus(inuc)%sep_e(1) )
      nucleus(inuc)%ML_mode(l_radiation)%gsf(1)%sr = 1.0d0
      f = ML_f(inuc, l_radiation, 7.0d0)

      ratio = f1/f

      nucleus(inuc)%ML_mode(l_radiation)%gsf(1)%sr = ratio/(0.0588*xA**(0.878))
   end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Higher electric multipoles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(e_l_max >= 2)then           !  Include electric quadrupole
      if(nucleus(inuc)%EL_mode(2)%default)then
         nucleus(inuc)%EL_mode(2)%num_gsf = 1
         nucleus(inuc)%EL_mode(2)%gsf(1)%er = 63.0d0/xA**(1.0d0/3.0d0)
         nucleus(inuc)%EL_mode(2)%gsf(1)%gr = 6.11d0 - 0.012d0*xA
         nucleus(inuc)%EL_mode(2)%gsf(1)%sr =                                                &
               0.00014d0*xZ**2*nucleus(inuc)%EL_mode(2)%gsf(1)%er/                           &
                 (xA**(1.0d0/3.0d0)*nucleus(inuc)%EL_mode(2)%gsf(1)%gr) 
      end if
   end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------    Electric multipoles above E2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do l_radiation = 3, e_l_max
      if(nucleus(inuc)%EL_mode(l_radiation)%default)then
         nucleus(inuc)%EL_mode(l_radiation)%num_gsf = 1
         nucleus(inuc)%EL_mode(l_radiation)%gsf(1)%er =                                      &
              nucleus(inuc)%EL_mode(l_radiation-1)%gsf(1)%er
         nucleus(inuc)%EL_mode(l_radiation)%gsf(1)%gr =                                      &
              nucleus(inuc)%EL_mode(l_radiation-1)%gsf(1)%gr
         nucleus(inuc)%EL_mode(l_radiation)%gsf(1)%sr =                                      &
              8.0d-4*nucleus(inuc)%EL_mode(l_radiation-1)%gsf(1)%sr
      end if
   end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------    Magnetic multipoles above M1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do l_radiation = 2, m_l_max
      if(nucleus(inuc)%ML_mode(l_radiation)%default)then
         nucleus(inuc)%ML_mode(l_radiation)%num_gsf = 1
         nucleus(inuc)%ML_mode(l_radiation)%gsf(1)%er =                                      &
              nucleus(inuc)%ML_mode(l_radiation-1)%gsf(1)%er
         nucleus(inuc)%ML_mode(l_radiation)%gsf(1)%gr =                                      &
              nucleus(inuc)%ML_mode(l_radiation-1)%gsf(1)%gr
         nucleus(inuc)%ML_mode(l_radiation)%gsf(1)%sr =                                      &
              8.0d-4*nucleus(inuc)%ML_mode(l_radiation-1)%gsf(1)%sr
      end if
   end do

   return
end subroutine EM_str_param
