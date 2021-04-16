!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to set up E&M strength function parameters for
!    each nucleus
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
subroutine EM_str_param
   use variable_kinds
   use options
   use constants
   use nodeinfo
   use nuclei
   implicit none
!-------------------------------------------------------------------------
   integer(kind=4) :: i, k
   integer(kind=4) :: num_res
   integer(kind=4) :: l_radiation
   real(kind=8) :: xZ, xA
   real(kind=8) :: f, f1, ratio
!--------------------External functions--------------------------------
   real(kind=8) :: EL_f
   real(kind=8) :: ML_f
!-------------------------------------------------------------------------
   do i = 1, num_comp
      nucleus(i)%lmax_E = e_l_max                     !   maximum Electro-magnetic multipole 
!----------------------------------------------------------------------
      xA = real(nucleus(i)%A,kind=8)
      xZ = real(nucleus(i)%Z,kind=8)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------   Electric dipole strength function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      num_res = 0
      do k = 1, 3
         if(nucleus(i)%sr_E1(k) > 1.0d-6)num_res = num_res +1
      end do
      nucleus(i)%num_res = num_res

      l_radiation = 1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Magnetic dipole
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(m_l_max >= 1)then
         nucleus(i)%lmax_M = m_l_max
         allocate(nucleus(i)%er_M(1:nucleus(i)%lmax_M))
         allocate(nucleus(i)%gr_M(1:nucleus(i)%lmax_M))
         allocate(nucleus(i)%sr_M(1:nucleus(i)%lmax_M))
         l_radiation = 1
         nucleus(i)%er_M(1) = 41.0d0/xa**(1.0d0/3.0d0)
         nucleus(i)%gr_M(1) = 4.0d0

!-----------   Get M1 value for f at 7 MeV

         f1 = EL_f(i, 1, 7.0d0, nucleus(i)%sep_e(1) )

!----------    Now M1 value for f at 7 MeV
 
         nucleus(i)%sr_M(1) = 1.0d0
         f = ML_f(i, l_radiation, 7.0d0)

         ratio=f1/f

         nucleus(i)%sr_M(1)=ratio/(0.0588*xA**(0.878))

      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Higher electric multipoles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(e_l_max > 1)then
         allocate(nucleus(i)%er_E(2:e_l_max))
         allocate(nucleus(i)%gr_E(2:e_l_max))
         allocate(nucleus(i)%sr_E(2:e_l_max))
      end if
      if(e_l_max >= 2)then           !  Include electric quadrupole
         nucleus(i)%er_E(2) = 63.0d0/xA**(1.0d0/3.0d0)
         nucleus(i)%gr_E(2) = 6.11-0.012*xA
         nucleus(i)%sr_E(2) = 0.00014d0*xZ**2*nucleus(i)%er_E(2)/        &
                              (xA**(1.0d0/3.0d0)*nucleus(i)%gr_E(2))

      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------    Electric multipoles above E2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do l_radiation = 3, e_l_max
         nucleus(i)%er_E(l_radiation) = nucleus(i)%er_E(l_radiation-1)
         nucleus(i)%gr_E(l_radiation) = nucleus(i)%gr_E(l_radiation-1)
         nucleus(i)%sr_E(l_radiation) = 8.0d-4*nucleus(i)%sr_E(l_radiation-1)
      end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------    Magnetic multipoles above M1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do l_radiation = 2, m_l_max
         nucleus(i)%er_M(l_radiation) = nucleus(i)%er_M(l_radiation-1)
         nucleus(i)%gr_M(l_radiation) = nucleus(i)%gr_M(l_radiation-1)
         nucleus(i)%sr_M(l_radiation) = 8.0d-4*nucleus(i)%sr_M(l_radiation-1)
      end do
   end do

   return
end subroutine EM_str_param