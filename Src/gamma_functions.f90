!
!*******************************************************************************
!
subroutine gdr_param(data_path, len_path, inuc)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine reads in and sets up data for E1 transitions
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
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
!*******************************************************************************
!
   use variable_kinds
   use options
   use nuclei
   implicit none
   character(len=200), intent(in) :: data_path      ! path where data files are kept
   integer(kind=4), intent(in) :: len_path             ! length of data_path
   integer(kind=4), intent(in) :: inuc
!--------------------------------------------------------------------------
   integer(kind=4) :: iz, ia
   integer(kind=4) :: i
   integer(kind=4) :: eof
   logical :: found
   integer(kind=4) :: izz,iaa
   real(kind=8) :: xA, xZ
   character(len=2) :: symbb
   real(kind=8) :: er1, gr1, sr1, er2, gr2, sr2
   integer(kind=4) :: L
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   iz = nucleus(inuc)%Z
   ia = nucleus(inuc)%A

   L = 1

   open(unit=52,                                                          &
        file=data_path(1:len_path)//'gdr-parameters-exp.dat',             &
        status='old')
   do i = 1, 4
      read(52,*)
   end do
   found = .false.
   eof = 0
   do while(eof == 0)
      read(52,'(2(1x,i3),1x,a2,6x,6(1x,f6.2))',iostat=eof)                &
           izz,iaa,symbb,er1,sr1,gr1,er2,sr2,gr2
      if(iz == izz .and. ia == iaa)then
         nucleus(inuc)%EL_mode(L)%num_gsf = 1
         nucleus(inuc)%EL_mode(L)%gsf(1)%er = er1
         nucleus(inuc)%EL_mode(L)%gsf(1)%gr = gr1
         nucleus(inuc)%EL_mode(L)%gsf(1)%sr = sr1
         if(sr2 > 1.0d-6)then
            nucleus(inuc)%EL_mode(L)%num_gsf = 2
            nucleus(inuc)%EL_mode(L)%gsf(2)%er = er2
            nucleus(inuc)%EL_mode(L)%gsf(2)%gr = gr2
            nucleus(inuc)%EL_mode(L)%gsf(2)%sr = sr2
         end if
         found = .true.
         exit
      end if
      if(iz == izz .and. iaa == 0)then     !   use natural GDR parameters
         nucleus(inuc)%EL_mode(L)%num_gsf = 1
         nucleus(inuc)%EL_mode(L)%gsf(1)%er = er1
         nucleus(inuc)%EL_mode(L)%gsf(1)%gr = gr1
         nucleus(inuc)%EL_mode(L)%gsf(1)%sr = sr1
         if(sr2 > 1.0d-6)then
            nucleus(inuc)%EL_mode(L)%num_gsf = 2
            nucleus(inuc)%EL_mode(L)%gsf(2)%er = er2
            nucleus(inuc)%EL_mode(L)%gsf(2)%gr = gr2
            nucleus(inuc)%EL_mode(L)%gsf(2)%sr = sr2
         end if
         found = .true.
         exit
      end if
   end do

   close(unit=52)

   if(found)return

   nucleus(inuc)%EL_mode(L)%num_gsf = 1
   xA = real(iA,kind=8)
   xZ = real(iZ,kind=8)
   nucleus(inuc)%EL_mode(L)%gsf(1)%er = 31.2d0*xA**(-1.d0/3.d0) + 20.6d0*xA**(-1.d0/6.d0)
   nucleus(inuc)%EL_mode(L)%gsf(1)%gr = 0.026d0*nucleus(inuc)%EL_mode(L)%gsf(1)%er**1.91d0
   nucleus(inuc)%EL_mode(L)%gsf(1)%sr = 144.0d0*real((ia-iz)*iz,kind=8)/                   &
                                        (xA*3.141593*nucleus(inuc)%EL_mode(L)%gsf(1)%gr)
   return
end subroutine gdr_param
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function E1_f_mode(e_gam, T, e1_model, er, gr, sr)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E1 strength function, f, as defiend by e1_model
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) gamma_t
!        real(kind=8) M_Lorentzian
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
!*******************************************************************************
!
   use variable_kinds
   use constants
   implicit none
   real(kind=8), intent(in) :: e_gam, T
   integer(kind=4), intent(in) :: e1_model
   real(kind=8), intent(in) :: er, gr, sr
   real(kind=8) gam, str
!-----   External functions   ---------------------------------------
   real(kind=8) gamma_t
   real(kind=8) M_Lorentzian
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   E1_f_mode = 0.0d0
   if(er < 1.0d-5 .or. gr < 1.0d-5)return
   if(e1_model > 1)then
      gam = gamma_t(e_gam,er,gr,T)
      if(e1_model == 3)then
         if(gam > 1.5*gr) gam = 1.5*gr
         if(gam < 0.25*gr) gam = 0.25*gr
      end if
   else
      gam = gr
   end if

   str = M_Lorentzian(e_gam,er,gam)
   str = e_gam*str
   if(e1_model > 1)str = str + 0.7d0*gr*(2.0d0*pi*T)**2/er**5      !  from RIPL handbook, note er**5
   E1_f_mode = sr*gr*str
   return
end function E1_f_mode
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function gamma_t(e_gam, er, gr, T)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E1 width as a function of energy and 
!    the parameter T
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
!
!     MPI routines:
!
!        None
!
!
!   Dependencies:
!
!     None
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
!*******************************************************************************
!
   use variable_kinds
   use constants
   implicit none
   real(kind=8), intent(in) :: e_gam,er,gr
   real(kind=8), intent(in) :: T
!----------------------------------------------------------------------
   real(kind=8) T_fac
   T_fac = two_pi*T
   gamma_t = gr*(e_gam**2+T_fac**2)/er**2
   return
end function gamma_t
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function Lorentzian(e_g,er,gr)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes a Lorentzian function
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
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
!*******************************************************************************
!
   use variable_kinds
   implicit none
   real(kind=8), intent(in) :: e_g,er,gr
!----------------------------------------------------------------------
   Lorentzian = e_g*gr**2/((e_g**2-er**2)**2+(e_g*gr)**2)
   return
end function Lorentzian
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function M_Lorentzian(e_g,er,gr)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes a Lorentzian function divided by the 
!    gamma energy eg
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
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
!*******************************************************************************
!
   use variable_kinds
   implicit none
   real(kind=8), intent(in) :: e_g,er,gr
!--------------------------------------------------------------------------
   M_Lorentzian = gr/((e_g**2-er**2)**2+(e_g*gr)**2)
   return
end function M_Lorentzian
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function Txl(egamma, L)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the multiplicative factor to convert the 
!    strength function to the transmission coefficient
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
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
!*******************************************************************************
!
   use variable_kinds
   use constants
   implicit none
   real(kind=8),intent(in) :: egamma
   integer(kind=4),intent(in) :: L
!------------------------------------------------------------------------------
   Txl = two_pi*egamma**(2*L+1)
   return
end function Txl
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function Kxl(L)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes k factor needed for E&M radiation
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
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
!*******************************************************************************
!
   use variable_kinds
   use constants
   implicit none
   integer(kind=4),intent(in) :: L
!------------------------------------------------------------------------------
   Kxl = real(2*L+1,kind=8)*(pi*hbar_c)**2*fmsq_eq_mbarn
   Kxl = 1.0d0/Kxl
   return
end function Kxl
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function temperature(U,ap)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes effective temperature with given excitation energy U
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
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
!*******************************************************************************
!
   use variable_kinds
   implicit none
   real(kind=8),intent(in) :: U, ap
!------------------------------------------------------------------------------
   temperature=0.0
   if(U > 0.0)temperature = sqrt(U/ap)
   return
end function temperature
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!real(kind=8) function Lfac(egam,T)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes effective temperature with given excitation energy U
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
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
!*******************************************************************************
!
!   use variable_kinds
!   implicit none
!   real(kind=8) egam, T
!!------------------------------------------------------------------------------
!   if(egam < 1.0d-6)then
!      Lfac = 0.0d0
!   else
!      Lfac = 1.0d0
!      if(T > 0.0)Lfac = 1.0/(1.0-exp(-egam/T))
!   end if
!   return
!end function Lfac
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function E1_f(i_f, e_gamma, energy)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E1 strength function, f
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: Ux
!        real(kind=8) :: aparam_U
!        real(kind=8) :: temperature
!        real(kind=8) :: E1_f_mode
!        real(kind=8) :: Kxl
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f
   real(kind=8),intent(in) :: e_gamma
   real(kind=8),intent(in) :: energy
!-----------------------------------------------------------
   integer(kind=4) :: kk
   integer(kind=4) :: l
   real(kind=8) :: U
   real(kind=8) :: ap
   real(kind=8) :: T
   real(kind=8) :: f
   real(kind=8) :: e_temp
!------------------------------------------------------------
!-----   External Functions
!------------------------------------------------------------
   real(kind=8) :: Ux
   real(kind=8) :: aparam_U
   real(kind=8) :: temperature
   real(kind=8) :: E1_f_mode
   real(kind=8) :: Kxl
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   E1_f = 0.0d0
   e_temp = energy
   if(.not. nucleus(i_f)%EL_mode(1)%gsf_read)then
      if(nucleus(i_f)%e1_model == 4)e_temp = nucleus(i_f)%sep_e(1)

      l = 1
      U = Ux(nucleus(i_f)%sep_e(1),nucleus(i_f)%level_param(3))
      ap = aparam_u(U,nucleus(i_f)%level_param(1),               &
                    nucleus(i_f)%level_param(4),                 &
                    nucleus(i_f)%level_param(5))
      U = Ux(e_temp - e_gamma, nucleus(i_f)%level_param(3))
      T = temperature(U,ap)

      f = 0.0

      do kk = 1, nucleus(i_f)%EL_mode(1)%num_gsf
         f = f + E1_f_mode(e_gamma,T,nucleus(i_f)%e1_model,      &
                           nucleus(i_f)%EL_mode(1)%gsf(kk)%er,   &
                           nucleus(i_f)%EL_mode(1)%gsf(kk)%gr,   &
                           nucleus(i_f)%EL_mode(1)%gsf(kk)%sr)
      end do
      E1_f = Kxl(l)*f
   end if

   return

end function E1_f
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function E1_f_component(i_f, k, e_gamma, energy)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns the kth component of the E1 strength function, f
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: Ux
!        real(kind=8) :: aparam_U
!        real(kind=8) :: temperature
!        real(kind=8) :: E1_f_mode
!        real(kind=8) :: Kxl
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, k
   real(kind=8),intent(in) :: e_gamma
   real(kind=8),intent(in) :: energy
!-----------------------------------------------------------
   integer(kind=4) :: l
   real(kind=8) :: U
   real(kind=8) :: ap
   real(kind=8) :: T
   real(kind=8) :: f
   real(kind=8) :: e_temp
!------------------------------------------------------------
!-----   External Functions
!------------------------------------------------------------
   real(kind=8) :: Ux
   real(kind=8) :: aparam_U
   real(kind=8) :: temperature
   real(kind=8) :: E1_f_mode
   real(kind=8) :: Kxl
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   E1_f_component = 0.0d0
   e_temp = energy
   if(.not. nucleus(i_f)%EL_mode(1)%gsf_read)then
      if(nucleus(i_f)%e1_model == 4)e_temp = nucleus(i_f)%sep_e(1)

      L = 1
      U = Ux(nucleus(i_f)%sep_e(1),nucleus(i_f)%level_param(3))
      ap = aparam_u(U,nucleus(i_f)%level_param(1),               &
                    nucleus(i_f)%level_param(4),                 &
                    nucleus(i_f)%level_param(5))
      U = Ux(e_temp - e_gamma, nucleus(i_f)%level_param(3))
      T = temperature(U,ap)

      f = 0.0

      E1_f_component = Kxl(L)*E1_f_mode(e_gamma,T,nucleus(i_f)%e1_model,    &
                                        nucleus(i_f)%EL_mode(1)%gsf(k)%er,  &
                                        nucleus(i_f)%EL_mode(1)%gsf(k)%gr,  &
                                        nucleus(i_f)%EL_mode(1)%gsf(k)%sr)
   end if

   return

end function E1_f_component
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function EL_f(i_f, L, e_gamma, energy)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E(L) strength function, f
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: E1_f_component
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, L
   real(kind=8),intent(in) :: e_gamma, energy
!------------------------------------------------------------------------------
   integer(kind=4) :: k
!-----   External functions      ----------------------------------------------
   real(kind=8) :: EL_f_component
!------------------------------------------------------------------------------

   EL_f = 0.0d0

   do k = 1, nucleus(i_f)%EL_mode(L)%num_gsf
      EL_f = EL_f + EL_f_component(i_f, L, k, e_gamma, energy)
   end do

   return

end function EL_f
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function EL_f_component(i_f, L, k, e_gamma, energy)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns the kth component of the E(L) strength function, f
!
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: E1_f_component
!        real(kind=8) :: Lorentzian
!        real(kind=8) :: Kxl
!        real(kind=8) :: interpolate_grid
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, L, k
   real(kind=8),intent(in) :: e_gamma, energy
!------------------------------------------------------------------------------
   real(kind=8) :: er, gr, sr
   real(kind=8) :: factor
!-----   External functions      ---------------------------------------
   real(kind=8) :: E1_f_component
   real(kind=8) :: Lorentzian
   real(kind=8) :: Kxl
   real(kind=8) :: interpolate_grid
!------------------------------------------------------------------------------

   EL_f_component = 0.0d0
   if(L == 1)then
      if(.not. nucleus(i_f)%EL_mode(L)%gsf_read)then
         EL_f_component = E1_f_component(i_f, k, e_gamma, energy)
      else
         EL_f_component = interpolate_grid(e_gamma, nucleus(i_f)%EL_mode(L)%gsf(k)%num_data, &
                                           nucleus(i_f)%EL_mode(L)%gsf(k)%e_gsf_r,           &
                                           nucleus(i_f)%EL_mode(L)%gsf(k)%gsf_r)*            &
                          nucleus(i_f)%EL_mode(L)%gsf(k)%gsf_norm
      end if
   else
      if(.not. nucleus(i_f)%EL_mode(L)%gsf_read)then
         factor = 0.0d0
         if(e_gamma > 0.0d0)factor = e_gamma**(2-2*L)
         er = nucleus(i_f)%EL_mode(L)%gsf(k)%er
         gr = nucleus(i_f)%EL_mode(L)%gsf(k)%gr
         sr = nucleus(i_f)%EL_mode(L)%gsf(k)%sr
         EL_f_component = Kxl(L)*factor*Lorentzian(e_gamma,er,gr)*sr
      else
         EL_f_component = interpolate_grid(e_gamma, nucleus(i_f)%EL_mode(L)%gsf(k)%num_data, &
                                           nucleus(i_f)%EL_mode(L)%gsf(k)%e_gsf_r,           &
                                           nucleus(i_f)%EL_mode(L)%gsf(k)%gsf_r)*            &
                          nucleus(i_f)%EL_mode(L)%gsf(k)%gsf_norm

      end if
   end if
   return

end function EL_f_component
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function EL_trans(i_f, L , e_gamma, energy)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E(L) transmission coefficient
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: EL_f
!        real(kind=8) :: Txl
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, L
   real(kind=8),intent(in) :: e_gamma, energy
!-----   External functions      -----------------------------------------------
   real(kind=8) :: EL_f
   real(kind=8) :: Txl
!-------------------------------------------------------------------------------
   
   EL_trans = Txl(e_gamma, L)*EL_f(i_f, L, e_gamma, energy)

   return

end function EL_trans
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function ML_f(i_f, L , e_gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the M(L) strength coefficient
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: ML_f_component
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, L
   real(kind=8),intent(in) :: e_gamma
!------------------------------------------------------------------------------
   integer(kind=4) :: k
!-----   External functions      ----------------------------------------------
   real(kind=8) :: ML_f_component
!------------------------------------------------------------------------------

   ML_f = 0.0d0
   do k = 1, nucleus(i_f)%ML_mode(L)%num_gsf
      ML_f = ML_f + ML_f_component(i_f, L, k, e_gamma)
   end do

   return

end function ML_f
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function ML_f_component(i_f, L, k, e_gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns the kth component of the M(L) strength function
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: Lorentzian
!        real(kind=8) :: Kxl
!        real(kind=8) :: interpolate_grid
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, L, k
   real(kind=8),intent(in) :: e_gamma
!------------------------------------------------------------------------------
   real(kind=8) :: er, gr, sr
   real(kind=8) :: factor
!-----   External functions      ---------------------------------------
   real(kind=8) :: Lorentzian
   real(kind=8) :: Kxl
   real(kind=8) :: interpolate_grid

   ML_f_component = 0.0d0
   if(.not. nucleus(i_f)%ML_mode(L)%gsf_read)then
      factor = 0.0d0
      if(e_gamma > 0.0d0)factor = e_gamma**(2-2*L)
      er = nucleus(i_f)%ML_mode(L)%gsf(k)%er
      gr = nucleus(i_f)%ML_mode(L)%gsf(k)%gr
      sr = nucleus(i_f)%ML_mode(L)%gsf(k)%sr
      ML_f_component = Kxl(L)*Lorentzian(e_gamma,er,gr)*sr
   else
      ML_f_component = interpolate_grid(e_gamma,nucleus(i_f)%ML_mode(L)%gsf(k)%num_data, &
                                        nucleus(i_f)%ML_mode(L)%gsf(k)%e_gsf_r,          &
                                        nucleus(i_f)%ML_mode(L)%gsf(k)%gsf_r)*           &
                       nucleus(i_f)%ML_mode(L)%gsf(k)%gsf_norm
   end if

   return
end function ML_f_component
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function ML_trans(i_f, L , e_gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the M(L) transmission coefficient
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: ML_f
!        real(kind=8) :: Txl
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, L
   real(kind=8),intent(in) :: e_gamma
!-------    External functions    ---------------------------------------------
   real(kind=8) :: ML_f
   real(kind=8) :: Txl
!------------------------------------------------------------------------------

   ML_trans = Txl(e_gamma, L)*ML_f(i_f, L, e_gamma)

   return

end function ML_trans
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function EL_absorption(icomp, L, egamma, ex)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E(L) photoabsorption cross section
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: EL_f
!        real(kind=8) :: Kxl
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: icomp
   integer(kind=4),intent(in) :: L
   real(kind=8),intent(in) :: egamma
   real(kind=8),intent(in) :: ex
!------------------------------------------------------------------------------
   real(kind=8) :: factor
!-----   External functions      ---------------------------------------
   real(kind=8) :: EL_f
   real(kind=8) :: Kxl
!------------------------------------------------------------------------------

   factor = egamma**(2*L-1)
   EL_absorption = 0.001d0*factor*EL_f(icomp, L, egamma, ex)/Kxl(L)

   return

end function EL_absorption
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function ML_absorption(icomp, L, egamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the M(L) photoabsorption cross section
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: ML_f
!        real(kind=8) :: Kxl
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: icomp
   integer(kind=4),intent(in) :: L
   real(kind=8),intent(in) :: egamma
!------------------------------------------------------------------------------
   real(kind=8) :: factor
!-----   External functions          ------------------------------------------
   real(kind=8) :: ML_f
   real(kind=8) :: Kxl
!------------------------------------------------------------------------------

   factor = egamma**(2*L-1)

   ML_absorption = 0.001d0*factor*ML_f(icomp, L, egamma)/Kxl(L)

   return

end function ML_absorption

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function photo_absorption(icomp, lmax_E, lmax_M, egamma, ex)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the M(L) photoabsorption cross section
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: EL_absorption
!        real(kind=8) :: ML_absorption
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
!*******************************************************************************
!
   use variable_kinds
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: icomp
   integer(kind=4),intent(in) :: lmax_E
   integer(kind=4),intent(in) :: lmax_M
   real(kind=8),intent(in) :: egamma
   real(kind=8),intent(in) :: ex
!------------------------------------------------------------------------------
   integer(kind=4) :: L
   real(kind=8) :: sum
!------   External functions            ---------------------------------------
   real(kind=8) :: EL_absorption
   real(kind=8) :: ML_absorption

   sum = 0.0d0
   do L = 1, lmax_E
      sum = sum + EL_absorption(icomp, L, egamma, ex)
   end do
   do L = 1, lmax_M
      sum = sum + ML_absorption(icomp, L, egamma)
   end do
   
   photo_absorption = sum

   return

end function photo_absorption
!
!*****************************************************************************80
!
real(kind=8) function interpolate_grid(x_in, num, x, y)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This function interpolates between two points on an exponential grid
!    here the x-grid was constructed by x(i) = x(1)*factor**(i-1)
!    routine finds a quadratic "fit" along the grid, namely y = a + b*x + c*x**2
!    for three points along the grid.
!    The function will do an interpolation that is either linear in y (itype == 0)
!    or log in y (itype == 1).
!
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!
!     Subroutines:
!
!        exit_YAHFC
!
!     External functions:
!
!         det3_3
!
!     MPI routines:
!
!        MPI_Abort   ----    via exit_YAHFC
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
   use nodeinfo
   implicit none
   real(kind=8), intent(in) :: x_in
   integer(kind=4), intent(in) :: num
   real(kind=8), intent(in) :: x(num), y(num)
!------------------------------------------------------
   real(kind=8) :: mat(3,3), mat0(3,3), vecy(3), a(0:2)
   integer(kind=4) :: i,j,k
   integer(kind=4) im1, i0, ip1
   real(kind=8) :: xm1, x0, xp1, xx
   real(kind=8) :: ym1, y0, yp1, yy
   real(kind=8) :: denom,numer
   real(kind=8) :: diff1, diff2
!-----   External functions      ---------------------------------------
   real(kind=8) :: det3_3
!-----------------------------------------------------------------------

  if(x_in > x(num))then
     if(iproc == 0)then
        write(6,*)num
        write(6,*)x_in,x(num)
        write(6,*)'E_in larger than input grid'
     end if
     call exit_YAHFC(101)
  end if

  if(abs(x_in - x(1)) <= 1.0d-6)then                  !  it is exactly on a grid point
     interpolate_grid = y(1)
     return
  end if

  i0 = 0
  if(x_in < x(1))then
     i = 1
  else
     do i = 1, num - 1                                  !  start from the bottom of the grid
        if(abs(x_in - x(i)) <= 1.0d-6)then                  !  it is exactly on a grid point
           interpolate_grid = y(i)
           return
        elseif(x_in > x(i) .and. x_in < x(i+1))then          ! grid(i) <= e <= grid(i+1)
!-----    Sandwiched between two points. Now which one is it closest too?
!-----    use the closest point on the grid for the center of the quadratic fit
           diff1 = abs(x(i) - x_in)
           diff2 = abs(x(i+1) - x_in)
           if(diff1 < diff2)then
              i0 = i
           else
              i0 = i + 1
           end if
           exit
        end if
     end do
  end if
  if(i == 1)then                   ! at the start, so we need to use i0 = 2
     im1 = 1
     i0 = 2
     ip1 = 3
  elseif(i == num)then             !  at the end, so we use i0 = num - 1
     ip1 = num
     i0 = num - 1
     im1 = num - 2
  else                             !  in the middle
     im1 = i0 - 1
     ip1 = i0 + 1
  end if
  if(ip1 >= num)then
     ip1 = num
     i0 = ip1 - 1
     im1 = i0 - 1
  end if
  xm1 = x(im1)
  x0 = x(i0)
  xp1 = x(ip1)
  ym1 = y(im1)
  y0 = y(i0)
  yp1 = y(ip1)
  mat0(1,1) = 1.0d0
  mat0(2,1) = 1.0d0
  mat0(3,1) = 1.0d0
  mat0(1,2) = xm1
  mat0(2,2) = x0
  mat0(3,2) = xp1
  mat0(1,3) = xm1**2
  mat0(2,3) = x0**2
  mat0(3,3) = xp1**2
  vecy(1) = ym1
  vecy(2) = y0
  vecy(3) = yp1

  denom = det3_3(mat0)

  xx = x_in
  yy = 0.0d0
  do k = 1, 3
     do i = 1, 3
        do j = 1, 3
           mat(i,j) = mat0(i,j)
        end do
     end do
     do i = 1, 3
        mat(i,k) = vecy(i)
     end do
     numer = det3_3(mat)
     a(k-1) = numer/denom
     yy = yy + a(k-1)*xx**(k-1)
  end do
  interpolate_grid = yy

  return
end function interpolate_grid

