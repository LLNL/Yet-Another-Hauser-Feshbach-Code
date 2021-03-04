!
!*******************************************************************************
!
subroutine gdr_param(data_path, len_path, iz, ia, err, grr, srr)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine reads in and sets up data for E1 transitions
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
   character(len=200), intent(in) :: data_path      ! path where data files are kept
   integer(kind=4), intent(in) :: len_path             ! length of data_path
   integer(kind=4), intent(in) :: iz, ia
   real(kind=8), intent(out) :: err(3), grr(3), srr(3)
!--------------------------------------------------------------------------
   integer(kind=4) :: i
   integer(kind=4) :: eof
   logical :: found
   integer(kind=4) :: izz,iaa
   real(kind=8) :: xA, xZ
   character(len=2) :: symbb
   real(kind=8) :: er1, gr1, sr1, er2, gr2, sr2
   integer(kind=4) :: iread, jread
!   real(kind=8) :: anat(100)
!   data anat/                                                                         &
!      1.0d0,  4.0d0,  6.9d0,  9.0d0, 10.8d0, 12.0d0, 14.0d0, 16.0d0, 19.0d0, 20.2d0,  &
!     23.0d0, 24.3d0, 27.0d0, 28.1d0, 31.0d0, 32.1d0, 35.5d0, 39.9d0, 39.1d0, 40.1d0,  &
!     45.0d0, 47.9d0, 50.9d0, 52.0d0, 54.9d0, 55.8d0, 58.9d0, 58.7d0, 63.5d0, 65.4d0,  &
!     69.7d0, 72.6d0, 74.9d0, 79.0d0, 79.9d0, 83.8d0, 85.5d0, 87.6d0, 88.9d0, 91.2d0,  &
!     92.9d0, 95.9d0, 98.0d0,101.1d0,102.9d0,106.4d0,107.8d0,112.4d0,114.8d0,118.7d0,  &
!    121.8d0,127.6d0,126.9d0,131.3d0,132.9d0,137.3d0,138.9d0,140.1d0,140.9d0,144.2d0,  &
!    145.0d0,150.4d0,152.0d0,157.3d0,158.9d0,162.5d0,164.9d0,167.3d0,168.9d0,173.0d0,  &
!    175.0d0,189.5d0,180.9d0,183.8d0,186.2d0,190.2d0,192.2d0,195.1d0,197.0d0,200.6d0,  &
!    204.4d0,207.2d0,209.0d0,209.0d0,210.0d0,222.0d0,223.0d0,226.0d0,227.0d0,232.0d0,  &
!    231.0d0,238.0d0,237.0d0,244.0d0,243.0d0,247.0d0,247.0d0,251.0d0,252.0d0,257.0d0/
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   jread  =  -1
   do i = 1, 3
      err(i) = 0.0d0
      grr(i) = 0.0d0
      srr(i) = 0.0d0
   end do
   open(unit=52,                                               &
        file=data_path(1:len_path)//'gdr-parameters-exp.dat',  &
        status='old')
   do i = 1, 4
      read(52,*)
   end do
   iread = 0
   found = .false.
   eof = 0
   do while(eof == 0)
      read(52,'(2(1x,i3),1x,a2,6x,6(1x,f6.2))',iostat=eof)             &
           izz,iaa,symbb,er1,sr1,gr1,er2,sr2,gr2
      iread = iread + 1
      if(iz == izz .and. ia == iaa)then
         err(1) = er1
         grr(1) = gr1
         srr(1) = sr1
         err(2) = er2
         grr(2) = gr2
         srr(2) = sr2
         found = .true.
         exit
      end if
      if(iz == izz .and. iaa == 0)then     !   use natural GDR parameters
         err(1) = er1
         grr(1) = gr1
         srr(1) = sr1
         err(2) = er2
         grr(2) = gr2
         srr(2) = sr2
         found = .true.
         exit
      end if
   end do

   close(unit=52)

   if(found)return

   xA = real(iA,kind=8)
   xZ = real(iZ,kind=8)
   err(1) = 31.2d0*xA**(-1.d0/3.d0)+20.6d0*xA**(-1.d0/6.d0)
   grr(1) = 0.026d0*err(1)**1.91d0
   srr(1) = 144.0d0*real((ia-iz)*iz,kind=8)/(xA*3.141593*grr(1))
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
   use constants
   implicit none
   real(kind=8), intent(in) :: e_gam,T
   integer(kind=4), intent(in) :: e1_model
   real(kind=8), intent(in) :: er, gr, sr
!---------------------------------------------------------------------
   real(kind=8) gam, str
   real(kind=8) gamma_t
   real(kind=8) M_Lorentzian
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   E1_f_mode = 0.0d0
   if(er < 1.0d-5.or.gr < 1.0d-5)return
   if(e1_model > 1)then
      gam = gamma_t(e_gam,er,gr,T)
      if(e1_model == 3)then
         if(gam > 1.5*gr) gam = 1.5*gr
         if(gam < 0.25*gr) gam = 0.25*gr
      end if
   else
      gam = gr
   end if
   str = M_lorentzian(e_gam,er,gam)
   str = e_gam*str
   if(e1_model > 1)str = str + 0.7*gr*(2.0*pi*T)**2/er**5      !  from RIPL handbook
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
   real(kind=8), intent(in) :: e_g,er,gr
!----------------------------------------------------------------------
   lorentzian = e_g*gr**2/((e_g**2-er**2)**2+(e_g*gr)**2)
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
   real(kind=8), intent(in) :: e_g,er,gr
!--------------------------------------------------------------------------
   M_Lorentzian = gr/((e_g**2-er**2)**2+(e_g*gr)**2)
   return
end function M_Lorentzian
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function txl(egamma, l_radiation)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the multiplicative factor to convert the 
!    strength function to the transmission coefficient
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
   use constants
   implicit none
   real(kind=8),intent(in) :: egamma
   integer(kind=4),intent(in) :: l_radiation
!------------------------------------------------------------------------------
   txl = two_pi*egamma**(2*l_radiation+1)
   return
end function txl
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function kxl(l_radiation)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes k factor needed for E&M radiation
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
   use constants
   implicit none
   integer(kind=4),intent(in) :: l_radiation
!------------------------------------------------------------------------------
   kxl=dfloat(2*l_radiation+1)*(pi*hbar_c)**2*fmsq_eq_mbarn
   kxl=1.0d0/kxl
   return
end function kxl
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
   real(kind=8),intent(in) :: U, ap
!------------------------------------------------------------------------------
   temperature=0.0
   if(U > 0.0)temperature = sqrt(U/ap)
   return
end function temperature
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function Lfac(egam,T)
   use variable_kinds
   implicit none
   real(kind=8) egam, T
!------------------------------------------------------------------------------
   if(egam < 1.0d-6)then
      Lfac = 0.0d0
   else
      Lfac = 1.0d0
      if(T > 0.0)Lfac = 1.0/(1.0-exp(-egam/T))
   end if
   return
end function Lfac
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function E1_trans(i_f, energy, e_gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes E1 transmission coefficient
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
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f
   real(kind=8),intent(in) :: energy
   real(kind=8),intent(in) :: e_gamma
!-----------------------------------------------------------
   integer(kind=4) :: l
   real(kind=8) :: U
   real(kind=8) :: ap
   real(kind=8) :: T
!------------------------------------------------------------
!-----   External Functions
!------------------------------------------------------------
   real(kind=8) :: Ux
   real(kind=8) :: aparam_U
   real(kind=8) :: temperature
   real(kind=8) :: E1_f
   real(kind=8) :: Kxl
   real(kind=8) :: Txl   

   l = 1
   U = Ux(nucleus(i_f)%sep_e(1),nucleus(i_f)%level_param(3))
   ap = aparam_u(U,nucleus(i_f)%level_param(1),               &
               nucleus(i_f)%level_param(4),                   &
               nucleus(i_f)%level_param(5))
   U = Ux(energy-e_gamma,nucleus(i_f)%level_param(3))
   T = temperature(U,ap)

   E1_trans = Txl(e_gamma,l)*Kxl(l)*E1_f(i_f, energy, e_gamma)

   return

end function E1_trans
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function E1_f(i_f,energy,e_gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E1 strength function, f
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
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f
   real(kind=8),intent(in) :: energy
   real(kind=8),intent(in) :: e_gamma
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
   e_temp = energy
   if(nucleus(i_f)%e1_model == 4)e_temp = nucleus(i_f)%sep_e(1)

   l = 1
   U = Ux(nucleus(i_f)%sep_e(1),nucleus(i_f)%level_param(3))
   ap = aparam_u(U,nucleus(i_f)%level_param(1),               &
                 nucleus(i_f)%level_param(4),                 &
                 nucleus(i_f)%level_param(5))
   U = Ux(e_temp - e_gamma, nucleus(i_f)%level_param(3))
   T = temperature(U,ap)

   f = 0.0

   do kk = 1, nucleus(i_f)%num_res
      f = f + E1_f_mode(e_gamma,T,nucleus(i_f)%e1_model,   &
                        nucleus(i_f)%er_E1(kk),            &
                        nucleus(i_f)%gr_E1(kk),            &
                        nucleus(i_f)%sr_E1(kk))
   end do
   E1_f = Kxl(l)*f

   return

end function E1_f
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function Kopecky_trans(i_f,energy,e_gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E1 transmision coefficient with the Kopecky
!    model for the strength function
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
   use nuclei
   implicit none
   integer(kind=4),intent(in):: i_f
   real(kind=8),intent(in) :: energy
   real(kind=8),intent(in) :: e_gamma
!-----------------------------------------------------------
   integer(kind=4) :: l
!------------------------------------------------------------
!-----   External Functions
!------------------------------------------------------------
   real(kind=8) :: Txl   
   real(kind=8) :: Kopecky_f
!------------------------------------------------------------
!-----   Start Calculation


   l = 1
   Kopecky_trans = Txl(e_gamma,l)*Kopecky_f(i_f, energy, e_gamma)
   return
end function Kopecky_trans
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function Kopecky_f(i_f,energy,e_gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E1 strength function with the Kopecky model
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
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f
   real(kind=8),intent(in) :: energy
   real(kind=8),intent(in) :: e_gamma
!-----------------------------------------------------------
   integer(kind=4) :: kk
   integer(kind=4) :: l
   real(kind=8) :: U
   real(kind=8) :: ap
   real(kind=8) :: T
   real(kind=8) :: f
!------------------------------------------------------------
!-----   External Functions
!------------------------------------------------------------
   real(kind=8) :: Ux
   real(kind=8) :: aparam_U
   real(kind=8) :: temperature
   real(kind=8) :: E1_f_mode
   real(kind=8) :: Kxl
!------------------------------------------------------------
!-----   Start Calculation
   l = 1
   U = Ux(nucleus(i_f)%sep_e(1),nucleus(i_f)%level_param(3))
   ap = aparam_u(U,nucleus(i_f)%level_param(1),               &
               nucleus(i_f)%level_param(4),                   &
               nucleus(i_f)%level_param(5))
   U = Ux(energy-e_gamma,nucleus(i_f)%level_param(3))
   T = temperature(U,ap)
   f = 0.0
   do kk = 1, nucleus(i_f)%num_res
      f = f + E1_f_mode(e_gamma,T,nucleus(i_f)%e1_model,       &
                        nucleus(i_f)%er_E1(kk),                &
                        nucleus(i_f)%gr_E1(kk),                &
                        nucleus(i_f)%sr_E1(kk))
   end do
   f = Kxl(l)*f
   Kopecky_f = f
   return
end function Kopecky_f
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function EL_f(i_f, l_radiation , e_gamma, energy)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E(L) strength function, f
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
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, l_radiation
   real(kind=8),intent(in) :: e_gamma, energy
!------------------------------------------------------------------------------
   real(kind=8) :: E1_f, Lorentzian, Kxl

   if(l_radiation == 1)then
      EL_f = E1_f(i_f,energy,e_gamma)
   else
      EL_f = Kxl(l_radiation)*                                                         &
           Lorentzian(e_gamma,nucleus(i_f)%er_E(l_radiation),                        &
                      nucleus(i_f)%gr_E(l_radiation))*nucleus(i_f)%sr_E(l_radiation)
   end if
   return

end function EL_f
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function EL_trans(i_f, l_radiation , e_gamma, energy)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the E(L) transmission coefficient
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
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, l_radiation
   real(kind=8),intent(in) :: e_gamma, energy
!------------------------------------------------------------------------------
   real(kind=8) :: E1_f, EL_f, Txl
   
   if(l_radiation == 1)then
      EL_trans = Txl(e_gamma, l_radiation)*E1_f(i_f, energy, e_gamma)
   else
      EL_trans = Txl(e_gamma, l_radiation)*EL_f(i_f, l_radiation, e_gamma, energy)
   end if
   return

end function EL_trans
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function ML_f(i_f, l_radiation , e_gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the M(L) strength coefficient
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
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, l_radiation
   real(kind=8),intent(in) :: e_gamma
!------------------------------------------------------------------------------
   real(kind=8) :: Lorentzian, Kxl

   ML_f = Kxl(l_radiation)*                                                       &
        Lorentzian(e_gamma,nucleus(i_f)%er_M(l_radiation),                        &
                   nucleus(i_f)%gr_M(l_radiation))*nucleus(i_f)%sr_M(l_radiation)

   return

end function ML_f
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real (kind=8) function ML_trans(i_f, l_radiation , e_gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the M(L) transmission coefficient
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
   use nuclei
   implicit none
   integer(kind=4),intent(in) :: i_f, l_radiation
   real(kind=8),intent(in) :: e_gamma
!------------------------------------------------------------------------------
   real(kind=8) :: ML_f, Txl

   ML_trans = Txl(e_gamma, l_radiation)*ML_f(i_f, l_radiation, e_gamma)

   return

end function ML_trans


