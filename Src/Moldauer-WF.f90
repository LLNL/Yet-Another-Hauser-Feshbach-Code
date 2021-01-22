!
!*******************************************************************************
!
subroutine Moldauer_WF(icomp,                                       &
                       k_a,l_a,xj_a,istate_a,xI_a,trans_a,          &
                       k_b,l_b,xj_b,istate_b,xI_b,                  &
                       ibin,ip,xI,energy,                           &
                       exp_gamma,HF_den,F_trans,CHnorm,WF)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine computes the Moldauer width-fluctuation correction
!
!  Reference:
!
!  S. Hilaire, Ch. Lagrange, and A. J. Koning, Ann. Phys. 306, 209 (2003)
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
   use options
   use nuclei
   use particles_def
   use constants 
   use nodeinfo
   use Gauss_integration
   implicit none
   integer(kind=4), intent(in) :: icomp
   integer(kind=4), intent(in) :: k_a, l_a, istate_a
   real(kind=8), intent(in) :: xj_a, xI_a, trans_a
   integer(kind=4), intent(in) :: k_b, l_b, istate_b
   real(kind=8), intent(in) :: xj_b, xI_b
   integer(kind=4), intent(in) :: ibin, ip
   real(kind=8), intent(in) :: xI
   real(kind=8), intent(in) :: energy, exp_gamma, HF_den, F_trans(4)
   real(kind=8), intent(in) :: CHnorm
   real(kind=8), intent(out) :: WF
!-------------------------------------------------------------------------+
   integer(kind=4) j_max
   integer(kind=4) nbin
   integer(kind=4) :: ix
   real(kind=8) :: x,dx
   real(kind=8) :: xnu_a
   real(kind=8) :: Product,Product_p,Product_g
   real(kind=8) :: weight(0:1)
   real(kind=8) :: factor,elastic

   real(kind=8) :: xstop

   integer(kind=4) :: ndata
   parameter (ndata = 10)
   real(kind=8) :: xx(ndata), yy(ndata), xxxx, yyxx
   real(kind=8) :: afit
   save xx
   data xx/ 0.2d0,  0.4d0, 0.6d0, 0.8d0, 1.0d0, 1.2d0, 1.4d0, 1.6d0, 1.8d0, 2.0d0/
   

   logical :: verbose

   real(kind=8) :: xnu
!-------------------------------------------------------------------------+
!------                                                                   +
!--------------------   Start subroutine                                  +
!------                                                                   +
!-------------------------------------------------------------------------+
   verbose=.false.

   weight(0)=2.0d0/3.0d0
   weight(1)=4.0d0/3.0d0
   xstop = 100.0
   dx=0.001
   j_max=nucleus(icomp)%j_max
   nbin=nucleus(icomp)%nbin

   factor = 1.0d0
   elastic = 0.0d0
   if(k_a == k_b .and. istate_a == istate_b .and. l_a == l_b .and.       &
      (abs(xj_a - xj_b) < 1.0d-3) .and. (abs(xI_a - xI_b) < 1.0d-3) )       &
       elastic = 1.0d0
   xnu_a = xnu(trans_a,HF_den)
   factor = factor + elastic*2.0d0/xnu_a

   xxxx = 0.0d0
   yyxx = 0.0d0
   do ix = 1, ndata
      x = xx(ix)
      call Moldauer_product(icomp,                                       &
                            k_a,l_a,xj_a,istate_a,                       &
                            k_b,l_b,xj_b,istate_b,xI_b,                  &
                            ibin,ip,xI,energy,                           &
                            HF_den,F_trans,x,CHnorm,Product_p)
      Product_g = exp(-exp_gamma*x)
      Product = Product_p*Product_g*factor
      yy(ix) = log(Product)
      xxxx = xxxx + xx(ix)**2
      yyxx = yyxx + yy(ix)*xx(ix)
   end do

   afit = abs(yyxx/xxxx)

   WF = 0.0d0
   do ix = 1, n_glag
      x = x_glag(ix)/afit
      call Moldauer_product(icomp,                                       &
                            k_a,l_a,xj_a,istate_a,                       &
                            k_b,l_b,xj_b,istate_b,xI_b,                  &
                            ibin,ip,xI,energy,                           &
                            HF_den,F_trans,x,CHnorm,Product_p)
      Product_g = exp(-exp_gamma*x)
      Product = Product_p*Product_g*factor
      WF = WF + Product/exp(-afit*x)*w_glag(ix)
   end do
   WF = WF/afit



return
end subroutine Moldauer_WF
!
!*******************************************************************************
!
real(kind=8) function xnu(trans_eff,HF_den)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the number of degrees of freedom nu
!
!  Reference:
!
!  S. Hilaire, Ch. lagrange, and A. J. Koning, Ann. Phys. 306, 209 (2003)
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
   real(kind=8) trans_eff,HF_den
!-----------    Calculation 
   xnu = 1.78d0 + (trans_eff**1.212d0 - 0.78d0)*exp(-0.228d0*HF_den)
   return
end function xnu
!
!*******************************************************************************
!
subroutine Moldauer_product(icomp,                               &
                            k_a,l_a,xj_a,istate_a,               &
                            k_b,l_b,xj_b,istate_b,xI_b,          &
                            ibin,ip,xI,energy,                   &
                            HF_den,F_trans,x,CHnorm,Product)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine computes the Moldauer product 
!
!  Reference:
!
!  S. Hilaire, Ch. lagrange, and A. J. Koning, Ann. Phys. 306, 209 (2003)
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
   use options
   use nuclei
   use particles_def
   use constants 
   use nodeinfo
   implicit none
   integer(kind=4), intent(in) :: icomp
   integer(kind=4), intent(in) :: k_a,l_a,istate_a
   real(kind=8), intent(in) :: xj_a
   integer(kind=4), intent(in) :: k_b,l_b,istate_b
   real(kind=8), intent(in) :: xj_b,xI_b
   integer(kind=4), intent(in) :: ibin,ip
   real(kind=8), intent(in) :: xI
   real(kind=8), intent(in) :: energy, HF_den, F_trans(4)
   real(kind=8), intent(in) :: x
   real(kind=8), intent(in) :: CHnorm
   real(kind=8), intent(out) :: Product
!-------------------------------------------------------------------------+
   integer(kind=4) :: num, num_bin
   integer(kind=4) :: isc,iss
   integer(kind=4) :: num_discrete
   integer(kind=4) j
   integer(kind=4) j_max
   integer(kind=4) l_c_min,l_c_max
   integer(kind=4) nbin, num_j
   integer(kind=4) if1, i_c, n_c, ip_c
   integer(kind=4) k_c,l_c,EM_c
   real(kind=8) :: xj_c,xj_c_min,xj_c_max,xj_c_min1,xj_c_max1
   real(kind=8) :: xI_c,xI_c_min,xI_c_max
   integer(kind=4) :: Ix_c,Ix_c_min,Ix_c_max
   real(kind=8) :: p_spin
   real(kind=8) :: E_c, E_f
   real(kind=8) :: trans 
   real(kind=8) :: par,cpar2
   real(kind=8) :: sum_Tg,sum_n,sum_p
   integer(kind=4) :: nptsx
   real(kind=8) :: trans_eff
   real(kind=8) :: xnu_c
   real(kind=8) :: exponent
   real(kind=8) :: N_eff
   real(kind=8) :: P_f
   real(kind=8) :: de
   logical converged
   real(kind=8) :: K_vib, K_rot
   real(kind=8) :: jfac, pfac, rho, tt, T_f, Ef
   integer(kind=4) :: ib
   real(kind=8) :: rho_Fm, apu, sig2
   real(kind=8) :: aa, bb, cc
   real(kind=8) :: F_Barrier, F_Barrier_hbw
   real(kind=8) :: Max_J
   real(kind=8) :: e1, b, bbb, mode
   real(kind=8) :: E0, T, E11,ecut
   integer(kind=4) :: j_min

!-------------------------------------------------------------------------+
!------     Function declarations
   real(kind=8) :: tco_interpolate
   real(kind=8) :: xnu
   real(kind=8) :: exp_1,exp_2
   real(kind=8) :: spin_fac
   real(kind=8) :: parity_fac
   real(kind=8) :: HW_trans
!-------------------------------------------------------------------------+
!------                                                                   +
!--------------------   Start subroutine                                  +
!------                                                                   +
!-------------------------------------------------------------------------+
   sum_n = 0.0d0
   sum_p = 0.0d0
   sum_Tg = 0.0d0
   nptsx = 1000
   j_max = nucleus(icomp)%j_max
   nbin = nucleus(icomp)%nbin

   Product = 0.0d0
   num = 0
   num_bin = 0
   par = 2.0*real(ip)-1.0
   do if1=1,nucleus(icomp)%num_decay                       !  loop over nuclei in the decay chain
      i_c = nucleus(icomp)%decay_to(if1)
      k_c = nucleus(icomp)%decay_particle(if1)
      if(k_c == 0)cycle                                    ! particle n,p,d,t,h,a  ! Skip photons - treated woth exp(-exp_gamma*x)
      if(energy < nucleus(icomp)%sep_e(k_c))cycle           !   not eneough energy to decay - cycle out
      xj_c_max1 = real(nucleus(i_c)%j_max,kind=8) + nucleus(i_c)%jshift
      EM_c = 0
!--------------------------   particle decay to continuous level bins
      p_spin = particle(k_c)%spin
      do n_c = nucleus(i_c)%nbin - ibin,1, -1                !  loop over final excitation energies
         e_f = energy - nucleus(icomp)%sep_e(k_c) -                                       &
               nucleus(i_c)%e_grid(n_c)
         if(e_f <= 0.0d0)cycle
         do l_c = 0, particle(k_c)%lmax                      !  loop over l-partial wave
            cpar2 = par*particle(k_c)%par*(-1.0d0)**l_c      !  parity for channel c
            ip_c = nint((cpar2 + 1.0d0)/2.0d0)                     !  parity index for channel c
            xj_c_min = abs(real(l_c,kind=8) - p_spin)
            xj_c_max = real(l_c,kind=8) + p_spin
            isc = nint(xj_c_max - xj_c_min)
            do iss = 0, isc
               trans = tco_interpolate(e_f,particle(k_c)%nume,                            &
                                       particle(k_c)%e_grid,                              &
                                       particle(k_c)%trans_read(1,iss,l_c))
               if(trans < trans_p_cut)cycle
               xj_c = xj_c_min + real(iss,kind=8)
               xI_c_min = abs(xI - xj_c)
               xI_c_max = xI + xj_c
               Ix_c_min = max(nint(xI_c_min - nucleus(i_c)%jshift),0)
               Ix_c_max = min(nint(xI_c_max - nucleus(i_c)%jshift),nucleus(i_c)%j_max)
               do Ix_c = Ix_c_min, Ix_c_max                    !  loop over final j
                  xI_c = real(Ix_c,kind=8) + nucleus(i_c)%jshift
                  N_eff = nucleus(i_c)%bins(Ix_c,ip_c,n_c)%rho*   &
                          nucleus(i_c)%delta_e(n_c)
                  trans_eff = trans*N_eff
                  if(trans_eff/CHnorm < prob_cut)cycle
                  exp_2 = 0.0d0
                  if(k_b == k_c .and. l_b == l_c .and. istate_b == -n_c .and.             &
                     abs(xj_b - xj_c) < 1.0d-3 .and. abs(xI_c - xI_b) < 1.0d-3)           &
                     exp_2 = exp_2 + 1.0d0
                  xnu_c = xnu(trans,HF_den)
                  exponent = -0.5d0*xnu_c*N_eff - exp_2
                  Product = Product +                                                     &
                    log(1.0d0 + 2.0d0*trans*x/(xnu_c*HF_den))*exponent
               end do
            end do
         end do
      end do
!---------------------------   particle decay to discrete states
      num_discrete = nucleus(i_c)%ncut
      if(All_gammas)num_discrete = nucleus(i_c)%num_discrete
      do n_c = 1, num_discrete
         E_c = energy - nucleus(icomp)%sep_e(k_c) - nucleus(i_c)%state(n_c)%energy
         if(E_c < 0.0)cycle
         xI_c = nucleus(i_c)%state(n_c)%spin
         xj_c_min = abs(xI - xI_c)
         xj_c_max = xI + xI_c
         num_j = nint(xj_c_max - xj_c_min)
         do j = 0, num_j
            xj_c = real(j,kind=8) + xj_c_min

            l_c_min = nint(abs(xj_c - p_spin))
            l_c_max = min(particle(k_c)%lmax, nint(xj_c + p_spin))
            cpar2 = particle(k_c)%par*nucleus(i_c)%state(n_c)%parity
            ip_c = nint((cpar2+1.0d0)/2.0d0)
            if(ip == ip_c)then                            !   parities are the same, l=even
               if(iand(l_c_min,1) == 1)l_c_min=l_c_min+1                  !   odd l_c_min, add 1 to make it even
               if(iand(l_c_max,1) == 1)l_c_max=l_c_max-1                  !   odd l_c_max, subtract 1 to make it even
            else                                          !   parities are different, l=odd
               if(iand(l_c_min,1) == 0)l_c_min=l_c_min+1                  !   even l_c_min, add 1 to make it even
               if(iand(l_c_max,1) == 0)l_c_max=l_c_max-1                  !   even l_c_max, subtract 1 to make it even
            end if
            do l_c = l_c_min,l_c_max,2
               xj_c_min1 = abs(real(l_c,kind=8) - p_spin)
               iss = nint(xj_c - xj_c_min1)
               if(iss < 0 .or. iss > nint(2*p_spin))cycle
               trans = tco_interpolate(E_c,particle(k_c)%nume,                            &
                                      particle(k_c)%e_grid,                               &
                                      particle(k_c)%trans_read(1,iss,l_c))  
               if(trans < trans_p_cut)cycle
               if(trans/CHnorm < prob_cut)cycle
               exp_1 = 0.0d0
               exp_2 = 0.0d0
               if(k_a == k_c .and. l_a == l_c .and.                                       &
                  abs(xj_a - xj_c) < 1.0d-3 .and. istate_a == n_c )exp_1=1.0d0
               if(k_b == k_c .and. l_b == l_c .and.                                       &
                  abs(xj_b - xj_c) < 1.0d-3 .and. istate_b == n_c )exp_2=1.0d0
               xnu_c = xnu(trans,HF_den)
               exponent = -0.5d0*xnu_c - exp_1 - exp_2
               Product = Product +                                                        &
                  log(1.0d0 + 2.0d0*trans*x/(xnu_c*HF_den))*exponent
            end do
         end do
      end do
   end do

!-------------------------------------------------------------------*
!---------    Now fission, new approach                             *
!-------------------------------------------------------------------*
!------   Talys approach is a bit odd for double-humped fission barriers
!------   as the total transmission coefficient is T = T_1*T_2/(T_1+T_2)
!------   where each T_i is integrated over the fission level density
!------   So, a single prescription as they describe is hard to do except 
!------   for a single barrier
!------   Thus, I convert to an effective single barrier for the first barrier
!------   Define P_F= T_2/(T_1+T_2)  then
!------   T = integral tt(E)*P_F*rho(E)dE with integral over transition states 
!------   above barrier 1. This becomes more like particle decay to the continuum
!------   We treat this integral as separate channels in the Moldauer expression.
   if(nucleus(icomp)%fission) then
      if(F_trans(4)/CHnorm > prob_cut)then
         P_f = 1.0d0
         ib = 1
         if(nucleus(icomp)%F_n_barr == 1)then
            ib = 1
            P_f = 1.0d0
         else if(nucleus(icomp)%F_n_barr == 2)then
            ib = 1
            P_f = F_trans(2)/(F_trans(1)+F_trans(2))
            if(nucleus(icomp)%F_Barrier(1)%barrier < nucleus(icomp)%F_Barrier(2)%barrier)then
               ib = 2
               P_f = F_trans(1)/(F_trans(1)+F_trans(2))
           end if
         else if(nucleus(icomp)%F_n_barr == 3)then
            ib = 1
            P_f = (F_trans(2) + F_trans(3))/(F_trans(1) + F_trans(2) + F_trans(3))
!         P_f = trans(2)*trans(3)/(trans(1)*trans(2)+trans(1)*trans(3)+trans(2)*trans(3))
         end if
         Max_j = nucleus(icomp)%F_barrier(ib)%Max_J
         aa = nucleus(icomp)%F_Barrier(ib)%barrier_damp(1)
         bb = nucleus(icomp)%F_Barrier(ib)%barrier_damp(2)
         cc = nucleus(icomp)%F_Barrier(ib)%barrier_damp(3)
!-----  Sum over discrete transition states
         do j = 1, nucleus(icomp)%F_Barrier(1)%num_discrete
            par = -1.0d0
            if(ip == 1) par = 1.0d0
            if(nucleus(icomp)%F_barrier(ib)%state_j(j) == xI .and.                        &
               nucleus(icomp)%F_barrier(ib)%state_pi(j) == par)then
               F_Barrier = nucleus(icomp)%F_Barrier(ib)%barrier
               F_Barrier = F_Barrier*aa*exp(-cc**2*(energy-bb)**2)
               if(Max_J > 0.0d0)then
                   b = 1.0d0/(Max_J*(Max_J+1.0d0))
                  if(xj_a <= Max_J)then
                     F_Barrier = F_Barrier*(1.0d0 - b*xj_a*(xj_a+1.0d0))
                  else
                     F_Barrier = 0.1d0
                  end if
               end if
               F_Barrier = max(F_Barrier,1.0d0)
               F_Barrier_hbw = nucleus(icomp)%F_Barrier(ib)%hbw
               Ef = nucleus(icomp)%F_barrier(ib)%state_e(j)
               tt = HW_trans(energy, Ef, F_Barrier, F_Barrier_hbw)
               xnu_c = xnu(tt,HF_den)
               trans = tt*P_f
               exponent = -0.5d0*xnu_c
               Product = Product + log(1.0d0 + 2.0d0*trans*x/(xnu_c*HF_den))*exponent
            end if
         end do

!-----   Integral over continuous energy bins for transition states
!-----   12/16/17 WEO ----   Had to replace to conform to level denbsity usage
!-----   in subroutine Fission_transmission. Namely we are no longer using 
!-----   subroutine Fiss_lev. All subroutines now use rhoe found in level-density.f
         ecut=nucleus(icomp)%F_Barrier(ib)%ecut
         E0 = nucleus(icomp)%F_barrier(ib)%level_param(15)
         T = nucleus(icomp)%F_barrier(ib)%level_param(14)
         if(E0 < 0.0d0)then
            E11 = T*log(1.0d0-exp(E0/T))
         else
            E11 = -5.0d0
         end if
         j_min = nint(E11/de)
         if(ecut > 0.0d0)then
            j = 0
         else
            j = j_min
         end if
         T_f = 0.0d0
         de=0.01
         j = 0
         converged = .false.
         do while (.not. converged)
            Ef = real(j,kind=8)*de + nucleus(icomp)%F_Barrier(ib)%ecut  
            call rhoe(Ef,nucleus(icomp)%F_barrier(ib)%level_param,                        &
                      nucleus(icomp)%F_barrier(ib)%vib_enh,                               &
                      nucleus(icomp)%F_barrier(ib)%rot_enh,                               &
                      nucleus(icomp)%A,rho_Fm,apu,sig2,K_vib,K_rot)
!
            mode = nucleus(icomp)%F_barrier(ib)%level_param(16)
            e1 = nucleus(icomp)%F_barrier(ib)%level_param(17)
            bbb = nucleus(icomp)%F_barrier(ib)%level_param(18)
            jfac = spin_fac(xI,sig2)
            pfac = parity_fac(Ef,xI,ip,mode,e1,bbb)
            rho = rho_FM*jfac*pfac
            F_Barrier = nucleus(icomp)%F_Barrier(ib)%barrier
            F_Barrier = F_Barrier*aa*exp(-cc**2*(energy-bb)**2)
            if(Max_J > 0.0d0)then
                b = 1.0d0/(Max_J*(Max_J+1.0d0))
               if(xj_a <= Max_J)then
                  F_Barrier = F_Barrier*(1.0d0 - b*xj_a*(xj_a+1.0d0))
               else
                  F_Barrier = 0.1d0
               end if
            end if
            F_Barrier = max(F_Barrier,1.0d0)
            F_Barrier_hbw = nucleus(icomp)%F_Barrier(ib)%hbw
            tt = HW_trans(energy, Ef, F_Barrier, F_Barrier_hbw)
            trans = tt*P_f
            N_eff = rho*de
            xnu_c = xnu(trans,HF_den)
            exponent = -0.5d0*xnu_c*N_eff
            Product = Product +                                                           &
                      log(1.0d0 + 2.0d0*trans*x/(xnu_c*HF_den))*exponent
            tt = tt*rho*de
            if(T_f > 0.0 .and. tt/T_f < 1.0d-6)converged = .true.
            if(tt < 1.0d-7)converged = .true.
            T_f = T_f + tt
            j = j + 1
         end do
      end if
   end if





   Product=exp(Product)
!write(6,*)'Product=',Product
!write(6,*)'sum_Tg=',sum_Tg
!write(6,*)'sum_n= ',sum_n
!write(6,*)'sum_p=',sum_p
!write(6,*)'HF_den ',HF_den
return
end subroutine Moldauer_product

