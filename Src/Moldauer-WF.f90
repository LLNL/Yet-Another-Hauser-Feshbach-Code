!
!*******************************************************************************
!
subroutine Moldauer_WF(icomp,                                       &
                       k_a,l_a,xj_a,istate_a,xI_a,trans_a,          &
                       k_b,l_b,xj_b,istate_b,xI_b,                  &
                       ip,xI,energy,                                &
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
!
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!        variable_kinds
!        options
!        nuclei
!        particles_def
!        constants 
!        Gauss_integration
!
!     Subroutines:
!
!        Moldauer_product
!
!     External functions:
!
!        real(kind=8) :: xnu
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
   use Gauss_integration
   implicit none
   integer(kind=4), intent(in) :: icomp
   integer(kind=4), intent(in) :: k_a, l_a, istate_a
   real(kind=8), intent(in) :: xj_a, xI_a, trans_a
   integer(kind=4), intent(in) :: k_b, l_b, istate_b
   real(kind=8), intent(in) :: xj_b, xI_b
   integer(kind=4), intent(in) :: ip
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
   real(kind=8) :: factor, elastic
   integer(kind=4) :: g_index

   real(kind=8) :: xstop

   integer(kind=4) :: ndata
!   parameter (ndata = 10)
   parameter (ndata = 4)
   real(kind=8) :: xx(ndata), yy(ndata), xxxx, yyxx
   real(kind=8) :: afit
   save xx
!   data xx/ 0.2d0,  0.4d0, 0.6d0, 0.8d0, 1.0d0, 1.2d0, 1.4d0, 1.6d0, 1.8d0, 2.0d0/
!   data xx/ 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0/
   data xx/ 1.0d0, 2.0d0, 3.0d0, 4.0d0/
!   data xx/0.5d0, 1.0d0, 1.5d0, 2.0d0, 3.0d0, 3.50d0, 4.0d0, 4.50d0, 5.0d0, 5.5d0/
!-------------     External Functions   
   real(kind=8) :: xnu
!-------------------------------------------------------------------------+
!------                                                                   +
!--------------------   Start subroutine                                  +
!------                                                                   +
!-------------------------------------------------------------------------+

   weight(0)=2.0d0/3.0d0
   weight(1)=4.0d0/3.0d0
   xstop = 100.0
   dx=0.001
   j_max=nucleus(icomp)%j_max
   nbin=nucleus(icomp)%nbin

   factor = 1.0d0
   elastic = 0.0d0
   if(k_a == k_b .and. istate_a == istate_b .and. l_a == l_b .and.       &
      (abs(xj_a - xj_b) < 1.0d-3) .and. (abs(xI_a - xI_b) < 1.0d-3) )    &
       elastic = 1.0d0
   xnu_a = xnu(trans_a,HF_den)
   factor = factor + elastic*2.0d0/xnu_a

   xxxx = 0.0d0
   yyxx = 0.0d0
!   write(40,*)'exp_gamma ',exp_gamma,' HF_den ',HF_den

   do ix = 1, ndata
      x = xx(ix)
      call Moldauer_product(icomp,                                       &
                            k_a,l_a,xj_a,istate_a,                       &
                            k_b,l_b,xj_b,istate_b,xI_b,                  &
                            ip,xI,energy,                                &
                            exp_gamma,HF_den,F_trans,x,CHnorm,Product_p)
      Product_g = exp(-exp_gamma*x)
      Product = Product_p*Product_g*factor
      yy(ix) = log(Product)
      xxxx = xxxx + xx(ix)**2
      yyxx = yyxx + yy(ix)*xx(ix)
   end do

   afit = abs(yyxx/xxxx)

!   write(40,'(3(1x,e15.7))')yyxx,xxxx,afit

!   flush(40)

!  write(6,*)'afit ',yyxx,xxxx,afit

   if(HF_den <= 2.0d0)then
      g_index = 1
   elseif(HF_den > 2.0d0 .and. HF_den <= 5.0d0)then
      g_index = 2
   elseif(HF_den > 5.0d0 .and. HF_den <= 10.0d0)then
      g_index = 3
   elseif(HF_den > 10.0d0)then
      g_index = 4
   end if


   WF = 0.0d0
   do ix = 1, gauss_laguerre(g_index)%num
      x = gauss_laguerre(g_index)%nodes(ix)/afit
      call Moldauer_product(icomp,                                       &
                            k_a,l_a,xj_a,istate_a,                       &
                            k_b,l_b,xj_b,istate_b,xI_b,                  &
                            ip,xI,energy,                                &
                            exp_gamma, HF_den,F_trans,x,CHnorm,Product_p)
      Product_g = exp(-exp_gamma*x)
      Product = Product_p*Product_g*factor
      WF = WF + Product/exp(-afit*x)*gauss_laguerre(g_index)%weights(ix)
!      write(40,'(1x,f10.5,1x,e15.7,1x,f10.5,4(1x,e15.7))')x, w_glag(ix), factor,  &
!           Product_g, Product_p, Product, Product/exp(-afit*x)
   end do
   WF = WF/afit


!   WF = 0.0d0
!   do ix = 1, n_glag
!      x = x_glag(ix)/afit
!      call Moldauer_product(icomp,                                       &
!                            k_a,l_a,xj_a,istate_a,                       &
!                            k_b,l_b,xj_b,istate_b,xI_b,                  &
!                            ip,xI,energy,                                &
!                            exp_gamma, HF_den,F_trans,x,CHnorm,Product_p)
!      Product_g = exp(-exp_gamma*x)
!      Product = Product_p*Product_g*factor
!      WF = WF + Product/exp(-afit*x)*w_glag(ix)
!      write(40,'(1x,f10.5,1x,e15.7,1x,f10.5,4(1x,e15.7))')x, w_glag(ix), factor,  &
!           Product_g, Product_p, Product, Product/exp(-afit*x)
!   end do
!   WF = WF/afit

!   write(6,*)exp_gamma, exp_gamma*HF_den
!
!   write(40,'(3(1x,i4),2(1x,f4.1))')k_a,l_a,istate_a,xj_a,xI 
!   write(40,'(3(1x,i4),2(1x,f4.1))')k_b,l_b,istate_b,xj_b,xI_b 

!   write(40,'(''HFden = '',1xf15.5,'' WF = '',1x,f10.6)')HF_den, WF

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
                            ip,xI,energy,                        &
                            exp_gamma,HF_den,F_trans,x,CHnorm,Product)
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
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!        variable_kinds
!        options
!        nuclei
!        particles_def
!        constants
!        Channel_info
!
!     Subroutines:
!
!        rho_J_par_e
!
!     External functions:
!
!       real(kind=8) :: xnu
!       real(kind=8) :: exp_1,exp_2
!       real(kind=8) :: HW_trans
!       logical :: real8_equal
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
   use options
   use nuclei
   use particles_def
   use constants
   use Channel_info
   implicit none
   integer(kind=4), intent(in) :: icomp
   integer(kind=4), intent(in) :: k_a,l_a,istate_a
   real(kind=8), intent(in) :: xj_a
   integer(kind=4), intent(in) :: k_b,l_b,istate_b
   real(kind=8), intent(in) :: xj_b,xI_b
   integer(kind=4), intent(in) :: ip
   real(kind=8), intent(in) :: xI
   real(kind=8), intent(in) :: energy, exp_gamma, HF_den, F_trans(4)
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
   real(kind=8) :: par
   integer(kind=4) :: cpar2
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
   real(kind=8) :: rho, tt, T_f, Ef
   integer(kind=4) :: ib
   real(kind=8) :: apu, sig2
   real(kind=8) :: aa, bb, cc
   real(kind=8) :: F_Barrier, F_Barrier_hbw
   real(kind=8) :: Max_J
   real(kind=8) :: b
   real(kind=8) :: E0, T, E11, ecut, emin
   real(kind=8) :: xZ_i, xA_i, xZ_part, xA_part
   real(kind=8) :: Coulomb_Barrier(6)
   real(kind=8) :: HFcheck, exp_g
   logical :: elastic
!-------------------------------------------------------------------------+
!------   External  Functions
   real(kind=8) :: xnu
   real(kind=8) :: exp_1,exp_2
   real(kind=8) :: HW_trans
   logical :: real8_equal
!-------------------------------------------------------------------------+
!------                                                                   +
!--------------------   Start subroutine                                  +
!------                                                                   +
!-------------------------------------------------------------------------+
   Coulomb_barrier(1:6) = 0.0d0
   if(Apply_Coulomb_Barrier)then
      do k_c = 1, 6
         xZ_part = real(particle(k_c)%Z,kind=8)
         xA_part = real(particle(k_c)%A,kind=8)
         xZ_i = real(nucleus(icomp)%Z,kind=8)
         xA_i = real(nucleus(icomp)%A,kind=8)
         Coulomb_Barrier(k_c) = 0.2d0*e_sq*(xZ_i-xZ_part)*xZ_part/               &
            (1.2d0*((xA_i-xA_part)**(1.0d0/3.0d0) + xA_part**(1.0d0/3.0d0)))
      end do
   end if


   sum_n = 0.0d0
   sum_p = 0.0d0
   sum_Tg = 0.0d0
   nptsx = 1000
   j_max = nucleus(icomp)%j_max
   nbin = nucleus(icomp)%nbin

   Product = 0.0d0
   HFcheck = 0.0d0
   num = 0
   num_bin = 0
   par = 2.0*real(ip)-1.0
   do if1 = 1, nucleus(icomp)%num_decay                       !  loop over nuclei in the decay chain
      i_c = nucleus(icomp)%decay_to(if1)
      k_c = nucleus(icomp)%decay_particle(if1)
      if(k_c == 0)cycle                                    ! particle n,p,d,t,h,a  ! Skip photons - treated woth exp(-exp_gamma*x)
      if(energy < nucleus(icomp)%sep_e(k_c))cycle           !   not eneough energy to decay - cycle out
      xj_c_max1 = real(nucleus(i_c)%j_max,kind=8) + nucleus(i_c)%jshift
      EM_c = 0
!--------------------------   particle decay to continuous level bins
      p_spin = particle(k_c)%spin
      isc = nint(2.0d0*p_spin)
      do n_c = 1, nucleus(i_c)%nbin                          !  loop over final excitation energies
         e_f = energy - nucleus(icomp)%sep_e(k_c) -                                       &
               nucleus(i_c)%e_grid(n_c)
         if(e_f - Coulomb_Barrier(k_c) <= 1.0d-6)exit
         do l_c = 0, particle(k_c)%lmax                      !  loop over l-partial wave
            cpar2 = nint(par*particle(k_c)%par*(-1.0d0)**l_c)      !  parity for channel c
            ip_c = (cpar2 + 1)/2                     !  parity index for channel c
            xj_c_min = real(l_c) - p_spin
            do iss = 0, isc
               xj_c = xj_c_min + real(iss,kind=8)
               if(xj_c < 0.0d0)cycle
!               trans = tco_interpolate(e_f,particle(k_c)%nume,                            &
!                                       particle(k_c)%e_grid,                              &
!                                       particle(k_c)%trans_read(1,iss,l_c))
               trans = particle(k_c)%trans_bin(iss,l_c,n_c)
               if(trans < trans_p_cut)cycle
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
                  HFcheck = HFcheck + trans_eff
                  exp_2 = 0.0d0
                  if(k_b == k_c .and. l_b == l_c .and. abs(xI_b - xI_c) < 1.0d-3 .and.    &
                     abs(xj_b - xj_c) < 1.0d-3 .and. istate_b == -n_c )exp_2 = 1.0d0
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
      if(all_discrete_states)num_discrete = nucleus(i_c)%num_discrete
      do n_c = 1, num_discrete
         E_c = energy - nucleus(icomp)%sep_e(k_c) - nucleus(i_c)%state(n_c)%energy
         if(E_c - Coulomb_Barrier(k_c) < 1.0d-6)exit
         xI_c = nucleus(i_c)%state(n_c)%spin
         xj_c_min = abs(xI_c - xI)
         xj_c_max = xI_c + xI
         num_j = nint(xj_c_max - xj_c_min)
         do j = 0, num_j
            xj_c = real(j,kind=8) + xj_c_min

            l_c_min = nint(abs(xj_c - p_spin))
            l_c_max = min(particle(k_c)%lmax, nint(xj_c + p_spin))
            cpar2 = nint(particle(k_c)%par*nucleus(i_c)%state(n_c)%parity)
            ip_c = (cpar2+1)/2
            if(ip == ip_c)then                            !   parities are the same, l=even
               if(iand(l_c_min,1) == 1)l_c_min = l_c_min + 1                  !   odd l_c_min, add 1 to make it even
               if(iand(l_c_max,1) == 1)l_c_max = l_c_max - 1                  !   odd l_c_max, subtract 1 to make it even
            else                                          !   parities are different, l=odd
               if(iand(l_c_min,1) == 0)l_c_min = l_c_min + 1                  !   even l_c_min, add 1 to make it even
               if(iand(l_c_max,1) == 0)l_c_max = l_c_max - 1                  !   even l_c_max, subtract 1 to make it even
            end if
            do l_c = l_c_min, l_c_max, 2
               xj_c_min1 = (real(l_c,kind=8) - p_spin)
               iss = nint(xj_c - xj_c_min1)
               if(iss < 0 .or. iss > nint(2*p_spin))cycle
!               trans = tco_interpolate(E_c,particle(k_c)%nume,                            &
!                                      particle(k_c)%e_grid,                               &
!                                      particle(k_c)%trans_read(1,iss,l_c))  
               trans = particle(k_c)%trans_discrete(iss,l_c,n_c)
               elastic = .false. 
               if(trans < trans_p_cut)cycle
               if(trans/CHnorm < prob_cut)cycle
               HFcheck = HFcheck + trans
               exp_1 = 0.0d0
               exp_2 = 0.0d0
               if(k_a == k_c .and. l_a == l_c .and.                                       &
                  abs(xj_a - xj_c) < 1.0d-3 .and. istate_a == n_c )exp_1 = 1.0d0
               if(k_b == k_c .and. l_b == l_c .and.                                       &
                  abs(xj_b - xj_c) < 1.0d-3 .and. istate_b == n_c )exp_2 = 1.0d0
               xnu_c = xnu(trans,HF_den)
               exponent = -0.5d0*xnu_c - exp_1 - exp_2
               Product = Product +                                                        &
                  log(1.0d0 + 2.0d0*trans*x/(xnu_c*HF_den))*exponent
!   if(k_a == 1)then
!   write(41,*)x,i_c,n_c,l_c,xj_c,exp_1,exp_2
!   write(41,*)x,trans,xnu_c,exponent,log(1.0d0 + 2.0d0*trans*x/(xnu_c*HF_den))*exponent,   &
!              Product 
!   end if
            end do
         end do
      end do
   end do

   exp_g = exp_gamma*HF_den
!   write(41,*)'Check ',HFcheck, exp_gamma*HF_den, HF_den, (HFcheck+exp_gamma*HF_den)/HF_den
!   write(6,*)'Check ',HFcheck, exp_gamma*HF_den, HF_den, (HFcheck+F_trans(4)+exp_gamma*HF_den)/HF_den

!-------------------------------------------------------------------*
!---------    Now fission, new approach                             *
!-------------------------------------------------------------------*
!------   The Moldauer width fluctuation is a bit complicated for fission, especially 
!------   for multiple barriers.
!------   As for two barriers, the total transmission coefficient is T = T_1*T_2/(T_1+T_2)
!------   where each T_i is integrated over the fission level density
!------   So, a single prescription is hard to do except for a single barrier.
!------   In general, one of the T_i is the larger, define T_L the larger, and T_s
!------   the smaller. Note that if T_L >> T_s, then T ~ T_s. One way to 
!------   proceed, which is similar to Hilaire, except they pick an arbitrary
!------   break up of T with Sum_k T_k. Here, each transmission coefficient T_i is
!------   summed over the transition states m, T_i = Sum_m T_i(m). Thus, one can break
!------   up the fission transmission coefficient as
!------   T = Sum_m T_s(m)*P_f, with
!------   P_f = T_L/(T_L+T_s). 
!------   One can treat each sum over transition states in an analagous manner
!------   in the Moldauer expression for each channel. This is similar to what Hilaire et al.
!------   does, but makes use of each term that gets summed in the transmission coefficient.
!------   Frankly, it isn't obvious what the best treatment for two or more barriers is. For 
!------   a single barrier, treating each transition state is analagous to the various channels
!------   for the other decays. For the most part, each transition state is a "channel" to fission.

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
         end if
         Max_j = nucleus(icomp)%F_barrier(ib)%Max_J
         aa = nucleus(icomp)%F_Barrier(ib)%barrier_damp(1)
         bb = nucleus(icomp)%F_Barrier(ib)%barrier_damp(2)
         cc = nucleus(icomp)%F_Barrier(ib)%barrier_damp(3)
!-----  Sum over discrete transition states
         do j = 1, nucleus(icomp)%F_Barrier(ib)%num_discrete
            par = -1.0d0
            if(ip == 1) par = 1.0d0
            if(real8_equal(nucleus(icomp)%F_barrier(ib)%state_j(j),xI) .and.                        &
               real8_equal(nucleus(icomp)%F_barrier(ib)%state_pi(j),par))then
               F_Barrier = nucleus(icomp)%F_Barrier(ib)%barrier
               F_Barrier = F_Barrier*aa*exp(-((energy-bb)/cc)**2)
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
               if(k_c == -1)exponent = exponent - 1.0d0
               Product = Product + log(1.0d0 + 2.0d0*trans*x/(xnu_c*HF_den))*exponent
            end if
         end do

!-----   Integral over continuous energy bins for transition states
!-----   12/16/17 WEO ----   Had to replace to conform to level denbsity usage
!-----   in subroutine Fission_transmission. Namely we are no longer using 
!-----   subroutine Fiss_lev. All subroutines now use rhoe found in level-density.f
         de = 0.01
         ecut = nucleus(icomp)%F_Barrier(ib)%ecut
         E0 = nucleus(icomp)%F_barrier(ib)%level_param(15)
         T = nucleus(icomp)%F_barrier(ib)%level_param(14)
         if(E0 < 0.0d0)then
            E11 = T*log(1.0d0-exp(E0/T))
         else
            E11 = -5.0d0
         end if
         E11 = 0.0d0

         T_f = 0.0d0

         emin = 0.0d0
         if(nucleus(icomp)%F_Barrier(ib)%num_discrete > 0)  &
            emin = nucleus(icomp)%F_Barrier(ib)%ecut

!----   k_c = -1 signals fission for this event

         converged = .false.
         Ef = emin - 0.5d0*de
         do while (.not. converged)
            Ef = Ef + de

            call rho_J_par_e(Ef,xI, ip, nucleus(icomp)%F_barrier(ib)%level_param,         &
                             nucleus(icomp)%F_barrier(ib)%vib_enh,                        &
                             nucleus(icomp)%F_barrier(ib)%rot_enh,                        &
                             nucleus(icomp)%A, rho , apu, sig2, K_vib, K_rot)

            F_Barrier = nucleus(icomp)%F_Barrier(ib)%barrier
            F_Barrier = F_Barrier*aa*exp(-((energy-bb)/cc)**2)
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
            if(k_c == -1)exponent = exponent - 1.0d0
            Product = Product +                                                           &
                      log(1.0d0 + 2.0d0*trans*x/(xnu_c*HF_den))*exponent
            tt = tt*rho*de
            if(T_f > 0.0 .and. tt/T_f < 1.0d-5)converged = .true.
            if(Ef > energy .and. tt < 1.0d-6)converged = .true.
            T_f = T_f + tt
         end do
      end if
   end if

   Product = exp(Product)
   return

end subroutine Moldauer_product

