!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine HF_primary_decay_setup(e_in,iproj,itarget,icomp,                  &
                                  istate,energy,absorption_cs)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine sets up the decay of the first compound nucleus
!    First step only - Populates the initial compound nucleus
!    with Optical model transmission coefficients, remove
!    population due to pre-equilibrium emission, and finally
!    apply effects of width fluctuations
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        Channel_info
!        particles_def
!        constants
!        nodeinfo
!
!     Subroutines:
!
!        sum_HFden
!        Fission_transmission
!        Moldauer_WF
!        pack_data
!
!     External functions:
!
!        real(kind=8) :: tco_interpolate
!        real(kind=8) :: jhat
!        real(kind=8) :: EL_trans
!        real(kind=8) :: ML_trans
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
   use Channel_info
   use particles_def
   use constants
   use nodeinfo
   implicit none
!------------------------------------    Argument data
   real(kind=8), intent(in) :: e_in
   integer(kind=4), intent(in) :: iproj,itarget,istate
   integer(kind=4), intent(in) :: icomp
   real(kind=8), intent(in) :: energy
   real(kind=8), intent(in) :: absorption_cs
!------------------------------------    Internal Data
   integer(kind=4) :: i, j, k, l, ii
   integer(kind=4) :: j_max
   integer(kind=4) :: lmin, lmax, le_min, lm_min
   integer(kind=4) :: nbin
   real(kind=8) :: xI_i,xI_i_min,xI_i_max
   real(kind=8) :: xI_f,xI_f_min,xI_f_max,xI_f_max1
   integer(kind=4) :: Ix_i,Ix_i_min,Ix_i_max
   integer(kind=4) :: Ix_f,Ix_f_min,Ix_f_max
   real(kind=8) :: xj_f, xj_f_min, xj_f_max, xj_f_min1
   real(kind=8) :: xj_i
   integer(kind=4) :: ipt
   integer(kind=4) :: ip_i, ip_f
   integer(kind=4) :: if1, i_f, n_f
   integer(kind=4) :: ifi
   integer(kind=4) :: num_j
   real(kind=8) :: p_spin
   real(kind=8) :: e_f,e_gamma
   real(kind=8) :: trans
   real(kind=8) :: trans_eff
   real(kind=8) :: F_trans(4)
   real(kind=8) :: exp_gamma
   real(kind=8) :: WF,WF_def
   real(kind=8) :: N_eff
   integer(kind=4) :: num_discrete
   integer(kind=4) :: nnn
!-------------------------------------------------------------------------

   real(kind=8) :: sum_cs

   real(kind=8) :: xZ_i, xA_i, xZ_part, xA_part
   real(kind=8) :: Coulomb_Barrier(6)

   integer(kind=4) :: num
   integer(kind=4) :: itemp, idb

   real(kind=8) cs,cs_fac
   real(kind=8) :: mass_i, mass_t, mass_rel, e_rel
   real(kind=8) :: momentum, wave_number
   real(kind=8) :: spin_proj, spin_target
   real(kind=8) sp1,sp2
   integer(kind=4) :: isp, is_i, is_f
   integer(kind=4) :: isp_f_max, isp_max
   integer(kind=4) nume
   integer(kind=4) :: l_proj
   integer(kind=4) :: EM_proj, EM_k
   integer(kind=4) cpar, cpar2, par, par_f
   real(kind=8) :: tcoef
   real(kind=8) :: HFden, CHnorm, CHnorm2, prob_sum, prob_norm
   integer(kind=4) :: WF_calc
   real(kind=8) :: WF_avg, WF_min, WF_max
   real(kind=8) :: chan_cs
   integer(kind=4) :: WF_num
!   real(kind=8) :: xnu_a
!-------------------------------------------------------------------------+
!--------------------------   External Function declarations
   real(kind=8) :: tco_interpolate
   real(kind=8) :: jhat
   real(kind=8) :: EL_trans, ML_trans
!-rem   real(kind=8) :: xnu
!-------------------------------------------------------------------------+
!------                                                                   +
!--------------------   Start subroutine                                  +
!------                                                                   +
!-------------------------------------------------------------------------+

   nume = particle(iproj)%nume
   j_max = nucleus(icomp)%j_max
   nbin = nucleus(icomp)%nbin
   e_rel = e_in*real(nucleus(itarget)%A,kind=8)/                           &
           real(nucleus(itarget)%A+projectile%A,kind=8)
   spin_target = nucleus(itarget)%state(istate)%spin
   spin_proj = particle(iproj)%spin
   isp_max = nint(2.0d0*spin_proj)
   EM_proj = 0
   mass_i = particle(iproj)%Mass
   mass_t = nucleus(itarget)%Mass + nucleus(itarget)%state(istate)%energy
   mass_rel = mass_i*mass_t/(mass_t+mass_i)
   momentum = dsqrt(2.0d0*e_rel*mass_rel)
   wave_number = momentum/hbar_c
   cs = pi/wave_number**2*fmsq_eq_barn
   sp1 = abs(spin_target-spin_proj)
   sp2 = spin_target+spin_proj
   cpar = nint(particle(iproj)%par*nucleus(itarget)%state(istate)%parity)
   ipt = nint((nucleus(itarget)%state(istate)%parity + 1.0d0)/2.0d0)

   Coulomb_barrier(1:6) = 0.0d0
   if(Apply_Coulomb_Barrier)then
      do k = 1, 6
         xZ_part = real(particle(k)%Z,kind=8)
         xA_part = real(particle(k)%A,kind=8)
         xZ_i = real(nucleus(icomp)%Z,kind=8)
         xA_i = real(nucleus(icomp)%A,kind=8)
         Coulomb_Barrier(k) = 0.2d0*e_sq*(xZ_i-xZ_part)*xZ_part/               &
            (1.2d0*((xA_i-xA_part)**(1.0d0/3.0d0) + xA_part**(1.0d0/3.0d0)))
      end do
   end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Zero out Arrays used to cache transmission coefficients              +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do if1 = 1, nucleus(1)%num_decay                                      !  loop over nuclei in the decay chain
      i_f = nucleus(1)%decay_to(if1)
      k = nucleus(1)%decay_particle(if1)
      isp_f_max = nint(2.0d0*particle(k)%spin)
      if(k == 0)cycle
      particle(k)%trans_bin(0:isp_f_max,0:particle(k)%lmax,1:nucleus(i_f)%nbin) = 0.0d0
      particle(k)%trans_discrete(0:isp_f_max,0:particle(k)%lmax,1:nucleus(i_f)%num_discrete) = 0.0d0
   end do


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------      Loop over possible entrance channel spins                +
!-------      (projectile+target)                                      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   isp = nint(2.0d0*spin_proj)

   if(.not.allocated(Channel))                                                     &
      allocate(Channel(0:particle(iproj)%lmax,0:isp,0:nucleus(icomp)%j_max))

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   First check to see if HF Denominator > 70. If so, turn     +
!-------   width fluctuations off.                                    +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   WF_calc = 0

   if(print_me)then
         write(6,*)'****************************************************************************'
      if(WF_model == 0)then
         write(6,*)'*  Width fluctuation corrections are not being applied - option not active *'
         write(6,*)'*  If desired use option "wf_model 1" (default)                            *'
      else
         write(6,*)'*  Width fluctuation corrections are being applied                         *'
      end if
         write(6,*)'****************************************************************************'
    end if

   sum_cs = 0.0d0
   do l_proj = 0, particle(iproj)%lmax
      par = cpar*(-1)**l_proj
      ip_i = (par+1)/2
      xj_i = real(l_proj,kind=8) - spin_proj
      do is_i = 0, isp_max
         xj_i = xj_i + real(is_i,kind=8)
         if(xj_i < 0.0d0)cycle
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------      Loop over possible entrance orbital angiular momenta     +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         xI_i_min = abs(xj_i - spin_target)
         xI_i_max = xj_i + spin_target
         Ix_i_min = max(nint(xI_i_min-nucleus(icomp)%jshift),0)
         Ix_i_max = min(nint(xI_i_max-nucleus(icomp)%jshift),nucleus(icomp)%j_max)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------      Loop over Compound nucleus angular momenta               +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do Ix_i = Ix_i_min, Ix_i_max

            xI_i = real(Ix_i,kind=8) + nucleus(icomp)%jshift

            if(.not.allocated(Channel(l_proj,is_i,Ix_i)%Channel_decay))                         &
               allocate(Channel(l_proj,is_i,Ix_i)%Channel_decay(nucleus(icomp)%num_decay+1))
            if(.not.allocated(Channel(l_proj,is_i,Ix_i)%Channel_prob))                          &
               allocate(Channel(l_proj,is_i,Ix_i)%Channel_prob(1:nucleus(icomp)%num_decay+1))
            Channel(l_proj,is_i,Ix_i)%Channel_prob(1:nucleus(icomp)%num_decay+1) = 0.0d0
!-rem            if(.not.allocated(Channel(l_proj,is_i,Ix_i)%Channel_trans))                         &
!-rem               allocate(Channel(l_proj,is_i,Ix_i)%Channel_trans(1:nucleus(icomp)%num_decay+1))
            if(.not.allocated(Channel(l_proj,is_i,Ix_i)%decay_to))                              &
               allocate(Channel(l_proj,is_i,Ix_i)%decay_to(1:nucleus(icomp)%num_decay+1))
            if(.not.allocated(Channel(l_proj,is_i,Ix_i)%decay_particle))                        &
               allocate(Channel(l_proj,is_i,Ix_i)%decay_particle(1:nucleus(icomp)%num_decay+1))
!-rem            Channel(l_proj,is_i,Ix_i)%Channel_trans(1:nucleus(icomp)%num_decay+1) = 0.0d0
            Channel(l_proj,is_i,Ix_i)%Channel_decay(1:nucleus(icomp)%num_decay+1)%num_decay = 0 
            Channel(l_proj,is_i,Ix_i)%num_decay = 0
            Channel(l_proj,is_i,Ix_i)%decay_to(1:nucleus(icomp)%num_decay+1) = 0
            Channel(l_proj,is_i,Ix_i)%decay_particle(1:nucleus(icomp)%num_decay+1) = 0

            WF = 1.0d0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   This part is to just check if this entrance l-value can be    +
!-----   skipped because it is too small.                              +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            cs_fac = jhat(xI_i)/(jhat(spin_target)*jhat(spin_proj))
            tcoef = tco_interpolate(e_rel,particle(iproj)%nume,                           &
                                    particle(iproj)%e_grid,                               &
                                    particle(iproj)%trans_read(1,is_i,l_proj))

            chan_cs = tcoef*cs*cs_fac

            sum_cs = sum_cs + chan_cs
            if(chan_cs/absorption_cs <= 1.0d-7)cycle
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   First compute the Hauser-Feshbach denominator
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            call sum_HFden(icomp, xI_i, par, energy, HFden, exp_gamma)
!            Channel(l_proj,is_i,Ix_i)%Channel_HFden = HFden
!            Channel(l_proj,is_i,Ix_i)%Channel_exp_gamma = exp_gamma

!            WF_calc = 0
!            if(WF_model > 0 .and. HFden < 100.0d0)WF_calc = 1
!            if(WF_model > 0 )WF_calc = 1

            CHnorm2 = 0.0d0
            exp_gamma = 0.0d0
            HFden = 0.0d0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Now loop over all nuclei this compound nucleus decays     +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            do if1 = 1, nucleus(icomp)%num_decay                                      !  loop over nuclei in the decay chain
               prob_sum = 0.0d0
               i_f = nucleus(icomp)%decay_to(if1)
               k = nucleus(icomp)%decay_particle(if1)
               if(energy < nucleus(icomp)%sep_e(k))cycle                              !   not eneough energy to decay - cycle out
               xI_f_max1 = real(nucleus(i_f)%j_max,kind=8) + nucleus(i_f)%jshift
               if(k > 0)then                                                          ! k > 0 - particle n,p,d,t,h,a decay
                  EM_k = 0
!--------------------------   particle decay to continuous level bins
                  p_spin = particle(k)%spin
                  isp_f_max = nint(2.0d0*p_spin)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------     Start loop over final states, first loop over final        +
!------     energy bins                                                +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  do n_f = 1, nucleus(i_f)%nbin                        !  loop over final excitation energies
                     e_f = energy - nucleus(icomp)%sep_e(k)-                &
                           nucleus(i_f)%e_grid(n_f)
                     if(e_f - Coulomb_Barrier(k) <= 1.0d-6)exit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    Loop over orbital angular momenta of emitted particle       +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                     do l = 0, particle(k)%lmax                                       !  loop over l-partial wave
                        cpar2 = par*(-1)**l                                           !  Parity of nculeus and emitted part
                        par_f = cpar2*nint(particle(k)%par)                           !  Parity of final nucleus
                        ip_f = (par_f + 1)/2
                        xj_f_min = real(l,kind=8) - p_spin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Now Loop over possible exit particle angular momenta j       +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        do is_f = 0, isp_f_max
                           xj_f = xj_f_min + real(is_f,kind=8)
                           if(xj_f < 0.0d0)cycle
                           trans = tco_interpolate(e_f,particle(k)%nume,             &
                                                   particle(k)%e_grid,               &
                                                   particle(k)%trans_read(1,is_f,l))
                           particle(k)%trans_bin(is_f,l,n_f) = trans
                           if(trans < trans_p_cut)cycle
                           xI_f_min = abs(xI_i - xj_f)
                           xI_f_max = xI_i + xj_f
                           Ix_f_min = max(nint(xI_f_min-nucleus(i_f)%jshift),0)                    !  min j-index
                           Ix_f_max = min(nint(xI_f_max-nucleus(i_f)%jshift),nucleus(i_f)%j_max)   !  max j-index                           
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Loop over angular momenta in the residual nucleus            +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                           do Ix_f = Ix_f_min, Ix_f_max                                            !  loop over final j
                              xI_f = real(Ix_f,kind=8) + xI_f_min
                              N_eff = nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*          &
                                      nucleus(i_f)%delta_e(n_f)
                              trans_eff = trans*N_eff
                              CHnorm2 = CHnorm2 + trans_eff
                           end do
                        end do
                     end do
                  end do
!---------------------------   particle decay to discrete states
                  num_discrete = nucleus(i_f)%ncut
                  if(all_discrete_states)num_discrete = nucleus(i_f)%num_discrete
                  do n_f = 1, num_discrete
                     e_f = energy-nucleus(icomp)%sep_e(k)-                                &
                               nucleus(i_f)%state(n_f)%energy
                     if(e_f - Coulomb_Barrier(k) <= 1.0d-6)exit
                     xI_f = nucleus(i_f)%state(n_f)%spin
                     Ix_f = nint(xI_f - nucleus(i_f)%jshift)
                     xj_f_min = abs(xI_f - xI_i)
                     xj_f_max = xI_f + xI_i
                     num_j = nint(xj_f_max - xj_f_min)
                     do j = 0, num_j
                        xj_f = real(j,kind=8) + xj_f_min
                        lmin = nint(abs(xj_f - p_spin))
                        lmax = min(particle(k)%lmax, nint(xj_f + p_spin))
                        cpar2 = nint(particle(k)%par*nucleus(i_f)%state(n_f)%parity)
                        ip_f = (cpar2+1)/2
                        if(ip_i == ip_f)then                                               !   parities are the same, l=even
                           if(iand(lmin,1) == 1)lmin = lmin + 1                            !   odd lmin, add 1 to make it even
                           if(iand(lmax,1) == 1)lmax = lmax - 1                            !   odd lmax, subtract 1 to make it even
                        else                                                               !   parities are different, l=odd
                           if(iand(lmin,1) == 0)lmin = lmin + 1                            !   even lmin, add 1 to make it even
                           if(iand(lmax,1) == 0)lmax = lmax - 1                            !   even lmax, subtract 1 to make it even
                        end if
                        do l = lmin, lmax, 2
                           xj_f_min1 = (real(l,kind=8) - p_spin)
                           is_f = nint(xj_f - xj_f_min1)
                           if(is_f < 0 .or. is_f > nint(2*p_spin))cycle
                           trans = tco_interpolate(e_f,particle(k)%nume,                  &
                                                   particle(k)%e_grid,                    &
                                                   particle(k)%trans_read(1,is_f,l))
                           particle(k)%trans_discrete(is_f,l,n_f) = trans
                           if(trans < trans_p_cut)cycle
                           CHnorm2 = CHnorm2 + trans
                        end do
                     end do
                  end do
               else                                              !  photons
!---------------------------   gamma decay to continuous level bins
                  do n_f = 1, nucleus(i_f)%nbin                                     !  loop over final excitation energies
                     e_gamma = energy - nucleus(i_f)%e_grid(n_f)
                     if(e_gamma <= 1.0d-6)exit
                     EM_k = 1
                     do l = 1,nucleus(i_f)%lmax_E                                              !  loop over EL decays
                        trans = EL_trans(i_f, l, e_gamma, energy)
                        if(trans < trans_e_cut)cycle
                        ip_f = iand((ip_i + l),1)                                              !  parity of final state
                        xI_f_min = abs(xI_i - real(l,kind=8))                                  !  min final spin
                        xI_f_max = min(xI_f_max1,xI_i + real(l,kind=8))                        !  max final spin
                        Ix_f_min = max(nint(xI_f_min-nucleus(i_f)%jshift),0)                   !  min j-index
                        Ix_f_max = min(nint(xI_f_max-nucleus(i_f)%jshift),nucleus(i_f)%j_max)  !  max j-index                           
                        do Ix_f = Ix_f_min,Ix_f_max                                            !  loop over final j
                           xI_f = real(Ix_f,kind=8)+nucleus(i_f)%jshift
                           if(xI_i < 1.0d-5.and. xI_f <= 1.0d-5)cycle                          !  O -> 0 not allowed
                           trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*         &
                                             nucleus(i_f)%delta_e(n_f)
                           CHnorm2 = CHnorm2 + trans_eff
                           exp_gamma = exp_gamma + trans_eff
                        end do
                     end do
!---------------------------   Now Magnetic decay 
                     EM_k = 2
                     do l = 1,nucleus(i_f)%lmax_M                                              !  loop over ML decays
                        trans = ML_trans(i_f, l, e_gamma)
                        if(trans < trans_e_cut)cycle
                        ip_f = iand((ip_i + l - 1),1)                                          !  parity of final state
                        xI_f_min = abs(xI_i - real(l,kind=8))                                       !  min final spin
                        xI_f_max = min(xI_f_max1, xI_i + real(l,kind=8))                            !  max final spin
                        Ix_f_min = max(nint(xI_f_min-nucleus(i_f)%jshift),0)                   !  min j-index
                        Ix_f_max = min(nint(xI_f_max-nucleus(i_f)%jshift),nucleus(i_f)%j_max)  !  max j-index                           
                        do Ix_f = Ix_f_min, Ix_f_max                                           !  loop over final j
                           xI_f = real(Ix_f,kind=8) + nucleus(i_f)%jshift
                           if(xI_i < 1.0d-5.and. xI_f <= 1.0d-5)cycle                          !  O -> 0 not allowed
                           trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*         &
                                             nucleus(i_f)%delta_e(n_f)
                           CHnorm2 = CHnorm2 + trans_eff
                           exp_gamma = exp_gamma + trans_eff
                        end do
                     end do
                  end do
!---------------------------   gamma decay to discrete states
                  num_discrete = nucleus(i_f)%ncut
                  if(all_discrete_states)num_discrete = nucleus(i_f)%num_discrete
                  do n_f = 1, num_discrete
                     e_f = energy - nucleus(i_f)%state(n_f)%energy
                     e_gamma = e_f
                     if(e_gamma <= 1.0d-6)cycle
                     xI_f = nucleus(i_f)%state(n_f)%spin
                     if(xI_i < 1.0d-5.and. xI_f <= 1.0d-5)cycle                !  O -> 0 not allowed
                     lmin = max(1,int(abs(xI_f - xI_i)))                       !   can't have L=0
                     ip_f = iabs(nint((nucleus(i_f)%state(n_f)%parity+1)/2.))
                     if(ip_i == ip_f)then                                      !  parity the same; even L for E, odd L for M
                        if(iand(lmin,1) == 0)then
                           le_min = lmin
                           lm_min = lmin + 1
                        else
                           le_min = lmin + 1
                           lm_min = lmin
                        end if
                     else                                                     !  parity the same; odd L for E, even L for M             
                        if(iand(lmin,1) == 0)then
                           le_min = lmin + 1
                           lm_min = lmin
                        else
                           le_min = lmin
                           lm_min = lmin + 1
                        end if
                     end if
                     EM_k = 1
                     do l = le_min, nucleus(i_f)%lmax_E, 2
                        trans = EL_trans(i_f, l, e_gamma, energy)
                        if(trans < trans_e_cut)cycle
                        CHnorm2 = CHnorm2 + trans
                        exp_gamma = exp_gamma + trans
                     end do
                     EM_k = 2
                     do l = lm_min, nucleus(i_f)%lmax_M, 2
                        trans = ML_trans(i_f, l, e_gamma)
                        if(trans < trans_e_cut)cycle
                        CHnorm2 = CHnorm2 + trans
                        exp_gamma = exp_gamma + trans
                     end do
                  end do
               end if
            end do                        ! ----  Finish if1 loop
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Fission decay channel                      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            F_trans(1:4) = 0.0d0
            if(nucleus(icomp)%fission)then
               call Fission_transmission(icomp,energy,xI_i,ip_i,F_trans)
               CHnorm2 = CHnorm2 + F_trans(4)
            end if

            HFden = CHnorm2
            exp_gamma = exp_gamma/CHnorm2

            WF_calc = 0
            if(WF_model == 1 .and. HFden < mold_cutoff)WF_calc = 1
            if(WF_model == 1 .and. chan_cs/absorption_cs <= 1.0d-4)WF_calc = 0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------      Run again to count and to implement other cuts on prob   +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----  Control for calculating width fluctuation corrections.         +
!-----  If HFden > 50-100, WF ~ 1.0                                    +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            WF_def = 1.0d0

   WF_avg = 0.0d0
   WF_min = 1.0d6
   WF_max = -1.0d0
   WF_num = 0

!   write(6,*)'CHnorm2 = ',CHnorm2, sum_cs

            HFden = 0.0d0
            exp_gamma = 0.0d0
            do ii = 1, 2
               if(ii == 1)CHnorm = 0.0d0
               if(ii == 2)then
                  exp_gamma = exp_gamma/HFden
                  Channel(l_proj,is_i,Ix_i)%Channel_HFden = HFden
                  Channel(l_proj,is_i,Ix_i)%Channel_exp_gamma = exp_gamma
               end if
               ifi = 1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Loop over all nuclei this compound nucleus decays to      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               do if1 = 1,nucleus(icomp)%num_decay                                 !  loop over nuclei in the decay chain
                  num = 0
                  prob_sum = 0.0d0
                  if(ii == 2) prob_norm = Channel(l_proj,is_i,iX_i)%Channel_prob(ifi)
                  i_f = nucleus(icomp)%decay_to(if1)
                  k = nucleus(icomp)%decay_particle(if1)
                  if(energy < nucleus(icomp)%sep_e(k))cycle                        !   not eneough energy to decay - cycle out
                  xI_f_max1 = real(nucleus(i_f)%j_max,kind=8) + nucleus(i_f)%jshift
                  if(k > 0)then                                                    ! k > 0 - particle n,p,d,t,h,a decay
                     EM_k = 0
!--------------------------  particle decay to continuous level bins
                     p_spin = particle(k)%spin
                     isp_f_max = nint(2.0d0*p_spin)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------     Start loop over final states, first loop over final        +
!------     energy bins                                                +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                     do n_f = 1, nucleus(i_f)%nbin                                   !  loop over final excitation energies
                        e_f = energy-nucleus(icomp)%sep_e(k)-                     &
                              nucleus(i_f)%e_grid(n_f)
                        if(e_f - Coulomb_Barrier(k) <= 1.0d-6)exit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    Loop over orbital angular momenta of emitted particle       +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        do l = 0, particle(k)%lmax                                  !  loop over l-partial wave
                           cpar2 = par*(-1)**l                                      !  Parity of nculeus and emitted part
                           par_f = cpar2*nint(particle(k)%par)                      !  Parity of final nucleus
                           ip_f = (par_f + 1)/2
                           xj_f_min = real(l,kind=8) - p_spin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Now Loop over possible exit particle angular momenta j       +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                           do is_f = 0, isp_f_max
                              xj_f = xj_f_min + real(is_f,kind=8)
                              if(xj_f < 0.0d0)cycle
!                              trans = tco_interpolate(e_f,particle(k)%nume,             &
!                                                      particle(k)%e_grid,               &
!                                                      particle(k)%trans_read(1,is_f,l))
                              trans = particle(k)%trans_bin(is_f,l,n_f)
                              if(trans < trans_p_cut)cycle
                              xI_f_min = abs(xI_i - xj_f)
                              xI_f_max = xI_i + xj_f
                              Ix_f_min = max(nint(xI_f_min-nucleus(i_f)%jshift),0)                        !  min j-index
                              Ix_f_max = min(nint(xI_f_max-nucleus(i_f)%jshift),nucleus(i_f)%j_max)      !  max j-index                           
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Loop over angular momenta in the residual nucleus            +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                              do Ix_f = Ix_f_min, Ix_f_max                                     !  loop over final j
                                 xI_f = real(Ix_f,kind=8) + xI_f_min
                                 N_eff = nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*            &
                                         nucleus(i_f)%delta_e(n_f)
                                 trans_eff = trans*N_eff

                                 if(trans_eff/CHnorm2 < prob_cut)cycle
                                 WF = WF_def
                                 if(WF_model == 1 .and. WF_calc == 1 .and. ii == 2)then
                                    call  Moldauer_WF(icomp,                              &
                                                      iproj,l_proj,xj_i,istate,           &
                                                      spin_target,tcoef,                  &
                                                      k,l,xj_f,-n_f,xI_f,                 &
                                                      ip_i,xI_i,energy,                   &
                                                      exp_gamma,HFden,F_trans,CHnorm2,WF)
                                    WF_avg = WF_avg + WF
                                    WF_num = WF_num + 1
                                    if(WF < WF_min)WF_min = WF
                                    if(WF > WF_max)WF_max = WF
                                 end if
                                 num = num + 1
                                 if(ii == 1)then
                                    HFden = HFden + trans_eff
                                  else
                                    CHnorm = CHnorm + WF*trans_eff
                                    prob_sum = prob_sum + WF*trans_eff
!   write(81,'(''bin '',2(1x,i3,1x,f4.1,1x,f4.1),3(1x,e15.7))')l_proj,xj_i,xI_i,l,xj_f,xI_f, &
!             trans_eff/CHnorm2,trans_eff, prob_sum
                                    Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(num) = WF*trans_eff
                                    idb = 0
                                    call pack_data(Ix_f,ip_f,n_f,idb,l,is_f,itemp)
                                    Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_list(num) = itemp
                                 end if
                              end do
                           end do
                        end do
                     end do
!---------------------------   particle decay to discrete states
                     num_discrete = nucleus(i_f)%ncut
                     if(all_discrete_states)num_discrete = nucleus(i_f)%num_discrete
                     do n_f = 1, num_discrete
                        e_f = energy-nucleus(icomp)%sep_e(k)-                                 &
                              nucleus(i_f)%state(n_f)%energy
                        if(e_f - Coulomb_Barrier(k) <= 1.0d-6)exit
                        xI_f = nucleus(i_f)%state(n_f)%spin
                        Ix_f = nint(xI_f - nucleus(i_f)%jshift)
                        xj_f_min = abs(xI_f - xI_i)
                        xj_f_max = xI_f + xI_i
                        num_j = nint(xj_f_max - xj_f_min)
                        do j = 0, num_j
                           xj_f = real(j,kind=8) + xj_f_min
                           lmin = nint(abs(xj_f - p_spin))
                           lmax = min(particle(k)%lmax, nint(xj_f + p_spin))
                           cpar2 = nint(particle(k)%par*nucleus(i_f)%state(n_f)%parity)
                           ip_f = (cpar2+1)/2
                           if(ip_i == ip_f)then                                               !   parities are the same, l=even
                              if(iand(lmin,1) == 1)lmin = lmin + 1                            !   odd lmin, add 1 to make it even
                              if(iand(lmax,1) == 1)lmax = lmax - 1                            !   odd lmax, subtract 1 to make it even
                           else                                                               !   parities are different, l=odd
                              if(iand(lmin,1) == 0)lmin = lmin + 1                            !   even lmin, add 1 to make it even
                              if(iand(lmax,1) == 0)lmax = lmax - 1                            !   even lmax, subtract 1 to make it even
                           end if
                           do l = lmin, lmax, 2
                              xj_f_min1 = real(l,kind=8) - p_spin
                              is_f = nint(xj_f - xj_f_min1)
                              if(is_f < 0 .or. is_f > nint(2*p_spin))cycle
!                              trans = tco_interpolate(e_f,particle(k)%nume,                     &
!                                                      particle(k)%e_grid,                       &
!                                                      particle(k)%trans_read(1,is_f,l))  
                              trans = particle(k)%trans_discrete(is_f,l,n_f)
                              if(trans < trans_p_cut)cycle
                              if(trans/CHnorm2 < prob_cut)cycle
                              WF = WF_def
                              if(WF_model == 1 .and. WF_calc == 1 .and. ii == 2)then
                                 call  Moldauer_WF(icomp,                                       &
                                                   iproj,l_proj,xj_i,istate,spin_target,tcoef,  &
                                                   k,l,xj_f,n_f,xI_f,                           &
                                                   ip_i,xI_i,energy,                            &
                                                   exp_gamma,HFden,F_trans,CHnorm2,WF)
                                    WF_avg = WF_avg + WF
                                    WF_num = WF_num + 1
                                    if(WF < WF_min)WF_min = WF
                                    if(WF > WF_max)WF_max = WF
                              end if
                              num = num + 1
                              if(ii == 1)then
                                 HFden = HFden + trans
                              else
                                 CHnorm = CHnorm + WF*trans
                                 prob_sum = prob_sum + WF*trans
!   write(81,'(''discrete '',2(1x,i3,1x,f4.1,1x,f4.1),3(1x,e15.7))')l_proj,xj_i,xI_i,l,xj_f,xI_f, &
!             trans/CHnorm2,trans, prob_sum
                                 Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(num) = WF*trans
                                 idb = 1

                                 call pack_data(Ix_f, ip_f, n_f, idb, l, is_f, itemp)

                                 Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_list(num) = itemp
                              end if
                           end do
                        end do
                     end do

                  else                                              !  photons
!---------------------------   gamma decay to continuous level bins
                     EM_k = 1
                     WF = WF_def
                     n_f = 1
                     if(WF_model == 1 .and. WF_calc == 1 .and. ii == 2)then
                        l = 1
                        xj_f = 0.0
                        xI_f_min = xI_i
                        call Moldauer_WF(icomp,                                           &
                                         iproj,l_proj,xj_i,istate,spin_target,tcoef,      &
                                         0,l,xj_f,-n_f,xI_f_min,                          &
                                         ip_i,xI_i,energy,                                &
                                         exp_gamma,HFden,F_trans,CHnorm2,WF)
                        WF_avg = WF_avg + WF
                        WF_num = WF_num + 1
                        if(WF < WF_min)WF_min = WF
                        if(WF > WF_max)WF_max = WF
                     end if
                     do n_f = 1, nucleus(i_f)%nbin                                                   !  loop over final excitation energies
                        e_gamma = energy - nucleus(i_f)%e_grid(n_f)
                        if(e_gamma <= 1.0d-6)exit
                        EM_k = 1
                        do l = 1,nucleus(i_f)%lmax_E                                                !  loop over EL decays
                           trans = EL_trans(i_f, l, e_gamma, energy)
                           if(trans < trans_e_cut)cycle
                           ip_f = iand((ip_i + l),1)                                                !  parity of final state
                           xI_f_min = abs(xI_i-real(l,kind=8))                                           !  min final spin
                           xI_f_max = min(xI_f_max1,xI_i + real(l,kind=8))                               !  max final spin
                           Ix_f_min = max(nint(xI_f_min-nucleus(i_f)%jshift),0)                     !  min j-index
                           Ix_f_max = min(nint(xI_f_max-nucleus(i_f)%jshift),nucleus(i_f)%j_max)    !  max j-index                           
                           do Ix_f = Ix_f_min,Ix_f_max                                              !  loop over final j
                              xI_f = real(Ix_f,kind=8)+nucleus(i_f)%jshift
                              if(xI_i < 1.0d-5.and. xI_f <= 1.0d-5)cycle                            !  O -> 0 not allowed
                              trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*         &
                                                nucleus(i_f)%delta_e(n_f)
                              if(trans_eff/CHnorm2 < prob_cut)cycle
                              num = num + 1
                              if(ii == 1)then
                                 HFden = HFden + trans_eff
                                 exp_gamma = exp_gamma + trans_eff
                              else
                                 CHnorm = CHnorm + WF*trans_eff
                                 prob_sum = prob_sum + WF*trans_eff
                                 Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(num) = WF*trans_eff
                                 idb = 0
                                 is_f = 0
                                 call pack_data(Ix_f, ip_f, n_f, idb, l, is_f, itemp)
                                 Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_list(num) = itemp
                              end if
                           end do
                        end do
!---------------------------   Now Magnetic decay 
                        EM_k = 2
                        do l = 1,nucleus(i_f)%lmax_M                                                !  loop over ML decays
                           trans = ML_trans(i_f, l, e_gamma)
                           if(trans < trans_e_cut)cycle
                           ip_f = iand((ip_i + l - 1),1)                                            !  parity of final state
                           xI_f_min = abs(xI_i - real(l,kind=8))                                         !  min final spin
                           xI_f_max = min(xI_f_max1, xI_i + real(l,kind=8))                              !  max final spin
                           Ix_f_min = max(nint(xI_f_min-nucleus(i_f)%jshift),0)                     !  min j-index
                           Ix_f_max = min(nint(xI_f_max-nucleus(i_f)%jshift),nucleus(i_f)%j_max)    !  max j-index                           
                           do Ix_f = Ix_f_min, Ix_f_max                                             !  loop over final j
                              xI_f = real(Ix_f,kind=8) + nucleus(i_f)%jshift
                              if(xI_i < 1.0d-5.and. xI_f <= 1.0d-5)cycle                            !  O -> 0 not allowed
                              trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*         &
                                                nucleus(i_f)%delta_e(n_f)
                              if(trans_eff/CHnorm2 < prob_cut)cycle
                              num = num + 1
                              if(ii == 1)then
                                 HFden = HFden + trans_eff
                                 exp_gamma = exp_gamma + trans_eff
                              else
                                 CHnorm = CHnorm + WF*trans_eff
                                 prob_sum = prob_sum + WF*trans_eff
                                 Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(num) = WF*trans_eff
                                 idb = 0
                                 is_f = 1
                                 call pack_data(Ix_f, ip_f, n_f, idb, l, is_f, itemp)
                                 Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_list(num) = itemp
                              end if
                           end do
                        end do
                     end do
!---------------------------   gamma decay to discrete states
                     num_discrete = nucleus(i_f)%ncut
                     if(all_discrete_states)num_discrete = nucleus(i_f)%num_discrete
                     do n_f = 1, num_discrete
                        e_f = energy - nucleus(i_f)%state(n_f)%energy
                        e_gamma = e_f
                        if(e_gamma <= 1.0d-6)exit
                        xI_f = nucleus(i_f)%state(n_f)%spin
                        Ix_f = nint(xI_f - nucleus(i_f)%jshift)
                        if(xI_i < 1.0d-5.and. xI_f <= 1.0d-5)cycle                         !  O -> 0 not allowed
                        lmin = max(1,int(abs(xI_f - xI_i)))                                !   can't have L=0
                        ip_f = iabs(nint((nucleus(i_f)%state(n_f)%parity+1)/2.))
                        if(ip_i == ip_f)then                                               !  parity the same; even L for E, odd L for M
                           if(iand(lmin,1) == 0)then
                              le_min = lmin
                              lm_min = lmin + 1
                           else
                              le_min = lmin + 1
                              lm_min = lmin
                           end if
                        else                                                              !  parity the same; odd L for E, even L for M             
                           if(iand(lmin,1) == 0)then
                              le_min = lmin + 1
                              lm_min = lmin
                           else
                              le_min = lmin
                              lm_min = lmin + 1
                           end if
                        end if
                        EM_k = 1
                        do l = le_min,nucleus(i_f)%lmax_E,2
                           trans = EL_trans(i_f, l, e_gamma, energy)
                           if(trans < trans_e_cut)cycle
                           if(trans/CHnorm2 < prob_cut)cycle
                           num = num + 1
                           if(ii == 1)then
                              HFden = HFden + trans
                              exp_gamma = exp_gamma + trans
                           else
                              CHnorm = CHnorm + WF*trans
                              prob_sum = prob_sum + WF*trans
                              Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(num) = WF*trans
                              idb = 1
                              is_f = 0
                              call pack_data(Ix_f,ip_f,n_f,idb,l,is_f,itemp)
                              Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_list(num) = itemp
                           end if
                        end do
                        EM_k = 2
                        do l = lm_min,nucleus(i_f)%lmax_M,2
                           trans = ML_trans(i_f, l, e_gamma)
                           if(trans < trans_e_cut)cycle
                           if(trans/CHnorm2 < prob_cut)cycle
                           num = num + 1
                           if(ii == 1)then
                              HFden = HFden + trans
                              exp_gamma = exp_gamma + trans
                           else
                              CHnorm = CHnorm + WF*trans
                              prob_sum = prob_sum + WF*trans
                              Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(num) = WF*trans
                              idb = 1
                              is_f = 1
                              call pack_data(Ix_f,ip_f,n_f,idb,l,is_f,itemp)
                              Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_list(num) = itemp
                           end if
                        end do
                     end do
                  end if
                  Channel(l_proj,is_i,Ix_i)%Channel_prob(ifi) = prob_sum
!-rem                  Channel(l_proj,is_i,Ix_i)%Channel_trans(ifi) = prob_sum

                  if(num > 0)then
                     if(ii == 1)then
                        Channel(l_proj,is_i,Ix_i)%num_decay = Channel(l_proj,is_i,Ix_i)%num_decay + 1
                        Channel(l_proj,is_i,Ix_i)%decay_to(ifi) = i_f
                        Channel(l_proj,is_i,Ix_i)%decay_particle(ifi) = k
                        Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%num_decay = num
                        if(.not.allocated(Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob))   &
                           allocate(Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(1:num))
                        if(.not.allocated(Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_list))   &
                           allocate(Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_list(1:num))
                        Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(1:num) = 0.0d0
                        Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_list(1:num)= 0
                        Channel(l_proj,is_i,Ix_i)%Channel_prob(ifi) = prob_sum
                     else
                        do i = 1, num
                           Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(i) =               &
                                  Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(i)/prob_sum
                        end do
   
                        do i = 2, num
                           Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(i) =               &
                                  Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(i) +        &
                                  Channel(l_proj,is_i,Ix_i)%Channel_decay(ifi)%decay_prob(i-1)
                        end do
                     end if
                     ifi = ifi + 1
                  end if
               end do                        ! ----  Finish if1 loop
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Fission decay channel         ----------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               if(nucleus(icomp)%fission)then
                  call Fission_transmission(icomp,energy,xI_i,ip_i,F_trans)
                  WF = WF_def
                  n_f = 1
                  if(WF_model == 1 .and. WF_calc == 1 .and. ii == 2)then
                     l = 1
                     xj_f = 0.0
                     call Moldauer_WF(icomp,                                              &
                                      iproj,l_proj,xj_i,istate,spin_target,tcoef,         &
                                      -1,l,xj_f,-n_f,xI_f_min,                            &
                                      ip_i,xI_i,energy,                                   &
                                      exp_gamma,HFden,F_trans,CHnorm2,WF)
                     WF_avg = WF_avg + WF
                     WF_num = WF_num + 1
                     if(WF < WF_min)WF_min = WF
                     if(WF > WF_max)WF_max = WF
!    write(6,*)'4',WF
                  end if
                  if(F_trans(4)/CHnorm2 > prob_cut)then
                     if(ii == 1)then
                        Channel(l_proj,is_i,Ix_i)%num_decay =                             &
                           Channel(l_proj,is_i,Ix_i)%num_decay + 1
                        CHnorm = CHnorm + F_trans(4)*WF
                        HFden = HFden + F_trans(4)
                     else
!   write(81,'(''Fission '',(1x,i3,1x,f4.1,1x,f4.1),2(1x,e15.7))')l_proj,xj_i,xI_i,F_trans(4)/CHnorm2,F_trans(4)
                        nnn = Channel(l_proj,is_i,Ix_i)%num_decay
                        Channel(l_proj,is_i,Ix_i)%Channel_prob(nnn) = F_trans(4)*WF
!-rem                        Channel(l_proj,is_i,Ix_i)%Channel_trans(nnn) = F_trans(4)*WF
                        Channel(l_proj,is_i,Ix_i)%decay_to(nnn) = -1
                        Channel(l_proj,is_i,Ix_i)%decay_particle(nnn) = 7
                        Channel(l_proj,is_i,Ix_i)%Channel_decay(nnn)%num_decay = -1
                     end if
                  end if
               end if
!   if(WF_num > 0)then
!      WF_avg = WF_avg/real(WF_num,kind=8)
!      write(6,*)l_proj, xj_i,Ix_i
!      write(6,*)'Num computed = ',WF_num
!   end if
!---------------------------------------------------------------------------------------------------
            end do                           ! ----  Finish ii loop
!      if(WF_num > 0)WF_avg = WF_avg/real(WF_num,kind=8)
!      write(6,*)'HFden ', HFden,'num = ',WF_num,' WF_avg = ',WF_avg,' WF_min = ', WF_min,' WF_max = ', WF_max
!      write(60,*)HFden,WF_avg,WF_min, WF_max


            do if1 = 1, Channel(l_proj,is_i,Ix_i)%num_decay
               Channel(l_proj,is_i,Ix_i)%Channel_prob(if1) =                              &
                  Channel(l_proj,is_i,Ix_i)%Channel_prob(if1)/CHnorm
            end do
            do if1 = 2, Channel(l_proj,is_i,Ix_i)%num_decay
               Channel(l_proj,is_i,Ix_i)%Channel_prob(if1) =                              &
                   Channel(l_proj,is_i,Ix_i)%Channel_prob(if1) +                          &
                   Channel(l_proj,is_i,Ix_i)%Channel_prob(if1-1)
            end do
            do if1 = 1, Channel(l_proj,is_i,Ix_i)%num_decay
               Channel(l_proj,is_i,Ix_i)%Channel_prob(if1) =                              &
                   Channel(l_proj,is_i,Ix_i)%Channel_prob(if1)/                           &
                   Channel(l_proj,is_i,Ix_i)%Channel_prob(Channel(l_proj,is_i,Ix_i)%num_decay)
            end do     


         end do                           ! ----  Finish j loop
      end do                              ! ----  Finish is loop
   end do                                 ! ----  Finish l_proj loop


   return
end subroutine HF_primary_decay_setup
