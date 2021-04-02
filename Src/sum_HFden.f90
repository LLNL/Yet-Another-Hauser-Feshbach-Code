!
!*******************************************************************************
!
subroutine sum_HFden(icomp, xI_i, par, energy, HFden, exp_gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine computes and sums the Hauser-Feshbach denominator
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
   use Channel_info
   use particles_def
   use constants
   use nodeinfo
   implicit none
!-------------------------------------------------------------------
   integer(kind=4), intent(in) :: icomp
   real(kind=8), intent(in) :: xI_i
   integer(kind=4), intent(in) :: par
   real(kind=8), intent(in) :: energy
   real(kind=8), intent(out) :: HFden, exp_gamma

!-------------------------------------------------------------------
   integer(kind=4) :: if1
   integer(kind=4) :: i_f, k, n_f, ip_f, isp_f, is_f
   integer(kind=4) :: Ix_f_min, Ix_f_max, Ix_f
   integer(kind=4) :: l, lmin, lmax
   integer(kind=4) :: le_min
   integer(kind=4) :: lm_min
   integer(kind=4) :: ip_i
   integer(kind=4) :: j, num_j

   real(kind=8) :: xI_f_max1, xj_f_min, xj_f_max, xj_f
   real(kind=8) :: xI_f_min, xI_f_max, xI_f
   real(kind=8) :: xj_f_min1
   real(kind=8) :: e_f
   real(kind=8) :: p_spin
   real(kind=8) :: trans, N_eff, trans_eff
   real(kind=8) :: cpar2, par_f
   real(kind=8) :: e_gamma
   real(kind=8) :: F_trans(4)

   integer(kind=4) :: num_discrete

   real(kind=8) :: xZ_i, xA_i, xZ_part, xA_part
   real(kind=8) :: Coulomb_Barrier(6)

!--------------------------------------------------------------------
   real(kind=8) :: tco_interpolate
   real(kind=8) :: EL_trans
   real(kind=8) :: ML_trans
   
!-----------------------------------------------------------------------
   HFden = 0.0d0
   exp_gamma = 0.0d0
   ip_i = (par + 1)/2

   Coulomb_barrier(1:6) = 0.0d0
   if(Apply_Coulomb_Barrier)then
      do k = 1, 6
         xZ_part = real(particle(k)%Z,kind=8)
         xA_part = real(particle(k)%A,kind=8)
         xZ_i = real(nucleus(icomp)%Z,kind=8)
         xA_i = real(nucleus(icomp)%A,kind=8)
         Coulomb_Barrier(k) = 0.6d0*e_sq*(xZ_i-xZ_part)*xZ_part/               &
            (1.2d0*((xA_i-xA_part)**(1.0d0/3.0d0) + xA_part**(1.0d0/3.0d0)))
      end do
   end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Now loop over all nuclei this compound nucleus decays to  +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do if1 = 1, nucleus(icomp)%num_decay                                          !  loop over nuclei in the decay chain
      i_f = nucleus(icomp)%decay_to(if1)
      k = nucleus(icomp)%decay_particle(if1)
      if(energy < nucleus(icomp)%sep_e(k))cycle                                  !   not eneough energy to decay - cycle out
      xI_f_max1 = dfloat(nucleus(i_f)%j_max) + nucleus(i_f)%jshift
      if(k > 0)then                                                              ! k > 0 - particle n,p,d,t,h,a decay
!--------------------------   particle decay to continuous level bins
         p_spin = particle(k)%spin
         isp_f = nint(2.0d0*p_spin)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------     Start loop over final states, first loop over final        +
!------     energy bins                                                +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do n_f = 1, nucleus(i_f)%nbin                                          !  loop over final excitation energies
            e_f = energy - nucleus(icomp)%sep_e(k) -      &
                  nucleus(i_f)%e_grid(n_f)
            if(e_f - Coulomb_Barrier(k) <= 1.0d-6)exit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    Loop over orbital angular momenta of emitted particle       +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            do l = 0, particle(k)%lmax                                            !  loop over l-partial wave
               cpar2 = par*(-1)**l                                                !  Parity of nculeus and emitted part
               par_f = cpar2*particle(k)%par                                      !  Parity of final nucleus
               ip_f = nint((par_f + 1.0d0)/2.0d0)
               xj_f = real(l,kind=8) - p_spin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Now Loop over possible exit particle angular momenta j       +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               do is_f = 0, isp_f
                  xj_f = xj_f + dfloat(is_f)
                  if(xj_f < 0.0d0)cycle
                  trans = tco_interpolate(e_f,particle(k)%nume,                 &
                                          particle(k)%e_grid,                   &
                                          particle(k)%trans_read(1,is_f,l))
                  if(trans < trans_p_cut)cycle
                  xI_f_min = abs(xI_i - xj_f)
                  xI_f_max = xI_i + xj_f
                  Ix_f_min = max(nint(xI_f_min-nucleus(i_f)%jshift),0)                       !  min j-index
                  Ix_f_max = min(nint(xI_f_max-nucleus(i_f)%jshift),nucleus(i_f)%j_max)      !  max j-index                           
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Loop over angular momenta in the residual nucleus            +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  do Ix_f = Ix_f_min, Ix_f_max                                               !  loop over final j
                     xI_f = dfloat(Ix_f) + xI_f_min
                     N_eff = nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*          &
                             nucleus(i_f)%delta_e(n_f)
                     trans_eff = trans*N_eff
                     HFden = HFden + trans_eff
                  end do
               end do
            end do
         end do

!---------------------------   particle decay to discrete states
         num_discrete = nucleus(i_f)%ncut
         if(all_discrete_states)num_discrete = nucleus(i_f)%num_discrete
         do n_f = 1, num_discrete
            e_f = energy-nucleus(icomp)%sep_e(k) -          &
            nucleus(i_f)%state(n_f)%energy
            if(e_f - Coulomb_Barrier(k) <= 1.0d-6)exit
            xI_f = nucleus(i_f)%state(n_f)%spin
            Ix_f = nint(xI_f - nucleus(i_f)%jshift)
            xj_f_min = abs(xI_f - xI_i)
            xj_f_max = xI_f + xI_i
            num_j = nint(xj_f_max - xj_f_min)
            do j = 0, num_j
               xj_f = dfloat(j) + xj_f_min
               lmin = nint(abs(xj_f - p_spin))
               lmax = min(particle(k)%lmax, nint(xj_f + p_spin))
               cpar2 = particle(k)%par*nucleus(i_f)%state(n_f)%parity
               ip_f = nint((cpar2+1.0d0)/2.0d0)
               if(ip_i == ip_f)then                                 !   parities are the same, l=even
                  if(iand(lmin,1) == 1)lmin=lmin+1                  !   odd lmin, add 1 to make it even
                  if(iand(lmax,1) == 1)lmax=lmax-1                  !   odd lmax, subtract 1 to make it even
               else                                                 !   parities are different, l=odd
                  if(iand(lmin,1) == 0)lmin=lmin+1                  !   even lmin, add 1 to make it even
                  if(iand(lmax,1) == 0)lmax=lmax-1                  !   even lmax, subtract 1 to make it even
               end if
               do l = lmin,lmax,2
                  xj_f_min1 = (dfloat(l) - p_spin)
                  is_f = nint(xj_f - xj_f_min1)
                  if(is_f < 0 .or. is_f > nint(2*p_spin))cycle
                  trans = tco_interpolate(e_f,particle(k)%nume,             &
                                          particle(k)%e_grid,               &
                                          particle(k)%trans_read(1,is_f,l))

                  if(trans < trans_p_cut)cycle
                  HFden = HFden + trans
               end do
            end do
         end do

      else                                                                    !  photons
!---------------------------   gamma decay to continuous level bins
         do n_f = 1, nucleus(i_f)%nbin                                      !  loop over final excitation energies
            e_gamma = energy - nucleus(i_f)%e_grid(n_f)
            if(e_gamma <= 0.0d0)exit
            do l = 1,nucleus(i_f)%lmax_E                   !  loop over EL decays

               trans = EL_trans(i_f, l, e_gamma, energy)

               ip_f = iand((ip_i+l),1)                                                    !  parity of final state
               xI_f_min = abs(xI_i-dfloat(l))                                             !  min final spin
               xI_f_max = min(xI_f_max1,xI_i + dfloat(l))                                 !  max final spin
               Ix_f_min = max(nint(xI_f_min-nucleus(i_f)%jshift),0)                       !  min j-index
               Ix_f_max = min(nint(xI_f_max-nucleus(i_f)%jshift),nucleus(i_f)%j_max)      !  max j-index                           
               do Ix_f = Ix_f_min,Ix_f_max                                                !  loop over final j
                  xI_f = dfloat(Ix_f)+nucleus(i_f)%jshift
                  if(xI_i < 1.0d-5.and. xI_f <= 1.0d-5)cycle                              !  O -> 0 not allowed
                  trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*           &
                                    nucleus(i_f)%delta_e(n_f)
                  HFden = HFden + trans_eff
                  exp_gamma = exp_gamma + trans_eff
               end do
            end do
!---------------------------   Now Magnetic decay 
            do l = 1,nucleus(i_f)%lmax_M                                                 !  loop over ML decays
               trans = ML_trans(i_f, l, e_gamma)
               ip_f = iand((ip_i + l - 1),1)                                             !  parity of final state
               xI_f_min = abs(xI_i - dfloat(l))                                          !  min final spin
               xI_f_max = min(xI_f_max1, xI_i + dfloat(l))                               !  max final spin
               Ix_f_min = max(nint(xI_f_min-nucleus(i_f)%jshift),0)                      !  min j-index
               Ix_f_max = min(nint(xI_f_max-nucleus(i_f)%jshift),nucleus(i_f)%j_max)     !  max j-index                           
               do Ix_f = Ix_f_min, Ix_f_max                                              !  loop over final j
                  xI_f = dfloat(Ix_f) + nucleus(i_f)%jshift
                  if(xI_i < 1.0d-5.and. xI_f <= 1.0d-5)cycle                             !  O -> 0 not allowed
                  trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*         &
                                    nucleus(i_f)%delta_e(n_f)
                  HFden = HFden + trans_eff
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
            if(e_gamma <= 0.0d0)exit
            xI_f = nucleus(i_f)%state(n_f)%spin
            if(xI_i < 1.0d-5.and. xI_f <= 1.0d-5)cycle                          !  O -> 0 not allowed
            lmin = max(1,int(abs(xI_f - xI_i)))                                 !   can't have L=0
            ip_f = iabs(nint((nucleus(i_f)%state(n_f)%parity-1)/2.))
                     if(ip_i == ip_f)then                                       !  parity the same; even L for E, odd L for M
               if(iand(lmin,1) == 0)then
                  le_min = lmin
                  lm_min = lmin + 1
               else
                  le_min = lmin + 1
                  lm_min = lmin
               end if
            else                                                                !  parity the same; odd L for E, even L for M             
               if(iand(lmin,1) == 0)then
                  le_min = lmin + 1
                  lm_min = lmin
               else
                  le_min = lmin
                  lm_min = lmin + 1
               end if
            end if
            do l = le_min,nucleus(i_f)%lmax_E,2
               trans = EL_trans(i_f, l, e_gamma, energy)
               HFden = HFden + trans
               exp_gamma = exp_gamma + trans
            end do
            do l = lm_min,nucleus(i_f)%lmax_M,2
               trans = ML_trans(i_f, l, e_gamma)
               HFden = HFden + trans
               exp_gamma = exp_gamma + trans
            end do
         end do
      end if
   end do                        ! ----  Finish if1 loop
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Fission decay channel                      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(nucleus(icomp)%fission)then
      call Fission_transmission(icomp, energy, xI_i, ip_i, F_trans)
      HFden = HFden + F_trans(4)
   end if
   exp_gamma = exp_gamma/HFden

   return

end subroutine sum_HFden
