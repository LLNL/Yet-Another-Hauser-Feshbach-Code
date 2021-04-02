!
!*******************************************************************************
!
subroutine HF_denominator(icomp)
!
!*******************************************************************************
!
!  Discussion:
!
!   This subroutine precalculates Hauser-Feshbach denominators.       
!   The primary purpose is to speed up the Hauser-Feshbach decay      
!   loops by having the important information pre-stored for          
!   each initial (E,J,par) bin for each compound nucleus.             
!   Program calculates HF-denominators for each (j,ip,n) bin          
!   Then creates an array of decay probabilites for each possible     
!   final nucleus from the starting compound nucleus.                 
!   To accomplish this, the HF-denominator loops are cycled through   
!   three times. The first is to get the overal magnitude of the      
!   denominator. The second is to attempt to cull paths within the    
!   a decay chain so that the more probable paths exist. The third    
!   finalizes the process. Note that all data is sotred in the        
!   derived types                                                     
!      nucleus(icomp)%bins(j,ip,n)%nuke_decay(if1)%decay_prob(num)    
!   and are dynamically allocated to use only the required memory.    
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
   use channel_info
   use nodeinfo
   implicit none
   include 'mpif.h'
   integer(kind=4), intent(in) :: icomp
!-------------------------------------------------------------------------+
   integer(kind=4) j, k, l, nnn, nnnn
   integer(kind=4) :: Ix_i_max, Ix_i
   integer(kind=4) :: Ix_f_min, Ix_f_max, Ix_f
   real(kind=8) :: xI_i, xI_i_shift
   real(kind=8) :: xI_f_min, xI_f_max, xI_f, xI_f_max1
   integer(kind=4) :: lmin, lmax, le_min, lm_min
   integer(kind=4) :: nbin
   real(kind=8) :: xj_f, xj_f_min1, xj_f_min, xj_f_max
   integer(kind=4) n
   integer(kind=4) :: ip
   integer(kind=4) :: if1, i_f, n_f, ns_f, ip_f
   integer(kind=4) :: ifi
   integer(kind=4) :: iss
   real(kind=8) :: par, par_i, par_f
   real(kind=8) :: hf_den, hf_den2, hf_denp, hf_prob(0:7)
   real(kind=8) :: p_spin
   real(kind=8) :: energy, e_f, e_gamma
   real(kind=8) :: cpar2
   real(kind=8) :: trans, trans_eff
   real(kind=8) :: F_trans(4)

   real(kind=8) :: xZ_i, xA_i, xZ_part, xA_part
   real(kind=8) :: Coulomb_Barrier(6)

!--------------------------------------------------------------------------

   integer(kind=4) :: num_discrete
   integer(kind=4) :: num_j
   integer(kind=4) :: ii
   integer(kind=4) :: isp
   integer(kind=4) :: num_prob
   integer(kind=8) :: num_tot
   integer(kind=4) :: itemp, idb
   integer(kind=4) :: ndecay
   integer(kind=4) :: num_Ix_ip
   integer(kind=4) :: iloop
   integer(kind=4) :: num_data
   integer(kind=4) :: my_proc

   real(kind=8) :: prob, prob_sum, prob_norm

!-------------------------------------------------------------------------+
!------     Function declarations
   real(kind=8) :: tco_interpolate
   real(kind=8) :: EL_trans
   real(kind=8) :: ML_trans
!-------------------------------------------------------------------------+
!------                                                                   +
!--------------------   Start subroutine                                  +
!------                                                                   +
!-------------------------------------------------------------------------+

   num_tot = 0

   Ix_i_max = nucleus(icomp)%j_max
   xI_i_shift = nucleus(icomp)%jshift
   nbin = nucleus(icomp)%nbin

!----   Coulomb barrier cutoff for low energy charged particles

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

   prob_norm = 0.0d0
!----------------------------------------------------------------------------------+
!-----   Compute decay probabilities for each continuous energy bin                +
!-----   Cycle trhough each bin and make a list of possible decays that satisfy    +
!-----   cuts imposed. The loop needs to be done three times.                      +
!-----                                                                             +
!-----   1. Compute overall size of the Hauser-Feshbach denominator for all decays +
!-----   2. Check if decay probability is > prob_cut and count number and allocate +
!-----   3. Set up decay arrays                                                    +
!-----                                                                             +
!----------------------------------------------------------------------------------+

   ndecay = nucleus(icomp)%num_decay
!   do n = 1, nbin
!      do ip = 0, 1
!         do Ix_i = 0, Ix_i_max
!            if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob))                  &
!                allocate(nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(ndecay+1))
!            if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%decay_to))                 &
!                allocate(nucleus(icomp)%bins(Ix_i,ip,n)%decay_to(ndecay+1))
!            if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle))           &
!                allocate(nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle(ndecay+1))
!            if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans))                 &
!                allocate(nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(ndecay+3))      ! added to see effect of each barrier
!            if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay))               &
!                allocate(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ndecay+1))
!            nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(1:ndecay+1) = 0.0d0
!            nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(1:ndecay+3) = 0.0d0
!            nucleus(icomp)%bins(Ix_i,ip,n)%num_decay = 0
!            nucleus(icomp)%bins(Ix_i,ip,n)%decay_to(1:ndecay+1) = 0
!            nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle(1:ndecay+1) = 0
!         end do
!      end do
!   end do


   num_Ix_ip = 2*(Ix_i_max+1)-1
   do n = 1, nbin                  
      energy = nucleus(icomp)%e_grid(n)
      do iloop = iproc, num_Ix_ip, nproc
         Ix_i = mod(iloop,Ix_i_max+1)
         ip = iloop/(Ix_i_max+1)
         par = 2.0d0*real(ip) - 1.0d0
         par_i = par

         xI_i = real(Ix_i,kind=8)+xI_i_shift
         hf_prob(0:7) = 0.0d0
         hf_den2 = 0.0d0

         if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob))                  &
             allocate(nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(ndecay+1))
         if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%decay_to))                 &
             allocate(nucleus(icomp)%bins(Ix_i,ip,n)%decay_to(ndecay+1))
         if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle))           &
             allocate(nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle(ndecay+1))
         if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans))                 &
             allocate(nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(ndecay+3))      ! added to see effect of each barrier
         if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay))               &
             allocate(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ndecay+1))
         nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(1:ndecay+1) = 0.0d0
         nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(1:ndecay+3) = 0.0d0
         nucleus(icomp)%bins(Ix_i,ip,n)%num_decay = 0
         nucleus(icomp)%bins(Ix_i,ip,n)%decay_to(1:ndecay+1) = 0
         nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle(1:ndecay+1) = 0

         do if1 = 1, nucleus(icomp)%num_decay                                               !  loop over nuclei in the decay chain
            prob_sum = 0.0d0
            i_f = nucleus(icomp)%decay_to(if1)
            k = nucleus(icomp)%decay_particle(if1)
            if(energy < nucleus(icomp)%sep_e(k))cycle                                       !   not eneough energy to decay - cycle out
            xI_f_max1 = real(nucleus(i_f)%j_max,kind=8) + nucleus(i_f)%jshift
            if(k > 0)then                                                                   ! particle n,p,d,t,h,a
!--------------------------   particle decay
               p_spin = particle(k)%spin
               isp = nint(2.0d0*p_spin)
!
!---------------    Particle decay to discrete states
!
!               do ns_f = 1, nucleus(i_f)%num_discrete
               num_discrete = nucleus(i_f)%ncut
               if(all_discrete_states)num_discrete = nucleus(i_f)%num_discrete
               do ns_f = 1, num_discrete
                  e_f = energy - nucleus(icomp)%sep_e(k) -                             &
                             nucleus(i_f)%state(ns_f)%energy
                  if(e_f - Coulomb_Barrier(k) <= 1.0d-6)exit
                  xI_f = nucleus(i_f)%state(ns_f)%spin
                  Ix_f = nint(xI_f - nucleus(i_f)%jshift)
                  xj_f_min = abs(xI_f - xI_i)
                  xj_f_max = xI_f + xI_i
                  num_j = nint(xj_f_max - xj_f_min)
                  do j = 0, num_j
                     xj_f = real(j,kind=8) + xj_f_min
                     lmin = nint(abs(xj_f - p_spin))
                     lmax = min(particle(k)%lmax, nint(xj_f + p_spin))
                     cpar2 = particle(k)%par*nucleus(i_f)%state(ns_f)%parity
                     ip_f = nint((cpar2+1.0d0)/2.0d0)
                     if(ip == ip_f)then                                                    !   parities are the same, l=even
                        if(iand(lmin,1) == 1)lmin = lmin + 1                               !   odd lmin, add 1 to make it even
                        if(iand(lmax,1) == 1)lmax = lmax - 1                               !   odd lmax, subtract 1 to make it even
                     else                                                                  !   parities are different, l=odd
                        if(iand(lmin,1) == 0)lmin = lmin + 1                               !   even lmin, add 1 to make it even
                        if(iand(lmax,1) == 0)lmax = lmax - 1                               !   even lmax, subtract 1 to make it even
                     end if
                     do l = lmin,lmax,2
                        xj_f_min1 = real(l,kind=8) - p_spin
                        iss = nint(xj_f - xj_f_min1)
                        if(iss < 0 .or. iss > nint(2*p_spin))cycle
                        trans = tco_interpolate(e_f,particle(k)%nume,                    &
                                                particle(k)%e_grid,                      &
                                                particle(k)%trans_read(1,iss,l))  
                        if(trans < trans_p_cut)cycle
                        hf_den2 = hf_den2 + trans
                        prob_sum = prob_sum + trans
                      end do
                  end do
               end do

               do n_f = 1, nucleus(i_f)%nbin                                               !  loop over final excitation energies
                  e_f = energy - nucleus(icomp)%sep_e(k) - nucleus(i_f)%e_grid(n_f)
                  if(e_f - Coulomb_Barrier(k) <= 1.0d-5)exit
                  do l = 0, particle(k)%lmax                                               !  loop over l-partial wave
                     cpar2 = par*(-1.0d0)**l                                               !  Parity of nucleus and emitted part
                     par_f = cpar2*particle(k)%par                                         !  Parity of final nucleus
                     ip_f = nint((par_f + 1.0)/2.0)
                     xj_f = real(l,kind=8) - p_spin    
                     do iss = 0, isp                                                       !  loop over particle spins
                        xj_f = xj_f + real(iss,kind=8)
                        if(xj_f < 0.0d0)cycle
                        trans = tco_interpolate(e_f,particle(k)%nume,                    &
                                                particle(k)%e_grid,                      &
                                                particle(k)%trans_read(1,iss,l))
                        if(trans < trans_p_cut)cycle
                        xI_f_min = abs(xj_f - xI_i)
                        xI_f_max = xj_f + xI_i
                        Ix_f_min = nint(xI_f_min - nucleus(i_f)%jshift)                    !  min j-index
                        Ix_f_max = nint(xI_f_max - nucleus(i_f)%jshift)                    !  max j-index
                        Ix_f_max = min(Ix_f_max, nucleus(i_f)%j_max)
                        do Ix_f = Ix_f_min, Ix_f_max                                       !  loop over final j
                           xI_f = real(Ix_f,kind=8) + nucleus(i_f)%jshift
                           trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*       &
                                             nucleus(i_f)%delta_e(n_f)
                           hf_den2 = hf_den2 + trans_eff
                           prob_sum = prob_sum + trans_eff
                        end do
                     end do
                  end do
               end do
            else                                                                           !  photons
!---------------------------   Gamma decay
!-------   Start with discrete states below ecut
!               do ns_f = 1, nucleus(i_f)%num_discrete
               num_discrete = nucleus(i_f)%ncut
               if(all_discrete_states)num_discrete = nucleus(i_f)%num_discrete
               do ns_f = 1, num_discrete
                  e_f = energy - nucleus(i_f)%state(ns_f)%energy
                  e_gamma = e_f
                  if(e_gamma <= 1.0d-6)cycle
                  xI_f = nucleus(i_f)%state(ns_f)%spin
                  if(xI_i < 1.0d-3 .and. xI_f < 1.0d-3)cycle                               !  O -> 0 not allowed
                  lmin = max(1, nint(abs(xI_f - xI_i)))                                    !   can't have L=0
                  ip_f = iabs(nint((nucleus(i_f)%state(ns_f)%parity+1.)/2.))
                  if(ip == ip_f)then                                                       !  parity the same even L for E odd L for M
                     if(iand(lmin,1) == 0)then
                        le_min = lmin
                        lm_min = lmin + 1
                     else
                        le_min = lmin + 1
                        lm_min = lmin
                     end if
                  else                                                                     !  parity the same odd L for E even L for M             
                     if(iand(lmin,1) == 0)then
                        le_min = lmin + 1
                        lm_min = lmin
                     else
                        le_min = lmin
                        lm_min = lmin + 1
                     end if
                  end if
                  do l = le_min, nucleus(i_f)%lmax_E, 2
                     trans = EL_trans(i_f, l, e_gamma, energy)
                     if(trans < trans_e_cut)cycle
                     hf_den2 = hf_den2 + trans
                     prob_sum = prob_sum + trans
                  end do
                  do l = lm_min, nucleus(i_f)%lmax_M, 2
                     trans = ML_trans(i_f, l, e_gamma)
                     if(trans < trans_e_cut)cycle
                     hf_den2 = hf_den2 + trans
                     prob_sum = prob_sum + trans
                  end do
               end do
!--------------    Now continuous bins
               do n_f = 1, n                                                            !  loop over final excitation energies
                  e_gamma = energy - nucleus(i_f)%e_grid(n_f)
                  if(e_gamma <= 1.0d-6)cycle
                  e_gamma = max(e_gamma,nucleus(i_f)%delta_e(n_f)/10.0)
!---------------------------   Start with Electric decay 
                  do l = 1, nucleus(i_f)%lmax_E                                            !  loop over EL decays
                     trans = EL_trans(i_f, l, e_gamma, energy)
                     if(trans < trans_e_cut)cycle
                     ip_f = iand((ip+l),1)                                                 !  parity of final state
                     xI_f_min = abs(xI_i-real(l,kind=8))                                   !  min final spin
                     xI_f_max = min(xI_f_max1,xI_i+real(l,kind=8))                         !  max final spin
                     Ix_f_min = nint(xI_f_min - nucleus(i_f)%jshift)                       !  min j-index
                     Ix_f_max = nint(xI_f_max - nucleus(i_f)%jshift)                       !  max j-index
                     do Ix_f = Ix_f_min, Ix_f_max                                          !  loop over final j
                        xI_f = real(Ix_f,kind=8) + nucleus(i_f)%jshift
                        if(xI_i < 1.0d-3 .and. xI_f < 1.0d-3)cycle                         !  O -> 0 not allowed
                        trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*       &
                                          nucleus(i_f)%delta_e(n_f)
                        hf_den2 = hf_den2 + trans_eff
                        prob_sum = prob_sum + trans_eff
                     end do
                  end do
!---------------------------   Now Magnetic decay 
                  do l = 1, nucleus(i_f)%lmax_M                                            !  loop over ML decays
                     trans = ML_trans(i_f, l, e_gamma)  
                     if(trans < trans_e_cut)cycle
                     ip_f = iand((ip+l+1),1)                                               !  parity of final state
                     xI_f_min = abs(xI_i-real(l,kind=8))                                   !  min final spin
                     xI_f_max = min(xI_f_max1,xI_i+real(l,kind=8))                         !  max final spin
                     Ix_f_min = nint(xI_f_min-nucleus(i_f)%jshift)                         !  min j-index
                     Ix_f_max = nint(xI_f_max-nucleus(i_f)%jshift)                         !  max j-index
                     do Ix_f = Ix_f_min, Ix_f_max                                          !  loop over final j
                        xI_f = real(Ix_f,kind=8) + nucleus(i_f)%jshift
                        if(abs(xI_i) < 1.0d-3.and. abs(xI_f) < 1.0d-3)cycle                !  O -> 0 not allowed
                        trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*       &
                                          nucleus(i_f)%delta_e(n_f)
                        hf_den2 = hf_den2 + trans_eff
                        prob_sum = prob_sum + trans_eff
                     end do
                  end do
               end do
            end if
         end do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Decay probabilities to each channel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(nucleus(icomp)%fission)then
            call Fission_transmission(icomp,energy,xI_i,ip,F_trans)
            hf_den2 = hf_den2 + F_trans(4)
         end if

         if(hf_den2 <= min(trans_p_cut,trans_e_cut))then                                   !   cannot decay this bin
            nucleus(icomp)%bins(Ix_i,ip,n)%HF_den = 0.0d0
            do if1 = 1, nucleus(icomp)%num_decay
               nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(if1) = 0.0d0
               nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(if1)%num_decay = 0
            end do
            goto 10101
         end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------  Do again, and reduce possible paths if need be     ----------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(1:nucleus(icomp)%num_decay+1) = 0.0d0
         nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(1:nucleus(icomp)%num_decay)%num_decay = 0
         do ii = 1, 2                                        ! do twice, the first is to count
            if(ii == 1)hf_den = 0.0d0
            ifi = 1
            do if1 = 1,nucleus(icomp)%num_decay                                            !  loop over nuclei in the decay chain
               hf_denp = 0.0d0
               num_prob = 0
               prob_sum = 0.0d0
               if(ii == 2) prob_norm = nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(ifi)
               i_f = nucleus(icomp)%decay_to(if1)
               k = nucleus(icomp)%decay_particle(if1)
               if(energy < nucleus(icomp)%sep_e(k))cycle                                   !   not eneough energy to decay - cycle out
               xI_f_max1 = real(nucleus(i_f)%j_max,kind=8) + nucleus(i_f)%jshift
               if(k > 0)then                                                               ! particle n,p,d,t,h,a
!--------------------------   particle decay to continuous level bins
                  p_spin=particle(k)%spin
                  isp = nint(2.0d0*p_spin)
!
!------------   Start with discrete states below ecut
!
!                  do ns_f = 1, nucleus(i_f)%num_discrete
                  num_discrete = nucleus(i_f)%ncut
                  if(all_discrete_states)num_discrete = nucleus(i_f)%num_discrete
                  do ns_f = 1, num_discrete
                     e_f = energy - nucleus(icomp)%sep_e(k)-           &
                                    nucleus(i_f)%state(ns_f)%energy
                     if(e_f - Coulomb_Barrier(k) <= 1.0d-6)exit
                     xI_f = nucleus(i_f)%state(ns_f)%spin
                     Ix_f = nint(xI_f - nucleus(i_f)%jshift)
                     xj_f_min = abs(xI_f - xI_i)
                     xj_f_max = xI_f + xI_i
                     num_j = nint(xj_f_max - xj_f_min)
                     do j = 0, num_j
                        xj_f = real(j,kind=8) + xj_f_min
                        lmin = nint(abs(xj_f - p_spin))
                        lmax = min(particle(k)%lmax, nint(xj_f + p_spin))
                        cpar2 = particle(k)%par*nucleus(i_f)%state(ns_f)%parity
                        ip_f = nint((cpar2+1.0d0)/2.0d0)
                        if(ip == ip_f)then                                                 !   parities are the same, l=even
                           if(iand(lmin,1) == 1)lmin = lmin + 1                            !   odd lmin, add 1 to make it even
                           if(iand(lmax,1) == 1)lmax = lmax - 1                            !   odd lmax, subtract 1 to make it even
                        else                                                               !   parities are different, l=odd
                           if(iand(lmin,1) == 0)lmin = lmin + 1                            !   even lmin, add 1 to make it even
                           if(iand(lmax,1) == 0)lmax = lmax - 1                            !   even lmax, subtract 1 to make it even
                        end if
                        do l = lmin, lmax, 2
                           xj_f_min1 = real(l,kind=8) - p_spin
                           iss = nint(xj_f - xj_f_min1)
                           if(iss < 0 .or. iss > nint(2*p_spin))cycle
                           trans = tco_interpolate(e_f,particle(k)%nume,           &
                                                   particle(k)%e_grid,             &
                                                   particle(k)%trans_read(1,iss,l))  
                           if(trans < trans_p_cut)cycle
                           prob = trans/hf_den2
                           if(prob <= prob_cut) cycle
                           num_prob = num_prob + 1
                           prob_sum = prob_sum + trans
                           if(ii == 1)then
                              hf_den = hf_den + trans
                           else
                              nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob(num_prob) = &
                                  prob_sum/prob_norm
                              idb = 1
                              call pack_data(Ix_f, ip_f, ns_f, idb, l, iss, itemp)
                              nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list(num_prob) = itemp
                           end if
                        end do
                     end do
                  end do

                  do n_f = 1, nucleus(i_f)%nbin                                             !  loop over final excitation energies
                     e_f = energy-nucleus(icomp)%sep_e(k) - nucleus(i_f)%e_grid(n_f)
                     if(e_f - Coulomb_Barrier(k) <= 1.0d-6)exit
                     do l = 0, particle(k)%lmax                                            !  loop over l-partial wave
                        xj_f = real(l,kind=8) - p_spin
                        cpar2 = par*(-1.)**l                                               !  Parity of nucleus and emitted part
                        par_f = cpar2*particle(k)%par                                      !  Parity of final nucleus
                        ip_f = nint((par_f + 1.0)/2.0)
                        do iss = 0, isp                                                    !  loop over particle spins
                           xj_f = xj_f + real(iss,kind=8)
                           if(xj_f < 0.0d0)cycle
                           trans = tco_interpolate(e_f,particle(k)%nume,                 &
                                                   particle(k)%e_grid,                   &
                                                   particle(k)%trans_read(1,iss,l))  
                           if(trans < trans_p_cut)cycle
                           xI_f_min = abs(xj_f - xI_i)
                           xI_f_max = xj_f + xI_i
                           Ix_f_min = nint(xI_f_min - nucleus(i_f)%jshift)                 !  min j-index
                           Ix_f_max = nint(xI_f_max - nucleus(i_f)%jshift)                 !  max j-index
                           Ix_f_max = min(Ix_f_max,nucleus(i_f)%j_max)
                           do Ix_f = Ix_f_min, Ix_f_max                                    !  loop over final j
                              xI_f = real(Ix_f,kind=8) + nucleus(i_f)%jshift
                              trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*    &
                                                nucleus(i_f)%delta_e(n_f)

                              prob = trans_eff/hf_den2
                              if(prob <= prob_cut)cycle
                              num_prob = num_prob + 1
                              prob_sum = prob_sum + trans_eff
                              if(ii == 1)then
                                 hf_den = hf_den + trans_eff
                              else
                                 nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob(num_prob) = prob_sum/prob_norm
                                 idb = 0
                                 call pack_data(Ix_f,ip_f,n_f,idb,l,iss,itemp)
                                 nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list(num_prob) = itemp
                              end if
                           end do
                        end do
                     end do
                  end do
               else                                                                        !  photons
!---------------------------   Gamma decay
!
!----------   First check states below ecut
!
!                  do ns_f = 1, nucleus(i_f)%num_discrete
                  num_discrete = nucleus(i_f)%ncut
                  if(all_discrete_states)num_discrete = nucleus(i_f)%num_discrete
                  do ns_f = 1, num_discrete
                     e_f = energy - nucleus(i_f)%state(ns_f)%energy
                     e_gamma = e_f
                     if(e_gamma <= 1.0d-6)cycle
                     xI_f = nucleus(i_f)%state(ns_f)%spin
                     if(xI_i < 1.0d-3 .and. xI_f < 1.0d-3)cycle                            !  O -> 0 not allowed
                     Ix_f = nint(xI_f - nucleus(i_f)%jshift)
                     lmin = max(1, nint(abs(xI_f-xI_i)))                                   !   can't have L=0
                     ip_f=iabs(nint((nucleus(i_f)%state(ns_f)%parity+1.)/2.))
                     if(ip == ip_f)then                                                    !  parity the same even L for E odd L for M
                        if(iand(lmin,1) == 0)then
                           le_min = lmin
                           lm_min = lmin + 1
                        else
                           le_min = lmin + 1
                           lm_min = lmin
                        end if
                     else                                                                  !  parity the same odd L for E even L for M             
                        if(iand(lmin,1) == 0)then
                           le_min = lmin + 1
                           lm_min = lmin
                        else
                           le_min = lmin
                           lm_min = lmin + 1
                        end if
                     end if
                     do l = le_min, nucleus(i_f)%lmax_E, 2
                        trans = EL_trans(i_f, l, e_gamma, energy)
                        if(trans < trans_e_cut)cycle
                        prob = trans/hf_den2
                        if(prob <= prob_cut) cycle
                        num_prob = num_prob + 1
                        prob_sum = prob_sum + trans
                        if(ii == 1)then
                           hf_den = hf_den + trans
                        else
                           nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob(num_prob) = &
                               prob_sum/prob_norm
                           idb = 1
                           iss = 0
                           call pack_data(Ix_f, ip_f, ns_f, idb, l, iss, itemp)
                           nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list(num_prob) = itemp
                        end if
                     end do
                     do l = lm_min, nucleus(i_f)%lmax_M, 2
                        trans = ML_trans(i_f,l,e_gamma)
                        if(trans < trans_e_cut)cycle
                        prob = trans/hf_den2
                        if(prob <= prob_cut) cycle
                        num_prob = num_prob + 1
                        prob_sum = prob_sum + trans
                        if(ii == 1)then
                           hf_den = hf_den + trans
                        else
                           nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob(num_prob) = &
                               prob_sum/prob_norm
                           idb = 1
                           iss = 1
                           call pack_data(Ix_f,ip_f,ns_f,idb,l,iss,itemp)
                           nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list(num_prob) = itemp
                        end if
                     end do
                  end do

                  do n_f = 1, n, 1              !  loop over final excitation energies
                     e_gamma=energy-nucleus(i_f)%e_grid(n_f)
                     if(e_gamma <= 1.0d-6)cycle
                     e_gamma = max(e_gamma,nucleus(i_f)%delta_e(n_f)/10.0)
!---------------------------   Start with Electric decay 
                     do l = 1, nucleus(i_f)%lmax_E                                             !  loop over EL decays
                        trans = EL_trans(i_f, l, e_gamma, energy)
                        if(trans < trans_e_cut)cycle
                        ip_f=iand((ip+l),1)                                                !  parity of final state
                        xI_f_min=abs(xI_i-real(l,kind=8))                                  !  min final spin
                        xI_f_max=min(xI_f_max1,xI_i+real(l,kind=8))                        !  max final spin
                        Ix_f_min = nint(xI_f_min-nucleus(i_f)%jshift)                      !  min j-index
                        Ix_f_max = nint(xI_f_max-nucleus(i_f)%jshift)                      !  max j-index
                        do Ix_f = Ix_f_min, Ix_f_max                                       !  loop over final j
                           xI_f = real(Ix_f,kind=8) + nucleus(i_f)%jshift
                           if(xI_i < 1.0d-3 .and. xI_f < 1.0d-3)cycle                      !  O -> 0 not allowed
                           trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*     &
                                             nucleus(i_f)%delta_e(n_f)
                           prob = trans_eff/hf_den2

                           if(prob <= prob_cut) cycle
                           num_prob = num_prob + 1
                           prob_sum = prob_sum + trans_eff
                           if(ii == 1)then
                              hf_den = hf_den + trans_eff
                           else
                              nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob(num_prob) = prob_sum/prob_norm
                              idb = 0
                              iss = 0
                              call pack_data(Ix_f, ip_f, n_f, idb, l, iss, itemp)
                              nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list(num_prob) = itemp
                           end if
                        end do
                     end do
!---------------------------   Now Magnetic decay 
                     do l = 1, nucleus(i_f)%lmax_M                                         !  loop over ML decays
                        trans = ML_trans(i_f, l, e_gamma)
                        if(trans < trans_e_cut)cycle
                        ip_f = iand((ip+l+1),1)                                            !  parity of final state
                        xI_f_min = abs(xI_i-real(l,kind=8))                                !  min final spin
                        xI_f_max = min(xI_f_max1,xI_i+real(l,kind=8))                      !  max final spin
                        Ix_f_min = nint(xI_f_min-nucleus(i_f)%jshift)                      !  min j-index
                        Ix_f_max = nint(xI_f_max-nucleus(i_f)%jshift)                      !  max j-index
                        do Ix_f = Ix_f_min, Ix_f_max                                       !  loop over final j
                           xI_f = real(Ix_f,kind=8) + nucleus(i_f)%jshift
                           if(xI_i < 1.0d-3.and.xI_f < 1.0d-3)cycle                        !  O -> 0 not allowed
                           trans_eff = trans*nucleus(i_f)%bins(Ix_f,ip_f,n_f)%rho*       &
                                             nucleus(i_f)%delta_e(n_f)
                           prob = trans_eff/hf_den2
                           if(prob <= prob_cut) cycle
                           num_prob = num_prob + 1
                           prob_sum = prob_sum + trans_eff
                           if(ii == 1)then
                              hf_den = hf_den + trans_eff
                           else
                              nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob(num_prob) = prob_sum/prob_norm
                              idb = 0
                              iss = 1

                              call pack_data(Ix_f, ip_f, n_f, idb, l, iss, itemp)

                              nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list(num_prob) = itemp
                           end if
                        end do
                     end do
                  end do
               end if

               if(num_prob > 0)then
                  if(ii == 1)then
                     nucleus(icomp)%bins(Ix_i,ip,n)%num_decay =                        &
                        nucleus(icomp)%bins(Ix_i,ip,n)%num_decay + 1
                     nucleus(icomp)%bins(Ix_i,ip,n)%decay_to(ifi) = i_f
                     nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle(ifi) = k
                     nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%num_decay = num_prob
                     allocate(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob(num_prob))
                     allocate(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list(num_prob))
                     nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob(1:num_prob) = 0.0d0
                     nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list(1:num_prob) = 0
                     nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(ifi) = prob_sum
                     nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(if1) = prob_sum
                     num_tot = num_tot + num_prob
                  end if
                  ifi = ifi + 1
               end if
            end do                           !  Finish do if1 = 1, nucleus(icomp)%num_decay
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Fission decay channel         ----------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(nucleus(icomp)%fission)then
               call Fission_transmission(icomp,energy,xI_i,ip,F_trans)
               prob = F_trans(4)/hf_den2

               if(ii == 2)then
                  nnnn = nucleus(icomp)%num_decay + 1
                  nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(nnnn) = F_trans(4)
                  nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(nnnn+1) = F_trans(1)
                  nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(nnnn+2) = F_trans(2)
               end if

               if(prob > prob_cut)then
                  prob_sum = prob_sum + F_trans(4)
                  if(ii == 1)then
                     hf_den = hf_den + F_trans(4)
                     nucleus(icomp)%bins(Ix_i,ip,n)%num_decay =                        &
                        nucleus(icomp)%bins(Ix_i,ip,n)%num_decay + 1
                     nnn = nucleus(icomp)%bins(Ix_i,ip,n)%num_decay
                  else
                     nnn = nucleus(icomp)%bins(Ix_i,ip,n)%num_decay
                     nnnn = nucleus(icomp)%num_decay + 1
                     nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(nnn) = F_trans(4)
                     nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(nnnn) = F_trans(4)
                     nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(nnnn+1) = F_trans(1)
                     nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(nnnn+2) = F_trans(2)
                     nucleus(icomp)%bins(Ix_i,ip,n)%decay_to(nnn) = -1                       !  Another signal that it is fission
                     nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle(nnn) = 7                  !  Another signal that this a fission decay
                     nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(nnn)%num_decay = -1           !  Signal that this is a fission event
                  end if
               end if
            end if
!---------------------------------------------------------------------------------------------------
         end do                  !   Finished do ii = 1, 2 

         nucleus(icomp)%bins(Ix_i,ip,n)%HF_den = hf_den

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Clean up probabilities to get them to line up to 1.00000000000, so we can eliminate some
!----   potential traps
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
10101   continue


         nnn = nucleus(icomp)%bins(Ix_i,ip,n)%num_decay
         if(nnn > 0)then
            prob_sum = 0.0d0
            do if1 = 1, nnn
               prob_sum = prob_sum + nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(if1)
            end do

            if(nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(nnn) > 0.0d0)then
               do if1 = 1, nnn
                  nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(if1) =                           &
                     nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(if1)/prob_sum
               end do
            end if
         end if
         nnn = nucleus(icomp)%num_decay
         if(nucleus(icomp)%fission)nnn = nnn + 1

!---------------------------------------------------------------------------------------------
      end do                            !   Finish: ip = 0, 1 
   end do                               !   Finish: do n = 1,nbin 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Now that the decay probabilities have been computed on separate nodes
!----   They need to be sent to all the other nodes
!----   my_proc is the processor that worked on the (Ix_i,ip,n) block of data above
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(nproc > 1) then
      do n = 1, nbin                  
         do iloop = 0, num_Ix_ip
            my_proc = mod(iloop,nproc)
            Ix_i = mod(iloop,Ix_i_max+1)
            ip = iloop/(Ix_i_max+1)
!----   Sync the MPI processes to be ready to Bcast data
            call MPI_Barrier(icomm,ierr)
!----  %num_decay
            num_data = 1
            call MPI_BCAST(nucleus(icomp)%bins(Ix_i,ip,n)%num_decay, num_data, MPI_INTEGER,        &
                           my_proc, icomm, ierr)
!----  %HF_prob
            if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob))                  &
                allocate(nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob(ndecay+1))
            num_data = ndecay + 1
            call MPI_BCAST(nucleus(icomp)%bins(Ix_i,ip,n)%HF_prob, num_data, MPI_REAL8, my_proc,   &
                           icomm, ierr)
!----  %decay_to
            if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%decay_to))                 &
                allocate(nucleus(icomp)%bins(Ix_i,ip,n)%decay_to(ndecay+1))
            num_data = ndecay + 1
            call MPI_BCAST(nucleus(icomp)%bins(Ix_i,ip,n)%decay_to, num_data, MPI_INTEGER,         &
                           my_proc, icomm, ierr)
!----  %decay_particle
            if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle))           &
                allocate(nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle(ndecay+1))
            num_data = ndecay + 1
            call MPI_BCAST(nucleus(icomp)%bins(Ix_i,ip,n)%decay_particle, num_data, MPI_INTEGER,   &
                           my_proc, icomm, ierr)
!----  %HF_trans
            if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans))                 &
                allocate(nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans(ndecay+3))
            num_data = ndecay + 3
            call MPI_BCAST(nucleus(icomp)%bins(Ix_i,ip,n)%HF_trans, num_data, MPI_REAL8, my_proc,  &
                           icomm, ierr)
!----   %nuke_decay
            if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay))               &
                allocate(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ndecay+1))
!---   Sync up again to start filling nuke_decay data types
            call MPI_Barrier(icomm,ierr)
!----   loop over elements attached to %nuke_decay
            do ifi = 1, ndecay + 1
                num_data = 1
                call MPI_BCAST(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%num_decay, num_data, &
                               MPI_INTEGER, my_proc, icomm, ierr)
                num_data = nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%num_decay
                if(num_data >= 0)then
                   if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob))     &
                       allocate(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob(num_data))
                   call MPI_BCAST(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_prob, num_data,&
                                  MPI_REAL8, my_proc, icomm, ierr)
                   num_data = nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%num_decay
                   if(.not. allocated(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list))     &
                       allocate(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list(num_data))
                   call MPI_BCAST(nucleus(icomp)%bins(Ix_i,ip,n)%nuke_decay(ifi)%decay_list, num_data,&
                                  MPI_INTEGER, my_proc, icomm, ierr)
               end if
            end do
         end do
      end do
   end if

   call MPI_Barrier(icomm,ierr)
   num_data = 1
   call MPI_Allreduce(MPI_IN_PLACE, num_tot, num_data, MPI_INTEGER, MPI_SUM, icomm, ierr)
   if(print_me)then
      write(6,*)'Finished calculating the decay probabilities for ',nucleus(icomp)%Label
      write(6,*)'Total number = ', num_tot
   end if


return
end subroutine HF_denominator
