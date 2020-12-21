!
!*******************************************************************************
!
   subroutine MC_decay_bin(icomp_i, Ix_i, ip_i, nbin_i,                &
                           icomp_f, Ix_f, ip_f, nbin_f, idb,           &
                           n_dat,dim_part,num_part,part_data,          &
                           Ang_L_max, part_Ang_data,                   &
                           Boost_Lab, Boost_COM) 
!
!*******************************************************************************
!
!  Discussion:
!
!    idbThis subroutine Monte Carlo decays a continuous energy bin
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
   use useful_data
   implicit none
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Passed Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!-----------  Input  ------------------------------------------------
   integer(kind=4), intent(in) :: icomp_i, Ix_i, ip_i, nbin_i 
   integer(kind=4), intent(out) :: icomp_f, Ix_f, ip_f, nbin_f, idb
   integer(kind=4), intent(in) :: n_dat, dim_part
   integer(kind=4), intent(inout) :: num_part
   real(kind=8), intent(inout) :: part_data(n_dat,dim_part)
   integer(kind=4), intent(in) :: Ang_L_max
   real(kind=8), intent(out) :: part_Ang_data(0:Ang_L_max,dim_part)
   real(kind=8), intent(inout) :: Boost_Lab(0:3,0:3), Boost_COM(0:3,0:3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Internal Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   real(kind=8) :: prob
   integer(kind=4) :: if1, iif1
   integer(kind=4) idex
   integer(kind=4) :: itemp
   integer(kind=4) :: mask6, mask10
   integer(kind=4) :: i,j, k, l, iss
   real(kind=8) :: e_i, e_f, ex_i, ex_f
   real(kind=8) :: xI_i, xI_f, xj_f, xj_f_min, xip_f
   real(kind=8) :: shift


   real(kind=8) :: costhp, sinthp, theta_0, phi_0
   real(kind=8) :: T_1, T_2, mass_1, mass_2, theta, phi
   real(kind=8) :: T_L, theta_L, phi_L

   integer(kind=4) :: icomp_ff, Ix_FF, ip_ff, nbin_ff, idb_ff, l_ff, iss_ff
   integer(kind=4) :: kk
   real(kind=8) :: e_ff

   real(kind=8) :: tally_prob
   real(kind=8) :: xnnn
   integer(kind=4) :: nnn, nnn_map(100)

   real(kind=8) :: test_prob(15)
   integer(kind=4) :: test_idex(15), test_if1(15)

   real(kind=8), allocatable :: prob_if1(:)
   integer(kind=4), allocatable :: num_if1(:)
   integer(kind=4) :: num_i, ii, ic
   real(kind=8) :: prob_sum, prob_sum_i

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   External Functions
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   real(kind=8) :: random_64
!-----------------    Data stored for each event
!-----   part_data(1,num_part) = icomp_f
!-----   part_data(2,num_part) = k             !--------   Type of particle emitted
!-----   part_data(3,num_part) = xI_f          !--------   Spin of final state
!-----   part_data(4,num_part) = xip_f         !--------   parity of final state
!-----   part_data(5,num_part) = nbin_f        !--------   Bin or state number of final state
!-----   part_data(6,num_part) = idb           !--------   = 0 for continuous bin, = 1 for a discrete state
!-----   part_data(7,num_part) = l             !--------   Orbital angular momentum
!-----   part_data(8,num_part) = xj_f          !--------   xj for emitted particle
!-----   part_data(9,num_part) = KE            !--------   KE of emitted particle  in emission rest frame
!-----   part_data(10,num_part) = thetap       !--------   Angle theta in emission rest frame
!-----   part_data(11,num_part) = phi          !--------   Angle phi
!-----   part_data(12,num_part) = T_1          !--------   emitted kinetic energy in COM frame
!-----   part_data(13,num_part) = theta        !--------   theta angle in COM frame
!-----   part_data(14,num_part) = phi          !--------   phi angle in COM frame
!-----   part_data(15,num_part) = T_L          !--------   kinetic energy of emitted particle in Lab frame
!-----   part_data(16,num_part) = theta_L      !--------   theta angle in Lab frame
!-----   part_data(17,num_part) = phi_L        !--------   phi angle in Lab frame
!-----   part_data(18,num_part) = T_2          !--------   recoil kinetic energy in Lab frame
!-----   part_data(19,num_part) = Tally        !--------   tally probability
!-----   part_data(20,num_part) = ex_i         !--------   excitation energy of the system after particle emission
!-----   part_data(21,num_part) = icomp_i      !--------   compound nucleus that decayed



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Start Program
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

   e_i = nucleus(icomp_i)%e_grid(nbin_i)
   ex_i = e_i
   xI_i = real(iX_i) + nucleus(icomp_i)%jshift

   if(num_part >= 1)ex_i = part_data(20,num_part)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  First, check which nucleus it decays to
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   tries = 0

22   prob = random_64(iseed)                !  starting probability

   tries = tries + 1

   test_prob(tries) = prob

   if1 = 0
   icomp_f = icomp_i
   nbin_f = -1
   if(nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay == 0)return   ! no way to decay, hung up
!
   tally_prob = 1.0d0

   if(biased_sampling)then
      do if1 = 1, nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
         if(prob <= nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(if1))exit
      end do
   else
      xnnn = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
      if1 = int(random_64(iseed)*xnnn) + 1
      tally_prob = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(if1)
      if(if1 > 1) tally_prob = tally_prob - nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(if1-1)
      tally_prob = tally_prob*xnnn
   end if

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  set icomp_f, note if1 > num_decay ->  Fission 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if(nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1) == -10)then     !  Fission 
      icomp_f = icomp_i
      nbin_f = 0
      num_part = num_part + 1
      part_data(1,num_part) = real(icomp_f,kind=8)
      part_data(2,num_part) = -2.0d0        !  one signal this is a fission event
      part_data(3,num_part) = xI_i
      part_data(4,num_part) = 0.0d0
      part_data(5,num_part) = nbin_f      !  Other signal this is a fission event
      part_data(6,num_part) = 0.0d0
      part_data(7,num_part) = 0.0d0
      part_data(19,num_part) = tally_prob
      part_data(20,num_part) = -10.0d0
      part_data(21,num_part) = icomp_i
      return
   end if

   icomp_f = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_to(if1)
   if(icomp_f < 1) stop 'icomp_f < 1 after attempting to decay'



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  Next, find which state it decays to
!--------  Establish probability, then search for it in list
!--------  using bisection between upper and lower
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if(nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%num_decay < 1)then
      write(6,*)'Bin Decay'
      write(6,*)icomp_i,Ix_i,ip_i,nbin_i,if1
   end if


   prob = random_64(iseed)

   call find_prob(nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%num_decay,     &
                  nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%decay_prob,    &
                  prob, idex)

   test_idex(tries) = idex
   test_if1(tries) = if1
   k = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1)

   mask6 = 2**6 - 1
   mask10 = 2**9 - 1

   itemp = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%decay_list(idex)

   call unpack_data(Ix_f, ip_f, nbin_f, idb, l, iss, itemp)

   if(nbin_f == 0)then
      write(6,*)'Error nbin_f = 0 in MC_decay_bin'
      write(6,*)nucleus(icomp_i)%A, Ix_i, ip_i, nbin_i
      write(6,*)idb,l,k
      stop
   end if


   xI_f = real(Ix_f) + nucleus(icomp_f)%jshift
   xj_f_min = abs(dfloat(l) - particle(k)%spin) 
   xj_f = real(iss) + xj_f_min
   xip_f = 2*ip_f - 1

   if(idb == 0)then
      e_f = ex_i - nucleus(icomp_i)%sep_e(k) -           &
                   nucleus(icomp_f)%e_grid(nbin_f)
      if(e_f < -0.5d0*de)then
         write(6,*)'icomp_i = ',icomp_i,' icomp_f= ',icomp_f
         write(6,*)'k= ', k
         write(6,*)'nbin_i = ',nbin_i,' nbin_f = ',nbin_f
         write(6,*)'ex_f = ',ex_f
         write(6,*)'ex_i = ',ex_i
         write(6,*)'Separation energy = ',nucleus(icomp_i)%sep_e(k)
         write(6,*)'Initial Bin energy = ',nucleus(icomp_f)%e_grid(nbin_i)
         write(6,*)'Final Bin energy = ',nucleus(icomp_f)%e_grid(nbin_f)
         write(6,*)'num_part = ',num_part
         stop 'e_f < -de/2 in MC_decay_bin (1)'
      end if
      if(e_f < 0.0d0)then
         e_f = (e_f + 0.5d0*de)*random_64(iseed)
         ex_f = ex_i - e_f - nucleus(icomp_i)%sep_e(k)
      elseif(e_f <= de)then           !   decays to same bin
         e_f = e_f*random_64(iseed)
         ex_f = ex_i - e_f - nucleus(icomp_i)%sep_e(k)
      elseif(e_f > de)then
         shift = 0.5d0*de*(2.0d0*random_64(iseed) - 1.0d0)
         ex_f = nucleus(icomp_f)%e_grid(nbin_f) + shift
         e_f = ex_i - nucleus(icomp_i)%sep_e(k) - ex_f
      end if
   else  
      e_f = ex_i - nucleus(icomp_i)%sep_e(k) -         &
                   nucleus(icomp_f)%state(nbin_f)%energy
      ex_f = nucleus(icomp_f)%state(nbin_f)%energy
      if(e_f < 0.0d0)then                !  Decay forbidden due to energy conservation issue
         if(tries < 8)then              !  try a few more times
            goto 22              
         else                            !  Tried too many times force decay so we can finish
            re_baseline = re_baseline + 1
!------------   Print out to analyze situation -- caused by previous decay, with "jitter" dropping
!------------   the initial excitation energy below a discrete state and decay not allowed.
!------------   First, try again, up to 15 times. If failed, then recompute probabilities for this
!------------   special case.
!            write(6,*)
!            write(6,'(''*************************************************************************************'')')
!            write(6,'(''Warning MC attempted to decay with negative energy of e_f = '',1pe15.7,'' MeV'')')e_f
!            write(6,'(''Caused by jitter in emission spectrum and decay to discrete states near threshold'')')
!            write(6,'(''Recomputing decay probabilities to find new decay path'')')
!            write(6,'(''*************************************************************************************'')')
!            write(6,*)
            write(26,*)
            write(26,'(''*************************************************************************************'')')
            write(26,'(''Warning MC attempted to decay with negative energy of e_f = '',1pe15.7,'' MeV'')')e_f
            write(26,'(''Caused by jitter in emission spectrum and decay to discrete states near threshold'')')
            write(26,'(''Recomputing decay probabilities to find new decay path'')')
            write(26,'(''*************************************************************************************'')')
            write(26,'(''Current sample # '',i10)')nsamp
            write(26,*)
            write(26,*)'Biased sampling = ',biased_sampling
            write(26,'(''Decay from continuous bin to discrete state'')')
!            write(26,'(''Tried 15 times and could not make the decay work. Will force decay with e_f = 0.001 MeV'')')
!            write(26,'(''decay particle = '',i4)')k
!  80        format('Initial nucleus ',i4,'  J = ',f3.1,' Ix_i = ',i3,' ip_i = ',i3,' initial bin ',i4,   &
!                   ' energy = ',f8.4,' bin energy ',f8.4,' Separation energy = ',f8.4)
!            write(26,80)icomp_i,xI_i,Ix_i,ip_i,nbin_i,ex_i,nucleus(icomp_i)%e_grid(nbin_i),nucleus(icomp_i)%sep_e(k)
!            write(26,'(''Final nucleus '',i4,'' Final state# '',i4,'' energy = '',f8.4)')  &
!                 icomp_f,nbin_f,nucleus(icomp_f)%state(nbin_f)%energy
!            write(26,'(''Emitted energy = '',f8.4)')e_f
!           write(26,'(''Number of possible particle decays '',i10)')nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
!            write(26,*)'Last tried to decay to path #',if1,'Ran = ',prob
!            do i = 1, nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
!               write(26,'(i6,1x,f10.7)')i,nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(i)
!            end do
!            do i = 1, tries
!               write(26,'(''Rolled probabilities'',i6,1x,f10.7,'' if1 = '',i4,'' idex = '',i4)')             &
!                    i,test_prob(i),test_if1(i),test_idex(i)
!            end do
            do i = 1, nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
               icomp_ff = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_to(i)
!               write(26,'(''Decay path# '',i4,'' icomp_f = '',i6,'' Num decays = '',i6)')i,icomp_ff,         &
!                   nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(i)%num_decay
!               if(i == if1)write(26,'(''Last tried to decay to posiibility #'',i6)')idex
!               write(26,'('' Number of decays in this channel '',i10)')                                      &
!                    nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(i)%num_decay
               do j = 1, nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(i)%num_decay                !
                  itemp = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(i)%decay_list(j)
                  call unpack_data(Ix_ff, ip_ff, nbin_ff, idb_ff, l_ff, iss_ff, itemp)
!                  if(idb_ff == 0)write(26,'(''Possibile decay to bin  '',i6,1pe15.7,2i6,1x,0pf8.4)')j,       &
!                       nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(i)%decay_prob(j),                  &
!                       idb_ff,nbin_ff,nucleus(icomp_ff)%e_grid(nbin_ff)
!                  if(idb_ff == 1)write(26,'(''Possible decay to state '',i6,1pe15.7,2i6,1x,0pf8.4)')j,       &
!                       nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(i)%decay_prob(j),                  &
!                       idb_ff,nbin_ff,nucleus(icomp_ff)%state(nbin_ff)%energy
               end do
            end do
            flush(26)
            ex_f = nucleus(icomp_f)%state(nbin_f)%energy
            num_i = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
            allocate(prob_if1(num_i))
            allocate(num_if1(num_i))
            prob_sum_i = 0.0d0
            do ii = 1, num_i
               prob_sum_i = 0.0d0
               ic = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_to(ii)
               if(nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(ii) == -10)then     !  Fission 
                  prob_sum_i = 1.0d0
               else
                  kk = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(ii)
!      write(26,*)'ii = ',ii,' kk = ',kk,' icomp = ',ic
!      write(26,*)'num ',nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(ii)%num_decay
                  num_if1(ii) = 0
                  do j = 1, nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(ii)%num_decay
                     itemp = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(ii)%decay_list(j)

                     call unpack_data(Ix_ff, ip_ff, nbin_ff, idb_ff, l_ff, iss_ff, itemp)
                     e_ff = 0.0d0
                     if(idb_ff == 0)then
                        e_ff = ex_i - nucleus(icomp_i)%sep_e(kk) -         &
                               nucleus(ic)%e_grid(nbin_ff)
                     elseif(idb_ff == 1)then
                        e_ff = ex_i - nucleus(icomp_i)%sep_e(kk) -         &
                               nucleus(ic)%state(nbin_ff)%energy
                     end if
                     if(e_ff < 1.0d-6)exit
                     prob_sum_i = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(ii)%decay_prob(j)
!                  if(idb_ff == 0)write(26,'(''Bin   '',i6,1pe15.7,i6,1x,0pf8.4,1x,1pE15.7)')j,        &
!                       nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(ii)%decay_prob(j),          &
!                       nbin_ff,nucleus(ic)%e_grid(nbin_ff),prob_sum_i
!                  if(idb_ff == 1)write(26,'(''State '',i6,1pe15.7,i6,1x,0pf8.4,1x,1pE15.7)')j,        &
!                       nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(ii)%decay_prob(j),          &
!                       nbin_ff,nucleus(ic)%state(nbin_ff)%energy,prob_sum_i
                     num_if1(ii) = num_if1(ii) + 1
                  end do
               end if
               if(ii == 1)then
                  prob_if1(ii) = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(ii)*prob_sum_i
               else
                  prob_if1(ii) = (nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(ii)-    &
                                   nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(ii-1))*prob_sum_i
               end if
!               write(26,'('' Prob for decay path '',i4,'' = '',f10.7)')ii,prob_if1(ii)
            end do
!----    Check if no decay is possible, leading to the bin being hung up, i.e., computational isomer
            if(prob_sum_i < 1.0d-8)then
               if1 = 0
               icomp_f = icomp_i
               nbin_f = -1
               return
            end if

!----    Check if any decays have been excluded and make map for unbiased sampling
!----    This fixes a bug where a channel is now excluded, but with unbiased sampling
!----    all possibilities are equally probable, hecen to need to skip it
!----    nnnn is the number of decay paths with non-zero probability, i.e. has a decay possibility
!----    nnnn_map is the original index
            nnn = 0
            nnn_map(1:100) = 0
            do ii = 1, num_i
               if(prob_if1(ii) > 1.0d-10)then
                  nnn = nnn + 1
                  nnn_map(nnn) = ii
               end if
            end do

!----    have baseline probabilities excluding decays to states that are forbidden
!----    renormalize probabilities for particle decays

            prob_sum = prob_if1(1)
            do ii = 2, num_i
               prob_if1(ii) = prob_if1(ii) + prob_if1(ii-1)               
            end do
            do ii = 1, num_i
               prob_if1(ii) = prob_if1(ii)/prob_if1(num_i)
!    write(26,*)'new Probabilities',ii,prob_if1(ii), num_if1(ii)
            end do

            prob = random_64(iseed)                !  starting probability
!    write(26,*)'Finding new path with prob = ',prob
            if(biased_sampling)then
               do if1 = 1, num_i
                  if(prob <= prob_if1(if1))exit
               end do
            else
               xnnn = nnn
               iif1 = int(random_64(iseed)*xnnn) + 1
               if1 = nnn_map(iif1)
               tally_prob = prob_if1(if1)
               if(if1 > 1) tally_prob = tally_prob - prob_if1(if1-1)
               tally_prob = tally_prob*xnnn
            end if
!            write(26,*)'Decay path = ',if1,' k = ',nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1)
            if(nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1) == -10)then     !  Fission 
               icomp_f = icomp_i
               nbin_f = 0
               num_part = num_part + 1
               part_data(1,num_part) = real(icomp_f,kind=8)
               part_data(2,num_part) = -2.0d0        !  one signal this is a fission event
               part_data(3,num_part) = 0.0d0
               part_data(4,num_part) = 0.0d0
               part_data(5,num_part) = nbin_f      !  Other signal this is a fission event
               part_data(6,num_part) = 0.0d0
               part_data(7,num_part) = 0.0d0
               part_data(19,num_part) = tally_prob
               part_data(20,num_part) = -10.0d0
               part_data(21,num_part) = icomp_i
               return
            end if

            icomp_f = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_to(if1)
            if(icomp_f < 1) stop 'icomp_f < 1 after renormalization of probabilities'
            prob = random_64(iseed)
            k = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1)

!    write(26,*)'Search for decay. if1 = ',if1,num_if1(if1),prob
            call find_prob(num_if1(if1),                                                 &
                  nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%decay_prob,    &
                  prob, idex)
            itemp = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%decay_list(idex)

            call unpack_data(Ix_f, ip_f, nbin_f, idb, l, iss, itemp)


!            write(26,*)idex,itemp,Ix_f,ip_f,nbin_f,idb,l,iss
!                  if(idb == 0)write(26,'(''Bin   '',i6,1pe15.7,i6,1x,0pf8.4,1x,1pE15.7)')idex,        &
!                       nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%decay_prob(idex),      &
!                       nbin_ff,nucleus(icomp_f)%e_grid(nbin_f),prob_sum_i
!                  if(idb == 1)write(26,'(''State '',i6,1pe15.7,i6,1x,0pf8.4,1x,1pE15.7)')idex,        &
!                       nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%decay_prob(idex),      &
!                       nbin_ff,nucleus(icomp_f)%state(nbin_f)%energy,prob_sum_i
            

            if(nbin_f == 0)then
               write(6,*)'Error nbin_f = 0 in MC_decay_bin'
               write(6,*)nucleus(icomp_i)%A, Ix_i, ip_i, nbin_i
               write(6,*)idb,l,k
               stop
            end if
            if(idb == 0)then
               e_f = ex_i - nucleus(icomp_i)%sep_e(k) -           &
                     nucleus(icomp_f)%e_grid(nbin_f)
               if(e_f < -0.5d0*de)stop 'e_f < -de/2 in MC_decay_bin (2)'
               if(e_f < 0.0d0)then
                  e_f = (e_f + 0.5d0*de)*random_64(iseed)
                  ex_f = ex_i - e_f - nucleus(icomp_i)%sep_e(k)
               elseif(e_f <= de)then           !   decays to same bin
                   e_f = e_f*random_64(iseed)
                   ex_f = ex_i - e_f - nucleus(icomp_i)%sep_e(k)
               elseif(e_f > de)then
                  shift = 0.5d0*de*(2.0d0*random_64(iseed) - 1.0d0)
                  ex_f = nucleus(icomp_f)%e_grid(nbin_f) + shift
                  e_f = ex_i - nucleus(icomp_i)%sep_e(k) - ex_f
               end if
            else  
               e_f = ex_i - nucleus(icomp_i)%sep_e(k) -         &
                     nucleus(icomp_f)%state(nbin_f)%energy
               ex_f = nucleus(icomp_f)%state(nbin_f)%energy
               if(e_f < 0.0d0)then                !  Decay forbidden due to energy conservation issue
                  write(6,*)'Still does not work negative energy in MC_decay_HF'
                  stop
               end if
            end if
            deallocate(prob_if1)
            deallocate(num_if1)
         end if
      end if
   end if

   num_part = num_part +1
   if(num_part > dim_part) stop 'num_part > dim_part'

   part_data(1,num_part) = real(icomp_f,kind=8)
   part_data(2,num_part) = real(k,kind=8)
   part_data(3,num_part) = xI_f
   part_data(4,num_part) = xip_f
   part_data(5,num_part) = real(nbin_f,kind=8)
   part_data(6,num_part) = real(idb,kind=8)
   part_data(7,num_part) = real(l,kind=8)
   if( k == 0)then
      part_data(8,num_part) = real(iss,kind=8)
   else
      part_data(8,num_part) = xj_f
   end if

   costhp = 2.0d0*random_64(iseed) - 1.0d0
   if(abs(costhp) > 1.0d0)stop 'cos(theta) wrong in MC_decay_bin'
   sinthp = sqrt(1.0d0-costhp**2)
   theta_0 = acos(costhp)
   phi_0 = two_pi*random_64(iseed)

   mass_1 = particle(k)%Mass
   if(idb == 0)then
      mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%e_grid(nbin_f)
   else
      mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%state(nbin_f)%energy
   end if
   call Boost_frame(e_f, mass_1, mass_2, theta_0, phi_0,                   &
                    Boost_Lab, Boost_COM, T_1, theta, phi,                 &
                    T_2, T_L, theta_L, phi_L)

   part_data(9,num_part) = e_f
   part_data(10,num_part) = theta_0
   part_data(11,num_part) = phi_0
   part_data(12,num_part) = T_1
   part_data(13,num_part) = theta
   part_data(14,num_part) = phi
   part_data(15,num_part) = T_L
   part_data(16,num_part) = theta_L
   part_data(17,num_part) = phi_L
   part_data(18,num_part) = T_2
   part_data(19,num_part) = tally_prob
   part_data(20,num_part) = ex_f
   part_data(21,num_part) = icomp_i
   nucleus(icomp_f)%Kinetic_Energy = T_2
   if(.not.pop_calc)part_Ang_data(0,num_part) = 0.5d0
   return

   end subroutine MC_decay_bin
!
!*******************************************************************************
!
   subroutine MC_decay_state(icomp_i, istate_i,                       &
                             n_dat,dim_part,num_part,part_data,       &
                             Ang_L_max,part_Ang_data, ichan, in,      &
                             Boost_Lab, Boost_COM)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine Monte Carlo decays a discrete state
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
   use Channel_info
   implicit none
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Passed Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!-----------  Input  ------------------------------------------------
   integer(kind=4), intent(in) :: icomp_i, istate_i
!-----------  Output ------------------------------------------------
   integer(kind=4), intent(in) :: n_dat, dim_part
   integer(kind=4), intent(inout) :: num_part
   real(kind=8), intent(inout) :: part_data(n_dat,dim_part)
   integer(kind=4), intent(in) :: Ang_L_max
   real(kind=8), intent(out) :: part_Ang_data(0:Ang_L_max,dim_part)
   integer(kind=4), intent(in) :: ichan, in
   real(kind=8), intent(inout) :: Boost_Lab(0:3,0:3), Boost_COM(0:3,0:3)


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Internal Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   real(kind=8) :: prob
   integer(kind=4) :: j, istate, n_f, k, icomp_f
   integer(kind=4) :: idb, l, iss
   real(kind=8) :: check
   real(kind=8) :: xI_f
   real(kind=8) :: xj_f, xip_f, e_gamma


   real(kind=8) :: costhp, sinthp, theta_0, phi_0
   real(kind=8) :: mass_1, mass_2, T_1, T_2, theta, phi
   real(kind=8) :: T_L, theta_L, phi_L

   real(kind=8) :: tally_prob, tally_weight
   integer(kind=4) :: n

!   real(kind=8) :: mass_i, mass_f, mass
!   real(kind=8) :: xkpxc, xk_0xc, xKKp_0xc, gamma, theta_f, phi_f

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   External Functions
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   real(kind=8) :: random_64

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Start Program
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

   n_f = -1
   icomp_f = icomp_i
   istate = istate_i

   k = 0

!   write(6,*)'istate ', icomp_i, istate

 1 continue

   tally_prob = 1.0d0

   if(istate > nucleus(icomp_i)%num_discrete)then
       write(6,*)'Trying to decay to a state # greater than num_discrete'
       write(6,*)icomp_i,istate,nucleus(icomp_i)%num_discrete
!   stop
   end if
   if(istate == 1 .or. nucleus(icomp_i)%state(istate)%isomer) return

   prob = random_64(iseed)
!    write(46,'(1x,e30.16)')prob
!    flush(46)

   check = 0.0d0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Put in a check to see if the state has branches to decay to +
!-----   If not, force decay to ground state. This is an override to +
!-----   complete calculation just in case a state slips through     +
!-----   with no decay path. Note, isomers are trapped up above      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   if(nsamp == 3391)then
!      write(36,*)'Decay discrete state'
!      write(36,*)icomp_i,istate
!      write(36,*)nucleus(icomp_i)%state(istate)%nbranch
!      flush(36)
!   end if

   if(nucleus(icomp_i)%state(istate)%nbranch == 0)then
      n_f = 1
      idb = 1
      xj_f=nucleus(icomp_i)%state(n_f)%spin
      xip_f=nucleus(icomp_i)%state(n_f)%parity
      iss = 0
      l = 0
      e_gamma = nucleus(icomp_i)%state(istate)%energy - nucleus(icomp_i)%state(n_f)%energy

      num_part = num_part +1
      if(num_part > dim_part) stop 'num_part > dim_part'

      part_data(1,num_part) = real(icomp_i,kind=8)
      part_data(2,num_part) = real(k,kind=8)
      part_data(3,num_part) = xj_f
      part_data(4,num_part) = xip_f
      part_data(5,num_part) = real(n_f,kind=8)
      part_data(6,num_part) = real(idb,kind=8)
      part_data(7,num_part) = real(l,kind=8)
      part_data(8,num_part) = real(iss,kind=8)

      costhp = 2.0d0*random_64(iseed) - 1.0d0
      if(abs(costhp) > 1.0d0)stop 'cos(theta) wrong in MC_decay_state'
      sinthp = sqrt(1.0d0-costhp**2)
      theta_0 = acos(costhp)
      phi_0 = two_pi*random_64(iseed)


      mass_1 = particle(k)%Mass
      mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%state(n_f)%energy

      call Boost_frame(e_gamma, mass_1, mass_2, theta_0, phi_0,                &
                       Boost_Lab, Boost_COM, T_1, theta, phi,                  &
                       T_2, T_L, theta_L, phi_L)

      part_data(9,num_part) = e_gamma
      part_data(10,num_part) = theta_0
      part_data(11,num_part) = phi_0
      part_data(12,num_part) = T_1
      part_data(13,num_part) = theta
      part_data(14,num_part) = phi
      part_data(15,num_part) = T_L
      part_data(16,num_part) = theta_L
      part_data(17,num_part) = phi_L
      part_data(18,num_part) = T_2
      part_data(19,num_part) = tally_prob
      part_data(20,num_part) = nucleus(icomp_i)%state(n_f)%energy
      part_data(21,num_part) = icomp_i
      nucleus(icomp_f)%Kinetic_Energy = T_2


      if(.not.pop_calc)part_Ang_data(0,num_part) = 0.5e0
      return 
   end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Check which decay branch is followed ------------------------+
!----   This is the actual decay check       ------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   do j = 1, nucleus(icomp_i)%state(istate)%nbranch
      check = check + nucleus(icomp_i)%state(istate)%branch(j)
      n_f = nucleus(icomp_i)%state(istate)%ibranch(j)
      if(prob <= check) exit
   end do

   prob = random_64(iseed)
   if(prob <= nucleus(icomp_i)%state(istate)%p_gamma(j))then   !  Gamma decay
      k = 0
      if(track_gammas)then
         tally_weight = 1.0d0
         if(.not. biased_sampling)then
            do n = 1, num_part
               tally_weight = tally_weight*part_data(19,n)
            end do
         end if
         Exit_channel(ichan)%state(istate)%cs(j,in) =                           &  
                      Exit_channel(ichan)%state(istate)%cs(j,in) + tally_weight

      end if
   else
      k = -1                                                              !  Internal conversion 
   end if

   idb = 1
   xI_f=nucleus(icomp_i)%state(n_f)%spin
   xip_f=nucleus(icomp_i)%state(n_f)%parity
   iss = 0
   l = 0

   e_gamma = nucleus(icomp_i)%state(istate)%energy - nucleus(icomp_i)%state(n_f)%energy

   num_part = num_part + 1
   if(num_part > dim_part) stop 'num_part > dim_part'

   part_data(1,num_part) = real(icomp_i,kind=8)
   part_data(2,num_part) = real(k,kind=8)
   part_data(3,num_part) = xI_f
   part_data(4,num_part) = xip_f
   part_data(5,num_part) = real(n_f,kind=8)
   part_data(6,num_part) = real(idb,kind=8)
   part_data(7,num_part) = real(l,kind=8)
   part_data(8,num_part) = real(iss,kind=8)

   costhp = 2.0d0*random_64(iseed) - 1.0d0
   if(abs(costhp) > 1.0d0)stop 'cos(theta) wrong in MC_decay_state #2'
   sinthp = sqrt(1.0d0-costhp**2)
   theta_0 = acos(costhp)
   phi_0 = two_pi*random_64(iseed)

   mass_1 = particle(k)%Mass
   mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%state(n_f)%energy

   call Boost_frame(e_gamma, mass_1, mass_2, theta_0, phi_0,                &
                    Boost_Lab, Boost_COM, T_1, theta, phi,                  &
                    T_2, T_L, theta_L, phi_L)

   part_data(9,num_part) = e_gamma
   part_data(10,num_part) = theta_0
   part_data(11,num_part) = phi_0
   part_data(12,num_part) = T_1
   part_data(13,num_part) = theta
   part_data(14,num_part) = phi
   part_data(15,num_part) = T_L
   part_data(16,num_part) = theta_L
   part_data(17,num_part) = phi_L
   part_data(18,num_part) = T_2
   part_data(19,num_part) = tally_prob
   part_data(20,num_part) = nucleus(icomp_i)%state(n_f)%energy
   part_data(21,num_part) = icomp_i
   nucleus(icomp_f)%Kinetic_Energy = T_2

   if(.not.pop_calc)part_Ang_data(0,num_part) = 0.5e0         !  Isotropic emission

   istate = n_f

   goto 1

   return
   end subroutine MC_decay_state
!
!
   subroutine find_prob(num,prob_array,prob,ifind)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine search the array prob_array for the value prob
!    using bisection
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
   integer(kind=4), intent(in) :: num
   real(kind=8), intent(in) :: prob_array(num)
   real(kind=8), intent(in) :: prob
   integer(kind=4), intent(out) :: ifind
!----------   Internal Data
   integer(kind=4) lower, mid, upper
   logical found
   real(kind=8) :: p_norm

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Check if at extremes  -------------------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   p_norm = prob_array(num)
   if(p_norm <= 1.0d-20)then
       write(6,*)'p_norm too small in find_prob = ',p_norm
       stop
   end if
   upper = num
   lower = 1

   if(upper == lower)then
      ifind = upper
      return
   end if

   if(prob < prob_array(lower)/p_norm)then
     ifind = lower
     return
   end if

   if(prob > prob_array(upper-1)/p_norm .and. prob <= prob_array(upper)/p_norm)then
     ifind = upper
     return
   end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   It isn't at extremes, search by bisecting the array   -----+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   found = .false.

   do while(.not. found)
      mid = (upper - lower)/2 + lower
      if(mid == lower)then
         ifind = mid
         found = .true.
         return
      end if
      if(prob > prob_array(mid-1)/p_norm .and. prob <= prob_array(mid)/p_norm)then
         ifind = mid
         found = .true.
         return
      end if

      if(prob < prob_array(mid)/p_norm)then
         upper = mid
      else
         lower = mid
      end if
   end do

   return
   end
!
!*******************************************************************************
!
   subroutine find_prob_point(num,prob_array,prob,ifind)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine to find position of a probability within the pointer array
!    prob_array(ifind-1) <= prob <= prob_array(ifind)
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
   integer(kind=4), intent(in) :: num
   real(kind=8),pointer, intent(in) :: prob_array(:)
   real(kind=8), intent(in) :: prob
   integer(kind=4), intent(out) :: ifind
!----------   Internal Data
   integer(kind=4) lower, mid, upper
   logical found

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Check if at extremes  -------------------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   upper = num
   lower = 1

   if(upper == lower)then
      ifind = upper
      return
   end if

   if(prob < prob_array(lower))then
     ifind = lower
     return
   end if

   if(prob > prob_array(upper-1) .and. prob < prob_array(upper))then
     ifind = upper
     return
   end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   It isn't at extremes, search by bisecting the array   -----+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   found = .false.

   do while(.not. found)
      mid = (upper - lower)/2 + lower
      if(mid == lower)then
         ifind = mid
         return
      end if
      if(prob > prob_array(mid-1) .and. prob <= prob_array(mid))then
         ifind = mid
         return
      end if

      if(prob < prob_array(mid))then
         upper = mid
      else
         lower = mid
      end if
   end do

   return
   end
!
!*******************************************************************************
!
   subroutine MC_primary_decay(iproj,spin_target,                        &
                               l_i, is_i, Ix_i, e_i, icomp_i,            &
                               icomp_f, Ix_f, ip_f, nbin_f, idb,         &
                               n_dat, dim_part, num_part, part_data,     &
                               Ang_L_max, part_Ang_data,                 &
                               ixx_max, delta_x,                         &
                               Boost_Lab, Boost_COM)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine decays first compound state - done differently than others
!    due to width fluctuations                                            
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
   use useful_data
   implicit none
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Passed Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!-----------  Input  ------------------------------------------------
!   integer(kind=4) :: iseed
   integer(kind=4), intent(in) :: iproj
   real(kind=8), intent(in) :: spin_target
   integer(kind=4), intent(in) :: l_i, is_i, Ix_i, icomp_i
!-----------  Output ------------------------------------------------
   integer(kind=4), intent(out) :: icomp_f, Ix_f, ip_f, nbin_f, idb
   integer(kind=4), intent(in) :: n_dat, dim_part
   integer(kind=4), intent(inout) :: num_part
   real(kind=8), intent(inout) :: part_data(n_dat,dim_part)
   integer(kind=4), intent(in) :: Ang_L_max
   real(kind=8), intent(inout) :: part_Ang_data(0:Ang_L_max,dim_part)
   integer(kind=4), intent(in) :: ixx_max
   real(kind=8), intent(in) :: delta_x
   real(kind=8), intent(inout) :: Boost_Lab(0:3,0:3), Boost_COM(0:3,0:3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Internal Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   real(kind=8) :: prob
   integer(kind=4) :: if1
   integer(kind=4) idex
   integer(kind=4) :: itemp
   integer(kind=4) :: mask6, mask10
   integer(kind=4) :: k, l_f, iss
   integer(kind=4) :: i
   real(kind=8) :: e_i, e_f, ex_i, ex_f
   real(kind=8) :: xI_i, xI_f
   real(kind=8) :: xl_i, xl_f
   real(kind=8) :: xj_i_min, xj_i, xL_ang
   real(kind=8) :: xj_f_min, xj_f, xip_f
   real(kind=8) :: spin_proj, spin_eject, spin_target_4
   integer(kind=4) :: L_ang, max_L
   real(kind=8) :: factor
   real(kind=8) :: shift
   real(kind=8) :: x, x1, check, sum, ran

!   real(kind=8) :: t, costhp, sinthp, tanth, thetap, phi, mass_i, mass_f, mass
!   real(kind=8) :: xkpxc, xk_0xc, xKKp_0xc, gamma, theta_f, phi_f

   real(kind=8) :: theta_0, phi_0
   real(kind=8) :: T_1, T_2, mass_1, mass_2, theta, phi
   real(kind=8) :: T_L, theta_L, phi_L

   real(kind=8) :: tally_prob
   real(kind=8) :: xnnn
   real(kind=8) :: alf, bet

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   External Functions
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!   real(kind=8) :: ran3
   real(kind=8) :: random_64
   real(kind=8) :: racahr
!   real(kind=8) :: Legendre
   real(kind=8) :: poly

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Start Program
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

   tally_prob = 1.0d0

!   e_i = nucleus(icomp_i)%e_grid(nbin_i)
   ex_i = e_i
   xI_i = real(Ix_i,kind=8) + real(nucleus(icomp_i)%jshift,kind=8)
   spin_proj = real(particle(iproj)%spin,kind=8)
   xl_i = l_i
   xj_i_min = abs(xl_i-spin_proj)
   xj_i = xj_i_min + is_i
   spin_target_4 = real(spin_target,kind=8)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  First, check which nucleus it decays to
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   prob = random_64(iseed)                !  starting probability

   if1 = 0
   icomp_f = icomp_i
   nbin_f = -1
   if(Channel(l_i,is_i,Ix_i)%num_decay == 0)return   ! no way to decay, hung up

   if(biased_sampling)then
      do if1 = 1, Channel(l_i,is_i,Ix_i)%num_decay
         if(prob <= Channel(l_i,is_i,Ix_i)%Channel_prob(if1))exit
      end do
   else
      xnnn = Channel(l_i,is_i,Ix_i)%num_decay
      if1 = int(random_64(iseed)*xnnn) + 1
      tally_prob = Channel(l_i,is_i,Ix_i)%Channel_prob(if1)
      if(if1 > 1) tally_prob = tally_prob - Channel(l_i,is_i,Ix_i)%Channel_prob(if1-1)
      tally_prob = tally_prob*xnnn
!      write(44,'(''++++++++++++++++++++++++++++++++++++++++++++++++'')')
!      write(44,*)if1
!      write(44,'(5(1x,i10))')l_i,is_i,Ix_i,if1,Channel(l_i,is_i,Ix_i)%num_decay
!      write(44,'(10(1x,f20.15))')(Channel(l_i,is_i,Ix_i)%Channel_prob(i),i=1,Channel(l_i,is_i,Ix_i)%num_decay)
!      write(44,*)'tally_prob = ',tally_prob
   end if

   if(Channel(l_i,is_i,Ix_i)%decay_particle(if1) == -10)then      !!  Fission
      icomp_f = icomp_i
      nbin_f = 0
      num_part = num_part + 1
      part_data(1,num_part) = real(icomp_f,kind=8)
      part_data(2,num_part) = -2.0d0        !  one signal this is a fission event
      part_data(3,num_part) = xI_i
      part_data(4,num_part) = 0.0d0
      part_data(5,num_part) = real(nbin_f,kind=8)      !  Other signal this is a fission event
      part_data(6,num_part) = 0.0d0
      part_data(7,num_part) = 0.0d0
      part_data(19,num_part) = tally_prob
      part_data(20,num_part) = -10.0d0
      part_data(21,num_part) = icomp_i
      return
   end if

   icomp_f = Channel(l_i,is_i,Ix_i)%decay_to(if1)
   k = Channel(l_i,is_i,Ix_i)%decay_particle(if1)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  Next, find which state it decays to
!--------  Establish probability, then search for it in list
!--------  using bisection between upper and lower
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if(Channel(l_i,is_i,Ix_i)%Channel_decay(if1)%num_decay < 1)then
      write(6,*)'Primary Decay'
      write(6,*)l_i, is_i, Ix_i, if1
   end if


   prob = random_64(iseed)
   call find_prob(Channel(l_i,is_i,Ix_i)%Channel_decay(if1)%num_decay,     &
                  Channel(l_i,is_i,Ix_i)%Channel_decay(if1)%decay_prob,    &
                  prob, idex)




!   num_dec(k) = num_dec(k) + 1

   num_part = num_part +1
   if(num_part > dim_part) stop 'num_part > dim_part'

   mask6 = 2**6 - 1
   mask10 = 2**9 - 1

   itemp = Channel(l_i,is_i,Ix_i)%Channel_decay(if1)%decay_list(idex)

   call unpack_data(Ix_f, ip_f, nbin_f, idb, l_f, iss, itemp)

!   Ix_f = iand(itemp,mask6)
!   ip_f = iand(ishft(itemp,-6),1)
!   nbin_f = iand(ishft(itemp,-7),mask10)
!   idb = iand(ishft(itemp,-17),1)
!   l_f = iand(ishft(itemp,-18),mask6)
!   iss = iand(ishft(itemp,-24),mask6)

   xl_f = l_f
   xI_f = real(Ix_f,kind=8) + real(nucleus(icomp_f)%jshift,kind=8)
   xj_f_min = abs(real(l_f,kind=8) - real(particle(k)%spin,kind=8)) 
   xj_f = real(iss,kind=8) + xj_f_min
   xip_f = real(2*ip_f - 1,kind=8)

   if(idb == 0)then
      e_f = ex_i - nucleus(icomp_i)%sep_e(k) -           &
                   nucleus(icomp_f)%e_grid(nbin_f)
      if(e_f < -0.5d0*de)stop 'e_f < -de/2 in MC_primary_decay (3)'
      if(e_f < 0.0d0)then
         e_f = (e_f + 0.5d0*de)*random_64(iseed)
         ex_f = ex_i - e_f - nucleus(icomp_i)%sep_e(k)
      elseif(e_f <= de)then           !   decays to same bin
         e_f = e_f*random_64(iseed)
         ex_f = ex_i - e_f - nucleus(icomp_i)%sep_e(k)
      elseif(e_f > de)then
         shift = 0.5d0*de*(2.0d0*random_64(iseed) - 1.0d0)
         ex_f = nucleus(icomp_f)%e_grid(nbin_f) + shift
         e_f = ex_i - nucleus(icomp_i)%sep_e(k) - ex_f
      end if
   else  
      e_f = ex_i - nucleus(icomp_i)%sep_e(k) -         &
                   nucleus(icomp_f)%state(nbin_f)%energy
      ex_f = nucleus(icomp_f)%state(nbin_f)%energy
      if(e_f < 0.0d0)stop 'A problem arose with a decay to a discrete state with e_f < 0.0d0 in MC_primary_decay'
   end if


!   write(41,'(3(1x,i5),1x,f5.1,1x,i5,1x,i3,3(1x,i5),2(1x,f5.1),i5)')                  &
!      icomp_i, l_i, iX_i, xI_i, nbin_i, k, icomp_f, l_f, Ix_f, xI_f, xip_f, nbin_f


   part_data(1,num_part) = real(icomp_f,kind=8)
   part_data(2,num_part) = real(k,kind=8)
   part_data(3,num_part) = xI_f
   part_data(4,num_part) = xip_f
!   part_data(4,num_part) = real(2*ip_f-1,kind=8)
   part_data(5,num_part) = real(nbin_f,kind=8)
   part_data(6,num_part) = real(idb,kind=8)
   part_data(7,num_part) = real(l_f,kind=8)
   if( k == 0)then
      part_data(8,num_part) = real(iss,kind=8)
   else
      part_data(8,num_part) = xj_f
   end if

   part_Ang_data(0:0:Ang_L_max,num_part) = 0.0d0

   part_Ang_data(0,num_part) = 0.5d0
   ran = random_64(iseed)

   if(k == 0)then
      x = 2.0d0*ran - 1.0d0           !  Assume gammas to be isotropic (for now).
   else                               !  particles - compute
      max_L = min(2*l_i, 2*l_f, Ang_L_max)
      spin_eject = real(particle(k)%spin,kind=8)
      factor = (-1.0d0)**(xI_f-spin_eject-spin_target_4+spin_proj)*      &
               (2.0d0*xI_i+1.0d0)*                                       &
               (2.0d0*xj_i+1.0d0)*(2.0d0*xl_i+1.0d0)*                    &
               (2.0d0*xj_f+1.0d0)*(2.0d0*xl_f+1.0d0)*0.5d0
      xl_f = l_f
      do L_ang = 2, max_L, 2
         xL_ang = L_ang
         part_Ang_data(L_ang,num_part) =                             &
            factor*                                                  &
            clb_l(L_ang, l_i)*clb_l(L_ang,l_f)*                      &
            racahr(xI_i,xj_i,xI_i,xj_i,spin_target_4,xl_Ang)*        &
            racahr(xj_i,xj_i,xl_i,xl_i,xl_Ang,spin_proj)*            &
            racahr(xI_i,xj_f,xI_i,xj_f,xI_f,xl_Ang)*                 &
            racahr(xj_f,xj_f,xl_f,xl_f,xl_Ang,spin_eject)
      end do

      alf = 0.0d0
      bet = 0.0d0
      check = 0.0d0
      x = -1.0d0
      do i = 1, ixx_max
         x = real(i,kind=8)*delta_x - 1.0d0
         x1 = x - delta_x
         sum = 0.5d0*(poly(0,1,alf,bet,x) + poly(0,1,alf,bet,x1))
         do L_ang = 2, max_L,2
            sum = sum + part_Ang_data(L_ang,num_part)*                          &
                       (poly(L_ang,1,alf,bet,x) + poly(L_ang,1,alf,bet,x1))
         end do
         check = check + 0.5d0*sum*delta_x
         x = x - random_64(iseed)*delta_x*0.9999999d0
         if(check >= ran)exit
      end do
      if(abs(x) > 1.0d0)then
         write(6,*) 'cos(theta) =',x,' wrong in primary decay'
         stop
      end if
      if(abs(x) < -1.0d0)then
         write(6,*) 'cos(theta) =',x,' wrong in primary decay'
         stop
      end if
   end if


   theta_0 = acos(x)

   phi_0 = two_pi*random_64(iseed)

   mass_1 = particle(k)%Mass
   if(idb == 0)then
      mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%e_grid(nbin_f)
   else
      mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%state(nbin_f)%energy
   end if

!   write(6,*)'theta,phi', e_f, theta_0, phi_0
   part_data(9,num_part) = e_f
   part_data(10,num_part) = theta_0
   part_data(11,num_part) = phi_0

!   if(k == 2 .and. T_1 >= 3.40001d0 .and. T_1 <= 3.6d0)then
!       write(61,*)'before ',e_f, x
!   end if
!   if(k == 2 .and. T_1 >= 3.50001d0 .and. T_1 <= 3.6d0)then
!       write(62,*)'before ',e_f, x
!   end if



   call Boost_frame(e_f, mass_1, mass_2, theta_0, phi_0,                &
                    Boost_Lab, Boost_COM, T_1, theta, phi,              &
                    T_2, T_L, theta_L, phi_L)

   part_data(12,num_part) = T_1
   part_data(13,num_part) = theta
   part_data(14,num_part) = phi
   part_data(15,num_part) = T_L
   part_data(16,num_part) = theta_L
   part_data(17,num_part) = phi_L
   part_data(18,num_part) = T_2
   part_data(19,num_part) = tally_prob
   part_data(20,num_part) = ex_f
   part_data(21,num_part) = icomp_i
   nucleus(icomp_f)%Kinetic_Energy = T_2


!   if(k == 2 .and. T_L >= 3.40001d0 .and. T_L <= 3.5d0)then
!       write(61,*)e_f, T_1, T_L, rrr, x, cos(theta), cos(theta_L)
!   end if
!   if(k == 2 .and. T_L >= 3.50001d0 .and. T_L <= 3.6d0)then
!       write(62,*)e_f, T_1, T_L, rrr, x, cos(theta), cos(theta_L)
!   end if


   return

   end subroutine MC_primary_decay
!
!*******************************************************************************
!
subroutine force_decay(icomp_i, nbin_i,                              &
                       icomp_f, Ix_f, ip_f, nbin_f, idb,             &
                       n_dat, dim_part, num_part, part_data,         &
                       Ang_L_max, part_Ang_data,                     &
                       Boost_Lab, Boost_COM)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine forces a decay from a continuous energy bin where
!    the probability to decay is too small, and is hung up. That is it
!    has no way to decay to a discrete state. This is a so-called
!    comoutaitonal isomer                                           
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
   use useful_data
   implicit none
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Passed Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!-----------  Input  ------------------------------------------------
!   integer(kind=4) :: iseed
   integer(kind=4), intent(in) :: icomp_i, nbin_i 
!-----------  Output ------------------------------------------------
   integer(kind=4), intent(out) :: icomp_f, Ix_f, ip_f, nbin_f, idb
   integer(kind=4), intent(in) :: n_dat, dim_part
   integer(kind=4), intent(inout) :: num_part
   real(kind=8), intent(inout) :: part_data(n_dat,dim_part)
   integer(kind=4), intent(in) :: Ang_L_max
   real(kind=8), intent(inout) :: part_Ang_data(0:Ang_L_max,dim_part)
   real(kind=8), intent(inout) :: Boost_Lab(0:3,0:3), Boost_COM(0:3,0:3)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Internal Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   integer(kind=4) :: k, iss
   real(kind=8) :: e_i, e_f, ex_i, ex_f
   real(kind=8) :: xnstate

   real(kind=8) :: costhp, sinthp, theta_0, phi_0
   real(kind=8) :: T_1, T_2, mass_1, mass_2, theta, phi
   real(kind=8) :: T_L, theta_L, phi_L

   real(kind=8) :: tally_prob

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   External Functions
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   real(kind=8) :: random_64
!---------------------------------------------------------------------------

   tally_prob = 1.0d0

   icomp_f = icomp_i

   e_i = nucleus(icomp_f)%e_grid(nbin_i)
   ex_i = e_i
   if(num_part >= 1)ex_i = part_data(20,num_part)

   num_part = num_part + 1

   k = 0                                           !  gamma decay
!--------------   find final state randomly   -----------------------
   xnstate = real(nucleus(icomp_f)%ncut,kind=8)

   nbin_f = int(xnstate*random_64(iseed)) + 1
!--------------------------------------------------------------------
   ex_f = nucleus(icomp_f)%state(nbin_f)%energy

   e_f = ex_i - ex_f

   Ix_f = nint(nucleus(icomp_f)%state(nbin_f)%spin)
   ip_f=iabs(nint((nucleus(icomp_f)%state(nbin_f)%parity+1.)/2.))
   idb = 1
   iss = 0
   
   part_data(1,num_part) = real(icomp_f,kind=8)
   part_data(2,num_part) = real(k,kind=8)
   part_data(3,num_part) = real(Ix_f,kind=8)
   part_data(4,num_part) = real(ip_f,kind=8)
   part_data(5,num_part) = real(nbin_f,kind=8)
   part_data(6,num_part) = real(idb,kind=8)
   part_data(7,num_part) = -1.0d0
   part_data(8,num_part) = real(iss,kind=8)
!   part_data(9,num_part) = t

   costhp = 2.0d0*random_64(iseed) - 1.0d0
   sinthp = sqrt(1.0d0-costhp**2)
   theta_0 = acos(costhp)
   phi_0 = two_pi*random_64(iseed)

   mass_1 = particle(k)%Mass
   if(idb == 0)then
      mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%e_grid(nbin_f)
   else
      mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%state(nbin_f)%energy
   end if

   call Boost_frame(e_f, mass_1, mass_2, theta_0, phi_0,                &
                    Boost_Lab, Boost_COM, T_1, theta, phi,              &
                    T_2, T_L, theta_L, phi_L)

   part_data(9,num_part) = e_f
   part_data(10,num_part) = theta_0
   part_data(11,num_part) = phi_0
   part_data(12,num_part) = T_1
   part_data(13,num_part) = theta
   part_data(14,num_part) = phi
   part_data(15,num_part) = T_L
   part_data(16,num_part) = theta_L
   part_data(17,num_part) = phi_L
   part_data(18,num_part) = T_2
   part_data(19,num_part) = tally_prob
   part_data(20,num_part) = ex_f
   part_data(21,num_part) = icomp_i
   nucleus(icomp_f)%Kinetic_Energy = T_2

   if(.not.pop_calc)part_Ang_data(0,num_part) = 0.5e0

end subroutine force_decay
!
!*******************************************************************************
!
subroutine Boost_frame(e_f, mass_1, mass_2, theta_0, phi_0,                   &
                       Boost_Lab, Boost_COM, T_1, theta, phi,                 &
                       T_2, T_L, theta_L, phi_L)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates and apply relativistic boost to the 
!    current decay frame. Decays occur in rest frame of compound nucleus, 
!    therefore, boost to Center of momentum or Lab frame
!    returns T_1 = kinetic emergy of emitted fragment
!    theta - angle measured relative to z-axis, or initial beam direction
!    phi, azimuthal angle, which is not tracked                                          
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
!-------------    Data passed in
   real(kind=8), intent(in) :: e_f, mass_1, mass_2, theta_0, phi_0
!-------------    Data passed out
   real(kind=8), intent(inout) :: Boost_lab(0:3,0:3), Boost_COM(0:3,0:3)       !   Boost matrices are rewritten
   real(kind=8), intent(inout) :: T_1, theta, phi, T_2, T_L, theta_L, phi_L
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------------    Inteneral Data
   real(kind=8) :: E_T, EE, pp
   real(kind=8) :: p_1(0:3), P_2(0:3), Lor(0:3,0:3), Temp(0:3,0:3), v_2(1:3)
   real(kind=8) :: ptemp(0:3)
   real(kind=8) :: cos_theta, sin_theta, cos_phi, sin_phi
   real(kind=8) :: beta, gamma 
   real(kind=8) :: e_f1
   integer(kind=4) :: i, j, k
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------------     Start Calculation      -----------------

!---------------------------------    e_f is energy released in decay

   e_f1 = max(e_f, 1.0d-5)

   E_T = e_f1 + mass_1 + mass_2                                        !  Total energy
   EE = sqrt(E_T**2 - (mass_1**2 + mass_2**2))
   pp = sqrt((EE**4 - 4.0d0*mass_1**2*mass_2**2)/(4.0d0*E_T**2))       !  momentum


   ptemp(0) = sqrt(pp**2 + mass_1**2)                             !  momentum four-vector of emitted particle
   ptemp(1) = pp*cos(theta_0)
   ptemp(2) = pp*sin(theta_0)*cos(phi_0)
   ptemp(3) = pp*sin(theta_0)*sin(phi_0)

   p_2(0) = sqrt(pp**2 + mass_2**2)                               !  momentum four-vector of residual nucleus
   p_2(1) = -ptemp(1)                                             !  momentum conservation - opposite direction
   p_2(2) = -ptemp(2)
   p_2(3) = -ptemp(3)
!------------------------------------------------------   Compute quantities in COM and Lab frame
!------------------------------------------------------   Lorentz transform momentum to COM frame
   do i = 0, 3
      p_1(i) = 0.0d0
      do j = 0, 3
         p_1(i) = p_1(i) + Boost_COM(i,j)*ptemp(j)
      end do
!      write(6,*)'p_1',i,p_1(i)
   end do
!---------------------------   Kinetic energies in COM frame
   T_1 = p_1(0) - mass_1                                 !   Kinetic energy of emitted particle
   T_2 = p_2(0) - mass_2                                 !   Kinetic energy of emitted particle
   pp = 0.0d0
   do i = 1, 3
      pp = pp + p_1(i)**2                                !  magnitude of vector momentum
   end do
   pp = sqrt(pp)
   cos_theta = 0.0d0
   sin_theta = 1.0d0
   phi = 0.0d0
   cos_phi = cos(phi)
   sin_phi = sin(phi)
!--------------------------   Compute theta and phi in COM
   cos_theta = p_1(1)/pp
   sin_theta = sqrt(1.0d0 - cos_theta**2)
   theta = acos(cos_theta)

   if(sin_theta > 0.0d0)then
      cos_phi = p_1(2)/(sin_theta*pp)
      sin_phi = p_1(3)/(sin_theta*pp)
      phi = acos(cos_phi)
      if(sin_phi < 0.0)phi = 2.0d0*pi - phi
   end if

!---------------------------------------------------------     Lorentz transformation to the Lab frame
   do i =0, 3
      ptemp(i) = p_1(i)
   end do
   do i = 0, 3
      p_1(i) = 0.0d0
      do j = 0, 3
         p_1(i) = p_1(i) + Boost_Lab(i,j)*ptemp(j)
      end do
   end do
!---------------------------   Kinetic energies in Lab frame
   T_L = p_1(0) - mass_1                                 !   Kinetic energy of emitted particle
   pp = 0.0d0
   do i = 1, 3
      pp = pp + p_1(i)**2                                !  magnitude of vector momentum
   end do
   pp = sqrt(pp)
   cos_theta = 0.0d0
   sin_theta = 1.0d0
   phi_L = 0.0d0
!--------------------------   Compute theta and phi in Lab
   cos_theta = p_1(1)/pp
   sin_theta = sqrt(1.0d0 - cos_theta**2)
   theta_L = acos(cos_theta)
   if(sin_theta > 0.0d0)then
      cos_phi = p_1(2)/(sin_theta*pp)
      sin_phi = p_1(3)/(sin_theta*pp)
      phi_L = acos(cos_phi)
      if(sin_phi < 0.0)phi_L = 2.0d0*pi - phi_L
   end if

!------------------------------------    Update Lorentz transformation to COM frame for next decay
!------------------------------------    Velocity of residual nucleus in units of c, i.e., V_2(i) = beta(i)
   do i = 1, 3
      v_2(i) = p_2(i)/mass_2
   end do
   pp = 0.0d0
   do i = 1, 3
      pp = pp + p_2(i)**2
   end do
   pp = sqrt(pp)
   beta = pp/mass_2
   gamma = 1.0d0/sqrt(1.0d0 - beta**2)
!---------------------------    Lorentz transformation for the residual nucleus
   Lor(0,0) = gamma
   do i = 1, 3
      Lor(0,i) = -gamma*v_2(i)
      Lor(i,0) = Lor(0,i)
   end do
   do i =1, 3
      do j = i, 3
         if(i == j)then
            Lor(i,j) = 1.0d0 + (gamma-1.0d0)*v_2(i)**2/beta**2
         else
            Lor(i,j) = (gamma - 1.0d0)*V_2(i)*v_2(j)/beta**2
            Lor(j,i) = Lor(i,j)
         end if
      end do
   end do
!--------------------------   Update Boost for next decay Boost_COM = Boost_COM*Lor
   do i = 0, 3
      do j = 0, 3
         Temp(i,j) = 0.0d0
         do k = 0, 3
            Temp(i,j) = Temp(i,j) + Boost_COM(i,k)*Lor(k,j)
         end do
      end do
   end do
   do i = 0, 3
      do j = 0, 3
         Boost_COM(i,j) = Temp(i,j)
      end do
   end do
   return
end subroutine Boost_frame
   
