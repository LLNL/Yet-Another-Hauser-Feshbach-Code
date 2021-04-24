!
!*******************************************************************************
!
   subroutine MC_decay_bin(icomp_i, Ix_i, ip_i, nbin_i,                &
                           icomp_f, Ix_f, ip_f, nbin_f, idb,           &
                           n_dat, dim_part, num_part_type, part_fact,  &
                           num_part, part_data,                        &
                           Ang_L_max, part_Ang_data,                   &
                           num_theta, extra_angle_data)
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
   integer(kind=4), intent(inout) :: num_part_type(0:6)
   real(kind=8), intent(inout) :: part_fact(0:7)
   integer(kind=4), intent(in) :: n_dat, dim_part
   integer(kind=4), intent(inout) :: num_part
   real(kind=8), intent(inout) :: part_data(n_dat,dim_part)
   integer(kind=4), intent(in) :: Ang_L_max
   real(kind=8), intent(out) :: part_Ang_data(0:Ang_L_max,dim_part)
   integer(kind=4), intent(in) :: num_theta
   real(kind=8), intent(inout) :: extra_angle_data(3*num_theta,dim_part)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Internal Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   real(kind=8) :: prob
   integer(kind=4) idex
   integer(kind=4) :: itemp
   integer(kind=4) :: mask6, mask10
   integer(kind=4) :: k, l, iss
   real(kind=8) :: e_i, e_f, ex_i, ex_f
   real(kind=8) :: xI_i, xI_f, xj_f, xj_f_min, xip_f, xip_i

   integer(kind=4) :: if1

   real(kind=8) :: base_prob, check_prob, tally_norm


   real(kind=8) :: costhp, theta_0, phi_0
   integer(kind=4) :: nang

   real(kind=8) :: tally_prob
   real(kind=8) :: xnnn

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   External Functions
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   real(kind=8) :: random_64
   real(kind=8) :: random_32
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
!-----   part_data(22,num_part) = xI_i         !--------   Spin of initial state
!-----   part_data(23,num_part) = xip_i        !--------   parity of final state
!-----   part_data(24,num_part) = nbin_i       !--------   Bin or state number of initial state

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Start Program
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

   e_i = nucleus(icomp_i)%e_grid(nbin_i)
   ex_i = e_i
   xI_i = real(iX_i,kind=8) + nucleus(icomp_i)%jshift
   xip_i = 2.0d0*real(ip_i,kind=8) - 1.0d0
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  First, check which nucleus it decays to
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   prob = random_64(iseed_64)                !  starting probability
   prob = random_32(iseed_32)                !  starting probability

   if1 = 0
   icomp_f = icomp_i
   nbin_f = 0
   if(nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay == 0)return   ! no way to decay, hung up
!
   tally_prob = 1.0d0

!   part_fact(0:7) = 1.0d0
!   num_part_type(0:6) = 0
!   do i = 1, num_part
!      k = nint(part_data(2,i))
!      if(k > 0 .and. k <= 6)then
!         num_part_type(k) = num_part_type(k) + 1
!         if(num_part_type(k) >= max_particle(k))part_fact(k) = 0.0d0
!      end if
!   end do

   xnnn = 0.0d0
   tally_norm = 0.0d0
   do if1 = 1, nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
      k = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1)
      tally_norm = tally_norm + nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(if1)*part_fact(k)
      xnnn = xnnn + part_fact(k)
   end do
   tally_norm = 1.0d0/tally_norm

   if(nint(xnnn) == 0)return   ! no way to decay, hung up


   if(biased_sampling)then
      check_prob = 0.0d0
      do if1 = 1, nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
         k = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1)
         check_prob = check_prob + nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(if1)*part_fact(k)
         if(prob <= check_prob)exit
      end do
   else
      base_prob = 1.0d0/real(xnnn,kind=8)
      check_prob = 0.0d0
      do if1 = 1, nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
         k = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1)
         check_prob = check_prob + base_prob*part_fact(k)
!         check_prob = check_prob + base_prob*part_fact(k)*tally_norm
         if(prob <= check_prob)exit
      end do
      tally_prob = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(if1)*xnnn*tally_norm
   end if
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  set icomp_f, note if1 > num_decay ->  Fission 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if(nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1) == 7)then     !  Fission 
      icomp_f = icomp_i
      nbin_f = -nbin_i
      num_part = num_part + 1
      part_data(1,num_part) = real(icomp_f,kind=8)
      part_data(2,num_part) = 7.0d0        !  one signal this is a fission event
      part_data(3,num_part) = xI_i
      part_data(4,num_part) = 0.0d0
      part_data(5,num_part) = real(nbin_f,kind=8) 
      part_data(6,num_part) = 0.0d0
      part_data(7,num_part) = 0.0d0
      part_data(19,num_part) = tally_prob
      part_data(20,num_part) = -10.0d0     !  Other signal this is a fission event
      part_data(21,num_part) = icomp_i
      part_data(22,num_part) = xI_i
      part_data(23,num_part) = xip_i
      part_data(24,num_part) = real(nbin_i,kind=8)
      return
   end if

   icomp_f = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_to(if1)

   if(icomp_f < 1)then
! write(6,*)icomp_i, Ix_i,ip_i,nbin_i
! write(6,*)'num_decay = ',nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
      write(6,*)xnnn, base_prob
      do k = 1, 6
         write(6,*)num_part_type(k),max_particle(k)
         write(6,*)k,part_fact(k)
      end do
      check_prob = 0.0d0
      do if1 = 1, nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%num_decay
         k = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1)
         check_prob = check_prob + base_prob*part_fact(k)*tally_norm
      write(6,*)'if1 = ',if1,k,check_prob,prob
         if(prob <= check_prob)exit
      end do
      tally_prob = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%HF_prob(if1)/base_prob
      write(6,*)tally_prob

   end if

   if(icomp_f < 1)then
      write(6,*)'icomp_f < 1 after attempting to decay'
      write(6,*)'iproc = ',iproc
      call MPI_Abort(icomm,201,ierr)
   end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  Next, find which state it decays to
!--------  Establish probability, then search for it in list
!--------  using bisection between upper and lower
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if(nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%num_decay < 1)then
      write(6,*)'Bin Decay'
      write(6,*)icomp_i,Ix_i,ip_i,nbin_i,if1
   end if


   prob = random_64(iseed_64)

   call find_prob(nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%num_decay,     &
                  nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%decay_prob,    &
                  prob, idex)

   k = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%decay_particle(if1)

   if(k >= 0 .and. k <= 6)then
      num_part_type(k) = num_part_type(k) + 1
      if(k > 0 .and. num_part_type(k) >= max_particle(k))part_fact(k) = 0.0d0
   end if

   mask6 = 2**6 - 1
   mask10 = 2**9 - 1

   itemp = nucleus(icomp_i)%bins(Ix_i,ip_i,nbin_i)%nuke_decay(if1)%decay_list(idex)

   call unpack_data(Ix_f, ip_f, nbin_f, idb, l, iss, itemp)

   if(nbin_f == 0)then
      write(6,*)'Error nbin_f = 0 in MC_decay_bin'
      write(6,*)nucleus(icomp_i)%A, Ix_i, ip_i, nbin_i
      write(6,*)idb,l,k
      write(6,*)'iproc = ',iproc
      call MPI_Abort(icomm,201,ierr)
   end if


   xI_f = real(Ix_f) + nucleus(icomp_f)%jshift
   xj_f_min = abs(dfloat(l) - particle(k)%spin) 
   xj_f = real(iss) + xj_f_min
   xip_f = 2.0d0*real(ip_f,kind=8) - 1

   if(idb == 0)then
      e_f = ex_i - nucleus(icomp_i)%sep_e(k) -           &
                   nucleus(icomp_f)%e_grid(nbin_f)
      ex_f = nucleus(icomp_f)%e_grid(nbin_f)
      if(e_f < -0.5d0*de)then
         write(6,*)'icomp_i = ',icomp_i,' icomp_f= ',icomp_f
         write(6,*)'k= ', k
         write(6,*)'nbin_i = ',nbin_i,' nbin_f = ',nbin_f
         write(6,*)'e_f = ',e_f
         write(6,*)'ex_i = ',ex_i
         write(6,*)'Separation energy = ', nucleus(icomp_i)%sep_e(k)
         write(6,*)'Initial Bin energy = ', nucleus(icomp_i)%e_grid(nbin_i)
         write(6,*)'Final Bin energy = ', nucleus(icomp_f)%e_grid(nbin_f)
         write(6,*)'Energy bin width = ', de
         write(6,*)'num_part = ', num_part
         write(6,*)'e_f < -de/2 in MC_decay_bin (1)'
         write(6,*)'iproc = ',iproc
         call MPI_Abort(icomm,201,ierr)
      end if
   else  
      e_f = ex_i - nucleus(icomp_i)%sep_e(k) -         &
                   nucleus(icomp_f)%state(nbin_f)%energy
      ex_f = nucleus(icomp_f)%state(nbin_f)%energy

      if(e_f < 0.d0*de)then
         write(6,*)'icomp_i = ',icomp_i,' icomp_f= ',icomp_f
         write(6,*)'k= ', k
         write(6,*)'nbin_i = ',nbin_i,' nbin_f = ',nbin_f
         write(6,*)'e_f = ',e_f
         write(6,*)'ex_i = ',ex_i
         write(6,*)'Separation energy = ', nucleus(icomp_i)%sep_e(k)
         write(6,*)'Initial Bin energy = ', nucleus(icomp_i)%e_grid(nbin_i)
         write(6,*)'Final State energy = ', nucleus(icomp_f)%state(nbin_f)%energy
         write(6,*)'Energy bin width = ', de
         write(6,*)'num_part = ', num_part
         write(6,*)'e_f < 0.0 in MC_decay_bin (1)'
         write(6,*)'iproc = ',iproc
         call MPI_Abort(icomm,201,ierr)
      end if
   end if

   num_part = num_part +1
   if(num_part > dim_part)then
      write(6,*)'num_part > dim_part'
      write(6,*)'iproc = ',iproc
      call MPI_Abort(icomm,201,ierr)
   end if

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

   theta_0 = 0.0d0
   phi_0 = 0.0d0

   if(.not. xs_only)then
      do nang = 1, num_theta
!         costhp = 2.0d0*random_64(iseed_64) - 1.0d0
         costhp = 2.0d0*random_32(iseed_32) - 1.0d0
         if(abs(costhp) > 1.0d0)then
            write(6,*)'cos(theta) wrong in MC_decay_bin'
            call MPI_Abort(icomm, 101, ierr)
         end if
         extra_angle_data(nang,num_part) = acos(costhp)
      end do
      theta_0 = extra_angle_data(1,num_part)
      phi_0 = two_pi*random_32(iseed_32)
!     phi_0 = two_pi*random_64(iseed_64)
   end if

   part_data(9,num_part)  = e_f
   part_data(10,num_part) = theta_0
   part_data(11,num_part) = phi_0
!   part_data(12,num_part) = T_1
!   part_data(13,num_part) = theta
!   part_data(14,num_part) = phi
!   part_data(15,num_part) = T_L
!   part_data(16,num_part) = theta_L
!   part_data(17,num_part) = phi_L
!   part_data(18,num_part) = T_2
   part_data(19,num_part) = tally_prob
   part_data(20,num_part) = ex_f
   part_data(21,num_part) = icomp_i
   part_data(22,num_part) = xI_i
   part_data(23,num_part) = xip_i
   part_data(24,num_part) = real(nbin_i,kind=8)
!   nucleus(icomp_f)%Kinetic_Energy = T_2
   if(.not.pop_calc)part_Ang_data(0,num_part) = 0.5d0
   return

   end subroutine MC_decay_bin
!
!*******************************************************************************
!
   subroutine MC_decay_state(icomp_i, istate_i,                       &
                             n_dat,dim_part,num_part,part_data,       &
                             Ang_L_max,part_Ang_data,                 &
                             num_theta, extra_angle_data, ichan, in)
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
   integer(kind=4), intent(in) :: num_theta
   real(kind=8), intent(inout) :: extra_angle_data(3*num_theta,dim_part)
   integer(kind=4), intent(in) :: ichan, in

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Internal Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   real(kind=8) :: prob
   integer(kind=4) :: j, istate, n_f, k, icomp_f
   integer(kind=4) :: idb, l, iss
   real(kind=8) :: check
   real(kind=8) :: xI_f
   real(kind=8) :: xj_f, xip_f, e_gamma


   real(kind=8) :: costhp, theta_0, phi_0
   integer(kind=4) :: nang

   real(kind=8) :: tally_prob, tally_weight
   integer(kind=4) :: n

!   real(kind=8) :: mass_i, mass_f, mass
!   real(kind=8) :: xkpxc, xk_0xc, xKKp_0xc, gamma, theta_f, phi_f

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   External Functions
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!   real(kind=8) :: random_64
   real(kind=8) :: random_32

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

   prob = random_32(iseed_32)
!   prob = random_64(iseed_64)
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
      xj_f = nucleus(icomp_i)%state(n_f)%spin
      xip_f = nucleus(icomp_i)%state(n_f)%parity
      iss = 0
      l = 0
      e_gamma = nucleus(icomp_i)%state(istate)%energy - nucleus(icomp_i)%state(n_f)%energy

      num_part = num_part +1
      if(num_part > dim_part)then
         write(6,*)'num_part > dim_part'
         write(6,*)'iproc = ',iproc
         call MPI_Abort(icomm,201,ierr)
      end if

      part_data(1,num_part) = real(icomp_i,kind=8)
      part_data(2,num_part) = real(k,kind=8)
      part_data(3,num_part) = xj_f
      part_data(4,num_part) = xip_f
      part_data(5,num_part) = real(n_f,kind=8)
      part_data(6,num_part) = real(idb,kind=8)
      part_data(7,num_part) = real(l,kind=8)
      part_data(8,num_part) = real(iss,kind=8)

      theta_0 = 0.0d0
      phi_0 = 0.0d0

      if(.not. xs_only)then
         do nang = 1, num_theta
!            costhp = 2.0d0*random_64(iseed_64) - 1.0d0
            costhp = 2.0d0*random_32(iseed_32) - 1.0d0
            if(abs(costhp) > 1.0d0)then
               write(6,*) 'cos(theta) wrong in MC_decay_state'
               call MPI_Abort(icomm, 101, ierr)
            end if
            extra_angle_data(nang,num_part) = acos(costhp)
         end do
         theta_0 = extra_angle_data(1,num_part)
         phi_0 = two_pi*random_32(iseed_32)
!         phi_0 = two_pi*random_64(iseed_64)
      end if

      part_data(9,num_part) = e_gamma
      part_data(10,num_part) = theta_0
      part_data(11,num_part) = phi_0
!      part_data(12,num_part) = T_1
!      part_data(13,num_part) = theta
!      part_data(14,num_part) = phi
!      part_data(15,num_part) = T_L
!      part_data(16,num_part) = theta_L
!      part_data(17,num_part) = phi_L
!      part_data(18,num_part) = T_2
      part_data(19,num_part) = tally_prob
      part_data(20,num_part) = nucleus(icomp_i)%state(n_f)%energy
      part_data(21,num_part) = icomp_i
      part_data(22,num_part) = nucleus(icomp_i)%state(istate)%spin
      part_data(23,num_part) = nucleus(icomp_i)%state(istate)%parity
      part_data(24,num_part) = real(istate,kind=8)
!      nucleus(icomp_f)%Kinetic_Energy = T_2


      if(.not.pop_calc)part_Ang_data(0,num_part) = 0.5d0
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

!   prob = random_64(iseed_64)
   prob = random_32(iseed_32)
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
   if(num_part > dim_part)then
      write(6,*)'num_part > dim_part'
      write(6,*)'iproc = ',iproc
      call MPI_Abort(icomm,201,ierr)
   end if

   part_data(1,num_part) = real(icomp_i,kind=8)
   part_data(2,num_part) = real(k,kind=8)
   part_data(3,num_part) = xI_f
   part_data(4,num_part) = xip_f
   part_data(5,num_part) = real(n_f,kind=8)
   part_data(6,num_part) = real(idb,kind=8)
   part_data(7,num_part) = real(l,kind=8)
   part_data(8,num_part) = real(iss,kind=8)

!   costhp = 2.0d0*random_64(iseed_64) - 1.0d0
   costhp = 2.0d0*random_32(iseed_32) - 1.0d0
   if(abs(costhp) > 1.0d0)then
      write(6,*)'cos(theta) wrong in MC_decay_state #2'
      write(6,*)'iproc = ',iproc
      call MPI_Abort(icomm,201,ierr)
   end if

   theta_0 = acos(costhp)
   phi_0 = two_pi*random_32(iseed_32)
!   phi_0 = two_pi*random_64(iseed_64)


   part_data(9,num_part) = e_gamma
   part_data(10,num_part) = theta_0
   part_data(11,num_part) = phi_0
!   part_data(12,num_part) = T_1
!   part_data(13,num_part) = theta
!   part_data(14,num_part) = phi
!   part_data(15,num_part) = T_L
!   part_data(16,num_part) = theta_L
!   part_data(17,num_part) = phi_L
!   part_data(18,num_part) = T_2
   part_data(19,num_part) = tally_prob
   part_data(20,num_part) = nucleus(icomp_i)%state(n_f)%energy
   part_data(21,num_part) = icomp_i
   part_data(22,num_part) = nucleus(icomp_i)%state(istate)%spin
   part_data(23,num_part) = nucleus(icomp_i)%state(istate)%parity
   part_data(24,num_part) = real(istate,kind=8)
!   nucleus(icomp_f)%Kinetic_Energy = T_2

   if(.not.pop_calc)part_Ang_data(0,num_part) = 0.5d0         !  Isotropic emission

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
   use nodeinfo
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
!  write(6,*)num,p_norm
!   if(p_norm <= 1.0d-20)then
!      write(6,*)'p_norm too small in find_prob = ',p_norm
!      write(6,*)'iproc = ',iproc
!      call MPI_Abort(icomm,201,ierr)
!   end if
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
   use nodeinfo
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
   subroutine MC_primary_decay(iproj,spin_target,                          &
                               l_i, is_i, Ix_i, e_i, icomp_i,              &
                               icomp_f, Ix_f, ip_f, nbin_f, idb,           &
                               n_dat, dim_part, num_part_type, part_fact,  &
                               num_part, part_data,                        &
                               Ang_L_max, part_Ang_data,                   &
                               ixx_max, delta_x, Leg_poly,                 &
                               num_theta, extra_angle_data)
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
   integer(kind=4), intent(in) :: iproj
   real(kind=8), intent(in) :: spin_target
   integer(kind=4), intent(in) :: l_i, is_i, Ix_i, icomp_i
!-----------  Output ------------------------------------------------
   integer(kind=4), intent(out) :: icomp_f, Ix_f, ip_f, nbin_f, idb
   integer(kind=4), intent(in) :: n_dat, dim_part
   integer(kind=4), intent(inout) :: num_part_type(0:6)
   real(kind=8), intent(inout) :: part_fact(0:7)
   integer(kind=4), intent(inout) :: num_part
   real(kind=8), intent(inout) :: part_data(n_dat,dim_part)
   integer(kind=4), intent(in) :: Ang_L_max
   real(kind=8), intent(inout) :: part_Ang_data(0:Ang_L_max,dim_part)
   integer(kind=4), intent(in) :: ixx_max
   real(kind=8), intent(in) :: delta_x
   real(kind=8), intent(in) :: Leg_poly(0:Ang_L_max,0:ixx_max)
   integer(kind=4), intent(in) :: num_theta
   real(kind=8), intent(inout) :: extra_angle_data(3*num_theta,dim_part)
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
   real(kind=8) :: x, sum, ran
!   real(kind=8) :: x, x1, check, sum, ran
!   real(kind=8) :: check
   real(kind=8) :: theta_0, phi_0
   integer(kind=4) :: nang

   real(kind=8) :: tally_prob
   real(kind=8) :: xnnn
!   real(kind=8) :: alf, bet

   real(kind=8) :: ang_prob(0:ixx_max)
   real(kind=8) :: pnorm
   real(kind=8) :: shift

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   External Functions
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   real(kind=8) :: random_32
   real(kind=8) :: random_64
   real(kind=8) :: racahr
!   real(kind=8) :: poly

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Start Program
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

   tally_prob = 1.0d0

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
!   prob = random_64(iseed_64)                !  starting probability
   prob = random_32(iseed_32)                !  starting probability

   if1 = 0
   icomp_f = icomp_i
   nbin_f = 0
   if(Channel(l_i,is_i,Ix_i)%num_decay == 0)return   ! no way to decay, hung up

   if(biased_sampling)then
      do if1 = 1, Channel(l_i,is_i,Ix_i)%num_decay
         if(prob <= Channel(l_i,is_i,Ix_i)%Channel_prob(if1))exit
      end do
   else
      xnnn = Channel(l_i,is_i,Ix_i)%num_decay
      ran = random_32(iseed_32)
      if1 = min(int(ran*xnnn)+1,Channel(l_i,is_i,Ix_i)%num_decay)  !  just in case ran = 1.0
!      if(if1 > Channel(l_i,is_i,Ix_i)%num_decay)then
!         write(6,*)'xnnn, ran, if1 ',xnnn, ran, if1
!         write(6,*)'l_i, is_i, Ix_i ',l_i,is_i,Ix_i
!      end if
!      if1 = int(random_64(iseed_64)*xnnn) + 1
      tally_prob = Channel(l_i,is_i,Ix_i)%Channel_prob(if1)
      if(if1 > 1) tally_prob = tally_prob - Channel(l_i,is_i,Ix_i)%Channel_prob(if1-1)
      tally_prob = tally_prob*xnnn
   end if

   if(Channel(l_i,is_i,Ix_i)%decay_particle(if1) == 7)then      !!  Fission
      icomp_f = icomp_i
      nbin_f = -999
      num_part = num_part + 1
      part_data(1,num_part) = real(icomp_f,kind=8)
      part_data(2,num_part) = 7.0d0        !  one signal this is a fission event
      part_data(3,num_part) = xI_i
      part_data(4,num_part) = 0.0d0
      part_data(5,num_part) = real(nbin_f,kind=8)
      part_data(6,num_part) = 0.0d0
      part_data(7,num_part) = 0.0d0
      part_data(19,num_part) = tally_prob
      part_data(20,num_part) = -10.0d0     !  Other signal this is a fission event
      part_data(21,num_part) = icomp_i
      part_data(22,num_part) = xI_i
      part_data(23,num_part) = 0.0d0
      part_data(24,num_part) = -1.0d0
      return
   end if

   icomp_f = Channel(l_i,is_i,Ix_i)%decay_to(if1)
   k = Channel(l_i,is_i,Ix_i)%decay_particle(if1)

   if(k >= 0 .and. k <= 6)then
      num_part_type(k) = num_part_type(k) + 1
      if(k > 0 .and. num_part_type(k) >= max_particle(k))part_fact(k) = 0.0d0
   end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  Next, find which state it decays to
!--------  Establish probability, then search for it in list
!--------  using bisection between upper and lower
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if(Channel(l_i,is_i,Ix_i)%Channel_decay(if1)%num_decay < 1)then
      write(6,*)'Primary Decay'
      write(6,*)l_i, is_i, Ix_i, if1
   end if


   prob = random_64(iseed_64)
   call find_prob(Channel(l_i,is_i,Ix_i)%Channel_decay(if1)%num_decay,     &
                  Channel(l_i,is_i,Ix_i)%Channel_decay(if1)%decay_prob,    &
                  prob, idex)

   num_part = num_part +1
   if(num_part > dim_part)then
      write(6,*)'num_part > dim_part'
      write(6,*)'iproc = ',iproc
      call MPI_Abort(icomm,201,ierr)
   end if

   mask6 = 2**6 - 1
   mask10 = 2**9 - 1

   itemp = Channel(l_i,is_i,Ix_i)%Channel_decay(if1)%decay_list(idex)

   call unpack_data(Ix_f, ip_f, nbin_f, idb, l_f, iss, itemp)

!   if(Ix_f > max_J_allowed)then
!       write(6,*)Ix_f
!       write(6,*)icomp_f
!       write(6,*)idb
!       write(6,*)nbin_f
!       write(6,*)nucleus(icomp_f)%state(nbin_f)%energy
!       write(6,*)nucleus(icomp_f)%state(nbin_f)%spin
!   end if

   xl_f = l_f
   xI_f = real(Ix_f,kind=8) + real(nucleus(icomp_f)%jshift,kind=8)
   xj_f_min = abs(real(l_f,kind=8) - real(particle(k)%spin,kind=8)) 
   xj_f = real(iss,kind=8) + xj_f_min
   xip_f = real(2*ip_f - 1,kind=8)

   if(idb == 0)then
      e_f = ex_i - nucleus(icomp_i)%sep_e(k) -                           &
                   nucleus(icomp_f)%e_grid(nbin_f)
      ex_f = nucleus(icomp_f)%e_grid(nbin_f)
      if(e_f < 0.0d0)then
         write(6,*)'problem with primary decay: e_f < 0 in MC_primary_decay'
         write(6,*)'iproc = ',iproc
         call MPI_Abort(icomm,201,ierr)
      end if
   else  
      e_f = ex_i - nucleus(icomp_i)%sep_e(k) -                           &
                   nucleus(icomp_f)%state(nbin_f)%energy
      ex_f = nucleus(icomp_f)%state(nbin_f)%energy
      if(e_f < 0.0d0)then
         write(6,*)'A problem arose with a decay to a discrete state with e_f < 0.0d0 in MC_primary_decay'
         write(6,*)'iproc = ',iproc
         call MPI_Abort(icomm,201,ierr)
      end if
   end if

   part_data(1,num_part) = real(icomp_f,kind=8)
   part_data(2,num_part) = real(k,kind=8)
   part_data(3,num_part) = xI_f
   part_data(4,num_part) = xip_f
   part_data(5,num_part) = real(nbin_f,kind=8)
   part_data(6,num_part) = real(idb,kind=8)
   part_data(7,num_part) = real(l_f,kind=8)
   if(k == 0)then
      part_data(8,num_part) = real(iss,kind=8)
   else
      part_data(8,num_part) = xj_f
   end if

   part_Ang_data(0:0:Ang_L_max,num_part) = 0.0d0

   part_Ang_data(0,num_part) = 0.5d0
!   ran = random_64(iseed_64)
!   ran = random_32(iseed_32)

   max_L = 0

   if(k /= 0)then
      max_L = min(2*l_i, 2*l_f, Ang_L_max)
      spin_eject = real(particle(k)%spin,kind=8)
      factor = (-1.0d0)**(xI_f-spin_eject-spin_target_4+spin_proj)*      &
               (2.0d0*xI_i+1.0d0)*                                       &
               (2.0d0*xj_i+1.0d0)*(2.0d0*xl_i+1.0d0)*                    &
               (2.0d0*xj_f+1.0d0)*(2.0d0*xl_f+1.0d0)*0.5d0
      xl_f = l_f
      do L_ang = 2, max_L, 2
         part_Ang_data(L_ang,num_part) = 1.0d-3
         xL_ang = L_ang
         part_Ang_data(L_ang,num_part) =                                 &
            factor*                                                      &
            clb_l(L_ang,l_i)*clb_l(L_ang,l_f)*                           &
            racahr(xI_i,xj_i,xI_i,xj_i,spin_target_4,xL_Ang)*            &
            racahr(xj_i,xj_i,xl_i,xl_i,xL_Ang,spin_proj)*                &
            racahr(xI_i,xj_f,xI_i,xj_f,xI_f,xL_Ang)*                     &
            racahr(xj_f,xj_f,xl_f,xl_f,xL_Ang,spin_eject)
      end do
      ang_prob(0:ixx_max) = 0.0d0
      do i = 1, ixx_max
         sum = 0.5d0*(Leg_poly(0,i) + Leg_poly(0,i-1))
         do L_ang = 2, max_L, 2
            sum = sum + part_Ang_data(L_ang,num_part)*                   &
                       (Leg_poly(L_ang,i) + Leg_poly(L_ang,i-1))
         end do
         ang_prob(i) = ang_prob(i-1) + 0.5d0*sum*delta_x
      end do
   end if
   pnorm = ang_prob(ixx_max)
   ang_prob = ang_prob/pnorm

   theta_0 = 0.0d0
   phi_0 = 0.0d0

   if(.not. xs_only)then
      do nang = 1, num_theta
         ran = random_32(iseed_32)
         if(k == 0)then
            x = 2.0d0*ran - 1.0d0
         else
            call find_prob(ixx_max, ang_prob(1), ran, i)
!            if(abs(x) > 1.0d0)then
!               write(6,*) 'cos(theta) =',x,' wrong in primary decay'
!               call MPI_Abort(icomm,101,ierr)
!            end if
!            if(abs(x) < -1.0d0)then
!               write(6,*) 'cos(theta) =',x,' wrong in primary decay'
!               call MPI_Abort(icomm,101,ierr)
!            end if
            shift = random_32(iseed_32)*delta_x*0.9999999d0
            x = real(i,kind=8)*delta_x - 1.0d0 - shift
         end if
         extra_angle_data(nang,num_part) = acos(x)
      end do
      theta_0 = extra_angle_data(1,num_part)
      phi_0 = two_pi*random_32(iseed_32)
!      phi_0 = two_pi*random_64(iseed_64)
   end if

   part_data(9,num_part) = e_f
   part_data(10,num_part) = theta_0
   part_data(11,num_part) = phi_0


!   part_data(12,num_part) = T_1
!   part_data(13,num_part) = theta
!   part_data(14,num_part) = phi
!   part_data(15,num_part) = T_L
!   part_data(16,num_part) = theta_L
!   part_data(17,num_part) = phi_L
!   part_data(18,num_part) = T_2
   part_data(19,num_part) = tally_prob
   part_data(20,num_part) = ex_f
   part_data(21,num_part) = icomp_i
   part_data(22,num_part) = xI_i
   part_data(23,num_part) = 0.0d0
   part_data(24,num_part) = -1.0d0
!   nucleus(icomp_f)%Kinetic_Energy = T_2


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
                       num_theta, extra_angle_data)
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
   integer(kind=4), intent(in) :: icomp_i, nbin_i 
!-----------  Output ------------------------------------------------
   integer(kind=4), intent(out) :: icomp_f, Ix_f, ip_f, nbin_f, idb
   integer(kind=4), intent(in) :: n_dat, dim_part
   integer(kind=4), intent(inout) :: num_part
   real(kind=8), intent(inout) :: part_data(n_dat,dim_part)
   integer(kind=4), intent(in) :: Ang_L_max
   real(kind=8), intent(inout) :: part_Ang_data(0:Ang_L_max,dim_part)
   integer(kind=4), intent(in) :: num_theta
   real(kind=8), intent(inout) :: extra_angle_data(3*num_theta,dim_part)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Internal Data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
   integer(kind=4) :: k, iss
   real(kind=8) :: e_i, e_f, ex_i, ex_f
   real(kind=8) :: xnstate

   real(kind=8) :: costhp, theta_0, phi_0
   integer(kind=4) :: nang

   real(kind=8) :: tally_prob

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   External Functions
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!   real(kind=8) :: random_64
   real(kind=8) :: random_32
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

!   nbin_f = int(xnstate*random_64(iseed_64)) + 1
   nbin_f = int(xnstate*random_32(iseed_32)) + 1
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

   theta_0 = 0.0d0
   phi_0 = 0.0d0

   if(.not. xs_only)then
      do nang = 1, num_theta
!         costhp = 2.0d0*random_64(iseed_64) - 1.0d0
         costhp = 2.0d0*random_32(iseed_32) - 1.0d0
         extra_angle_data(nang,num_part) = acos(costhp)
      end do
      theta_0 = extra_angle_data(1,num_part)
      phi_0 = two_pi*random_32(iseed_32)
!      phi_0 = two_pi*random_64(iseed_64)
   end if

   part_data(9,num_part) = e_f
   part_data(10,num_part) = theta_0
   part_data(11,num_part) = phi_0
!   part_data(12,num_part) = T_1
!   part_data(13,num_part) = theta
!   part_data(14,num_part) = phi
!   part_data(15,num_part) = T_L
!   part_data(16,num_part) = theta_L
!   part_data(17,num_part) = phi_L
!   part_data(18,num_part) = T_2
   part_data(19,num_part) = tally_prob
   part_data(20,num_part) = ex_f
   part_data(21,num_part) = icomp_i
   part_data(22,num_part) = -1.0d0
   part_data(23,num_part) = 0.0d0
   part_data(24,num_part) = real(nbin_i,kind=8)
!   nucleus(icomp_f)%Kinetic_Energy = T_2

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
!    This subroutine calculates and applies relativistic boost to the 
!    current decay frame. Decays occur in rest frame of compound nucleus, 
!    therefore, boost to Center of momentum or Lab frame
!    returns T_1 = kinetic emergy of emitted fragment
!    theta - angle measured relative to z-axis, or initial beam direction
!    phi, azimuthal angle, which is not tracked                                          
!
!
!    mass_1:    Mass of emitted particle
!    mass_2:    Mass of residual nucleus
!    Boost_Lab: Relativisitic Boost matrix to Lab frame
!    Boost_COM: Relativisitic Boost matrix to COM frame
!    e_f:       Kinetic Energy of emitted particle in rest frame of the compound nucleus
!    theta_0:   Theta of emitted particle in rest frame
!    phi_0:     Azimutal angle of emitted particle in rest frame
!    T_1:       Kinetic Energy of the emitted particle in COM frame
!    theta_1:   Theta of emitted particle in COM frame
!    phi_1:     Azimutal angle of emitted particle in COM frame
!    T_2:       Recoil kinetic energy of the residual nucleus
!    T_L:       Kinetic energy of the emitted particle in the Lab frame
!    theta_L:   Theta of emitted particle in Lab frame
!    phi_L:     Azimutal angle of emitted particle in Lab frame
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
   real(kind=8) :: vtemp, gamma_m1
   real(kind=8) :: ptemp(0:3)
   real(kind=8) :: cos_theta, sin_theta, cos_phi, sin_phi
   real(kind=8) :: beta, gamma 
   real(kind=8) :: e_f1
   integer(kind=4) :: i, j
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------------     Start Calculation      -----------------

!---------------------------------    e_f is energy released in decay

   e_f1 = max(e_f, 1.0d-5)

   E_T = e_f1 + mass_1 + mass_2                                   !  Total energy
   EE = sqrt(E_T**2 - (mass_1**2 + mass_2**2))
   pp = sqrt((EE**4 - 4.0d0*mass_1**2*mass_2**2)/(4.0d0*E_T**2))  !  momentum


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
   p_1 = matmul(Boost_COM,ptemp)

!---------------------------   Kinetic energies in COM frame
   T_1 = p_1(0) - mass_1                                     !   Kinetic energy of emitted particle
   T_2 = p_2(0) - mass_2                                     !   Kinetic energy of residual nucleus
   pp = 0.0d0
   do i = 1, 3
      pp = pp + p_1(i)*p_1(i)                                !  magnitude of vector momentum
   end do
   pp = sqrt(pp)
   cos_theta = 0.0d0
   sin_theta = 1.0d0
   phi = 0.0d0
   cos_phi = 1.0d0
   sin_phi = 0.0d0
!--------------------------   Compute theta and phi in COM
   cos_theta = p_1(1)/pp
   sin_theta = sqrt(1.0d0 - cos_theta**2)
   theta = acos(cos_theta)

   if(sin_theta > 0.0d0)then
      cos_phi = p_1(2)/(sin_theta*pp)
      sin_phi = p_1(3)/(sin_theta*pp)
      phi = acos(cos_phi)
      if(sin_phi < 0.0d0)phi = two_pi - phi
   end if

!---------------------------------------------------------     Lorentz transformation to the Lab frame
   ptemp = p_1
   p_1 = matmul(Boost_Lab,ptemp)
!
   T_L = p_1(0) - mass_1                                 !   Kinetic energy of emitted particle
   pp = 0.0d0
   do i = 1, 3
      pp = pp + p_1(i)*p_1(i)                                !  magnitude of vector momentum
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
      if(sin_phi < 0.0d0)phi_L = two_pi - phi_L
   end if

!------------------------------------    Update Lorentz transformation to COM frame for next decay
!------------------------------------    Velocity of residual nucleus in units of c, i.e., V_2(i) = beta(i)
   pp = 0.0d0
   do i = 1, 3
      v_2(i) = p_2(i)/mass_2
      pp = pp + p_2(i)*p_2(i)
   end do
   pp = sqrt(pp)
   beta = pp/mass_2
   gamma = 1.0d0/sqrt(1.0d0 - beta**2)
!---------------------------    Lorentz transformation for the residual nucleus
   Lor = 0.0d0
   Lor(0,0) = gamma
   do i = 1, 3
      Lor(i,0) = -gamma*v_2(i)
      Lor(0,i) = Lor(i,0)
      Lor(i,i) = 1.0d0
   end do
!   v_2 = v_2/beta
   gamma_m1 = (gamma - 1.0d0)/beta**2
   do j = 1, 3
      vtemp = gamma_m1*v_2(j)
      do i = 1, 3
         Lor(i,j) = Lor(i,j) + v_2(i)*vtemp
!         Lor(i,j) = Lor(i,j) + (gamma - 1.0d0)*V_2(i)*vtemp
      end do
   end do
!--------------------------   Update Boost for next decay Boost_COM = Boost_COM*Lor
   Temp = matmul(Boost_COM,Lor)
   Boost_COM = Temp

   return
end subroutine Boost_frame
   
