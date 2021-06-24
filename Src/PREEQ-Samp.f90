!
!*******************************************************************************
!
subroutine PREEQ_sample(iproj, in, itarget, istate, e_in, ex_tot,      &
                        l_i, is_i, Ix_i, ip_i, icomp_i, icomp_f,       &
                        Ix_f, l_f, ip_f, nbin_f, idb,                  &
                        n_dat, dim_part, num_part_type, part_fact,     &
                        num_part, part_data,                           &
                        num_theta, extra_angle_data)
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a subroutine to Monte Carlo sample the pre-equilibrium decay 
!    of an energy bin
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
!        options
!        nodeinfo
!        print_control
!        useful_data
!        nuclei
!        particles_def
!        pre_equilibrium_no_1
!
!     Subroutines:
!
!        PREEQ_Angular
!
!     External functions:
!
!       integer(kind=4) :: preeq_l
!       integer(kind=4) :: find_ibin
!       logical :: real8_equal
!
!     MPI routines:
!
!        MPI_Abort
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
   use options
   use nodeinfo
   use print_control
   use useful_data
   use nuclei
   use particles_def
   use pre_equilibrium_no_1
   implicit none
!--------------------------------    Argument data
   integer(kind=4), intent(in) :: iproj, in, itarget, istate
   real(kind=8), intent(in) :: e_in, ex_tot
   integer(kind=4), intent(in) :: l_i, is_i, Ix_i, ip_i, icomp_i
   integer(kind=4), intent(out) :: icomp_f, Ix_f, l_f, ip_f, nbin_f, idb
   integer(kind=4), intent(in) :: n_dat, dim_part
   integer(kind=4), intent(inout) :: num_part_type(0:6)
   real(kind=8), intent(inout) :: part_fact(0:7)
   integer(kind=4), intent(out) :: num_part
   real(kind=8), intent(out) :: part_data(n_dat,dim_part)
   integer(kind=4), intent(in) :: num_theta
   real(kind=8), intent(inout) :: extra_angle_data(3*num_theta,dim_part)
!--------------------------------    Internal data
   real(kind=8) :: preeq_cs, preeq_cs_k
   integer(kind=4) :: i, k, kk, m
   integer(kind=4) :: l, l_min, l_max
   real(kind=8) :: xj_f, xI_i, xI_f
   real(kind=8) :: prob, ran, frac, rho_sum
   real(kind=8) :: energy, e_max, ex_final, e_cut
   integer(kind=4) :: nbin_end
   real(kind=8) :: mass_p, mass_e
   integer(kind=4) :: A_p, A_t
   real(kind=8) :: xj_i, xj_f_min, xj_f_max, xI_f_min, xI_f_max
   integer(kind=4) :: ixI_f_min, ixI_f_max
   real(kind=8) :: cpar, cpar2, par, par_f
   real(kind=8) :: ex_min_bin
   real(kind=8) :: x_Ang
   integer(kind=4) :: nang

   real(kind=8) :: e_state
   integer(kind=4) :: closest
   real(kind=8) :: e_closest
   real(kind=8) :: e_diff
   real(kind=8) :: xI_state
   integer(kind=4) :: ip_state

   real(kind=8) :: base_prob, tally_norm


   real(kind=8) :: dee

   real(kind=8) :: theta_0, phi_0

   real(kind=8) :: prob_part(0:6)
   real(kind=8) :: sum_prob
   real(kind=8) :: tally_prob

!--------------------------------    External functions
   real(kind=8) :: random_32
   integer(kind=4) :: preeq_l
   integer(kind=4) :: find_ibin
   logical :: real8_equal
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------------------------    Start Calculation

!----
!----  General procedure is to first select an emitting particle type
!----  usually selected based on total cross section (with unbiased
!----  sampling, each particle is nearly equally probable). Then, select
!----  and energy based on the exciton emission probabilities. Next, 
!----  find an emission angle based on Kalbach-Mann statistics. From
!----  this and the emission energy, estimate the final angular momentum
!----  under the assumption that the react occurs near the nuclear surface.
!----  The angular momentum is then coupled to the particle spin and the
!----  initial angular momentum to find a range of valied angular momenta.
!----  the angular momentum of the final state is then selected randomly
!----  weighted by the density of states for each angular momentum.
!----  For lighter nuclei with lower density of states and high-energy 
!----  pre-equilibrium particle emission, and high initial angular momentum,
!----  there is a chance that the level density is zero for all angular
!----  momentum bins. In this rare case, start at min angular momentum and 
!----  reduce spin and select the first bin with non-zero level density.
!----  If for some reason this fails, EXTREMELY unlikely, reinitiate 
!----  pre-equilibrium search with goto 22
!----  

 22  continue

     tally_prob =1.0d0

     cpar = particle(iproj)%par*nucleus(itarget)%state(istate)%parity
     par = cpar*(-1.0d0)**l_i
     xj_i = abs(real(l_i,kind=8) - particle(iproj)%spin) + real(is_i,kind=8) 

     mass_p = particle(iproj)%mass
     A_p = particle(iproj)%A
     A_t = nucleus(itarget)%A

     xI_i = dfloat(Ix_i) + nucleus(1)%jshift
     par = 2*ip_i - 1
     ran = random_32(iseed_32)
     preeq_cs = nucleus(icomp_i)%PREEQ_cs(in)
     k = -1

     prob = 0.0d0

     if(biased_sampling)then
        prob = 0.0d0
        do kk = 1, nucleus(icomp_i)%num_decay                       !   Loop over particle allowed to decay from this nucleus
           k = nucleus(icomp_i)%decay_particle(kk)
           icomp_f = nucleus(icomp_i)%decay_to(kk)
           prob = prob + nucleus(icomp_i)%PREEQ_part_cs(kk,in)/preeq_cs
           if(ran < prob)exit
        end do
        tally_prob = 1.0d0
     else
        prob_part(0:6) = 0.0d0
        sum_prob = 0.0d0
        tally_norm = 0.0d0
        do kk = 1, nucleus(icomp_i)%num_decay
           k = nucleus(icomp_i)%decay_particle(kk)
           if(nucleus(icomp_i)%PREEQ_part_cs(kk,in)/preeq_cs >= 1.0d-6)then
              prob_part(k) = 1.0d0
              tally_norm = tally_norm + nucleus(icomp_i)%PREEQ_part_cs(kk,in)
              sum_prob = sum_prob + 1.0d0
           end if
        end do
        base_prob = 1.0d0/sum_prob
        tally_norm = 1.0d0/tally_norm
        prob = 0.0d0
        do kk = 1, nucleus(icomp_i)%num_decay                       !   Loop over particle allowed to decay from this nucleus
           if(nucleus(icomp_i)%PREEQ_part_cs(kk,in)/preeq_cs < 1.0-6)cycle
           k = nucleus(icomp_i)%decay_particle(kk)
           icomp_f = nucleus(icomp_i)%decay_to(kk)
           prob = prob + prob_part(k)*base_prob
           if(ran < prob)exit
        end do
        tally_prob = nucleus(icomp_i)%PREEQ_part_cs(kk,in)*tally_norm*sum_prob
     end if

     if(k == -1)stop 'k not set properly in PREEQ_Samp'

     ran = random_32(iseed_32)
     preeq_cs_k = nucleus(icomp_i)%PREEQ_part_cs(kk,in)
     e_max = max(ex_tot - nucleus(icomp_i)%sep_e(k),0.0d0)
     nbin_end = min(int(e_max/de),nucleus(icomp_i)%nbin_part)
!------------------------------------------   Energy of emitted particle
     ex_final = -1.0d0
     prob = 0.0d0
     do m = 0, nbin_end
        energy = dfloat(m)*de
        ex_final = ex_tot - energy - nucleus(icomp_i)%sep_e(k)
        prob = prob + nucleus(icomp_i)%PREEQ_part_spectrum(kk,m)*de
        if(ran <= prob)exit
     end do

     e_cut = nucleus(icomp_f)%level_param(7)

     num_part = 1
!
!------   Find angular distribution quantities for pre-equilirum emission
!------   Finds angle based Kalbach-Mann statistics for pre-equilibrium
!------   emission
!
     theta_0 = 0.0d0

     do nang = 1, num_theta
        if(k > 0)then
           call PREEQ_Angular(icomp_i, icomp_f, iproj, e_in, k, energy, x_Ang)
           if(abs(x_Ang) > 1.0d0)then
              write(6,*)'cos(theta) wrong in PREEQ_sample'
              write(6,*)'iproc = ',iproc
#if(USE_MPI==1)
              call MPI_Abort(icomm,101,ierr)
#endif
           end if
           extra_angle_data(nang,num_part) = acos(x_Ang)
        else               !   make photons isotropic for now
           extra_angle_data(nang,num_part) = two_pi*random_32(iseed_32)
        end if
     end do
     theta_0 = extra_angle_data(1,num_part)
     phi_0 = two_pi*random_32(iseed_32)

!-----------------------------------------  Find energy bin, or discrete state
     ex_min_bin = nucleus(icomp_f)%e_grid(1) - 0.5d0*nucleus(icomp_f)%delta_e(1)

!------
!------   Decay to a continuous energy bin, defined by nbin_f
!------
     xj_f = -1.0d0
     xI_f = -1.0d0
     xj_f_min = -1.0d0
     xj_f_max = -1.0d0
     xI_f_min = -1.0d0
     xI_f_max = -1.0d0
     ixI_f_min = -1
     ixI_f_max = -1

     if(ex_final > e_cut .and. ex_final > ex_min_bin)then     !  energy bin
        nbin_f = 0
        nbin_f = find_ibin(ex_final,icomp_f)
!-------
        if(nbin_f <= 0) nbin_f = 1
        idb = 0

        ex_final = nucleus(icomp_f)%e_grid(nbin_f)
        energy = ex_tot - ex_final - nucleus(icomp_i)%sep_e(k)
        if(energy <= 0.0d0)nbin_f = nbin_f -1
        ex_final = nucleus(icomp_f)%e_grid(nbin_f)
        energy = ex_tot - ex_final - nucleus(icomp_i)%sep_e(k)
!------
!------   Estimate of angular momentum for emitted particle based on the 
!------   emitted angle and assuming the reaction occurs near the surface
!------   l_f is outgoing orbital angular momentum
!------
        if(k > 0)then
           mass_e = particle(k)%mass
           l_f = preeq_l(l_i, mass_p, E_in, mass_e, energy, theta_0, A_p, A_t)
           if(l_f > particle(k)%lmax)then
              if(iproc == 0)then
                 write(6,*)'*******  Warning l_f > particle%lmax in PREEQ_samp  ******'
                 write(6,*)'*******  l_f set = particle%lmax                    ******'
                 write(6,*)particle(k)%name
                 write(6,*)energy, theta_0
                 write(6,*)l_i, l_f, particle(k)%lmax
              end if
              l_f = particle(k)%lmax
           end if
           cpar2 = par*(-1.0d0)**l_f                        !  Parity of nucleus and emitted particle
           par_f = cpar2*particle(k)%par                 !  Parity of final nucleus
           ip_f = nint((par_f + 1.0d0)/2.0d0)
!------
!------   With estimate of outgoing orbital angular momentum, couple to spin
!------   and initial angular momentum to find possible final angular momenta
!------   Then choose final angular momentum based on density of states within
!------   the min and max ranges
!------
           xj_f_min = abs(real(l_f,kind=8) - particle(k)%spin)
           xj_f_max = real(l_f,kind=8) + particle(k)%spin
           xI_f_min = min(abs(xI_i - xj_f_min),abs(xI_i - xj_f_max))
           xI_f_max = xI_i + xj_f_max
           ixI_f_min = nint(xI_f_min - nucleus(icomp_f)%jshift)
           ixI_f_max = min(nint(xI_f_max - nucleus(icomp_f)%jshift), nucleus(icomp_f)%j_max)
        elseif(k == 0)then
           par_f = -par
           ip_f = nint((par_f + 1.0d0)/2.0d0)
           xI_f_min = abs(xI_i - 1.0d0)
           xI_f_max = xI_i + 1.0d0
           ixI_f_min = nint(xI_f_min - nucleus(icomp_f)%jshift)
           ixI_f_max = min(nint(xI_f_max - nucleus(icomp_f)%jshift), nucleus(icomp_f)%j_max)
        end if
!
!---- Distribute J and parity according to level density
!
        rho_sum = 0.0d0
        do Ix_f = ixI_f_min, ixI_f_max
           rho_sum = rho_sum + nucleus(icomp_f)%bins(Ix_f, ip_f, nbin_f)%rho
        end do
!
!---- If rho_sum > 0.0, is normal situation, choose angular momentum
!---- randomly with probability defined by the density of states
!----
        if(rho_sum >= 1.0d-6)then
           ran = random_32(iseed_32)
           prob = 0.0d0
           do Ix_f = ixI_f_min, ixI_f_max
              prob = prob + nucleus(icomp_f)%bins(Ix_f, ip_f, nbin_f)%rho/rho_sum
              if(ran <= prob)exit
           end do
        else
!---- If rho_sum = 0.0, then this decay can't happen due to no states
!---- being in the angular momentum window. Generally rare, but possible
!---- for lighter nuclei with lower level densities, high emission energy,
!---- and higher initial angular momentum. In this case, reduce spin to
!---- first angular momentum bin with non-zero level density. If, for some
!---- reason it still can't be distributed, start over.
           do Ix_f = iXI_f_min, 0, -1
              if(nucleus(icomp_f)%bins(Ix_f, ip_f, nbin_f)%rho > 1.0d-6)exit
           end do
           if(Ix_f < 0)goto 22    !---  It went past min value; no bin allowed, start over 
        end if
        xI_f = Ix_f + nucleus(icomp_f)%jshift
        xj_f = xj_f_min
     else                         ! Discrete state
!----- Trap that ex_final is actually greater than e_cut
        if(ex_final > e_cut)ex_final = e_cut - nucleus(icomp_f)%delta_e(1)*random_32(iseed_32)
        dee = 0.5d0*de
        idb = 1
        nbin_f = 0
        if(k > 0)then
           frac = 0
           do i = 1, nucleus(icomp_f)%num_discrete     ! start with 1st excited state - inelastic
              e_state = nucleus(icomp_f)%state(i)%energy
              if(e_state + dee >= ex_final .and. e_state - dee <= ex_final)frac = frac + 1
           end do
           if(frac > 0.0d0)then
              frac=1.0d0/frac
              ran = random_32(iseed_32)
              prob = 0.0d0
              do i = 1, nucleus(icomp_f)%num_discrete
                 e_state = nucleus(icomp_f)%state(i)%energy
                 if(e_state + dee >= ex_final .and. e_state - dee <= ex_final)then
                    prob = prob + frac
                    if(ran <= prob)then
                       nbin_f = i
                       exit
                    end if
                 end if
               end do
           else            ! special case where there are no states, Force into closest level below ex-final
              do i = nucleus(icomp_f)%num_discrete,1,-1
                 if(nucleus(icomp_f)%state(i)%energy < ex_final)then
                    nbin_f = i
                    exit
                 end if
              end do
!------------   Lastly, it is possible due to roundoff error that the energy could fall just below the 
!------------   tolerance for the ground-state, returning nbin_f = 0, which will cause a seg fault
!------------   force decay to ground state
              if(nbin_f == 0)nbin_f = 2                !   lowest state has to be first excited state
           end if

           ex_final = nucleus(icomp_f)%state(nbin_f)%energy
           energy = ex_tot - ex_final - nucleus(icomp_i)%sep_e(k)

           xI_f = nucleus(icomp_f)%state(nbin_f)%spin
           ip_f = nint((nucleus(icomp_f)%state(nbin_f)%parity+1.0d0)/2.0d0)
           xj_f_min = abs(xI_i - xI_f)
           xj_f_max = xI_i + xI_f
           l_min = nint(min(abs(xj_f_min - particle(k)%spin),abs(xj_f_max - particle(k)%spin)))
           l_max = nint(xj_f_max + particle(k)%spin)
           do l = l_min, l_max
              par_f = (-1.0d0)**l*nucleus(icomp_f)%state(nbin_f)%parity*particle(k)%par
!              if(par_f == par)exit
              if(real8_equal(par_f,par))exit
           end do
           l_f = l
           if(l_f == 0)then
              xj_f = particle(k)%spin
           else
              if(random_32(iseed_32) <= 0.5)then
                 xj_f = real(l_f) + particle(k)%spin
              else
              xj_f = real(l_f) - particle(k)%spin
              end if
           end if
        elseif(k == 0)then                       !   E1 decay, more restrictive, so force decay to closest discrete state
           xI_f_min = abs(xI_i - 1.0d0)
           xI_f_max = xI_i + 1.0d0
           par_f = -par
           ip_f = nint((par_f + 1.0d0)/2.0d0)
           e_closest = 99999.0d0
           closest = -1
           do i = 1, nucleus(icomp_f)%num_discrete
              e_diff = abs(ex_final - nucleus(icomp_f)%state(i)%energy)
              xI_state = nucleus(icomp_f)%state(i)%spin
              ip_state = nint((nucleus(icomp_f)%state(i)%parity + 1.0d0)/2.0d0)
              if((xI_state >= xI_f_min .and. xI_state <= xI_f_max) .and. ip_f == ip_state .and.     &
                  e_diff < e_closest)then
                  closest = i
                  e_closest = e_diff
              end if
           end do
           if(closest > 0)then
              nbin_f = closest
              xI_f = nucleus(icomp_f)%state(closest)%spin
              ex_final = nucleus(icomp_f)%state(closest)%energy
              energy = ex_tot - ex_final - nucleus(icomp_i)%sep_e(k)
           else
              goto 22                                         !  It didn't work, so start over
           end if
        end if
     end if
!
!----   Store date for decay in part_data, specifying all aspects of this decay
!
     part_data(1,num_part) = real(icomp_f,kind=8)
     part_data(2,num_part) = real(k,kind=8)
     part_data(3,num_part) = xI_f
     part_data(4,num_part) = real(2*ip_f-1,kind=8)
     part_data(5,num_part) = real(nbin_f,kind=8)
     part_data(6,num_part) = real(idb,kind=8)
     part_data(7,num_part) = real(l_f,kind=8)
     part_data(9,num_part) = energy

     if(k == 0)then
        part_data(8,num_part) = 0.0d0     !  set to zero for now because no photons in preeq
     else
        part_data(8,num_part) = xj_f
     end if


     part_data(10,num_part) = theta_0
     part_data(11,num_part) = phi_0
!     part_data(12,num_part) = T_1
!     part_data(13,num_part) = theta
!     part_data(14,num_part) = phi
!     part_data(15,num_part) = T_L
!     part_data(16,num_part) = theta_L
!     part_data(17,num_part) = phi_L
!     part_data(18,num_part) = T_2
     part_data(19,num_part) = tally_prob
     part_data(20,num_part) = ex_final
     part_data(21,num_part) = icomp_i
     part_data(22,num_part) = nucleus(itarget)%state(istate)%spin
     part_data(23,num_part) = nucleus(itarget)%state(istate)%parity
     part_data(24,num_part) = istate
!     nucleus(icomp_f)%Kinetic_Energy = T_2
!
!----    Before leaving update number of particles for this decay type
!
     if(k >= 0 .and. k <= 6)then
        num_part_type(k) = num_part_type(k) + 1
        if(k > 0 .and. num_part_type(k) >= max_particle(k))part_fact(k) = 0.0d0
     end if

   return
end subroutine PREEQ_Sample
!
!*******************************************************************************
!
subroutine PREEQ_Angular(icomp_C, icomp_B, iproj, E_A, ieject, E_B, x_Ang)
!
!*******************************************************************************
!
!  Discussion:
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options 
!        nuclei
!        particles_def
!        constants
!        nodeinfo
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: expdev
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
   use particles_def
   use constants
   use nodeinfo
   implicit none
!-------------------   Argument data
   integer(kind=4), intent(in) :: icomp_C, icomp_B
   integer(kind=4), intent(in) :: iproj
   real(kind=8), intent(in) :: E_A
   integer(kind=4), intent(in) :: ieject
   real(kind=8), intent(in) :: E_B
   real(kind=8), intent(out) :: x_Ang
!-------------------   Internal data
   real(kind=8) :: a, E1, E3, eap, ebp, Sa, Sb
   real(kind=8) :: Mb, ma
   real(kind=8) :: A_C, A_B, Z_C, Z_B, N_C, N_B, Ia, Ib

   real(kind=8) :: ex1, ex2, ex4
!------------   External functions    -------------------------------
   real(kind=8) :: expdev
!-------------------   Run program    -------------------------------

   Ma = 0.0d0
   if(iproj < 6)Ma = 1.0d0
   if(ieject == 1)then
      Mb = 0.5d0
   elseif(ieject > 1 .and. ieject < 6)then
      Mb = 1.0d0
   else
      Mb = 2.0d0
   end if
   Ia = 0.0d0
   if(iproj == 3)Ia = 2.225d0
   if(iproj == 4)Ia = 8.482d0
   if(iproj == 5)Ia = 7.718d0
   if(iproj == 6)Ia = 28.296d0
   Ib = 0.0d0
   if(ieject == 3)Ib = 2.225d0
   if(ieject == 4)Ib = 8.482d0
   if(ieject == 5)Ib = 7.718d0
   if(ieject == 6)Ib = 28.296d0


   A_C = nucleus(icomp_C)%A
   Z_C = nucleus(icomp_C)%Z
   N_C = A_C - Z_C
   A_B = nucleus(icomp_B)%A
   Z_B = nucleus(icomp_B)%Z
   N_B = A_B - Z_B

   ex1 = 1.0d0/3.0d0
   ex2 = 2.0d0/3.0d0
   ex4 = 4.0d0/3.0d0
   Sa = 15.68d0*(A_C-A_B) - 28.07*((N_C-Z_C)**2/A_C - (N_B-Z_B)**2/A_B) -    &
        18.56d0*(A_C**ex2 - A_B**ex2) +                                      &
        28.07d0*((N_C-Z_C)**2/A_C**ex4 - (N_B-Z_B)**2/A_B**ex4) -            &
        0.717d0*(Z_C**2/A_C**ex1 - Z_B**2/A_B**ex1) +                        &
        1.211d0*(Z_C**2/A_C - Z_B**2/A_B)
   Sb = Sa
   Sa = Sa - Ia
   Sb = Sb - Ib

   eap = E_a + Sa
   ebp = E_b + Sb

   E1 = min(eap, 130.0d0)
   E3 = min(eap, 41.0d0)

   a = 0.04d0*E1*ebp/eap + 1.8d-6*(E1*ebp/eap)**3 +      &
       6.7d-7*Ma*mb*(E3*ebp/eap)**4

   x_Ang = expdev(iseed_32,a)

   return
end subroutine PREEQ_Angular
!
!*******************************************************************************
!
real(kind=8) function exp_leg_int(n,a)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates the integral of exp(-a*x)*Leg(n,x)
!    Perfomed using Guass-Legendre integration with n_leg
!    points. Weights and abbscissas are calculated in main 
!    routine with a call to gauleg. The points x_gleg, and
!    w_gleg are stored in the module Gauss_integration.
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
!        Gauss_integration
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        poly
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
   use Gauss_integration
   implicit none
!------------   Argument declarations
   integer(kind=4), intent(in) :: n
   real(kind=8), intent(in) :: a
!------------  Internal data
   integer(kind=4) :: m
   real(kind=8) :: sum
!------------  External functions
   real(kind=8) :: poly
!------------  Start calculation
!-----   Legendre(n,x) polynomials are compouted with function poly(n,kind,alf,bet,x)
!-----   with kind=1, alf=bet=0.0
   exp_leg_int = 0.0d0
   alf = 0.0d0
   bet = 0.0d0
   sum = 0.0d0
   do m = 1 , n_gleg
      sum = sum + w_gleg(m)*exp(a*x_gleg(m))*poly(n,1,alf,bet,x_gleg(m))
   end do
   sum = sum*(0.5d0*a/sinh(a))
   exp_leg_int = sum
   return
end function exp_leg_int
!
!*******************************************************************************
!
integer(kind=4) function preeq_l(l_p, mass_p, E_p, mass_e, E_e, theta_e, A_p, A_t)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function uses a model to compute the orbital angular momentum l
!    of an emitted particle given is energy, mass and emitted angle.
!    It assumes momentum conservation and that the reaction occurs on 
!    the nuclear surface.
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
   use constants
   implicit none
!--------------------   Argument declarations
   integer(kind=4), intent(in) :: l_p
   real(kind=8), intent(in) :: mass_p, E_p, mass_e, E_e, theta_e
   integer(kind=4), intent(in) :: A_p, A_t
!--------------------   Internal data
   real(kind=8) :: phi
   real(kind=8) :: b_p, k_p, b_e, k_e
   real(kind=8) :: R
   real(kind=8) :: ex1
!--------------------   Calculation
   if(E_p < 1.0d-2)then
      preeq_l = 0
      return
   end if
   ex1 = 1.0d0/3.0d0
   R = 1.2*(real(A_p,kind=8)**ex1 + real(A_t,kind=8)**ex1)
   k_p = sqrt(2.0d0*mass_p*E_p)/hbar_c
   k_e = sqrt(2.0d0*mass_e*E_e)/hbar_c
   b_p = (real(l_p,kind=8)+0.5d0)/k_p
   if(b_p < R)then
      phi = asin(b_p/R)
      if(theta_e < pi/2.0d0)then
         b_e = sqrt(R**2-b_p**2)*sin(theta_e) + b_p*cos(theta_e)
      else
         b_e = abs(sqrt(R**2-b_p**2)*sin(theta_e) + b_p*cos(theta_e))
      end if
   else
      b_e = R*abs(cos(theta_e))
   end if

   preeq_l = nint(b_e*k_e)

   return
end function preeq_l
!
!*******************************************************************************
!
real(kind=8) function expdev(iseed_32,a)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function generates random variables distributed
!    according to exp(a*x) on the interval (-1,1)
!    actually, properly normalized is a*exp(a*x)/(2*sinh(a))
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
!       real(kind=8) random_32
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
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   use variable_kinds
   implicit none
   integer(kind=int_32), intent(inout) :: iseed_32
   real(kind=8), intent(in) :: a
   real(kind=8) dum
!--------    External Functions    -------------------------------------------
   real(kind=8) random_32
   logical :: real8_equal
!-----------------------------------------------------------------------------
   expdev = 10.0d0
   do while(expdev > 1.0d0)
      dum = 1.0d0 - random_32(iseed_32)*(1.0d0-exp(-2.0d0*a))
      do while(real8_equal(dum,0.0d0))
         dum = 1.0d0 - random_32(iseed_32)*(1.0d0-exp(-2.0d0*a))
      end do
      expdev = -log(dum)/a
      expdev = -(expdev - 1.0d0)
   end do

   return
end function expdev
