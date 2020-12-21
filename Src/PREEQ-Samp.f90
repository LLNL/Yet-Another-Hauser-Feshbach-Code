!
!*******************************************************************************
!
subroutine PREEQ_sample(iproj, in, itarget, istate, e_in, ex_tot,      &
                        l_i, is_i, Ix_i, ip_i, icomp_i, icomp_f,       &
                        Ix_f, l_f, ip_f, nbin_f, idb,                  &
                        n_dat, dim_part, num_part, part_data,          &
                        Ang_L_max, part_Ang_data, x_Ang,               &
                        Boost_Lab, Boost_COM)
!
!*******************************************************************************
!
!  Discussion:
!
!    This is a subroutine to Monte Carlo sample the pre-equilibrium decay 
!    of an energy bin
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
   integer(kind=4), intent(out) :: num_part
   real(kind=8), intent(out) :: part_data(n_dat,dim_part)
   integer(kind=4), intent(in) :: Ang_L_max
   real(kind=8), intent(out) :: part_Ang_data(0:Ang_L_max,dim_part)
   real(kind=8), intent(out) :: x_Ang
   real(kind=8), intent(out) :: Boost_lab(0:3,0:3), Boost_COM(0:3,0:3)
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


   real(kind=8) :: dee

   real(kind=8) :: costhp, sinthp, theta_0, phi_0
   real(kind=8) :: T_1, T_2, mass_1, mass_2, theta, phi
   real(kind=8) :: T_L, theta_L, phi_L

   real(kind=8) :: prob_part(0:6)
   real(kind=8) :: sum_prob
   real(kind=8) :: tally_prob

   real(kind=8) :: shift

!--------------------------------    External functions
   real(kind=8) :: random_64
   integer(kind=4) :: preeq_l
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------------------------    Start Calculation

!----
!----  Entry point to start pre-equilibrium emission selection
!----  used in extremely rare case where a decay failed to find a 
!----  valid final state. 
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

     dee = de/2.0d0
     cpar = particle(iproj)%par*nucleus(itarget)%state(istate)%parity
     par = cpar*(-1.0d0)**l_i
     xj_i = abs(real(l_i,kind=8) - particle(iproj)%spin) + real(is_i,kind=8) 

     mass_p = particle(iproj)%mass
     A_p = particle(iproj)%A
     A_t = nucleus(itarget)%A

     xI_i = dfloat(Ix_i) + nucleus(1)%jshift
     par = 2*ip_i - 1
     ran = random_64(iseed)
     preeq_cs = nucleus(icomp_i)%PREEQ_cs(in)
     k = -1

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
        do kk = 1, nucleus(icomp_i)%num_decay
           k = nucleus(icomp_i)%decay_particle(kk)
           if(nucleus(icomp_i)%PREEQ_part_cs(kk,in)/preeq_cs > 0.0d0)then
              prob_part(k) = 1.0d0
           end if
           sum_prob = sum_prob + prob_part(k)
        end do
        prob = 0.0d0
        do kk = 1, nucleus(icomp_i)%num_decay                       !   Loop over particle allowed to decay from this nucleus
           k = nucleus(icomp_i)%decay_particle(kk)
           icomp_f = nucleus(icomp_i)%decay_to(kk)
           prob = prob + prob_part(k)/sum_prob
!    write(20,*)kk,k,prob_part(k),prob,ran
           if(ran < prob)exit
        end do
        tally_prob = (nucleus(icomp_i)%PREEQ_part_cs(kk,in)/preeq_cs)*(prob_part(k)/sum_prob)
     end if


     if(k == -1)stop 'k not set properly in PREEQ_Samp'
     ran = random_64(iseed)
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
!   write(20,*)kk,m,prob,ran
        if(ran <= prob)exit
     end do

     e_cut = nucleus(icomp_f)%level_param(7)


     num_part = 1
!
!------   Find angular distribution quantities for pre-equilirum emission
!------   Finds angle based Kalbach-Mann statistics for pre-equilibrium
!------   emission
!

     call PREEQ_Angular(icomp_i, icomp_f, iproj, e_in, k, energy,            &
                        dim_part, num_part, Ang_L_max, part_Ang_data, x_Ang)


!     x_ang = 2.0d0*random_64(iseed) - 1.0d0       !   test to make it isotropic
     if(abs(x_Ang) > 1.0d0)stop 'cos(theta) wrong in PREEQ_sample'
     theta_0 = acos(x_Ang)
!     write(20,*)x_Ang, theta_0

!-----------------------------------------  Find energy bin, or discrete state
     ex_min_bin = nucleus(icomp_f)%e_grid(1) - de/2.0d0

!------
!------   Decay to a continuos energy bin, defined by nbin_f
!------

     if(ex_final > e_cut .and. ex_final > ex_min_bin)then     !  energy bin
        nbin_f = int((ex_final - ex_min_bin)/de) + 1
!-------
        if(nbin_f <= 0) nbin_f = 1
        idb = 0

        if(energy <= de)then
           energy = de*random_64(iseed)
           ex_final = ex_tot - energy - nucleus(icomp_i)%sep_e(k)
        else           
           shift = de*(random_64(iseed) - 0.5d0)
           energy = energy + shift
           ex_final = ex_tot - energy - nucleus(icomp_i)%sep_e(k)
        end if
        mass_e = particle(k)%mass

!------
!------   Estimate of angular momentum for emitted particle based on the 
!------   emitted angle and assuming the reaction occurs near the surface
!------   l_f is outgoing orbital angular momentum
!------

        l_f = preeq_l(l_i, mass_p, E_in, mass_e, energy, theta_0, A_p, A_t)

        if(l_f > particle(k)%lmax)then
           l_f = particle(k)%lmax
           write(6,*)'*******  Warning l_f > particle%lmax in PREEQ_samp  ******'
           write(6,*)'*******  l_f set = particle%lmax                    ******'
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
           ran = random_64(iseed)
!           ran = ran3(iseed)
           prob = 0.0d0
           do Ix_f = ixI_f_min, ixI_f_max
              prob = prob + nucleus(icomp_f)%bins(Ix_f, ip_f, nbin_f)%rho/rho_sum
              if(ran <= prob)exit
           end do
        else
!
!---- If rho_sum = 0.0, then this decay can't happen due to no states
!---- being in the angular momentum window. Generlaly rare, but possible
!---- for lighter nuclei with lower level densities, high emission energy,
!---- and higher initial angular momentum. In this case, reduce spin to
!---- first angular momentum bin with non-zero level density. If, for some
!---- reason it still can't be distributed, start over.
!----
           do Ix_f = iXI_f_min, 0, -1
              if(nucleus(icomp_f)%bins(Ix_f, ip_f, nbin_f)%rho > 1.0d-6)exit
           end do
           if(Ix_f < 0)goto 22    !---  It went past min value; no bin allowed, start over 
        end if
        xI_f = Ix_f + nucleus(icomp_f)%jshift
        xj_f = xj_f_min

     else                         ! Discrete state
!----- Trap that ex_final is actually greater than e_cut
        if(ex_final > e_cut)ex_final = e_cut - de*random_64(iseed)

        idb = 1
        nbin_f = 0
        frac = 0
        do i = 1, nucleus(icomp_f)%num_discrete     ! start with 1st excited state - inelastic
           if(nucleus(icomp_f)%state(i)%energy + dee >= ex_final .and.           &
              nucleus(icomp_f)%state(i)%energy -dee <= ex_final)frac = frac + 1
        end do
        if(frac > 0.0d0)then
           frac=1.0d0/frac
           ran = random_64(iseed)
           prob = 0.0d0
           do i = 1, nucleus(icomp_f)%num_discrete 
              if(nucleus(icomp_f)%state(i)%energy + dee >= ex_final .and.        &
                 nucleus(icomp_f)%state(i)%energy - dee <= ex_final)then
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
!------------   Lastly it is possible due to roundoff error that the energy could fall just below the 
!------------   tolerance for the ground-state, returning nbin_f = 0, which will cause a seg fault
!------------   force decay to ground state
           if(nbin_f == 0) nbin_f = 2                !   lowest state has to be first excited state
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
           if(par_f == par)exit
        end do
        l_f = l
        if(l_f == 0)then
           xj_f = particle(k)%spin
        else
           if(random_64(iseed) <= 0.5)then
              xj_f = real(l_f) + particle(k)%spin
           else
              xj_f = real(l_f) - particle(k)%spin
           end if
        end if

     end if

!----   Store date for decay in part_data, specifying all aspects of this decay

     part_data(1,num_part) = real(icomp_f,kind=8)
     part_data(2,num_part) = real(k,kind=8)
     part_data(3,num_part) = xI_f
     part_data(4,num_part) = real(2*ip_f-1,kind=8)
     part_data(5,num_part) = real(nbin_f,kind=8)
     part_data(6,num_part) = real(idb,kind=8)
     part_data(7,num_part) = real(l_f,kind=8)
     part_data(9,num_part) = energy

     if( k == 0)then
        part_data(8,num_part) = 0.0d0     !  set to zero for now because no photons in preeq
     else
        part_data(8,num_part) = xj_f
     end if

     costhp = 2.0*random_64(iseed) - 1.0
     sinthp = sqrt(1.0-costhp**2)
     phi_0 = two_pi*random_64(iseed)

     mass_1 = particle(k)%Mass
     if(idb == 0)then
        mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%e_grid(nbin_f)
     else
        mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%state(nbin_f)%energy
     end if
!
!----   Compute energy, angles, all decay data in original COM frame and
!----   Lab frame. 

     call Boost_frame(energy, mass_1, mass_2, theta_0, phi_0,                   &
                      Boost_Lab, Boost_COM, T_1, theta, phi,                    &
                      T_2, T_L, theta_L, phi_L)

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
     part_data(20,num_part) = ex_final
     part_data(21,num_part) = icomp_i
     nucleus(icomp_f)%Kinetic_Energy = T_2

   return
end subroutine PREEQ_Sample
!
!*******************************************************************************
!
subroutine PREEQ_Angular(icomp_C, icomp_B, iproj, E_A, ieject, E_B,           &
                         dim_part, num_part, Ang_L_max, part_Ang_data, x_Ang)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine returns the angles for a particle emitted by
!    pre-equilibirum emission
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
!-------------------   Argument data
   integer(kind=4), intent(in) :: icomp_C, icomp_B
   integer(kind=4), intent(in) :: iproj
   real(kind=8), intent(in) :: E_A
   integer(kind=4), intent(in) :: ieject
   real(kind=8), intent(in) :: E_B
   integer(kind=4), intent(in) :: dim_part, num_part
   integer(kind=4), intent(in) :: Ang_L_max
   real(kind=8), intent(out) :: part_Ang_data(0:Ang_L_max,dim_part)
   real(kind=8), intent(out) :: x_Ang
!-------------------   Internal data
   real(kind=8) :: a, E1, E3, eap, ebp, Sa, Sb
   real(kind=8) :: Mb, ma
   real(kind=8) :: A_C, A_B, Z_C, Z_B, N_C, N_B, Ia, Ib

   real(kind=8) :: ex1, ex2, ex4

   integer(kind=4) :: numx
   parameter (numx = 200)

   integer(kind=4) :: L
   real(kind=8) :: xL
!------------   External functions
   real(kind=8) :: exp_leg_int
   real(kind=8) :: expdev
!-------------------   Run program


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

   part_Ang_data(0,num_part) = 0.5d0
   if(output_mode >= 1)then
      do L = 1, Ang_L_max
         xL = real(L,kind=8)
         part_Ang_data(L,num_part) = exp_leg_int(L,a)*(2.0d0*xL+1.0d0)/2.0d0
      end do
   end if

   x_Ang = expdev(iseed,a)

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
   use Gauss_integration
   implicit none
!------------   Argument declarations
   integer(kind=4), intent(in) :: n
   real(kind=8), intent(in) :: a
!------------  Internal data
   integer(kind=4) :: m
   real(kind=8) :: sum
!------------  External functions
!   real(kind=8) :: Legendre
   real(kind=8) :: poly
!------------  Start calculation
!-----   Legendre(n,x) polynomials are compouted with function poly(n,kind,alf,bet,x)
!-----   with kind=1, alf=bet=0.0
   exp_leg_int = 0.0d0
   alf = 0.0d0
   bet = 0.0d0
   sum = 0.0d0
   do m = 1 , n_gleg
!      sum = sum + w_gleg(m)*exp(a*x_gleg(m))*Legendre(n,x_gleg(m))
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
      b_e = b_p*abs(cos(theta_e))
   end if

   preeq_l = nint(b_e*k_e)
   return
end function preeq_l
!
!*******************************************************************************
!
real(kind=8) function expdev(iseed,a)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function generates random variables distributed
!    according to exp(a*x) on the interval (-1,1)
!    actually, properly normalized is a*exp(a*x)/(2*sinh(a))
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
   integer(kind=int_64), intent(inout) :: iseed
   real(kind=8), intent(in) :: a
!-------------------------
   real(kind=8) dum
   real(kind=8) random_64
!------------------------------
 1 dum = 1.0d0 - random_64(iseed)*(1.0d0-exp(-2.0d0*a))
   if(dum == 0.0) goto 1
   expdev = -log(dum)/a
   expdev = -(expdev - 1.0d0)
   if(abs(expdev) > 1.0d0)goto 1
   return
end function expdev
