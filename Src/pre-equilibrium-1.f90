!
!*******************************************************************************
!
subroutine pre_equilibrium_1(icomp, istate, in, E_inc,    &
                             Ex_tot, de, reac_cs)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine sets up the data for the pre-equilibirum model #1
!    the two-component exciton model. It precalculates information needed 
!    in the main code (pre-equilibirum probability) PREEQ-Samp.f90 (spectra),
!    which Monte Carlo samples the decay.
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        channel_info
!        particles_def
!        constants 
!        nodeinfo
!        pre_equilibrium_no_1
!
!     Subroutines:
!
!        int_trans_rate
!
!     External functions:
!
!        real(kind=8) :: compound_cs
!        real(kind=8) :: omega2
!        real(kind=8) :: Well
!        real(kind=8) :: Delta_pre
!        real(kind=8) :: Pauli
!        real(kind=8) :: Prob_func
!        real(kind=8) :: finite_well
!        real(kind=8) :: aparam_u
!        real(kind=8) :: EL_absorption
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
!  Reference:
!     
!    A. J. Koning and M.C. Dujvestijn, A global pre-equilibrium analysis 
!    from 7 to 200 MeV based on the optical model potential, 
!    Nucl. Phys. A 744. 15 (2004)
!
!*******************************************************************************
!
   use variable_kinds
   use options
   use nuclei
   use channel_info
   use particles_def
   use constants 
   use nodeinfo
   use pre_equilibrium_no_1
   implicit none
!---------   input variables   -----------------------------
   integer(kind=4), intent(in) :: icomp, istate, in
   real(kind=8), intent(in) :: E_inc, Ex_tot
   real(kind=8), intent(in) :: de, reac_cs
!-----------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: j, k, m, kk
   integer(kind=4) :: h
   integer(kind=4) :: pn, pp, hn, hp, pn_tot, pp_tot
   integer(kind=4) :: pn_res, pp_res, hn_res, hp_res
   integer(kind=4) :: pn_calc, pp_calc
   integer(kind=4) :: p_tot, h_tot, n_tot, n0_tot
   integer(kind=4) :: pn_min, pp_min
   integer(kind=4) :: z_k, n_k, Z, N, A, Ap, Zf, Nf, Af
   integer(kind=4) :: j_max
   integer(kind=4) :: ipar
   real(kind=8) :: pre_eq_cs(0:6)

   real(kind=8) :: xn
   real(kind=8) :: energy, U
   real(kind=8) :: sig_inv
   real(kind=8) :: sp1, sp2
   real(kind=8) :: spin_target, spin_proj
   real(kind=8) :: mass_i, mass_t, mass_rel
   real(kind=8) :: factor, factor1
   real(kind=8) :: Msq
   real(kind=8) :: xA
   real(kind=8) :: xAp
   real(kind=8) :: Msq_pp, Msq_pn, Msq_np, Msq_nn
   real(kind=8) :: sum,sum2
   real(kind=8) :: L_p_p_an, L_p_n_an, L_0_pn_an, L_0_np_an
   real(kind=8) :: A_Pauli, B_Pauli
   real(kind=8) :: denom
   real(kind=8) :: omdenom, om, om1, om2
   real(kind=8) :: Sep
   real(kind=8) :: gp,gn, g
   real(kind=8) :: gpp, gnn, gg
   real(kind=8) :: gpf, gnf, gf
   real(kind=8) :: gppf, gnnf, ggf
   real(kind=8) :: Delta

   real(kind=8) :: V1, V3, xK

   real(kind=8) :: stest,stest2

   integer(kind=4) :: nbin_end
   real(kind=8)    :: ex_final, e_max, e_bin
   real(kind=8)    :: xji
   integer(kind=4) :: h_max
   real(kind=8) :: direct, semi_direct, term

   integer(kind=4) :: ifinal
   real(kind=8) :: e_cut

   real(kind=8) :: shell, gamma

   integer(kind=4) :: ifile
   character(len=20) :: outfile

   real(kind=8) :: xZ_i, xA_i, xZ_part, xA_part
   real(kind=8) :: Coulomb_Barrier(0:6)
!--------   External Functions   ---------------------------
   real(kind=8) :: compound_cs
   real(kind=8) :: omega2
   real(kind=8) :: Well
   real(kind=8) :: Delta_pre
   real(kind=8) :: Pauli
   real(kind=8) :: Prob_func
   real(kind=8) :: finite_well
   real(kind=8) :: aparam_u
   real(kind=8) :: EL_absorption
!--------   Start Calculation    ---------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------    Calculate internal transition rates
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   Coulomb_barrier(0:6) = 0.0d0
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

   j_max = nucleus(icomp)%j_max
   iproj = projectile%particle_type
   Z = nucleus(icomp)%Z
   A = nucleus(icomp)%A
   xA = real(A,kind=8)
   Ap = projectile%A
   xAp = real(Ap,kind=8)
   N = A - Z
   if(Preeq_g_a)then
      g = nucleus(icomp)%a_Sn*pi**2/6.0d0
      Preeq_g_div = xA/g
   end if
   gnn = real(N,kind=8)/Preeq_g_div
   gpp = real(Z,kind=8)/Preeq_g_div
   gg = gpp + gnn


   pn_min = p0(1)
   pp_min = p0(2)
   n0_tot = p0(1) + p0(2)
   pn_calc = pn_max
   pp_calc = pp_max
 

   h_max = pn_max + pp_max

   tau(0:pn_max+1,0:pp_max+1) = 0.0d0
   taup(0:pn_max+1,0:pp_max+1) = 0.0d0
   Pre_Prob(0:pn_max+1,0:pp_max+1) = 0.0d0
   Gam_p(0:pn_max+1,0:pp_max+1) = 0.0d0
   Gam_n(0:pn_max+1,0:pp_max+1) = 0.0d0
   Gam_pp(0:pn_max+1,0:pp_max+1) = 0.0d0
   Gam_pn(0:pn_max+1,0:pp_max+1) = 0.0d0
   Gam_0_pn(0:pn_max+1,0:pp_max+1) = 0.0d0
   Gam_0_np(0:pn_max+1,0:pp_max+1) = 0.0d0
   Lamb_p_p(0:pn_max+1,0:pp_max+1) = 0.0d0
   Lamb_p_n(0:pn_max+1,0:pp_max+1) = 0.0d0
   Lamb_0_pn(0:pn_max+1,0:pp_max+1) = 0.0d0
   Lamb_0_np(0:pn_max+1,0:pp_max+1) = 0.0d0
   dWk(0:pn_max+1,0:pp_max+1,0:nucleus(icomp)%nbin_part,0:6) = 0.0d0
   Wk(0:pn_max+1,0:pp_max+1,0:6) = 0.0d0
   W(0:pn_max+1,0:pp_max+1) = 0.0d0
   outfile(1:3) = 'dWn'
   ifile=50
   shell = 0.0d0
   gamma = 0.0d0
   gp = aparam_u(Ex_tot,gpp,shell,gamma)
   gn = aparam_u(Ex_tot,gnn,shell,gamma)
   g = gp + gn
   V1 = Preeq_V1
   V3 = Preeq_V
   xK = Preeq_K
   if(v1 > V3) V1 = V3
   do h = 0, h_max
      if(preeq_fwell == 0 .or. preeq_fwell == 2)Vwell(h) = Well(h,A,E_inc,V1,V3,xK)
      if(preeq_fwell == 1)Vwell(h) = V3
      Vwell_g(h) = V3
   end do

   do pn = p0(1), pn_max
      do pp = p0(2), pp_max
         hn = pn - p0(1)
         hp = pp - p0(2)
         pn_tot = pn + hn
         pp_tot = pp + hp
         p_tot = pn + pp
         h_tot = hn + hp
         n_tot = p_tot + h_tot
         Lamb_p_p(pn,pp) = 0.0d0
         Lamb_p_n(pn,pp) = 0.0d0
         Lamb_0_pn(pn,pp) = 0.0d0
         Lamb_0_np(pn,pp) = 0.0d0

         Delta = Delta_pre(pn,hn,pp,hp,Z,A,gp,gn,Ex_tot)

         omdenom = omega2(pn,hn,pp,hp,Z,A,gp,gn,Ex_tot,Delta,H_max,Vwell)

         if(omdenom <= 1.0d-10)cycle
         xn = real(n_tot,kind=8)
         Msq = (M2_C1*xAp/xA**3)*(7.48d0*M2_C2+4.62d5/(Ex_tot/(xn*xAp)+10.7*M2_C3)**3)
         if(analytic_preeq)Msq=1.2*Msq
         Msq_pp = M2_Rpp*Msq
         Msq_pn = M2_Rpn*Msq
         Msq_nn = M2_Rnn*Msq
         Msq_np = M2_Rnp*Msq
         if(n_tot > 1 .and. .not. analytic_preeq)then
            call int_trans_rate(pn,hn,pp,hp,Z,A,gp,gn,Ex_tot,Delta,h_max,Vwell,      &
                                Msq_nn,Msq_pp,Msq_pn,Msq_np,                         &
                                Lamb_p_p(pn,pp),Lamb_p_n(pn,pp),                     &
                                Lamb_0_pn(pn,pp),Lamb_0_np(pn,pp))
         else
            Lamb_p_p(pn,pp) = 0.0d0
            Lamb_p_n(pn,pp) = 0.0d0
            Lamb_0_pn(pn,pp) = 0.0d0
            Lamb_0_np(pn,pp) = 0.0d0
            L_p_p_an = 0.0d0
            L_p_n_an = 0.0d0
            L_0_pn_an = 0.0d0
            L_0_np_an = 0.0d0

            L_p_p_an = 0.0d0
            L_p_n_an = 0.0d0

            Delta = Delta_pre(pn,hn,pp,hp,Z,A,gp,gn,Ex_tot)
            U = Ex_tot - Delta
            if(U <= 0.0d0)cycle



            L_p_p_an=2.0d0*pi*gp**2/(2.0d0*xn*(xn+1.0d0))*((Ex_tot-Pauli(pn,hn,pp+1,hp+1,gp,gn))**(n_tot+1)/                     & 
                                                           (Ex_tot-Pauli(pn,hn,pp,hp,gp,gn))**(n_tot-1))*                        &
                                                           (real(pp+hp,kind=8)*gp*Msq_pp+2.0d0*real(pn+hn,kind=8)*gn*Msq_pn)*    &
                                                           finite_well(p_tot+1,h_tot+1,Ex_tot,h_max,Vwell)
            L_p_n_an=2.0d0*pi*gn**2/(2.0d0*xn*(xn+1.0d0))*((Ex_tot-Pauli(pn+1,hn+1,pp,hp,gp,gn))**(n_tot+1)/                     &
                                                           (Ex_tot-Pauli(pn,hn,pp,hp,gp,gn))**(n_tot-1))*                        &
                                                           (real(pn+hn,kind=8)*gn*Msq_nn+2.0d0*real(pp+hp,kind=8)*gp*Msq_np)*    &
                                                           finite_well(p_tot+1,h_tot+1,Ex_tot,h_max,Vwell)
            A_Pauli = Pauli(pn,hn,pp,hp,gp,gn)
            B_Pauli = max(Pauli(pn,hn,pp,hp,gp,gn),Pauli(pn+1,hn+1,pp-1,hp-1,gp,gn))
            L_0_pn_an = 2.0d0*pi*Msq_pn*(real(pp*hp,kind=8)/xn)*gn**2*                                         &
                        finite_well(p_tot,h_tot,Ex_tot,h_max,Vwell)*                                           &
                        (2.0d0*(Ex_tot-B_Pauli)+xn*dabs(A_Pauli-Pauli(pn+1,hn+1,pp-1,hp-1,gp,gn)))*            &
                         ((Ex_tot-B_Pauli)/(Ex_tot-A_Pauli))**(n_tot-1)
            B_Pauli = max(Pauli(pn,hn,pp,hp,gp,gn),Pauli(pn-1,hn-1,pp+1,hp+1,gp,gn))
            L_0_np_an = 2.0d0*pi*Msq_np*(real(pn*hn,kind=8)/xn)*gp**2*                                         &
                        finite_well(p_tot,h_tot,Ex_tot,h_max,Vwell)*                                           &
                        (2.0d0*(Ex_tot-B_Pauli)+xn*dabs(A_Pauli-Pauli(pn-1,hn-1,pp+1,hp+1,gp,gn)))*            &
                        ((Ex_tot-B_Pauli)/(Ex_tot-A_Pauli))**(n_tot-1)



            Lamb_p_p(pn,pp) = L_p_p_an/hbar
            Lamb_p_n(pn,pp) = L_p_n_an/hbar
            Lamb_0_pn(pn,pp) = L_0_pn_an/hbar
            Lamb_0_np(pn,pp) = L_0_np_an/hbar
         end if

     
         do kk = 1, nucleus(icomp)%num_decay                       !   Loop over particle allowed to decay from this nucleus
            k = nucleus(icomp)%decay_particle(kk)
            if(k > 0 .and. pn == p0(1) .and. pp == p0(2))cycle        !  only photons as there are no hole states, and thus emission 
                                                                   ! from here is just elastic scattering
            ifinal = nucleus(icomp)%decay_to(kk)


            z_k = particle(k)%Z
            n_k = particle(k)%A-z_k
            if(pn < max(p0(1),n_k))cycle                                  !  not enough excitons - cycle
            if(pp < max(p0(2),z_k))cycle
            Sep = nucleus(icomp)%sep_e(k)
            Zf = nucleus(ifinal)%Z
            Af = nucleus(ifinal)%A
            Nf = Af - Zf
            gnnf = real(Nf,kind=8)/Preeq_g_div
            gppf = real(Zf,kind=8)/Preeq_g_div
            ggf = gnnf + gppf
            if(k == iproj)then
               spin_target = nucleus(ifinal)%state(istate)%spin
            else
               spin_target = nucleus(ifinal)%state(1)%spin
            end if
            spin_proj = particle(k)%spin
            sp1 = max((spin_target+spin_proj),abs(spin_target-spin_proj))
            sp2 = max(particle(k)%lmax+sp1,abs(particle(k)%lmax-sp1))
            j_max = nint(sp2-nucleus(1)%jshift)
            if(j_max > nucleus(ifinal)%j_max)j_max = nucleus(ifinal)%j_max
            mass_i = particle(k)%Mass
            e_max = ex_tot - nucleus(icomp)%sep_e(k)

            e_cut = nucleus(ifinal)%level_param(7)
            e_bin = e_max

            nbin_end = min(int(e_bin/de),nucleus(icomp)%nbin_part)

            factor1 = barn_eq_fmsq/((pi*hbar_c)**2*hbar)

            do m = 1, nbin_end                        !  loop over output energies
               energy = real(m,kind=8)*de
               if(energy - Coulomb_Barrier(k) <= 1.0d-6)cycle
               Ex_final = Ex_tot - energy - Sep
               Delta = Delta_pre(pn_res,hn,pp_res,hp,Zf,Af,gpf,gnf,Ex_final)
               gnf = aparam_u(Ex_final, gnnf, shell, gamma)
               gpf = aparam_u(Ex_final, gppf, shell, gamma)
               gf = gnf + gpf
               term = 0.0d0
               semi_direct = 0.0d0
               direct = 0.0d0
               om = 0.0d0
               om1 = 0.0d0
               om2 = 0.0d0
               if(k > 0)then
                  mass_t = nucleus(ifinal)%Mass + ex_final
                  mass_rel = mass_i*mass_t/(mass_t+mass_i)
                  sig_inv = 0.0d0
                  do ipar = 0, 1                                              !   parity of compound nucleus
                     do j = 0, j_max                                          !   loop over J values
                        xji = real(j,kind=8) + nucleus(icomp)%jshift
                        sig_inv = sig_inv + compound_cs(energy,ipar,xji,ifinal,1,k)
                     end do
                  end do
                  factor = factor1*(2.0d0*spin_proj+1.0d0)*mass_rel*energy*sig_inv
                  pn_res = pn - n_k
                  pp_res = pp - z_k
                  hn_res = hn
                  hp_res = hp
                  Delta = Delta_pre(pn_res,hn_res,pp_res,hp_res,Zf,Af,gpf,gnf,Ex_final)
                  om = omega2(pn_res,hn_res,pp_res,hp_res,Zf,Af,gpf,gnf,Ex_final,Delta,h_max,Vwell)
                  term = om
               else
!----    Inverse cross section, uses photo-absorption
!----    TALYS uses total, including all multipoles
!                  sig_inv = photo_absorption(ifinal, nucleus(ifinal)%lmax_E, nucleus(ifinal)%lmax_M,    &
!                                             energy, Ex_final)
!----   Assumption is that pre-equilibrium photon emission is just E1
!----   Since E1 i s
                  sig_inv = EL_absorption(ifinal, 1, energy, nucleus(icomp)%sep_e(k)) 
                  factor = factor1*energy**2*sig_inv
!----    Semi-direct component
                  semi_direct = 0.0d0
                  if(n_tot >= n0_tot + 2)then
                     pn_res = pn - 1
                     pp_res = pp
                     hn_res = hn - 1
                     hp_res = hp
                     Delta = Delta_pre(pn_res,hn_res,pp_res,hp_res,Zf,Af,gpf,gnf,Ex_final)
                     om1 = omega2(pn_res,hn_res,pp_res,hp_res,Zf,Af,gpf,gnf,Ex_final,Delta,h_max,Vwell_g)
                     pn_res = pn
                     pp_res = pp - 1
                     hn_res = hn
                     hp_res = hp - 1
                     Delta = Delta_pre(pn_res,hn_res,pp_res,hp_res,Zf,Af,gpf,gnf,Ex_final)
                     om2 = omega2(pn_res,hn_res,pp_res,hp_res,Zf,Af,gpf,gnf,Ex_final,Delta,h_max,Vwell_g)
                     semi_direct = 0.5d0*gf*energy*(om1 + om2)/(xn - 2.0d0 + gf*energy)
                  end if
!----   Direct component 
                  direct = 0.0d0
                  pn_res = pn
                  pp_res = pp
                  hn_res = hn
                  hp_res = hp
                  Delta = Delta_pre(pn_res,hn_res,pp_res,hp_res,Zf,Af,gpf,gnf,Ex_final)
                  om = omega2(pn_res,hn_res,pp_res,hp_res,Zf,Af,gpf,gnf,Ex_final,Delta,h_max,Vwell_g)
                  direct = xn*om/(xn + gf*energy)
                  term = Preeq_gam_fact*(direct + semi_direct)
               end if

               dWk(pn,pp,m,k) = factor*term/omdenom                            !   decay rate for each particle type and energy

               Wk(pn,pp,k) = Wk(pn,pp,k) + dWk(pn,pp,m,k)*de                          !   its integral
            end do

            W(pn,pp) = W(pn,pp) + Wk(pn,pp,k)

         end do
         
         denom = Lamb_p_p(pn,pp) + Lamb_p_n(pn,pp) +                       &
                 Lamb_0_pn(pn,pp) + Lamb_0_np(pn,pp) + W(pn,pp)
         if(denom > 0.0d0)then
            tau(pn,pp) = 1.0d0/denom
            Gam_p(pn,pp) = Lamb_p_p(pn,pp)*tau(pn,pp)
            Gam_n(pn,pp) = Lamb_p_n(pn,pp)*tau(pn,pp)
            Gam_0_pn(pn,pp) = Lamb_0_pn(pn,pp)*tau(pn,pp)
            Gam_0_np(pn,pp) = Lamb_0_np(pn,pp)*tau(pn,pp)
         else
            tau(pn,pp) = 1.0d20
            Gam_p(pn,pp) = 0.0d0
            Gam_n(pn,pp) = 0.0d0
            Gam_0_pn(pn,pp) = 0.0d0
            Gam_0_np(pn,pp) = 0.0d0
         end if
         denom = Lamb_p_p(pn,pp) + Lamb_p_n(pn,pp) + W(pn,pp)
         if(denom > 0.0d0)then
            taup(pn,pp)=1.0d0/denom
            Gam_pp(pn,pp) = Lamb_p_p(pn,pp)*taup(pn,pp)
            Gam_pn(pn,pp) = Lamb_p_n(pn,pp)*taup(pn,pp)
         else
            taup(pn,pp) = 1.0d20
            Gam_pp(pn,pp) = 0.0d0
            Gam_pn(pn,pp) = 0.0d0
         end if
      end do
   end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Surviving pre-equilibrium flux
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Pre_Prob(p0(1),p0(2)) = 1.0d0
   do pn = p0(1), pn_calc
      do pp = p0(2), pp_calc
         hn = pn - p0(1)
         if((pn == p0(1) .and. pp == p0(2)) .or. pp < p0(2) .or. pp > pp_max .or. pn > pn_max)cycle
         hp = pp - p0(2)
         Pre_Prob(pn,pp) = Prob_func(pn,pp,pn_max,pp_max,                 &
                                     Gam_n,Gam_p,Gam_pn,Gam_pp,           &
                                     Gam_0_pn,Gam_0_np,Pre_Prob)
      end do
   end do

   nucleus(icomp)%PREEQ_cs(in) = 0.0d0
   nucleus(icomp)%PREEQ_part_cs(1:nucleus(icomp)%num_decay,in) = 0.0d0
!
   stest2 = 0.0d0
   sum2 = 0.0d0
   pre_eq_cs(0:6) = 0.0d0
   stest = 0.0d0
   do kk = 1,nucleus(icomp)%num_decay                       !   Loop over particle allowed to decay from this nucleus
      k = nucleus(icomp)%decay_particle(kk)
      z_k = particle(k)%Z
      n_k = particle(k)%A - z_k
      ifinal = nucleus(icomp)%decay_to(kk)
      e_max = ex_tot - nucleus(icomp)%sep_e(k)
      e_cut = nucleus(ifinal)%level_param(7)
      e_bin = e_max
      nbin_end = min(int(e_bin/de),nucleus(icomp)%nbin_part)
      do m = 0, nbin_end                           !  loop over output energies
         energy = real(m,kind=8)*de
         if(energy - Coulomb_Barrier(k) <= 1.0d-6)cycle
         nucleus(icomp)%PREEQ_part_spectrum(kk,m)=0.0d0
         ex_final = ex_tot - energy - nucleus(icomp)%sep_e(k)
         pn_min = max(p0(1),n_k)
         pp_min = max(p0(2),z_k)
         sum = 0.0d0
         do pn = p0(1), pn_calc
            do pp = p0(2), pp_calc
               hn = pn - p0(1)
               hp = pp - p0(2)
               n_tot = pn + hn + pp + hp
               if(k > 0 .and. pn == p0(1) .and. pp == p0(2))cycle                       !  Particle emission only for n_tot >= 3
               sum = sum + dWk(pn,pp,m,k)*tau(pn,pp)*Pre_Prob(pn,pp)


           end do
         end do
         nucleus(icomp)%PREEQ_part_spectrum(kk,m) = reac_cs*sum
         pre_eq_cs(k) = pre_eq_cs(k) + nucleus(icomp)%PREEQ_part_spectrum(kk,m)*de
      end do

      sum = 0.0d0
      do m = 0, nbin_end
         sum = sum + nucleus(icomp)%PREEQ_part_spectrum(kk,m)*de
      end do
      if(sum > 0.0d0)then
        do m = 0, nbin_end
            nucleus(icomp)%PREEQ_part_spectrum(kk,m) =                            &
              nucleus(icomp)%PREEQ_part_spectrum(kk,m)/sum
         end do
      end if
      nucleus(icomp)%PREEQ_part_cs(kk,in) = pre_eq_cs(k)
      nucleus(icomp)%PREEQ_cs(in) = nucleus(icomp)%PREEQ_cs(in) + pre_eq_cs(k)
   end do

   return
end subroutine pre_equilibrium_1
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function Well(h,A,energy,V1,V3,K)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function calculates well depth
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
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
!----------  Input Data          -------------------------------------
   use variable_kinds
   use constants
   implicit none
   integer(kind=4), intent(in) :: h, A
   real(kind=8), intent(in) :: energy, V1, V3, K
!----------  Internal Data     ---------------------------------------
   real(kind=8) xA,V2
!----------  Start Calculation
   xA = real(A,kind=8)**(1.0d0/3.0d0)
   V2 = V3 - V1
   if(h <= 1)then
      Well = V1 + V2*energy**4/(energy**4 + (K/xA)**4)
   else
      Well = V3
   end if
   return
end function Well
!
!*******************************************************************************
!
real(kind=8) function omega2(pn,hn,pp,hp,Z,A,gp,gn,Ex,Delta,h_max,Vwell)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function calculates the exciton level density
!
!  Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: Pauli
!        real(kind=8) :: Finite_well
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
!----------  Input Data          -------------------------------------
   use variable_kinds
   use constants
   implicit none
   integer(kind=4), intent(in) :: pn,hn,pp,hp,Z,A
   real(kind=8), intent(in) :: gp,gn
   real(kind=8), intent(in) :: Ex, Delta
   integer(kind=4), intent(in) :: h_max
   real(kind=8), intent(in) :: Vwell(0:h_max)
!----------  Internal Data       -------------------------------------
   integer(kind=4) :: N
   integer(kind=4) :: pn_tot,pp_tot,p_tot,h_tot,n_tot
   real(kind=8) :: Ux,Uxx
   real(kind=8) :: res1,res2
   real(kind=8) :: APauli,fp,fn
!--------   External Functions   -------------------------------------
   real(kind=8) :: Pauli
   real(kind=8) :: Finite_well
!----------  Start Calculation   -------------------------------------
   N = A - Z
   omega2 = 0.0d0
   Ux = Ex - Delta
   if(Ux <= 0.0d0)return
   if(pn + hn + pp + hp == 0) return
   if(pn < 0 .or. hn < 0 .or. pp < 0 .or. hp < 0)return
   pn_tot = pn + hn
   pp_tot = pp + hp
   p_tot = pn + pp
   h_tot = hn + hp
   n_tot = pn + hn + pp + hp
   fp = real(pp*2+hp**2 + pp + hp,kind=8)/(4.0d0*gp)
   fn = real(pn*2+hn**2 + pn + hn,kind=8)/(4.0d0*gn)
   APauli = Pauli(pn,hn,pp,hp,gp,gn)
   if(APauli >= Ux)return
   Uxx = Ux - APauli                         !  Ux-pauli defines omega - it shouldn't be negative
   res1 = gn**pn_tot*gp**pp_tot/                                                        &
         (factorial(pn)*factorial(hn)*factorial(pp)*factorial(hp)*factorial(n_tot-1))*  &
         Uxx**(n_tot-1)
   res2 = Finite_well(p_tot,h_tot,Ux,h_max,Vwell)
   omega2 = res1*res2
   return
end function omega2
!
!*******************************************************************************
!
real(kind=8) function Delta_pre(pn,hn,pp,hp,Z,A,gp,gn,Ex)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function calculates the Pre-equilibrium Pairing term 
!
!  Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
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
   use options
   implicit none
!--------   Input data
   integer(kind=4), intent(in) :: pn,hn,pp,hp,Z,A
   real(kind=8), intent(in) :: gp,gn,Ex
!--------  Internal Data       -------------------------------------
   integer(kind=4) :: N
   integer(kind=4) :: n_tot
   real(kind=8) :: g
   real(kind=8) :: xA
   real(kind=8) :: Delta,Delta_ex,ratio
   real(kind=8) :: T_crit,n_crit,xn, xdelta
!--------   Start Calculation    -------------------------------------
   Delta_pre = 0.0d0
   if(preeq_pair_model < 0)return
   if(preeq_pair_model == 2)then
      Delta_pre = preeq_delta
      return
   end if
   g = gp + gn
   n_tot = pn + pp + hn + hp
   Delta = 0.0d0
   xA = real(A,kind=8)
   N = A - Z
   if(iand(Z,1) == 1 .and. iand(N,1) == 1)then
      Delta = 0.0d0
      Delta_pre = Delta
      return
   elseif(iand(Z,1) == 0 .and. iand(N,1) == 1)then
      Delta=12.0d0/dsqrt(xA)
   elseif(iand(Z,1) == 1 .and. iand(N,1) == 0)then
      Delta=12.0d0/dsqrt(xA)
   else
      Delta=24.0d0/dsqrt(xA)
   end if
   Delta_pre = Delta
   if(preeq_pair_model == 0)return      
   T_crit = 2.0d0*dsqrt(Delta/(0.25d0*g))/3.5d0
   n_crit = 2.0d0*g*T_crit*log(2.0)
   xn = real(n_tot,kind=8)/n_crit
   ratio = 0.716d0+2.44d0*xn**2.17d0
   xdelta = (0.996-1.76*xn**1.6/(Ex/Delta)**0.68)**2
   Delta_ex = 0.0d0
   if(xdelta <= 1.0d0)Delta_ex = Delta*xdelta
   Delta_pre = Delta - Delta_ex
   return
end function Delta_pre
!
!*******************************************************************************
!
real(kind=8) function Pauli(pn,hn,pp,hp,gp,gn)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function calculates the Pauli Correction Factor
!
!  Dependencies:
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
   implicit none
!--------   Input data
   integer(kind=4), intent(in) :: pn,hn,pp,hp
   real(kind=8), intent(in) :: gp,gn
!--------   Start Calculation    -------------------------------------
   Pauli = real(max(pp,hp),kind=8)**2/gp + real(max(pn,hn),kind=8)**2/gn -      &
           real(pp**2 + hp**2 + pp + hp,kind=8)/(4.0d0*gp) -               &
           real(pn**2 + hn**2 + pn + hn,kind=8)/(4.0d0*gn)
   return
end function Pauli
!
!*******************************************************************************
!
real(kind=8) function finite_well(p,h,Ex,h_max,Vwell)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function calculates the Finte Well function 
!
!  Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        constants
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
   use options
   use constants
   implicit none
   integer(kind=4), intent(in) :: p,h
   real(kind=8), intent(in) :: Ex
   integer(kind=4), intent(in) :: h_max
   real(kind=8), intent(in) :: Vwell(0:h_max)
!------   Internal Variables   ---------------------------------------
   integer(kind=4) :: i,n
   real(kind=8) :: sum
   real(kind=8) :: Ux
   real(kind=8) :: Step
   real(kind=8) :: V, VV, R
!------   Start Calculation    ---------------------------------------
!------   Orignal model uses step function, which causes a kink to occur in the 
!------   outgoing spectrum, especially if the well depth is lowered. Made a modification
!------   to smooth out transition. Note that factor basically adds and subtracts ((U-i*V)/U)**(n-1)
!------   to the finite_well factor for U-i*V >. Essentially:
!------   finite_well =    1.0                                      U - i*V <= 0
!------                    1.0 + Sum_i x(i)*((U-i*V)/U)**(n-1).     U - i*V > 0
!------   The leading term in omega is also ((U-i*V)/U)**(n-1), so this factor subtracts off
!------   the increase in omega, flattening it out, but with a sharp transition.
!------   Note that this same effect can be accomplished using 
!------   ((U-VV)/UU)**(n-1) where VV = U for U -i*V <= 0 and VV = i*V for U - i*V > 0.
!------   The overall effect can be smoothed by making a smooth transition from VV = U to VV = i*V
!------   near U = iV. Do this here by make a attaching the straight linear line VV = U to a circle of
!------   radius R and then continuing along the circle to the top, where U = iV.
!------
   finite_well = 1.0d0
   if(h < 1)return
   R = 2.0d0
   n = p + h
   V = Vwell(h)
   sum = 0.0d0
   do i = 1, h
      VV = real(i,kind=8)*V
      Ux = Ex - VV
      Step = 0.0d0
      if(Ux > 0.0d0)Step = 1.0d0
      sum = sum + (-1.0d0)**i*factorial(h)/(factorial(i)*factorial(h-i))*                   &
                  (Ux/Ex)**(n-1)*Step
   end do
   finite_well = 1.0d0 + sum
   return
end function finite_well
!
!*******************************************************************************
!
real(kind=8) function Prob_func(pn,pp,pn_max,pp_max,                      &    
                                Gam_n,Gam_p,Gam_pn,Gam_pp,                &
                                Gam_0_pn,Gam_0_np,Pre_Prob)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function calculates the pre-equilibrium flux surviving previous
!    emission. Solve recurssively. Note depends only on previous
!    probabilities with 1 less total exciton  
!
!  Dependencies:
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
   implicit none
!---------   input variables   -----------------------------
   integer(kind=4), intent(in) :: pn,pp,pn_max,pp_max
   real(kind=8), intent(in) :: Gam_n(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Gam_p(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Gam_pn(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Gam_pp(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Gam_0_pn(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Gam_0_np(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Pre_Prob(0:pn_max,0:pp_max)
!---------------------    Internal Data        -----------------------------
   real(kind=8) sum
!---------------------    Start Calculation    -----------------------------
   sum=0.0d0
   if(pp-1 >= 0)sum = sum + Pre_Prob(pn,pp-1)*Gam_p(pn,pp-1)
   if(pn-1 >= 0)sum = sum + Pre_Prob(pn-1,pp)*Gam_n(pn-1,pp)
   if(pp-2 >= 0 .and. pn+1 <= pn_max)sum = sum + (Pre_Prob(pn+1,pp-2)*Gam_pp(pn+1,pp-2) +  &
                                                  Pre_Prob(pn,pp-1)*Gam_pn(pn,pp-1))*      &
                                                  Gam_0_np(pn+1,pp-1)
   if(pn-2 >= 0 .and. pp+1 <= pp_max)sum = sum + (Pre_Prob(pn-1,pp)*Gam_pp(pn-1,pp) +      &
                                                  Pre_Prob(pn-2,pp+1)*Gam_pn(pn-2,pp+1))*  &
                                                  Gam_0_pn(pn-1,pp+1)
   Prob_func=sum
   return
end function Prob_func
!
!*******************************************************************************
!
subroutine int_trans_rate(pn,hn,pp,hp,Z,A,gp,gn,Ex,Delta,h_max,Vwell,          &
                          Msq_nn,Msq_pp,Msq_pn,Msq_np,                         &
                          Lamb_p_p,Lamb_p_n,Lamb_0_pn,Lamb_0_np)
!
!*******************************************************************************
!
!  Discussion:
!
!  This subroutine calculates the internal transition rate for one exciton
!  state to another
!
!  Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: omega2
!        real(kind=8) :: Pauli
!
!     MPI routines:
!
!        None
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
!----------  Input Data          -------------------------------------
   integer(kind=4), intent(in) :: pn,hn,pp,hp
   integer(kind=4), intent(in) :: Z,A
   real(kind=8), intent(in) :: gp,gn
   real(kind=8), intent(in) :: Ex, Delta
   integer(kind=4), intent(in) :: h_max
   real(kind=8), intent(in) :: Vwell(0:h_max)  
   real(kind=8), intent(in) :: Msq_nn,Msq_pp,Msq_pn,Msq_np
!---------------------------------------------------------------------
   real(kind=8), intent(out) :: Lamb_p_n
   real(kind=8), intent(out) :: Lamb_p_p
   real(kind=8), intent(out) :: Lamb_0_pn
   real(kind=8), intent(out) :: Lamb_0_np
!----------  Internal Data       -------------------------------------
   integer(kind=4) :: i
   integer(kind=4) :: npoints
   real(kind=8) :: Lp_pp_1,Lp_pp_2
   real(kind=8) :: Lh_pp_1,Lh_pp_2
   real(kind=8) :: Lp_np_1,Lp_np_2
   real(kind=8) :: Lh_np_1,Lh_np_2
   real(kind=8) :: Lp_nn_1,Lp_nn_2
   real(kind=8) :: Lh_nn_1,Lh_nn_2
   real(kind=8) :: Lp_pn_1,Lp_pn_2
   real(kind=8) :: Lh_pn_1,Lh_pn_2
   real(kind=8) :: L0_pn_1,L0_pn_2
   real(kind=8) :: L0_np_1,L0_np_2
   real(kind=8) :: de_p_pp,de_h_pp,de_p_np,de_h_np
   real(kind=8) :: de_p_nn,de_h_nn,de_p_pn,de_h_pn
   real(kind=8) :: de_0_pn,de_0_np
   real(kind=8) :: u_p_pp,u_h_pp,u_p_np,u_h_np
   real(kind=8) :: u_p_nn,u_h_nn,u_p_pn,u_h_pn
   real(kind=8) :: u_0_pn,u_0_np
   real(kind=8) :: lamb_p_pp,lamb_h_pp
   real(kind=8) :: lamb_p_np,lamb_h_np
   real(kind=8) :: lamb_p_nn,lamb_h_nn
   real(kind=8) :: lamb_p_pn,lamb_h_pn
   real(kind=8) :: l_0_pn,l_0_np
   real(kind=8) :: xl_1p_pp,xl_1h_pp,xl_1p_np,xl_1h_np,xl_1p1h_pn
   real(kind=8) :: xl_1p_nn,xl_1h_nn,xl_1p_pn,xl_1h_pn,xl_1p1h_np
   real(kind=8) :: om_denom
   real(kind=8) :: U
   real(kind=8) :: dee
   !----------  External Functions      ---------------------------------
   real(kind=8) :: omega2
   real(kind=8) :: Pauli
!-----------   Start Calculation    ----------------------------------
   lamb_p_pp=0.0d0
   lamb_h_pp=0.0d0
   lamb_p_np=0.0d0
   lamb_h_np=0.0d0
   lamb_p_nn=0.0d0
   lamb_h_nn=0.0d0
   lamb_p_pn=0.0d0
   lamb_h_pn=0.0d0
   l_0_pn=0.0d0
   l_0_np=0.0d0

   Lamb_p_p = 0.0d0
   Lamb_p_n = 0.0d0
   Lamb_0_pn = 0.0d0
   Lamb_0_np = 0.0d0

   U = Ex - Delta
   if(U <= 0.0d0)return
   dee = 0.01d0
   npoints = nint(U/dee)

   om_denom = omega2(pn,hn,pp,hp,Z,A,gp,gn,Ex,Delta,h_max,Vwell)

   de_p_pp = 0.0d0
   de_h_pp = 0.0d0
   de_p_np = 0.0d0
   de_h_np = 0.0d0
   de_p_nn = 0.0d0
   de_h_nn = 0.0d0
   de_p_pn = 0.0d0
   de_h_pn = 0.0d0
   de_0_pn = 0.0d0
   de_0_np = 0.0d0
!----   Set up endpoints for numerical integrations

!------------------------------------------------ Limits for Lambda^+_p
   Lp_pp_1 = max(Pauli(pn,hn,pp+1,hp+1,gp,gn) - Pauli(pn,hn,pp-1,hp,gp,gn),0.0d0)
   Lp_pp_2 = max(Ex - Pauli(pn,hn,pp-1,hp,gp,gn),Lp_pp_1)
   if(Lp_pp_2 > Lp_pp_1)de_p_pp = (Lp_pp_2-Lp_pp_1)/real(npoints+1,kind=8)

   Lh_pp_1 = max(Pauli(pn,hn,pp+1,hp+1,gp,gn) - Pauli(pn,hn,pp,hp-1,gp,gn),0.0d0)
   Lh_pp_2 = max(Ex - Pauli(pn,hn,pp,hp-1,gp,gn),Lh_pp_1)
   if(Lh_pp_2 > Lh_pp_1)de_h_pp = (Lh_pp_2-Lh_pp_1)/real(npoints+1,kind=8)

   Lp_np_1 = max(Pauli(pn+1,hn+1,pp,hp,gp,gn) - Pauli(pn-1,hn,pp,hp,gp,gn),0.0d0)
   Lp_np_2 = max(Ex - Pauli(pn-1,hn,pp,hp,gp,gn),Lp_np_1)
   if(Lp_np_2 > Lp_np_1)de_p_np = (Lp_np_2-Lp_np_1)/real(npoints+1,kind=8)

   Lh_np_1 = max(Pauli(pn,hn,pp+1,hp+1,gp,gn) - Pauli(pn,hn-1,pp,hp,gp,gn),0.0d0)
   Lh_np_2 = max(Ex - Pauli(pn,hn-1,pp,hp,gp,gn),Lh_np_1)
   if(Lh_np_2 > Lh_np_1)de_h_np = (Lh_np_2-Lh_np_1)/real(npoints+1,kind=8)
!------------------------------------------------ Limits for Lambda^+_p
   Lp_nn_1 = max(Pauli(pn+1,hn+1,pp,hp,gp,gn) - Pauli(pn-1,hn,pp,hp,gp,gn),0.0d0)
   Lp_nn_2 = max(Ex - Pauli(pn-1,hn,pp,hp,gp,gn),Lp_nn_1)
   if(Lp_nn_2 > Lp_nn_1)de_p_nn = (Lp_nn_2-Lp_nn_1)/real(npoints+1,kind=8)

   Lh_nn_1 = max(Pauli(pn+1,hn+1,pp,hp,gp,gn) - Pauli(pn,hn-1,pp,hp,gp,gn),0.0d0)
   Lh_nn_2 = max(Ex - Pauli(pn,hn-1,pp,hp,gp,gn),lh_pp_1)
   if(Lh_nn_2 > Lh_nn_1)de_h_nn = (Lh_nn_2-Lh_nn_1)/real(npoints+1,kind=8)

   Lp_pn_1 = max(Pauli(pn,hn,pp+1,hp+1,gp,gn) - Pauli(pn,hn,pp-1,hp,gp,gn),0.0d0)
   Lp_pn_2 = max(Ex - Pauli(pn,hn,pp-1,hp,gp,gn),Lp_pn_1)
   if(Lp_pn_2 > Lp_pn_1)de_p_pn = (Lp_pn_2-Lp_pn_1)/real(npoints+1,kind=8)

   Lh_pn_1 = max(Pauli(pn+1,hn+1,pp,hp,gp,gn) - Pauli(pn,hn,pp,hp-1,gp,gn),0.0d0)
   Lh_pn_2 = max(Ex - Pauli(pn,hn,pp,hp-1,gp,gn),Lh_pn_1)
   if(Lh_pn_2 > Lh_pn_1)de_h_pn = (Lh_pn_2-Lh_pn_1)/real(npoints+1,kind=8)
!------------------------------------------------ Limits for Lambda^0_pn
   L0_pn_1 = max(Pauli(pn,hn,pp,hp,gp,gn) - Pauli(pn,hn,pp-1,hp-1,gp,gn),0.0d0)
   L0_pn_2 = max(ex - Pauli(pn,hn,pp-1,hp-1,gp,gn),L0_pn_1)
   if(L0_pn_2 > L0_pn_1)de_0_pn = (L0_pn_2-L0_pn_1)/real(npoints+1,kind=8)
!------------------------------------------------ Limits for Lambda^0_np
   L0_np_1 = max(Pauli(pn,hn,pp,hp,gp,gn) - Pauli(pn-1,hn-1,pp,hp,gp,gn),0.0d0)
   L0_np_2 = max(Ex - Pauli(pn-1,hn-1,pp,hp,gp,gn),L0_np_1)
   if(L0_np_2 > L0_np_1)de_0_np = (L0_np_2-L0_np_1)/real(npoints+1,kind=8)
!------------------------------------------------  Initialize integration terms
   do i = 1, npoints
!------------------------------------------------Lambda^+_p
      u_p_pp = real(i-1,kind=8)*de_p_pp + Lp_pp_1
      xl_1p_pp = 2.0d0*pi*Msq_pp*omega2(0,0,2,1,Z,A,gp,gn,u_p_pp,Delta,h_max,Vwell)/hbar
      lamb_p_pp = lamb_p_pp + omega2(pn,hn,pp-1,hp,Z,A,gp,gn,Ex-u_p_pp,Delta,h_max,Vwell)*     &
                              omega2(0,0,1,0,Z,A,gp,gn,u_p_pp,Delta,h_max,Vwell)*              &
                              xl_1p_pp*de_p_pp
      u_h_pp = real(i-1,kind=8)*de_h_pp + Lh_pp_1
      xl_1h_pp = 2.0d0*pi*Msq_pp*omega2(0,0,1,2,Z,A,gp,gn,u_h_pp,Delta,h_max,Vwell)/hbar
      lamb_h_pp = lamb_h_pp + omega2(pn,hn,pp,hp-1,Z,A,gp,gn,Ex-u_h_pp,Delta,h_max,Vwell)*     &
                              omega2(0,0,0,1,Z,A,gp,gn,u_h_pp,Delta,h_max,Vwell)*              &
                              xl_1h_pp*de_h_pp
      u_p_np = real(i-1,kind=8)*de_p_np + Lp_np_1
      xl_1p_np = 2.0d0*pi*Msq_np*omega2(1,0,1,1,Z,A,gp,gn,u_p_np,Delta,h_max,Vwell)/hbar
      lamb_p_np = lamb_p_np + omega2(pn-1,hn,pp,hp,Z,A,gp,gn,Ex-u_p_np,Delta,h_max,Vwell)*     &
                              omega2(1,0,0,0,Z,A,gp,gn,u_p_np,Delta,h_max,Vwell)*              &
                              xl_1p_np*de_p_np
      u_h_np = real(i-1,kind=8)*de_h_np + Lh_np_1
      xl_1h_np = 2.0d0*pi*Msq_np*omega2(0,1,1,1,Z,A,gp,gn,u_h_np,Delta,h_max,Vwell)/hbar
      lamb_h_np = lamb_h_np + omega2(pn,hn-1,pp,hp,Z,A,gp,gn,Ex-u_h_np,Delta,h_max,Vwell)*     &
                              omega2(0,1,0,0,Z,A,gp,gn,u_h_np,Delta,h_max,Vwell)*              &
                              xl_1h_np*de_h_np
!------------------------------------------------Lambda^+_n
      u_p_nn = real(i-1,kind=8)*de_p_nn + Lp_nn_1
      xl_1p_nn = 2.0d0*pi*Msq_nn*omega2(2,1,0,0,Z,A,gp,gn,u_p_nn,Delta,h_max,Vwell)/hbar
      lamb_p_nn = lamb_p_nn + omega2(pn-1,hn,pp,hp,Z,A,gp,gn,Ex-u_p_nn,Delta,h_max,Vwell)*     &
                              omega2(1,0,0,0,Z,A,gp,gn,u_p_nn,Delta,h_max,Vwell)*              &
                              xl_1p_nn*de_p_nn
      u_h_nn = real(i-1,kind=8)*de_h_nn + Lh_nn_1
      xl_1h_nn = 2.0d0*pi*Msq_nn*omega2(1,2,0,0,Z,A,gp,gn,u_h_nn,Delta,h_max,Vwell)/hbar
      lamb_h_nn = lamb_h_nn + omega2(pn,hn-1,pp,hp,Z,A,gp,gn,Ex-u_h_nn,Delta,h_max,Vwell)*     &
                              omega2(0,1,0,0,Z,A,gp,gn,u_h_nn,Delta,h_max,Vwell)*              &
                              xl_1h_nn*de_h_nn
      u_p_pn = real(i-1,kind=8)*de_p_pn + Lp_pn_1
      xl_1p_pn = 2.0d0*pi*Msq_pn*omega2(1,1,1,0,Z,A,gp,gn,u_p_pn,Delta,h_max,Vwell)/hbar
      lamb_p_pn = lamb_p_pn + omega2(pn,hn,pp-1,hp,Z,A,gp,gn,Ex-u_p_pn,Delta,h_max,Vwell)*     &
                              omega2(0,0,1,0,Z,A,gp,gn,u_p_pn,Delta,h_max,Vwell)*              &
                              xl_1p_pn*de_p_pn
      u_h_pn = real(i-1,kind=8)*de_h_pn + Lh_pn_1
      xl_1h_pn = 2.0d0*pi*Msq_pn*omega2(1,1,0,1,Z,A,gp,gn,u_h_pn,Delta,h_max,Vwell)/hbar
      lamb_h_pn = lamb_h_pn + omega2(pn,hn,pp,hp-1,Z,A,gp,gn,Ex-u_h_pn,Delta,h_max,Vwell)*     &
                              omega2(0,0,0,1,Z,A,gp,gn,u_h_pn,Delta,h_max,Vwell)*              &
                              gp*xl_1h_pn*de_h_pn
!-------------------------------------------------Lambda^0_pn
      u_0_pn = real(i-1,kind=8)*de_0_pn + L0_pn_1
      xl_1p1h_pn = 2.0d0*pi*Msq_pn*omega2(1,1,0,0,Z,A,gp,gn,u_0_pn,Delta,h_max,Vwell)/hbar
      l_0_pn = l_0_pn + omega2(pn,hn,pp-1,hp-1,Z,A,gp,gn,Ex-u_0_pn,Delta,h_max,Vwell)*         &
                        omega2(0,0,1,1,Z,A,gp,gn,u_0_pn,Delta,h_max,Vwell)*                    &
                        xl_1p1h_pn*de_0_pn
!-------------------------------------------------Lambda^0_np
      u_0_np = real(i-1,kind=8)*de_0_np + L0_np_1
      xl_1p1h_np = 2.0d0*pi*Msq_np*omega2(0,0,1,1,Z,A,gp,gn,u_0_np,Delta,h_max,Vwell)/hbar
      l_0_np = l_0_np + omega2(pn-1,hn-1,pp,hp,Z,A,gp,gn,Ex-u_0_np,Delta,h_max,Vwell)*         &
                        omega2(1,1,0,0,Z,A,gp,gn,u_0_np,Delta,h_max,Vwell)*                    &
                        xl_1p1h_np*de_0_np

   end do

   Lamb_p_p = (lamb_p_pp + lamb_h_pp + lamb_p_np + lamb_h_np)/om_denom
   Lamb_p_n = (lamb_p_nn + lamb_h_nn + lamb_p_pn + lamb_h_pn)/om_denom
   Lamb_0_pn = l_0_pn/om_denom
   Lamb_0_np = l_0_np/om_denom

   return
end subroutine int_trans_rate
