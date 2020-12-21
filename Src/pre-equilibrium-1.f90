!
!*******************************************************************************
!
subroutine pre_equilibrium_1(icomp,istate,in,E_inc,    &
                             Ex_tot,de,reac_cs)
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
   integer(kind=4) :: i,j,k,m,kk
   integer(kind=4) :: h
   integer(kind=4) :: pn,pp,hn,hp,pn_tot,pp_tot
   integer(kind=4) :: pn_res,pp_res
   integer(kind=4) :: pn_calc, pp_calc
   integer(kind=4) :: p_tot,h_tot,n_tot
   integer(kind=4) :: pn_min,pp_min
   integer(kind=4) :: z_k,n_k,Z,N,A,Ap,Zf,Nf,Af
   integer(kind=4) :: j_max
   integer(kind=4) :: ipar
   real(kind=8) :: pre_eq_cs(0:6)

   real(kind=8) :: xn
   real(kind=8) :: energy, U
   real(kind=8) :: sig_inv
   real(kind=8) :: sp1,sp2
   real(kind=8) :: spin_target,spin_proj
   real(kind=8) :: mass_i,mass_t,mass_rel
   real(kind=8) :: factor,factor1
   real(kind=8) :: Msq
!   real(kind=8) :: R_pp,R_nn,R_pn,R_np
   real(kind=8) :: xA
   real(kind=8) :: xAp
   real(kind=8) :: ratio
   real(kind=8) :: Msq_pp,Msq_pn,Msq_np,Msq_nn
   real(kind=8) :: sum,sum2
   real(kind=8) :: L_p_p_an,L_p_n_an,L_0_pn_an,L_0_np_an
   real(kind=8) :: A_Pauli,B_Pauli
   real(kind=8) :: denom
   real(kind=8) :: omdenom,om
   real(kind=8) :: Sep
   real(kind=8) :: gp,gn,g
   real(kind=8) :: gpp,gnn,gg
   real(kind=8) :: gpf,gnf,gf
   real(kind=8) :: gppf,gnnf,ggf
   real(kind=8) :: Delta

   real(kind=8) :: V1, V3, xK
   real(kind=8), allocatable :: Vwell(:)

   real(kind=8) :: stest,stest2

   integer(kind=4) :: nbin_f,nbin_end
   real(kind=8)    :: ex_final,e_max,e_bin
   real(kind=8)    :: xji
   integer(kind=4) :: h_max

   real(kind=8)    :: frac
   integer(kind=4) :: ifinal
   real(kind=8) :: e_cut

   real(kind=8) :: shell, gamma

   integer(kind=4) :: ifile
   character(len=20) :: outfile

!--------   External Functions   ---------------------------
   real(kind=8) :: compound_cs
   real(kind=8) :: omega2
   real(kind=8) :: Well
   real(kind=8) :: Delta_pre
   real(kind=8) :: Pauli
   real(kind=8) :: Prob_func
   real(kind=8) :: finite_well
   real(kind=8) :: aparam_u
!--------   Start Calculation    ---------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------    Calculate internal transition rates
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   j_max=nucleus(icomp)%j_max
   iproj=projectile%particle_type
   Z=nucleus(icomp)%Z
   A=nucleus(icomp)%A
   xA=dfloat(A)
   Ap=projectile%A
   xAp=dfloat(Ap)
   N=A-Z
!   R_pp = 1.0d0
!   R_nn = 1.0d0
!   R_pn = 1.5d0
!   R_np = 1.5d0
   if(Preeq_g_a)then
      g = nucleus(icomp)%a_Sn*pi**2/6.0d0
      Preeq_g_div = xA/g
   end if
   gnn = dfloat(N)/Preeq_g_div
   gpp = dfloat(Z)/Preeq_g_div
   gg = gpp + gnn
   pn_min = p0(1)
   pp_min = p0(2)

   pn_calc = pn_max
   pp_calc = pp_max
   allocate(dWk(0:pn_max+1,0:pp_max+1,0:nucleus(1)%nbin_part,0:6))
   allocate(Wk(0:pn_max+1,0:pp_max+1,0:6))
   allocate(W(0:pn_max+1,0:pp_max+1))
   allocate(tau(0:pn_max+1,0:pp_max+1))
   allocate(taup(0:pn_max+1,0:pp_max+1))
   allocate(Pre_Prob(0:pn_max+1,0:pp_max+1))
   allocate(Gam_p(0:pn_max+1,0:pp_max+1))
   allocate(Gam_n(0:pn_max+1,0:pp_max+1))
   allocate(Gam_pp(0:pn_max+1,0:pp_max+1))
   allocate(Gam_pn(0:pn_max+1,0:pp_max+1))
   allocate(Gam_0_pn(0:pn_max+1,0:pp_max+1))
   allocate(Gam_0_np(0:pn_max+1,0:pp_max+1))
   allocate(Lamb_p_p(0:pn_max+1,0:pp_max+1))
   allocate(Lamb_p_n(0:pn_max+1,0:pp_max+1))
   allocate(Lamb_0_pn(0:pn_max+1,0:pp_max+1))
   allocate(Lamb_0_np(0:pn_max+1,0:pp_max+1))

   h_max = pn_max + pp_max
   allocate(Vwell(0:h_max))

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
!write(52,*)Ex_tot
!   do pn=pn_min,pn_max
!      do pp=pp_min,pp_max
!   shell = nucleus(icomp)%level_param(4)
!   gamma = nucleus(icomp)%level_param(5)
   shell = 0.0d0
   gamma = 0.0d0
   gp = aparam_u(Ex_tot,gpp,shell,gamma)
   gn = aparam_u(Ex_tot,gnn,shell,gamma)
   g = gp + gn
!   V1 = Preeq_V_part(iproj)
   V1 = Preeq_V1
   V3 = Preeq_V
   xK = Preeq_K
   if(v1 > V3) V1 = V3
   do h = 0, h_max
      if(preeq_fwell == 0 .or. preeq_fwell == 2)Vwell(h) = Well(h,A,E_inc,V1,V3,xK)
      if(preeq_fwell == 1)Vwell(h) = V3
   end do
!   do i = p0(1) + p0(2), max(pn_max,pp_max)
!      do pn = p0(1), min(i,pn_max)
   do pn = p0(1), pn_max
      do pp = p0(2), pp_max
!         pp = i - pn
         hn = pn - p0(1)
         hp = pp - p0(2)
         pn_tot = pn + hn
         pp_tot = pp + hp
         p_tot = pn + pp
         h_tot = hn + hp
         n_tot = p_tot + h_tot
!         if(preeq_fwell == 0 .or. preeq_fwell == 2)V_well=Well(h_tot,A,E_inc,V1,V3)
!         if(preeq_fwell == 1)V_well = V3
!         if(preeq_fwell == 1)V_well = 38.0d0
         Lamb_p_p(pn,pp) = 0.0d0
         Lamb_p_n(pn,pp) = 0.0d0
         Lamb_0_pn(pn,pp) = 0.0d0
         Lamb_0_np(pn,pp) = 0.0d0

         Delta = Delta_pre(pn,hn,pp,hp,Z,A,gp,gn,Ex_tot)
!         U = max(Ex_tot - Delta,1.0d-1)
!         U = Ex_tot - Delta
!         if(U <= 0.0d0)cycle

         omdenom = omega2(pn,hn,pp,hp,Z,A,gp,gn,Ex_tot,Delta,H_max,Vwell)
!         omdenom = omega2(pn,hn,pp,hp,Z,A,gp,gn,Ex_tot,Delta,E_inc,H_max,Vwell)
!         omdenom = omega(pn,hn,pp,hp,Z,A,gp,gn,Ex_tot,E_inc,V1,V3)

!    write(60,'(''ph'',4(1x,i5),1x,e15.7)')pn,pp,hn,hp,omdenom
         if(omdenom <= 1.0d-10)cycle
         xn = dfloat(n_tot)
         Msq = (M2_C1*xAp/xA**3)*(7.48d0*M2_C2+4.62d5/(Ex_tot/(xn*xAp)+10.7*M2_C3)**3)
         if(analytic_preeq)Msq=1.2*Msq
!   write(31,'(1x,f10.3,1x,e15.7)')Ex_tot,Msq
         Msq_pp = M2_Rpp*Msq
         Msq_pn = M2_Rpn*Msq
         Msq_nn = M2_Rnn*Msq
         Msq_np = M2_Rnp*Msq
         if(n_tot > 1 .and. .not. analytic_preeq)then
            call int_trans_rate(pn,hn,pp,hp,Z,A,gp,gn,Ex_tot,Delta,h_max,Vwell,      &
                                Msq_nn,Msq_pp,Msq_pn,Msq_np,                         &
                                Lamb_p_p(pn,pp),Lamb_p_n(pn,pp),                     &
                                Lamb_0_pn(pn,pp),Lamb_0_np(pn,pp))
!            call int_trans_rate(pn,hn,pp,hp,pn_max,pp_max,Z,A,gp,gn,Ex_tot,Delta,E_inc,h_max,Vwell,   &
!                                Msq_nn,Msq_pp,Msq_pn,Msq_np,                                          &
!                                Lamb_p_p(pn,pp),Lamb_p_n(pn,pp),                                      &
!                                Lamb_0_pn(pn,pp),Lamb_0_np(pn,pp))
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
!            U = max(Ex_tot - Delta,1.0d-1)
            U = Ex_tot - Delta
            if(U <= 0.0d0)cycle

!            write(6,*)V_well,finite_well(p_tot+1,h_tot+1,Ex_tot,V_well)


            L_p_p_an=2.0d0*pi*gp**2/(2.0d0*xn*(xn+1.0d0))*((Ex_tot-Pauli(pn,hn,pp+1,hp+1,gp,gn))**(n_tot+1)/           & 
                                                           (Ex_tot-Pauli(pn,hn,pp,hp,gp,gn))**(n_tot-1))*              &
                                                           (dfloat(pp+hp)*gp*Msq_pp+2.0d0*dfloat(pn+hn)*gn*Msq_pn)*    &
                                                           finite_well(p_tot+1,h_tot+1,Ex_tot,h_max,Vwell)
!                                                           finite_well(p_tot+1,h_tot+1,Ex_tot,A,E_inc,h_max,Vwell)
            L_p_n_an=2.0d0*pi*gn**2/(2.0d0*xn*(xn+1.0d0))*((Ex_tot-Pauli(pn+1,hn+1,pp,hp,gp,gn))**(n_tot+1)/           &
                                                           (Ex_tot-Pauli(pn,hn,pp,hp,gp,gn))**(n_tot-1))*              &
                                                           (dfloat(pn+hn)*gn*Msq_nn+2.0d0*dfloat(pp+hp)*gp*Msq_np)*    &
                                                           finite_well(p_tot+1,h_tot+1,Ex_tot,h_max,Vwell)
!                                                           finite_well(p_tot+1,h_tot+1,Ex_tot,A,E_inc,h_max,Vwell)
            A_Pauli = Pauli(pn,hn,pp,hp,gp,gn)
            B_Pauli = max(Pauli(pn,hn,pp,hp,gp,gn),Pauli(pn+1,hn+1,pp-1,hp-1,gp,gn))
            L_0_pn_an = 2.0d0*pi*Msq_pn*(dfloat(pp)*dfloat(hp)/xn)*gn**2*                                              &
                        finite_well(p_tot,h_tot,Ex_tot,h_max,Vwell)*                                           &
!                        finite_well(p_tot,h_tot,Ex_tot,A,E_inc,h_max,Vwell)*                                           &
                        (2.0d0*(Ex_tot-B_Pauli)+xn*dabs(A_Pauli-Pauli(pn+1,hn+1,pp-1,hp-1,gp,gn)))*                    &
                         ((Ex_tot-B_Pauli)/(Ex_tot-A_Pauli))**(n_tot-1)
            B_Pauli = max(Pauli(pn,hn,pp,hp,gp,gn),Pauli(pn-1,hn-1,pp+1,hp+1,gp,gn))
            L_0_np_an = 2.0d0*pi*Msq_np*(dfloat(pn)*dfloat(hn)/xn)*gp**2*                                              &
                        finite_well(p_tot,h_tot,Ex_tot,h_max,Vwell)*                                           &
!                        finite_well(p_tot,h_tot,Ex_tot,A,E_inc,h_max,Vwell)*                                           &
                        (2.0d0*(Ex_tot-B_Pauli)+xn*dabs(A_Pauli-Pauli(pn-1,hn-1,pp+1,hp+1,gp,gn)))*                    &
                        ((Ex_tot-B_Pauli)/(Ex_tot-A_Pauli))**(n_tot-1)


!            L_p_p_an = 0.0d0
!            if(U-Pauli(pn,hn,pp+1,hp+1,gp,gn) > 0.0d0 .and.           &
!               U-Pauli(pn,hn,pp,hp,gp,gn) > 0.0d0)then
!               L_p_p_an=2.0d0*pi*gp**2/(2.0d0*xn*(xn+1.0d0))*((U-Pauli(pn,hn,pp+1,hp+1,gp,gn))**(n_tot+1)/                & 
!                                                              (U-Pauli(pn,hn,pp,hp,gp,gn))**(n_tot-1))*                   &
!                                                              (dfloat(pp+hp)*gp*Msq_pp+2.0d0*dfloat(pn+hn)*gn*Msq_pn)*    &
!                                                              finite_well(p_tot+1,h_tot+1,Ex_tot,V_well)
!            end if
!            L_p_n_an = 0.0d0
!            if(Ex_tot-Pauli(pn+1,hn+1,pp,hp,gp,gn) > 0.0d0 .and.           &
!               Ex_tot-Pauli(pn,hn,pp,hp,gp,gn) > 0.0d0)then
!               L_p_n_an=2.0d0*pi*gn**2/(2.0d0*xn*(xn+1.0d0))*((U-Pauli(pn+1,hn+1,pp,hp,gp,gn))**(n_tot+1)/                &
!                                                              (U-Pauli(pn,hn,pp,hp,gp,gn))**(n_tot-1))*                   &
!                                                              (dfloat(pn+hn)*gn*Msq_nn+2.0d0*dfloat(pp+hp)*gp*Msq_np)*    &
!                                                              finite_well(p_tot+1,h_tot+1,Ex_tot,V_well)
!!            end if
!            A_Pauli = Pauli(pn,hn,pp,hp,gp,gn)
!            B_Pauli = max(Pauli(pn,hn,pp,hp,gp,gn),Pauli(pn+1,hn+1,pp-1,hp-1,gp,gn))
!            if(U-B_Pauli > 0.0d0 .and. U-A_Pauli > 0.0d0)                                                              &
!               L_0_pn_an=2.0d0*pi*Msq_pn*(dfloat(pp)*dfloat(hp)/xn)*gn**2*finite_well(p_tot,h_tot,Ex_tot,V_well)*      &
!                         (2.0d0*(U-B_Pauli)+xn*dabs(A_Pauli-Pauli(pn+1,hn+1,pp-1,hp-1,gp,gn)))*                        &
!                         ((U-B_Pauli)/(U-A_Pauli))**(n_tot-1)
!            B_Pauli = max(Pauli(pn,hn,pp,hp,gp,gn),Pauli(pn-1,hn-1,pp+1,hp+1,gp,gn))
!            if(U-B_Pauli > 0.0d0 .and. U-A_Pauli > 0.0d0)                                                              &
!               L_0_np_an=2.0d0*pi*Msq_np*(dfloat(pn)*dfloat(hn)/xn)*gp**2*finite_well(p_tot,h_tot,Ex_tot,V_well)*      &
!                         (2.0d0*(U-B_Pauli)+xn*dabs(A_Pauli-Pauli(pn-1,hn-1,pp+1,hp+1,gp,gn)))*                        &
!                         ((U-B_Pauli)/(U-A_Pauli))**(n_tot-1)



            Lamb_p_p(pn,pp) = L_p_p_an/hbar
            Lamb_p_n(pn,pp) = L_p_n_an/hbar
            Lamb_0_pn(pn,pp) = L_0_pn_an/hbar
            Lamb_0_np(pn,pp) = L_0_np_an/hbar
         end if

!write(50,'(4(1x,i2),4(1x,e15.7))')pp,hp,pn,hn,Lamb_p_p(pn,pp),Lamb_p_n(pn,pp),Lamb_0_pn(pn,pp),Lamb_0_np(pn,pp)
!   if(pn == 1 .and. pp == 0)then
!      write(32,'(1x,f10.3,4(1x,e15.7))')Ex_tot,Lamb_p_p(pn,pp),Lamb_p_n(pn,pp),     &
!                                               L_p_p_an,L_p_n_an
!      write(6,'(1x,f10.3,4(1x,e15.7))')Ex_tot,Lamb_p_p(pn,pp),Lamb_p_n(pn,pp),      &
!                                               L_p_p_an,L_p_n_an
!   end if
!   if(pn == 2 .and. pp == 0)then
!      write(32,'(1x,f10.3,4(1x,e15.7))')Ex_tot,Lamb_p_p(pn,pp),Lamb_p_n(pn,pp),     &
!                                               L_p_p_an,L_p_n_an
!   end if
!-------------------   Calculate dWk and W for each particle type
!    write(60,'(4(1x,i4),1xf10.4,1x,e15.7)')pn,pp,hn,hp,Ex_tot,omdenom
     
         do kk=1,nucleus(icomp)%num_decay                       !   Loop over particle allowed to decay from this nucleus
            k=nucleus(icomp)%decay_particle(kk)
            if(k == 0)cycle
            if(k > 0 .and. pn == p0(1) .and. pp == p0(2))cycle        !  only photons as there are no hole states, and thus emission 
                                                                   ! from here is just elastic scattering
            ifinal = nucleus(icomp)%decay_to(kk)
!            e_cut=nucleus(icomp)%state(nucleus(icomp)%num_discrete)%energy
            z_k=particle(k)%Z
            n_k=particle(k)%A-z_k
            if(pn < max(p0(1),n_k))cycle                                  !  not enough excitons - cycle
            if(pp < max(p0(2),z_k))cycle
            Sep=nucleus(icomp)%sep_e(k)
!            Zf=nucleus(icomp)%Z-z_k
!            Af=nucleus(icomp)%A-(z_k+n_k)
            Zf=nucleus(ifinal)%Z
            Af=nucleus(ifinal)%A
            Nf = Af - Zf
            gnnf = dfloat(Nf)/Preeq_g_div
            gppf = dfloat(Zf)/Preeq_g_div
            ggf = gnnf + gppf
!            shell = nucleus(ifinal)%level_param(4)
!            gamma = nucleus(ifinal)%level_param(5)
            spin_target = nucleus(ifinal)%state(istate)%spin
            spin_proj = particle(k)%spin
            sp1 = max((spin_target+spin_proj),abs(spin_target-spin_proj))
            sp2 = max(particle(k)%lmax+sp1,abs(particle(k)%lmax-sp1))
            j_max = nint(sp2-nucleus(1)%jshift)
            if(j_max > nucleus(ifinal)%j_max)j_max = nucleus(ifinal)%j_max
            mass_i = particle(k)%Mass
!            mass_i = (particle(k)%ME + particle(k)%A*mass_u)/mass_u
            e_max = ex_tot - nucleus(icomp)%sep_e(k)

            e_cut = nucleus(ifinal)%level_param(7)
            e_bin = e_max

            nbin_end = min(int(e_bin/de),nucleus(icomp)%nbin_part)

            do m = 1, nbin_end                        !  loop over output energies
               energy = dfloat(m)*de
               Ex_final = Ex_tot - energy - Sep
               if(energy > Ex_tot - nucleus(icomp)%sep_e(k)+de/2.0)exit         !  Check if there is enough energy to emit this particle
               Delta = Delta_pre(pn_res,hn,pp_res,hp,Zf,Af,gpf,gnf,Ex_final)
!               U = Ex_final - Delta
!               if(U <= 0.0d0)cycle
               gnf = aparam_u(Ex_final, gnnf, shell, gamma)
               gpf = aparam_u(Ex_final, gppf, shell, gamma)
               gf = gnf + gpf
               mass_t = nucleus(ifinal)%Mass + ex_final
!               mass_t=(nucleus(ifinal)%ME+nucleus(ifinal)%A*mass_u+ex_final)/mass_u
               mass_rel=mass_i*mass_t/(mass_t+mass_i)
!               mass_rel=mass_i*mass_t/(mass_t+mass_i)
!               mass_rel=mass_rel*mass_u
               sig_inv=0.0d0
               do ipar=0,1                                              !   parity of compound nucleus
                  do j=0,j_max                                          !   loop over J values
                     xji=dfloat(j)+nucleus(icomp)%jshift
                     sig_inv = sig_inv + compound_cs(energy,ipar,xji,ifinal,1,k)
                  end do
               end do
               factor1 = (2.0d0*spin_proj+1.0d0)/(pi*hbar_c)**2*mass_rel*barn_eq_fmsq/hbar
               factor = factor1*energy*sig_inv
               pn_res = pn - n_k
               pp_res = pp - z_k
!               U = Ex_final - Delta
               om = omega2(pn_res,hn,pp_res,hp,Zf,Af,gpf,gnf,Ex_final,Delta,h_max,Vwell)
!               om = omega2(pn_res,hn,pp_res,hp,Zf,Af,gpf,gnf,Ex_final,Delta,E_inc,h_max,Vwell)
!               om = omega(pn_res,hn,pp_res,hp,Zf,Af,gpf,gnf,Ex_final,E_inc,V1,V3)
!           om = omega(pn_res,hn,pp_res,hp,Zf,Af,gpf,gnf,Ex_tot - energy,V_well)

               ratio = om/omdenom
               dWk(pn,pp,m,k) = factor*ratio                                      !   decay rate for each particle type and energy
!     if(hn + hp == 2)write(12,'(3(1x,f8.4),5(1x,e15.7))')energy,Ex_final,Vwell(hn+hp),om,ratio,factor,factor*ratio

!      n_tot2 = pn_res + pp_res + hn + hp
!      write(60,'(4(1x,f10.3),7(1x,i5),2(1x,e15.7),6(1x,i5),4(1x,e15.7))')                       &
!               energy,Ex_tot,Ex_final,Delta,m,k,pn,hn,pp,hp,n_tot,Vwell(h_tot),omdenom,               &
!               pn_res,hn,pp_res,hp,h_tot,n_tot2,Vwell(h_tot),om,ratio,dWk(pn,pp,m,k)

!      write(60,'(4(1x,f10.3),9(1x,i5),4(1x,e15.7))')energy,Ex_tot,Ex_final,Delta,m,k,pn,pp,pn_res,pp_res,hn,hp,h_tot,      &
!               Vwell(h),dWK(pn,pp,m,k),ratio,om
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
!   do i = p0(1) + p0(2) + 1, max(pn_max,pp_max)
!      do pn = p0(1), min(i,pn_max)
   do pn=p0(1), pn_calc
      do pp=p0(2), pp_calc
         hn = pn - p0(1)
!         pp = i - pn
!         if(pn+pp > 6)cycle
         if((pn == p0(1) .and. pp == p0(2)) .or. pp < p0(2) .or. pp > pp_max .or. pn > pn_max)cycle
         hp = pp - p0(2)
!         Pre_Prob(pn,pp) = Prob_func(pn,hn,pp,hp,p0(1),p0(2),pn_max,pp_max,      &
!                                     Gam_n,Gam_p,Gam_pn,Gam_pp,                  &
!                                     Gam_0_pn,Gam_0_np,Pre_Prob)
         Pre_Prob(pn,pp) = Prob_func(pn,pp,pn_max,pp_max,                 &
                                     Gam_n,Gam_p,Gam_pn,Gam_pp,           &
                                     Gam_0_pn,Gam_0_np,Pre_Prob)
      end do
   end do
!   j=60
!   outfile(1:4)='Cont'
!   do i=p0(1)+p0(2),max(pn_max,pp_max)
!      do pn=p0(1),min(i,pn_max)
!   do pn=p0(1), pn_calc
!      do pp=p0(2), pp_calc
!         pp=i-pn
!         hn=pn-p0(1)
!         hp=pp-p0(2)
!         n_tot=pn+hn+pp+hp
!         if( pn == p0(1) .and. pp == p0(2))cycle                       !  Particle emission only for n_tot >= 3
         
!         write(outfile(5:5),'(i1)')pn
!         write(outfile(6:6),'(i1)')hn
!         write(outfile(7:7),'(i1)')pp
!         write(outfile(8:8),'(i1)')hp
!         j=j+1
!         open(unit=j,file=outfile(1:8)//'.dat',status='unknown')
!      end do
!   end do
!   j = j + 1
!   open(unit=j,file = outfile(1:4)//'-total.dat',status='unknown')

   nucleus(icomp)%PREEQ_cs(in)=0.0d0
   nucleus(icomp)%PREEQ_part_cs(1:nucleus(icomp)%num_decay,in)=0.0d0
!
   stest2=0.0d0
   sum2=0.0d0
   pre_eq_cs(0:6)=0.0d0
   stest=0.0d0
   do kk=1,nucleus(icomp)%num_decay                       !   Loop over particle allowed to decay from this nucleus
      k=nucleus(icomp)%decay_particle(kk)
      z_k=particle(k)%Z
      n_k=particle(k)%A-z_k
      ifinal=nucleus(icomp)%decay_to(kk)
!      e_cut=nucleus(ifinal)%state(nucleus(ifinal)%num_discrete)%energy
      if(k == 0)cycle                                     !  We don't do photons yet
      e_max = ex_tot-nucleus(icomp)%sep_e(k)
      e_cut = nucleus(ifinal)%level_param(7)
      e_bin = e_max
      nbin_end = min(int(e_bin/de),nucleus(icomp)%nbin_part)
!   write(6,*)ifinal,e_cut,nucleus(ifinal)%num_discrete
      do m = 0, nbin_end                           !  loop over output energies
         energy = dfloat(m)*de
         if(energy > ex_tot - nucleus(icomp)%sep_e(k)+de/2.0) exit
         ex_final = ex_tot - energy - nucleus(icomp)%sep_e(k)
         pn_min = max(p0(1),n_k)
         pp_min = max(p0(2),z_k)
         nucleus(icomp)%PREEQ_part_spectrum(kk,m)=0.0d0
         sum = 0.0d0
!         j = 60
!         do i = p0(1) + p0(2), max(pn_max,pp_max)
!            do pn = p0(1), min(i,pn_max)
         do pn = p0(1), pn_calc
            do pp = p0(2), pp_calc
!               pp = i - pn
               hn = pn - p0(1)
               hp = pp - p0(2)
               n_tot = pn + hn + pp + hp
               if(k > 0 .and. pn == p0(1) .and. pp == p0(2))cycle                       !  Particle emission only for n_tot >= 3
               sum = sum + dWk(pn,pp,m,k)*tau(pn,pp)*Pre_Prob(pn,pp)
!               j = j + 1
!               if(k ==1)write(j,'(1x,f10.5,4(1x,e15.7),5(1x,i3))')            &
!                  energy,reac_cs*dWk(pn,pp,m,k)*tau(pn,pp)*Pre_Prob(pn,pp),   &
!                  reac_cs*dWk(pn,pp,m,k),tau(pn,pp),Pre_Prob(pn,pp),pn,hn,pp,hp,m
           end do
         end do
!         j = j + 1
!         if(k == 1)write(j,'(1x,f10.5,1x,e15.7)')energy,reac_cs*sum
         if(ex_final > e_cut)then
            nbin_f = nint((ex_final-nucleus(ifinal)%e_grid(1))/de)
            nucleus(icomp)%PREEQ_part_spectrum(kk,m) = reac_cs*sum
            pre_eq_cs(k) = pre_eq_cs(k) + reac_cs*sum*de
         else
            frac = 0.0d0
            do i = 1, nucleus(ifinal)%num_discrete
!               if(i == istate)cycle
!          write(6,'(2i5,3(1x,f10.5))')                                       &
!            m,i,nucleus(ifinal)%state(i)%energy,ex_final,ex_final+de
               if(nucleus(ifinal)%state(i)%energy >= ex_final .and.                     &
                  nucleus(ifinal)%state(i)%energy <= ex_final + de)frac = frac + 1.0d0
            end do
!   write(6,*)'frac',frac,sum
            if(frac > 0.0d0)then
               frac = 1.0d0/frac
               do i = 1, nucleus(ifinal)%num_discrete
!                  if(i == istate)cycle
                  if(nucleus(ifinal)%state(i)%energy >= ex_final .and.                  &
                     nucleus(ifinal)%state(i)%energy < ex_final + de )then
                        nucleus(icomp)%PREEQ_part_spectrum(kk,m) =                      &
                           nucleus(icomp)%PREEQ_part_spectrum(kk,m) + reac_cs*sum*frac
                        pre_eq_cs(k) = pre_eq_cs(k) + reac_cs*sum*de*frac
                  end if
               end do
            end if
         end if
      end do

!   j=60
!   do i=p0(1)+p0(2),max(pn_max,pp_max)
!      do pn=p0(1),min(i,pn_max)
!         close(unit=j)
!      end do
!   end do

      sum = 0.0d0
      do m = 0, nbin_end
         sum = sum + nucleus(icomp)%PREEQ_part_spectrum(kk,m)*de
      end do
      if(sum > 0.0d0)then
!   prob = 0.0d0
!   write(6,*)'kk = ',kk
         do m = 0, nbin_end
            nucleus(icomp)%PREEQ_part_spectrum(kk,m) =                            &
              nucleus(icomp)%PREEQ_part_spectrum(kk,m)/sum
!   prob = prob + nucleus(icomp)%PREEQ_part_spectrum(kk,m)*de
!   write(60+kk,'(1x,f10.4,1x,1pe15.7)')dfloat(m)*de,prob
         end do
      end if
!      write(6,*)'preeq',kk,pre_eq_cs(k)
      nucleus(icomp)%PREEQ_part_cs(kk,in) = pre_eq_cs(k)
      nucleus(icomp)%PREEQ_cs(in) = nucleus(icomp)%PREEQ_cs(in) + pre_eq_cs(k)
   end do

!-------   Now we need to remove the pre-equilibrium cross section from the
!-------   initial compound-nucleus cross section. This is done by assuming
!-------   a J-distribution according to the level density, in pops(j,ip).
!-------   This is then removed from each entrance channel according to the
!-------   fractional amount of the J,ip cross section in that S, l channel.

!   itarget = target%icomp
!   spin_target=nucleus(itarget)%state(istate)%spin
!   spin_proj=particle(iproj)%spin
!   sp1=abs(spin_target-spin_proj)
!   sp2=spin_target+spin_proj
!   isp=int(sp2-sp1)
!   par=particle(iproj)%par*nucleus(itarget)%state(istate)%parity
!   stest2=0.0d0
!   do ip=0,1
!      cpar=2*ip-1
!      do j=0,nucleus(icomp)%j_max
!         xj = real(j,kind=8) + nucleus(icomp)%jshift
!         nucleus(icomp)%PREEQ_pop(j,ip)=              &
!              nucleus(icomp)%PREEQ_cs(in)*pops(j,ip)/de
!!         nucleus(icomp)%pop(j,ip,nbin_i)=nucleus(icomp)%pop(j,ip,nbin_i)-   &
!!               nucleus(icomp)%PREEQ_cs(in)*pops(j,ip)/de
!         preeq_pop = nucleus(icomp)%PREEQ_cs(in)*pops(j,ip)/de
!         do i=0,isp                                         !   loop over channel spins
!            S=dfloat(i)+sp1                                 !   Channel spin
!            lmin=abs(xj-S)
!            lmax=xj+S
!            lmin=min(lmin,particle(iproj)%lmax)        !   min and max l-values
!            lmax=min(lmax,particle(iproj)%lmax)
!            do l=lmin,lmax,1                                !   loop over angular momentum
!               if(par*(-1)**l.ne.cpar .or. target%pop_xjpi(j,ip) < 1.0d-15)cycle
!                  frac_ch = target%pop_channel(j,ip,i,l)/target%pop_xjpi(j,ip)!
!                  target%pop_channel(j,ip,i,l) = target%pop_channel(j,ip,i,l) - &
!                      preeq_pop*frac_ch
!!   write(28,'(2(1x,i4),1x,f4.1,1x,i4,3(1x,f10.6))')j,ip,S,l,preeq_pop,frac_ch,target%pop_channel(j,ip,i,l)
!            end do
!         end do
!         target%pop_xjpi(j,ip) = target%pop_xjpi(j,ip) - preeq_pop
!!         stest2=stest2+nucleus(icomp)%pop(j,ip,nbin_i)
!      end do
!   end do
!write(6,*)stest,stest2,stest+stest2
!write(6,*)stest,nucleus(icomp)%PREEQ_cs(in)
!write(6,*)'check',stest2+nucleus(icomp)%PREEQ_cs(in),reac_cs

!   write(6,*)'FINISHED PRREQ'

!   deallocate(pops)

   deallocate(dWk)
   deallocate(Wk)
   deallocate(W)
   deallocate(tau)
   deallocate(taup)
   deallocate(Pre_Prob)
   deallocate(Gam_p)
   deallocate(Gam_n)
   deallocate(Gam_pp)
   deallocate(Gam_pn)
   deallocate(Gam_0_pn)
   deallocate(Gam_0_np)
   deallocate(Lamb_p_p)
   deallocate(Lamb_p_n)
   deallocate(Lamb_0_pn)
   deallocate(Lamb_0_np)

   deallocate(Vwell)

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
!----------  Input Data          -------------------------------------
   use variable_kinds
   use constants
   implicit none
   integer(kind=4), intent(in) :: h, A
   real(kind=8), intent(in) :: energy, V1, V3, K
!----------  Internal Data     ---------------------------------------
   real(kind=8) xA,V2
!----------  Start Calculation
   xA=dfloat(A)**(1.0d0/3.0d0)
   V2 = V3 - V1
   if(h <= 1)then
      Well = V1 + V2*energy**4/(energy**4 + (K/xA)**4)
!      Well = V1 + V2*energy**4/(energy**4 + (245.0d0/xA)**4)
   else
      Well = V3
   end if
!   write(21,*)k,h,energy,Well
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
!----------  Input Data          -------------------------------------
   use variable_kinds
   use constants
   use pre_equilibrium_no_1
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
   fp = dfloat(pp*2+hp**2 + pp + hp)/(4.0d0*gp)
   fn = dfloat(pn*2+hn**2 + pn + hn)/(4.0d0*gn)
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
   xA = dfloat(A)
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
   xn = dfloat(n_tot)/n_crit
   ratio = 0.716d0+2.44d0*xn**2.17d0
   xdelta = (0.996-1.76*xn**1.6/(Ex/Delta)**0.68)**2
   Delta_ex = 0.0d0
   if(xdelta <= 1.0d0)Delta_ex = Delta*xdelta
   Delta_pre = Delta - Delta_ex
!   write(20,'(4(1x,f10.5))')Ex,Delta,Delta_ex,Delta_pre
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
!--------   Input data
   integer(kind=4), intent(in) :: pn,hn,pp,hp
   real(kind=8), intent(in) :: gp,gn
!--------   Start Calculation    -------------------------------------
   Pauli = dfloat(max(pp,hp))**2/gp + dfloat(max(pn,hn))**2/gn -      &
           dfloat(pp**2 + hp**2 + pp + hp)/(4.0d0*gp) -               &
           dfloat(pn**2 + hn**2 + pn + hn)/(4.0d0*gn)
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
   real(kind=8) :: V,VV,R
!------------   External Functions
!---unused   real(kind=8) :: Well
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
!      V0 = real(i,kind=8)*V
!      if(n > 1)then
!         VV = Ex
!         X0 = V0 + (sqrt(2.0d0) - 1.0d0)*R
!         if(Ex > (X0 - R/sqrt(2.0d0)) .and. Ex < X0 )then
!            VV = V0 - R + sqrt(R**2 - (Ex - V0 + (1.0d0 - sqrt(2.0d0))*R)**2)
!         elseif(Ex >= X0)then
!            VV = V0
!         end if
!         Ux = Ex - VV
!         sum = sum + (-1.0d0)**i*factorial(h)/(factorial(i)*factorial(h-i))*    &
!                     (Ux/Ex)**(n-1)
!      else
!         if((Ex - V0) > 0.0d0)sum = sum + (-1.0d0)**i*factorial(h)/(factorial(i)*factorial(h-i))
!      end if
   end do
!   finite_well = finite_well + sum
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
!---------   input variables   -----------------------------
   integer(kind=4), intent(in) :: pn,pp,pn_max,pp_max
   real(kind=8), intent(in) :: Gam_n(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Gam_p(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Gam_pn(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Gam_pp(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Gam_0_pn(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Gam_0_np(0:pn_max,0:pp_max)
   real(kind=8), intent(in) :: Pre_Prob(0:pn_max,0:pp_max)
!---------------------    Internal Data 
   real(kind=8) sum
!---------------------    Start Calculation
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
!    This subroutine calculates the internal exciton transition rates  
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
   if(Lp_pp_2 > Lp_pp_1)de_p_pp = (Lp_pp_2-Lp_pp_1)/dfloat(npoints+1)

   Lh_pp_1 = max(Pauli(pn,hn,pp+1,hp+1,gp,gn) - Pauli(pn,hn,pp,hp-1,gp,gn),0.0d0)
   Lh_pp_2 = max(Ex - Pauli(pn,hn,pp,hp-1,gp,gn),Lh_pp_1)
   if(Lh_pp_2 > Lh_pp_1)de_h_pp = (Lh_pp_2-Lh_pp_1)/dfloat(npoints+1)

   Lp_np_1 = max(Pauli(pn+1,hn+1,pp,hp,gp,gn) - Pauli(pn-1,hn,pp,hp,gp,gn),0.0d0)
   Lp_np_2 = max(Ex - Pauli(pn-1,hn,pp,hp,gp,gn),Lp_np_1)
   if(Lp_np_2 > Lp_np_1)de_p_np = (Lp_np_2-Lp_np_1)/dfloat(npoints+1)

   Lh_np_1 = max(Pauli(pn,hn,pp+1,hp+1,gp,gn) - Pauli(pn,hn-1,pp,hp,gp,gn),0.0d0)
   Lh_np_2 = max(Ex - Pauli(pn,hn-1,pp,hp,gp,gn),Lh_np_1)
   if(Lh_np_2 > Lh_np_1)de_h_np = (Lh_np_2-Lh_np_1)/dfloat(npoints+1)
!------------------------------------------------ Limits for Lambda^+_p
   Lp_nn_1 = max(Pauli(pn+1,hn+1,pp,hp,gp,gn) - Pauli(pn-1,hn,pp,hp,gp,gn),0.0d0)
   Lp_nn_2 = max(Ex - Pauli(pn-1,hn,pp,hp,gp,gn),Lp_nn_1)
   if(Lp_nn_2 > Lp_nn_1)de_p_nn = (Lp_nn_2-Lp_nn_1)/dfloat(npoints+1)

   Lh_nn_1 = max(Pauli(pn+1,hn+1,pp,hp,gp,gn) - Pauli(pn,hn-1,pp,hp,gp,gn),0.0d0)
   Lh_nn_2 = max(Ex - Pauli(pn,hn-1,pp,hp,gp,gn),lh_pp_1)
   if(Lh_nn_2 > Lh_nn_1)de_h_nn = (Lh_nn_2-Lh_nn_1)/dfloat(npoints+1)

   Lp_pn_1 = max(Pauli(pn,hn,pp+1,hp+1,gp,gn) - Pauli(pn,hn,pp-1,hp,gp,gn),0.0d0)
   Lp_pn_2 = max(Ex - Pauli(pn,hn,pp-1,hp,gp,gn),Lp_pn_1)
   if(Lp_pn_2 > Lp_pn_1)de_p_pn = (Lp_pn_2-Lp_pn_1)/dfloat(npoints+1)

   Lh_pn_1 = max(Pauli(pn+1,hn+1,pp,hp,gp,gn) - Pauli(pn,hn,pp,hp-1,gp,gn),0.0d0)
   Lh_pn_2 = max(Ex - Pauli(pn,hn,pp,hp-1,gp,gn),Lh_pn_1)
   if(Lh_pn_2 > Lh_pn_1)de_h_pn = (Lh_pn_2-Lh_pn_1)/dfloat(npoints+1)
!------------------------------------------------ Limits for Lambda^0_pn
   L0_pn_1 = max(Pauli(pn,hn,pp,hp,gp,gn) - Pauli(pn,hn,pp-1,hp-1,gp,gn),0.0d0)
   L0_pn_2 = max(ex - Pauli(pn,hn,pp-1,hp-1,gp,gn),L0_pn_1)
   if(L0_pn_2 > L0_pn_1)de_0_pn = (L0_pn_2-L0_pn_1)/dfloat(npoints+1)
!------------------------------------------------ Limits for Lambda^0_np
   L0_np_1 = max(Pauli(pn,hn,pp,hp,gp,gn) - Pauli(pn-1,hn-1,pp,hp,gp,gn),0.0d0)
   L0_np_2 = max(Ex - Pauli(pn-1,hn-1,pp,hp,gp,gn),L0_np_1)
   if(L0_np_2 > L0_np_1)de_0_np = (L0_np_2-L0_np_1)/dfloat(npoints+1)
!------------------------------------------------  Initialize integration terms
   do i = 1, npoints
!------------------------------------------------Lambda^+_p
      u_p_pp = dfloat(i-1)*de_p_pp + Lp_pp_1
      xl_1p_pp = 2.0d0*pi*Msq_pp*omega2(0,0,2,1,Z,A,gp,gn,u_p_pp,Delta,h_max,Vwell)/hbar
      lamb_p_pp = lamb_p_pp + omega2(pn,hn,pp-1,hp,Z,A,gp,gn,Ex-u_p_pp,Delta,h_max,Vwell)*     &
                              omega2(0,0,1,0,Z,A,gp,gn,u_p_pp,Delta,h_max,Vwell)*              &
                              xl_1p_pp*de_p_pp
      u_h_pp = dfloat(i-1)*de_h_pp + Lh_pp_1
      xl_1h_pp = 2.0d0*pi*Msq_pp*omega2(0,0,1,2,Z,A,gp,gn,u_h_pp,Delta,h_max,Vwell)/hbar
      lamb_h_pp = lamb_h_pp + omega2(pn,hn,pp,hp-1,Z,A,gp,gn,Ex-u_h_pp,Delta,h_max,Vwell)*     &
                              omega2(0,0,0,1,Z,A,gp,gn,u_h_pp,Delta,h_max,Vwell)*              &
                              xl_1h_pp*de_h_pp
      u_p_np = dfloat(i-1)*de_p_np + Lp_np_1
      xl_1p_np = 2.0d0*pi*Msq_np*omega2(1,0,1,1,Z,A,gp,gn,u_p_np,Delta,h_max,Vwell)/hbar
      lamb_p_np = lamb_p_np + omega2(pn-1,hn,pp,hp,Z,A,gp,gn,Ex-u_p_np,Delta,h_max,Vwell)*     &
                              omega2(1,0,0,0,Z,A,gp,gn,u_p_np,Delta,h_max,Vwell)*              &
                              xl_1p_np*de_p_np
      u_h_np = dfloat(i-1)*de_h_np + Lh_np_1
      xl_1h_np = 2.0d0*pi*Msq_np*omega2(0,1,1,1,Z,A,gp,gn,u_h_np,Delta,h_max,Vwell)/hbar
      lamb_h_np = lamb_h_np + omega2(pn,hn-1,pp,hp,Z,A,gp,gn,Ex-u_h_np,Delta,h_max,Vwell)*     &
                              omega2(0,1,0,0,Z,A,gp,gn,u_h_np,Delta,h_max,Vwell)*              &
                              xl_1h_np*de_h_np
!------------------------------------------------Lambda^+_n
      u_p_nn = dfloat(i-1)*de_p_nn + Lp_nn_1
      xl_1p_nn = 2.0d0*pi*Msq_nn*omega2(2,1,0,0,Z,A,gp,gn,u_p_nn,Delta,h_max,Vwell)/hbar
      lamb_p_nn = lamb_p_nn + omega2(pn-1,hn,pp,hp,Z,A,gp,gn,Ex-u_p_nn,Delta,h_max,Vwell)*     &
                              omega2(1,0,0,0,Z,A,gp,gn,u_p_nn,Delta,h_max,Vwell)*              &
                              xl_1p_nn*de_p_nn
      u_h_nn = dfloat(i-1)*de_h_nn + Lh_nn_1
      xl_1h_nn = 2.0d0*pi*Msq_nn*omega2(1,2,0,0,Z,A,gp,gn,u_h_nn,Delta,h_max,Vwell)/hbar
      lamb_h_nn = lamb_h_nn + omega2(pn,hn-1,pp,hp,Z,A,gp,gn,Ex-u_h_nn,Delta,h_max,Vwell)*     &
                              omega2(0,1,0,0,Z,A,gp,gn,u_h_nn,Delta,h_max,Vwell)*              &
                              xl_1h_nn*de_h_nn
      u_p_pn = dfloat(i-1)*de_p_pn + Lp_pn_1
      xl_1p_pn = 2.0d0*pi*Msq_pn*omega2(1,1,1,0,Z,A,gp,gn,u_p_pn,Delta,h_max,Vwell)/hbar
      lamb_p_pn = lamb_p_pn + omega2(pn,hn,pp-1,hp,Z,A,gp,gn,Ex-u_p_pn,Delta,h_max,Vwell)*     &
                              omega2(0,0,1,0,Z,A,gp,gn,u_p_pn,Delta,h_max,Vwell)*              &
                              xl_1p_pn*de_p_pn
      u_h_pn = dfloat(i-1)*de_h_pn + Lh_pn_1
      xl_1h_pn = 2.0d0*pi*Msq_pn*omega2(1,1,0,1,Z,A,gp,gn,u_h_pn,Delta,h_max,Vwell)/hbar
      lamb_h_pn = lamb_h_pn + omega2(pn,hn,pp,hp-1,Z,A,gp,gn,Ex-u_h_pn,Delta,h_max,Vwell)*     &
                              omega2(0,0,0,1,Z,A,gp,gn,u_h_pn,Delta,h_max,Vwell)*              &
                              gp*xl_1h_pn*de_h_pn
!-------------------------------------------------Lambda^0_pn
      u_0_pn = dfloat(i-1)*de_0_pn + L0_pn_1
      xl_1p1h_pn = 2.0d0*pi*Msq_pn*omega2(1,1,0,0,Z,A,gp,gn,u_0_pn,Delta,h_max,Vwell)/hbar
      l_0_pn = l_0_pn + omega2(pn,hn,pp-1,hp-1,Z,A,gp,gn,Ex-u_0_pn,Delta,h_max,Vwell)*         &
                        omega2(0,0,1,1,Z,A,gp,gn,u_0_pn,Delta,h_max,Vwell)*                    &
                        xl_1p1h_pn*de_0_pn
!-------------------------------------------------Lambda^0_np
      u_0_np = dfloat(i-1)*de_0_np + L0_np_1
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
