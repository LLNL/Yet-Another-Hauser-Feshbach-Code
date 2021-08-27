!
!*******************************************************************************
!
subroutine Fission_data(data_path,len_path,icomp)                
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine to set up fission data for nuclei in the calculation
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        nodeinfo
!        constants
!
!     Subroutines:
!
!       init_barrier_data
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
   use nuclei
   use nodeinfo
   use constants
   implicit none
   character(len=200), intent(in) :: data_path
   integer(kind=4), intent(in) :: len_path
   integer(kind=4), intent(in) :: icomp
!-----------------------------------------------------------
   integer(kind=4) :: ia
   real(kind=8) :: Ea, hbwa, Eb, hbwb, Delf

   character(len=2) :: symb, symma, symmb
   integer(kind=4) :: izz, iaa
   integer(kind=4) :: i

!-------------------------------------------------------------------

   ia = nucleus(icomp)%A

   nucleus(icomp)%fission = .false.

!----    First set up fission parameters based on data in 'Fission-barrier.dat'
!----    Then check for data in 'Fission-Parameters.txt', which contains
!----    user specified parameters, likely due to previous fitting
!----
!   open(unit=53,file=data_path(1:len_path)//'Fission-barrier.dat',status='old')
   open(unit=53,file=data_path(1:len_path)//'Fission/empirical-barriers.dat',status='old')
   do i = 1, 4
      read(53,*)
   end do
! 1 read(53,'(2(1x,i3),1x,a2,2(f8.3,2x,a2),f8.3)',end=2)    &
!              izz,iaa,symb,Ea,symma,Eb,symmb,Delf
 1 read(53,'(2(1x,i3),1x,a2,1x,2(3x,a2,2(1x,f7.2)),1x,f8.3)',end=2)    &
              izz, iaa, symb, symma, Ea, hbwa, symmb, Eb, hbwb, Delf

   if(izz == nucleus(icomp)%Z .and. iaa == nucleus(icomp)%A)then

      nucleus(icomp)%fission = .true.
      nucleus(icomp)%fiss_tran_states = use_tran_states    !   New default for transition states
!      nucleus(icomp)%fiss_tran_states = .true.    !   New default for transition states
      nucleus(icomp)%fiss_user_levels = .false.
      nucleus(icomp)%F_n_barr = 2
      if(Eb <= 0.001)then
         nucleus(icomp)%fission = .true.
         nucleus(icomp)%F_n_barr = 1
         if(.not.allocated(nucleus(icomp)%F_Barrier))allocate(nucleus(icomp)%F_Barrier(1))
         nucleus(icomp)%F_Barrier(1)%state_scale = 1.0d0
         nucleus(icomp)%F_barrier(1)%barrier = Ea
         nucleus(icomp)%F_Barrier(1)%fit_ematch = .true.

         if(symmb == 'S ')then
             nucleus(icomp)%F_barrier(1)%symmetry = 1
         elseif(symmb == 'GA')then
             nucleus(icomp)%F_barrier(1)%symmetry = 3
         elseif(symmb == 'MA')then
             nucleus(icomp)%F_barrier(1)%symmetry = 2
         end if

!----   Use a default for the width, file does not contain hbw for single barriers
         if(iand(iaa,1) == 0)then
            if(iand(izz,1) == 0)then
               nucleus(icomp)%F_Barrier(1)%hbw = 0.9d0
            else
               nucleus(icomp)%F_Barrier(1)%hbw = 0.6d0
            end if
         else
            nucleus(icomp)%F_Barrier(1)%hbw = 0.8d0
         end if
      else
         nucleus(icomp)%fission = .true.
         nucleus(icomp)%F_n_barr = 2
         if(.not.allocated(nucleus(icomp)%F_Barrier))allocate(nucleus(icomp)%F_Barrier(2))
         nucleus(icomp)%F_Barrier(1)%state_scale = 1.0d0
         nucleus(icomp)%F_Barrier(2)%state_scale = 1.0d0
         nucleus(icomp)%F_Barrier(1)%barrier = Ea
         nucleus(icomp)%F_Barrier(2)%barrier = Eb
         nucleus(icomp)%F_Barrier(1)%hbw = hbwa
         nucleus(icomp)%F_Barrier(2)%hbw = hbwb
         nucleus(icomp)%F_Barrier(1)%level_param(3) = Delf
         nucleus(icomp)%F_Barrier(1)%level_param(6) = 2.5 + 150./real(iA,kind=8) + Delf
         nucleus(icomp)%F_Barrier(2)%level_param(3) = Delf
         nucleus(icomp)%F_Barrier(2)%level_param(6) = 2.5 + 150./real(iA,kind=8) + Delf
         nucleus(icomp)%F_Barrier(1)%fit_ematch = .true.
         nucleus(icomp)%F_Barrier(2)%fit_ematch = .true.

         if(symma == 'S ')then
             nucleus(icomp)%F_barrier(1)%symmetry = 1
         elseif(symma == 'GA')then
             nucleus(icomp)%F_barrier(1)%symmetry = 3
         elseif(symma == 'MA')then
             nucleus(icomp)%F_barrier(1)%symmetry = 2
         end if
         if(symmb == 'S ')then
             nucleus(icomp)%F_barrier(2)%symmetry = 1
         elseif(symmb == 'GA')then
             nucleus(icomp)%F_barrier(2)%symmetry = 3
         elseif(symmb == 'MA')then
             nucleus(icomp)%F_barrier(2)%symmetry = 2
         end if
      end if

      call init_barrier_data(icomp)

   end if
   goto 1
 2 continue
   close(unit=53)

   return
end subroutine Fission_data
!
!*******************************************************************************
!
subroutine Fission_levels(data_path,len_path,icomp)                
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine sets up transition states above the barrier (if requested)
!    and fits to matching energy given transition states above the fission barrier.
!    If transition states above the barrier are not used, it sets up conditions
!    to compute the spin cut off parameter below the user defined matching energy.
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        nodeinfo
!
!     Subroutines:
!
!       find_T_E0
!       cumm_rho
!
!     External functions:
!
!       real(kind=8) :: sig_2_param
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
   use nodeinfo
   implicit none
   character(len=200), intent(in) :: data_path
   integer(kind=4), intent(in) :: len_path
   integer(kind=4), intent(in) :: icomp
!-----------------------------------------------------------
   integer(kind=4) :: iZ, iA, iN
   real(kind=8) :: spin, sg2cut, sum
   real(kind=8) :: chisq, chisq_min
   real(kind=8) :: Ecut, Ematch, Ematch_best, Ematch_min
   real(kind=8) :: A
   integer(kind=4) :: nfit
!   real(kind=8) :: sig2
!   real(kind=8) :: sig2F
!   real(kind=8) :: sig2_em
   real(kind=8) :: sig2_min
!   real(kind=8) :: deriv

   real(kind=8), allocatable :: elv(:)
   real(kind=8), allocatable :: cum_rho(:),weight(:)
   real(kind=8), allocatable :: cum_fit(:)

   integer(kind=4) i, j, num
!-----------------------------------------------------------
!   real(kind=8) :: sig2_param
!-----------------------------------------------------------
   iZ = nucleus(icomp)%Z
   iA = nucleus(icomp)%A
   iN = ia - iz
   A = real(ia,kind=8)

   num = 0

!-------------------------------------------------------------------
   do i = 1, nucleus(icomp)%F_n_barr

!----   Find Temperature and E0 at start -  later fit discrete
      call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,                             &
                     nucleus(icomp)%F_Barrier(i)%vib_enh,                                    &
                     nucleus(icomp)%F_Barrier(i)%rot_enh)
       
      if(nucleus(icomp)%fiss_tran_states)then

         if(.not. nucleus(icomp)%fiss_user_levels)call make_transition_states(data_path, len_path, icomp, i)


         num = nucleus(icomp)%F_Barrier(i)%num_discrete

         do j = 1, num
            nucleus(icomp)%F_Barrier(i)%state_e(j) = nucleus(icomp)%F_Barrier(i)%state_e(j)*  &
                                                     nucleus(icomp)%F_Barrier(i)%state_scale
         end do
         nucleus(icomp)%F_Barrier(i)%ecut = nucleus(icomp)%F_Barrier(i)%ecut*                 &
                                            nucleus(icomp)%F_Barrier(i)%state_scale

!-----   Spin cut off parameter at Ecut
!-----   Taken from Reffo based on maximum liklihood arguments

         sg2cut = (0.83d0*A**0.26d0)**2

         sg2cut=0.0d0
         sum = 0.0d0
         do j = 1, num
            spin = nucleus(icomp)%F_Barrier(i)%state_j(j)
            sg2cut = sg2cut + spin*(spin + 1.0d0)*(2.0d0*spin+1.0d0)
            sum = sum + (2.0d0*spin+1.0d0)
         end do
         sg2cut = sg2cut/(2.0d0*sum)      !   standard model,const. sig**2, integrate over J
!         sg2cut = sg2cut/(3.0d0*sum)     !   TALYS, and other arguments
         if(sg2cut < 1.0d-4)sg2cut = (0.83d0*real(nucleus(icomp)%A,kind=8)**0.26d0)**2

         nucleus(icomp)%F_Barrier(i)%level_param(9) = nucleus(icomp)%level_param(9)

         nucleus(icomp)%F_Barrier(i)%ecut = nucleus(icomp)%F_Barrier(i)%state_e(num) + 0.0001

         nucleus(icomp)%F_Barrier(i)%level_param(8) = sg2cut
         nucleus(icomp)%F_Barrier(i)%ecut = nucleus(icomp)%F_Barrier(i)%state_e(num) + 0.001
         nucleus(icomp)%F_Barrier(i)%level_param(7) = nucleus(icomp)%F_Barrier(i)%ecut
         nfit = nucleus(icomp)%F_Barrier(i)%num_discrete

         if(nucleus(icomp)%F_Barrier(i)%level_param(19) == 0)then
             nucleus(icomp)%F_Barrier(i)%level_param(12) =                        &
                  nucleus(icomp)%F_Barrier(i)%level_param(2)*A**(2.d0/3.d0)
         else
             nucleus(icomp)%F_Barrier(i)%level_param(12) =                        &
                  nucleus(icomp)%F_Barrier(i)%level_param(2)*A**(5.d0/3.d0)
         end if

         nfit = nucleus(icomp)%F_Barrier(i)%num_discrete

         if(nucleus(icomp)%F_Barrier(i)%fit_ematch)then 
            if(nfit > 10)then
               if(.not.allocated(elv))allocate(elv(nfit))                   !  allocate cumulative density array
               if(.not.allocated(cum_rho))allocate(cum_rho(nfit))           !  allocate cumulative density array
               if(.not.allocated(weight))allocate(weight(nfit))
               if(.not.allocated(cum_fit))allocate(cum_fit(nfit))           !  allocate cumulative density array
               do j = 1, nfit                                               !  Make cumulative rho array
                  elv(j) = nucleus(icomp)%F_Barrier(i)%state_e(j)
                  cum_rho(j)=real(j,kind=8)
                  weight(j)=1.0d0/cum_rho(j)                                 !  weight higher values more
!                  weight(j)=1.0d0/sqrt(cum_rho(j))                          !  assume error 1/sqrt
               end do

               Ematch_min = nucleus(icomp)%F_Barrier(i)%level_param(3) + 0.4d0
               Ematch = 8.0
               nucleus(icomp)%F_Barrier(i)%level_param(6) = ematch
                  nucleus(icomp)%F_Barrier(i)%level_param(6) = ematch
                  call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,            &
                                 nucleus(icomp)%F_Barrier(i)%vib_enh,                   &
                                 nucleus(icomp)%F_Barrier(i)%rot_enh)
                  call cumm_rho(nfit,elv,iA,nucleus(icomp)%F_Barrier(i)%level_param,    &
                                nucleus(icomp)%F_Barrier(i)%vib_enh,                    &
                                nucleus(icomp)%F_Barrier(i)%rot_enh,cum_fit)
               chisq = 0.0d0
               do j = 1, nfit
                  chisq = chisq + (cum_rho(j)-cum_fit(j))**2/weight(j)**2
               end do
               chisq = chisq/real(nfit-2,kind=8)
               chisq_min = chisq
               Ematch_best = Ematch
               do while(Ematch > Ematch_min)
                  Ematch = Ematch - 0.01d0
                  nucleus(icomp)%F_Barrier(i)%level_param(6) = ematch
                  call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,            &
                                 nucleus(icomp)%F_Barrier(i)%vib_enh,                   &
                                 nucleus(icomp)%F_Barrier(i)%rot_enh)
                  call cumm_rho(nfit,elv,iA,nucleus(icomp)%F_Barrier(i)%level_param,    &
                                nucleus(icomp)%F_Barrier(i)%vib_enh,                    &
                                nucleus(icomp)%F_Barrier(i)%rot_enh,cum_fit)
                  chisq = 0.0d0
                  do j = 1, nfit
                     chisq = chisq + (cum_rho(j)-cum_fit(j))**2/weight(j)**2
                  end do
                  chisq = chisq/real(nfit-2,kind=8)
                  if(chisq < chisq_min)then
                     chisq_min = chisq
                     Ematch_best = Ematch
                  end if
               end do
               Ematch = Ematch_best
               nucleus(icomp)%F_Barrier(i)%level_param(6) = ematch
               call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,           &
                              nucleus(icomp)%F_Barrier(i)%vib_enh,                  &
                              nucleus(icomp)%F_Barrier(i)%rot_enh)
                  call cumm_rho(nfit,elv,iA,nucleus(icomp)%F_Barrier(i)%level_param,    &
                                nucleus(icomp)%F_Barrier(i)%vib_enh,                    &
                                nucleus(icomp)%F_Barrier(i)%rot_enh,cum_fit)

               if(allocated(elv))deallocate(elv)
               if(allocated(cum_rho))deallocate(cum_rho)
               if(allocated(weight))deallocate(weight)
               if(allocated(cum_fit))deallocate(cum_fit)
            else
!--------------------------------------------------------------------------------
!------   Fail safe, shouldn't get here, but need to do something to proceed
!------   too few levels to fit ematch, and ematch has not been specified
!------   punt and set to ematch of the compound nucleus
!--------------------------------------------------------------------------------
               if(iproc == 0)then
                  write(6,*)'Too few discrete transition states to fit ematch'
                  write(6,*)'Either create more transition states, or do not use'
               end if
               call exit_YAHFC(601)
            end if
         else
            if(iproc == 0 .and. verbose)write(6,*)'Ematch set by user'
         end if
      else
!--------------------------------------------------------------------------------
!------   No discrete states above the barrier, just continuous level density
!------   need to set spin cut off parameter so that it is not too different from
!------   the compound nucleus (CN). Previously, it was just the analytic formula
!------   But, this led to a big difference in the spin distributions at low energy 
!------   between the compound nucleus discrete spectrum and the barriers
!------   Plan is to have a spin cut off at an effective ecut, as if there are 
!------   discrete states. Scale the CN value by ratio of spin cutoff at Ematch,
!------   should be specified by user, with that at the barrier. Scale the CN spin 
!------   cut off by this amount.
!--------------------------------------------------------------------------------

!         nucleus(i)%F_Barrier(i)%fit_ematch = .false.
!         Ecut = nucleus(icomp)%level_param(7)
!         Ematch = nucleus(icomp)%level_param(6)

!         sg2cut = nucleus(icomp)%level_param(8)
!         sig2_em = sig2_param(Ematch+1.0d-5,nucleus(icomp)%level_param,A)
!         deriv = (sig2_em - sg2cut)/(ematch - ecut)
!         sig2_min = min((0.83*A**0.26)**2,sg2cut)

!         Ematch = nucleus(icomp)%F_Barrier(i)%level_param(6)
!         sig2 = sig2_param(Ematch+1.0d-5,nucleus(icomp)%level_param,A)
!         sig2F = sig2_param(Ematch+1.0d-5,nucleus(icomp)%F_Barrier(i)%level_param,A)


!         Ecut = max(Ematch - (sig2F - sig2_min)/deriv,0.0d0)

!         deriv = (sig2F - sig2_min)/(Ematch - Ecut)

!         nucleus(icomp)%F_Barrier(i)%Ecut = ecut
!         nucleus(icomp)%F_Barrier(i)%level_param(7) = nucleus(icomp)%F_Barrier(i)%ecut
!         nucleus(icomp)%F_Barrier(i)%level_param(20) = deriv

!         nucleus(icomp)%F_Barrier(i)%level_param(8) = sig2_min

         Ematch = nucleus(icomp)%level_param(6)
         sig2_min = (0.83*real(nucleus(icomp)%A,kind=8)**0.26)**2
!         deriv = (sig2F - sig2_min)/(Ematch - Ecut)
         nucleus(icomp)%F_Barrier(i)%level_param(9) = nucleus(icomp)%level_param(9)
         Ecut = 0.0001
         nucleus(icomp)%F_Barrier(i)%ecut = Ecut
         nucleus(icomp)%F_Barrier(i)%level_param(8) = sig2_min
         nucleus(icomp)%F_Barrier(i)%level_param(7) = nucleus(icomp)%F_Barrier(i)%ecut
!         deriv = (sig2F - sig2_min)/(Ematch - Ecut)

         call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,            &
                        nucleus(icomp)%F_Barrier(i)%vib_enh,                   &
                        nucleus(icomp)%F_Barrier(i)%rot_enh)
      end if
   end do

   return

end subroutine Fission_levels
!
!*******************************************************************************
!
subroutine Fission_transmission(icomp,Ex,xji,ipar,F_trans)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine to compute the fission transmission coefficients at
!    excitation energy ex, angular moment xji, and parity ipar 
!    for nucleus icomp
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        nodeinfo
!        constants
!
!     Subroutines:
!
!        rho_J_par_e
!
!     External functions:
!
!      real(kind=8) :: HW_trans
!      logical :: real8_equal
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
   use nodeinfo
   use constants
   implicit none
   integer(kind=4), intent(in) :: icomp
   real(kind=8), intent(in) :: Ex, xji
   integer(kind=4), intent(in) :: ipar
   real(kind=8), intent(out) :: F_trans(4)
!-----------------------------------------------------
   real(kind=8) :: trans(3)
   real(kind=8) :: de
   real(kind=8) :: T_f
   integer(kind=4) :: i, j
   real(kind=8) :: Ef
   real(kind=8) :: apar
   real(kind=8) :: spin_cut
   real(kind=8) :: delta
   real(kind=8) :: ecut
   real(kind=8) :: sg2cut
   real(kind=8) :: sig_mod
   real(kind=8) :: ematch
   real(kind=8) :: emin
   real(kind=8) :: par
   real(kind=8) :: rho, sig2, apu
   real(kind=8) :: tt
   real(kind=8) :: K_vib, K_rot
!   real(kind=8) :: rho_FM
   integer(kind=4) :: ia
   real(kind=8) :: e_old
   logical converged
   real(kind=8) :: pmode, pe1, pbb
   integer(kind=4) :: iipar
   integer(kind=4) :: ji
   integer(kind=4) :: jj
   save e_old

   real(kind=8) :: E0, T, E1
!   real(kind=8) :: pfac, jfac
   real(kind=8) :: b, Max_J, F_Barrier, F_Barrier_hbw
   real(kind=8) :: aa, bb, cc
!--------   External Functions  ----------------------
   real(kind=8) :: HW_trans
!   real(kind=8) :: spin_fac
!   real(kind=8) :: parity_fac
!   logical :: real8_equal
!-----------------------------------------------------

   par = 1.0d0
   if(ipar == 0)par = -1.0d0

   de = 0.01
   ia = nucleus(icomp)%A
   do i = 1, nucleus(icomp)%F_n_barr
      apar = nucleus(icomp)%F_barrier(i)%level_param(1)
      spin_cut = nucleus(icomp)%F_barrier(i)%level_param(2)
      delta = nucleus(icomp)%F_barrier(i)%level_param(3)
      ecut = nucleus(icomp)%F_Barrier(i)%ecut
      sg2cut = nucleus(icomp)%F_barrier(i)%level_param(8)
      sig_mod = nucleus(icomp)%F_barrier(i)%level_param(9)
      ematch = nucleus(icomp)%F_barrier(i)%level_param(6)
      E0 = nucleus(icomp)%F_barrier(i)%level_param(15)
      T = nucleus(icomp)%F_barrier(i)%level_param(14)
      E1 = 0.0d0

      aa = nucleus(icomp)%F_Barrier(i)%barrier_damp(1)
      bb = nucleus(icomp)%F_Barrier(i)%barrier_damp(2)
      cc = nucleus(icomp)%F_Barrier(i)%barrier_damp(3)

      trans(i)=0.0d0


      Max_j = nucleus(icomp)%F_barrier(i)%Max_J
      F_Barrier = nucleus(icomp)%F_Barrier(i)%barrier
      F_Barrier = F_Barrier*aa*exp(-((Ex-bb)/cc)**2)
      if(Max_J > 0.0d0)then
         b = 1.0d0/(Max_J*(Max_J+1.0d0))
         if(xji <= Max_J)then
            F_Barrier = F_Barrier*(1.0d0 - b*xji*(xji+1.0d0))
         else
            F_Barrier = 0.1d0
         end if
      end if
      F_Barrier = max(F_Barrier,0.1d0)
      F_Barrier_hbw = nucleus(icomp)%F_Barrier(i)%hbw

      ji = nint(2.0d0*xji)

      do j = 1, nucleus(icomp)%F_Barrier(i)%num_discrete
         jj = nint(2.0d0*nucleus(icomp)%F_barrier(i)%state_j(j))
         iipar = (nint(nucleus(icomp)%F_barrier(i)%state_pi(j))+1)/2
         if(jj == ji .and. iipar == ipar)then
            Ef = nucleus(icomp)%F_barrier(i)%state_e(j)
            tt = HW_trans(Ex,Ef,F_Barrier,F_barrier_hbw)
            trans(i) = trans(i) + tt
         end if
      end do

      pmode = nucleus(icomp)%F_barrier(i)%level_param(16)
      pe1 = nucleus(icomp)%F_barrier(i)%level_param(17)
      pbb = nucleus(icomp)%F_barrier(i)%level_param(18)

      emin = 0.0d0
      if(nucleus(icomp)%F_Barrier(i)%num_discrete > 0)emin = nucleus(icomp)%F_Barrier(i)%ecut

      T_f = 0.0d0
      converged = .false.
      Ef = emin - 0.5d0*de
      do while (.not. converged)
         Ef = Ef + de

         call rho_J_par_e(Ef, xji, ipar,                              &
                          nucleus(icomp)%F_barrier(i)%level_param,    &
                          nucleus(icomp)%F_barrier(i)%vib_enh,        &
                          nucleus(icomp)%F_barrier(i)%rot_enh,        &
                          ia, rho, apu, sig2, K_vib, K_rot)
!         call rhoe(Ef,nucleus(icomp)%F_barrier(i)%level_param,                 &
!                   nucleus(icomp)%F_barrier(i)%vib_enh,                        &
!                   nucleus(icomp)%F_barrier(i)%rot_enh,                        &
!                   ia,rho_Fm,apu,sig2,K_vib,K_rot)
!
!         jfac = spin_fac(xji,sig2)
!         pfac = parity_fac(Ef,xji,ipar,pmode,pe1,pbb)
!         rho = rho_FM*jfac*pfac
         tt = HW_trans(Ex,Ef,F_Barrier,F_barrier_hbw)
         tt = tt*rho*de
         if(T_f > 0.0d0 .and. tt/T_f < 1.0d-5)converged = .true.
         if(Ef > (Ex - F_Barrier + F_barrier_hbw) .and. tt < 1.0d-6)converged = .true.
         T_f = T_f + tt
!         if(converged)write(6,*)tt,T_f,tt/T_f
!  if(Ef > Ex)write(82,*)Ex,Ef,rho,F_barrier,F_barrier_hbw,tt,T_f,tt/T_f
!  if(converged)write(82,*)'finsihed'
      end do


      trans(i) = trans(i) + T_f
      if(trans(i) < 1.0d-20)trans(i)=1.0d-20
      F_trans(i) = trans(i)
   end do

   if(nucleus(icomp)%F_n_barr == 1)then
      F_trans(4) = F_trans(1)
   else if(nucleus(icomp)%F_n_barr == 2)then
      F_trans(4) = F_trans(1)*F_trans(2)/(F_trans(1)+F_trans(2))
   else if(nucleus(icomp)%F_n_barr == 3)then
      F_trans(4) = F_trans(1)*(F_trans(2) + F_trans(3))/(F_trans(1) + F_trans(2) + F_trans(3))
   end if

   if(F_trans(4) <= 1.0d-20)F_trans(4) = 0.0d0
   e_old = ex
    return
end subroutine Fission_transmission
!
!*******************************************************************************
!
real(kind=8) function HW_trans(Ex,ei,Barrier,hbw)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the Hill-Wheeler transmission coefficient for 
!    fission barrier
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
   use variable_kinds
   use constants
   implicit none
   real(kind=8), intent(in) :: Ex, ei, Barrier, hbw
!-------------------------------------------------------------------
   real(kind=8) :: hbwx
   real(kind=8) :: exponent
!-------------------------------------------------------------------
   hbwx = max(hbw,1.0d-2)
   exponent = -2.0d0*pi*(Ex-Barrier-ei)/hbwx
   if(exponent <= -100.0d0)exponent = -100.0d0
   if(exponent >= 100.0d0)exponent = 100.0d0
   HW_trans = 1.0d0/(1.0d0 + exp(exponent))
end function HW_trans
!
!*******************************************************************************
!
subroutine init_barrier_data(i)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine initializes fission parameter info for nucleus i
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nodeinfo
!        options
!        useful_data
!        nuclei
!        particles_def
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
   use nodeinfo
   use options
   use useful_data
   use nuclei
   use particles_def
   use constants
   implicit none
   integer(kind=4), intent(in) :: i
!-----------------------------------------------------------------------------
   integer(kind=4) :: j
   real(kind=8) :: beta_2
   real(kind=8) :: sig2_perp
   real(kind=8) :: sig2_ax
!-----------------------------------------------------------------------------

   do j = 1, nucleus(i)%F_n_barr
      if(Fiss_Max_J > 0.0d0)then
         nucleus(i)%F_Barrier(j)%Max_J = Fiss_Max_J
         if(iand(nucleus(i)%A,1) == 1)nucleus(i)%F_Barrier(j)%Max_J =        &
                                      nucleus(i)%F_Barrier(j)%Max_J + 0.5d0
      else
         nucleus(i)%F_Barrier(j)%Max_J = 0.0d0
      end if 
!---- Don't reset barrier heights and widths if already set
      if(nucleus(i)%F_Barrier(j)%barrier < 1.0d-1)nucleus(i)%F_Barrier(j)%barrier = 6.0d0
      if(nucleus(i)%F_Barrier(j)%hbw < 1.0d-1)nucleus(i)%F_Barrier(j)%hbw = 0.6d0
      nucleus(i)%F_Barrier(j)%ecut = 0.0d0
!----   If barrier symmetry has not been defined, then make guess for barrier symmetries.
!----   If barrier symmetry is not defined for this call, it is because user has 
!----   reset nucleus(i)%F_n_barr and the data has been reallocated. It is expected that 
!----   the user will reset the barrier symmetry as well. This is initialize based on
!----   the behavior of many actinides.
      if(nucleus(i)%F_Barrier(j)%symmetry == 0)then
         if(nucleus(i)%F_n_barr == 1)then
            nucleus(i)%F_Barrier(j)%symmetry = 1               !  Axially symmetric
         else
            if(j == 1)then
               nucleus(i)%F_Barrier(j)%symmetry = 1                            !  Axially symmetric
               if(nucleus(i)%A >= 144)nucleus(i)%F_Barrier(j)%symmetry = 3     !  Axially asymmetric
            end if
            if(j == 2)nucleus(i)%F_Barrier(j)%symmetry = 2     !  Mass asymmetric
         end if
      end if

      beta_2 = 2.0d0*nucleus(i)%beta(2)
      if(j == 2)then
         beta_2 = 4.0d0*nucleus(i)%beta(2)
      end if

      nucleus(i)%F_Barrier(j)%level_param(19) = nucleus(i)%level_param(19)

      nucleus(i)%F_Barrier(j)%beta_2 = beta_2

      sig2_perp = 0.0d0
      sig2_ax = 0.0d0
      if(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 2 .or.            &
         nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 3)then
         sig2_perp = (1.0d0 + beta_2/3.0d0)
         sig2_ax = (1.0d0 - 2.0d0*beta_2/3.0d0)
      elseif(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 4 .or.            &
         nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 5)then
         sig2_perp = (1.0d0 + sqrt(5.0d0/(4.0d0*pi))*beta_2)
         sig2_ax = (1.0d0 - 0.5d0*sqrt(5.0d0/(4.0d0*pi))*beta_2)
      end if

      nucleus(i)%F_Barrier(j)%level_param(1) = nucleus(i)%level_param(1)
      nucleus(i)%F_Barrier(j)%level_param(2) = nucleus(i)%level_param(2)
      nucleus(i)%F_Barrier(j)%level_param(3) = nucleus(i)%level_param(3)
      if(j == 1)nucleus(i)%F_Barrier(j)%level_param(4) = 1.5d0
      if(j == 2)nucleus(i)%F_Barrier(j)%level_param(4) = 0.6d0
      nucleus(i)%F_Barrier(j)%level_param(5) = nucleus(i)%level_param(5)
      nucleus(i)%F_Barrier(j)%level_param(6) = nucleus(i)%level_param(6)
      nucleus(i)%F_Barrier(j)%level_param(7) = nucleus(i)%F_Barrier(j)%ecut
      nucleus(i)%F_Barrier(j)%level_param(8) = (0.83*real(nucleus(i)%A,kind=8)**0.26)**2

      nucleus(i)%F_Barrier(j)%level_param(9) = 1.0d0
      nucleus(i)%F_Barrier(j)%level_param(10) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(11) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(12) = nucleus(i)%level_param(12)
      nucleus(i)%F_Barrier(j)%level_param(13) = nucleus(i)%level_param(13)
      nucleus(i)%F_Barrier(j)%level_param(14) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(15) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(16) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(17) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(18) = 0.0d0

      nucleus(i)%F_barrier(j)%barrier_damp(1) = 1.0d0
      nucleus(i)%F_barrier(j)%barrier_damp(2) = 0.0d0
      nucleus(i)%F_barrier(j)%barrier_damp(3) = 1.0d6

      if(nucleus(i)%lev_option >= 2)then
         nucleus(i)%F_Barrier(j)%level_param(10) = real(nucleus(i)%F_Barrier(j)%symmetry,kind=8)

         nucleus(i)%F_Barrier(j)%level_param(11) = nucleus(i)%level_param(11)
            
         nucleus(i)%F_Barrier(j)%rot_enh(1) = 1.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(2) = 30.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(3) = 5.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(4) = sig2_perp
         nucleus(i)%F_Barrier(j)%rot_enh(5) = sig2_ax
         nucleus(i)%F_Barrier(j)%vib_enh(1) = 1.0d0
         nucleus(i)%F_Barrier(j)%vib_enh(2) = 30.0d0
         nucleus(i)%F_Barrier(j)%vib_enh(3) = 5.0d0
      else
         nucleus(i)%F_Barrier(j)%rot_enh(1) = 0.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(2) = 30.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(3) = 5.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(4) = sig2_perp
         nucleus(i)%F_Barrier(j)%rot_enh(5) = sig2_ax
         nucleus(i)%F_Barrier(j)%vib_enh(1) = 0.0d0
         nucleus(i)%F_Barrier(j)%vib_enh(2) = 30.0d0
         nucleus(i)%F_Barrier(j)%vib_enh(3) = 5.0d0
      end if
      nucleus(i)%F_Barrier(j)%num_discrete = 0
      nucleus(i)%F_Barrier(j)%ncut = 0
      nucleus(i)%F_Barrier(j)%state_scale = 1.0d0
      nucleus(i)%F_Barrier(j)%fit_ematch = .true.
   end do

   return
end subroutine init_barrier_data
!
!*******************************************************************************
!
subroutine make_transition_states(data_path, len_path, inuc, ibarr)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine generates the transition states above barrier ibarr
!    based on the bandheads store in the $YAHFC_DATA/Fission directory
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nodeinfo
!        options
!        useful_data
!        nuclei
!        particles_def
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
   use nuclei
   use nodeinfo
   use constants
   implicit none
   character(len=200), intent(in) :: data_path
   integer(kind=4), intent(in) :: len_path
   integer(kind=4), intent(in) :: inuc
   integer(kind=4), intent(in) :: ibarr
!-----------------------------------------------------------
   integer(kind=4) :: iZ, iA, iN
   integer(kind=4) :: izz, iaa, inn
   integer(kind=4) :: i, j, n
   integer(kind=4) :: isymm
   real(kind=8) :: A
   real(kind=8) :: Mom_Inertia
   real(kind=8) :: Mass
   real(kind=8) :: beta2
   real(kind=8) :: r0, R
   real(kind=8) :: xJ, xK, xJ_step, xJ_min
   real(kind=8) :: energy
   character(len=132) :: file_name
   integer(kind=4) :: ilast
   integer(kind=4) :: num_bands
   real(kind=8)  :: e_max
   real(kind=8) :: band_head(20)
   real(kind=8) :: band_K(20)
   real(kind=8) :: band_par(20)
   integer(kind=4) :: num_discrete
   logical :: finished
   character(len=132) :: buffer
   real(kind=8) :: temp

!----------------------------------------------------------------   
   r0 = 1.2d0
   iZ = nucleus(inuc)%Z
   iA = nucleus(inuc)%A
   iN = iA - iZ

   izz = iand(iZ,1)   
   iaa = iand(iA,1)
   inn = iand(iN,1)
   isymm = nucleus(inuc)%F_barrier(ibarr)%symmetry
!-----------------------------------------------------------------------------
!-----    Get file name for band head information
!-----------------------------------------------------------------------------
   if(izz == 1 .and. inn == 0)then
      if(ibarr == 1)file_name = 'Fission-bands-odd-even-1.dat'
      if(ibarr == 2)file_name = 'Fission-bands-odd-even-2.dat'
   elseif(izz == 0 .and. inn == 0)then
      if(ibarr == 1 .and. isymm == 1)file_name = 'Fission-bands-even-even-1-ax.dat'
      if(ibarr == 1 .and. isymm == 3)file_name = 'Fission-bands-even-even-1-nonax.dat'
      if(ibarr == 2)file_name = 'Fission-bands-even-even-2.dat'
   elseif(izz == 0 .and. inn == 1)then
      if(ibarr == 1)file_name = 'Fission-bands-even-odd-1.dat'
      if(ibarr == 2)file_name = 'Fission-bands-even-odd-2.dat'
   elseif(izz == 1 .and. inn == 1)then
      if(ibarr == 1)file_name = 'Fission-bands-odd-odd-1.dat'
      if(ibarr == 2)file_name = 'Fission-bands-odd-odd-2.dat'
   end if

   ilast = index(file_name,' ')-1

!-----------------------------------------------------------------------------
!-----    Read band head data
!-----------------------------------------------------------------------------
   open(unit=19, file = data_path(1:len_path)//'/Fission/'//file_name(1:ilast), status = 'old')
   finished = .false.
   num_bands = 0
   n = 0
   do while(.not. finished)
      read(19,'(a)')buffer
      if(buffer(1:1) == '#' .or. buffer(1:1) == '!')cycle
      if(num_bands == 0)then
         read(buffer,*)num_bands, e_max
      else
         n = n + 1
         read(buffer,*)band_head(n), band_K(n), band_par(n)
         if(n == num_bands)finished = .true.
      end if
   end do
   close(unit=19)
   temp = nucleus(inuc)%F_barrier(ibarr)%ecut
   if(temp > 0.0d0)e_max = max(temp,e_max)
   nucleus(inuc)%F_barrier(ibarr)%ecut = e_max
!-----------------------------------------------------------------------------
!-----    calculate moment of inertia
!-----------------------------------------------------------------------------
   A = real(iA,kind=8) 
   Mass = nucleus(inuc)%mass
   R = r0*A**(1.0d0/3.0d0)
   beta2 = nucleus(inuc)%F_barrier(ibarr)%beta_2
   Mom_inertia = (1.0d0 + beta2/3.0d0)*0.4d0*Mass*R**2
!   Mom_inertia = (1.0d0 + sqrt(5.0d0/(4.0d0*pi))*beta2)*0.4d0*Mass*R**2
!-----------------------------------------------------------------------------
!-----    First run to compute energies and count how many levels below ecut
!-----------------------------------------------------------------------------
   num_discrete = 0
   do n = 1, num_bands
      xK = band_k(n)
      xJ_min = xK
      xJ_step = 1.0d0
      if(nint(xK) == 0)then
         if(nint(band_par(n)) == 1)xJ_step = 2.0d0
         if(nint(band_par(n)) == -1)xJ_min = 1.0d0
      end if
      xJ = xJ_min - xJ_step
      energy = 0.0d0
      do while(energy < e_max)
         xJ = xJ + xJ_step
         energy = (xJ*(xJ+1.0d0)-xJ_min*(xJ_min+1.0d0))/(2.0d0*Mom_Inertia)*hbar_c**2 + band_head(n)
         if(energy > e_max) exit
         num_discrete = num_discrete + 1
      end do
   end do
!-----------------------------------------------------------------------------
!----   Finished creating all elvels below ecut
!----   allocate arrays and put in the arrays
!-----------------------------------------------------------------------------
   nucleus(inuc)%F_Barrier(ibarr)%num_discrete = num_discrete
   nucleus(inuc)%F_Barrier(ibarr)%ncut = num_discrete
   nucleus(inuc)%F_Barrier(ibarr)%ecut = e_max
   allocate(nucleus(inuc)%F_Barrier(ibarr)%state_e(num_discrete))
   allocate(nucleus(inuc)%F_Barrier(ibarr)%state_J(num_discrete))
   allocate(nucleus(inuc)%F_Barrier(ibarr)%state_pi(num_discrete))
   i = 0
   do n = 1, num_bands
      xK = band_k(n)
      xJ_min = xK
      xJ_step = 1.0d0
      if(nint(xK) == 0)then
         if(nint(band_par(n)) == 1)xJ_step = 2.0d0
         if(nint(band_par(n)) == -1)xJ_min = 1.0d0
      end if
      xJ = xJ_min - xJ_step
      energy = 0.0d0
      do while(energy < e_max)
         xJ = xJ + xJ_step
         energy = (xJ*(xJ+1.0d0)-xJ_min*(xJ_min+1.0d0))/(2.0d0*Mom_Inertia)*hbar_c**2 + band_head(n)
         if(energy > e_max) exit
         i = i + 1
         nucleus(inuc)%F_Barrier(ibarr)%state_e(i) = energy
         nucleus(inuc)%F_Barrier(ibarr)%state_J(i) = xJ
         nucleus(inuc)%F_Barrier(ibarr)%state_pi(i) = band_par(n)
      end do
   end do
!-----------------------------------------------------------------------------
!-----    Put into ascending order
!-----------------------------------------------------------------------------
   do i = 1, num_discrete - 1
      do j = i, num_discrete
         if(nucleus(inuc)%F_Barrier(ibarr)%state_e(j) < nucleus(inuc)%F_Barrier(ibarr)%state_e(i))then
            temp = nucleus(inuc)%F_Barrier(ibarr)%state_e(i)
            nucleus(inuc)%F_Barrier(ibarr)%state_e(i) = nucleus(inuc)%F_Barrier(ibarr)%state_e(j)
            nucleus(inuc)%F_Barrier(ibarr)%state_e(j) = temp
            temp = nucleus(inuc)%F_Barrier(ibarr)%state_J(i)
            nucleus(inuc)%F_Barrier(ibarr)%state_J(i) = nucleus(inuc)%F_Barrier(ibarr)%state_J(j)
            nucleus(inuc)%F_Barrier(ibarr)%state_J(j) = temp
            temp = nucleus(inuc)%F_Barrier(ibarr)%state_pi(i)
            nucleus(inuc)%F_Barrier(ibarr)%state_pi(i) = nucleus(inuc)%F_Barrier(ibarr)%state_pi(j)
            nucleus(inuc)%F_Barrier(ibarr)%state_pi(j) = temp
         end if
      end do
   end do

   return

end subroutine make_transition_states
