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
   open(unit=53,file=data_path(1:len_path)//'empirical-barriers.dat',status='old')
   do i = 1, 4
      read(53,*)
   end do
! 1 read(53,'(2(1x,i3),1x,a2,2(f8.3,2x,a2),f8.3)',end=2)    &
!              izz,iaa,symb,Ea,symma,Eb,symmb,Delf
 1 read(53,'(2(1x,i3),1x,a2,1x,2(3x,a2,2(1x,f7.2)),1x,f8.3)',end=2)    &
              izz, iaa, symb, symma, Ea, hbwa, symmb, Eb, hbwb, Delf

   if(izz == nucleus(icomp)%Z .and. iaa == nucleus(icomp)%A)then

      nucleus(icomp)%fission = .true.
      nucleus(icomp)%F_n_barr = 2
      if(Ea <= 0.001)then
         nucleus(icomp)%fission = .true.
         nucleus(icomp)%F_n_barr = 1
         if(.not.allocated(nucleus(icomp)%F_Barrier))allocate(nucleus(icomp)%F_Barrier(1))
         nucleus(icomp)%F_barrier(1)%barrier = Eb
         nucleus(icomp)%F_Barrier(1)%hbw = hbwa

         if(symmb == 'S ')then
             nucleus(icomp)%F_barrier(1)%symmetry = 1
         elseif(symmb == 'GA')then
             nucleus(icomp)%F_barrier(1)%symmetry = 2
         elseif(symmb == 'MA')then
             nucleus(icomp)%F_barrier(1)%symmetry = 3
         end if

!----   Use a default for the width, file does not contain hbw for single barriers
         if(iand(iaa,1) == 0)then
            if(iand(izz,1) == 0)then
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
         nucleus(icomp)%F_Barrier(1)%barrier = Ea
         nucleus(icomp)%F_Barrier(2)%barrier = Eb
         nucleus(icomp)%F_Barrier(1)%hbw = hbwa
         nucleus(icomp)%F_Barrier(2)%hbw = hbwb
         nucleus(icomp)%F_Barrier(1)%level_param(3) = Delf
         nucleus(icomp)%F_Barrier(1)%level_param(6) = 2.5 + 150./real(iA,kind=8) + Delf
         nucleus(icomp)%F_Barrier(2)%level_param(3) = Delf
         nucleus(icomp)%F_Barrier(2)%level_param(6) = 2.5 + 150./real(iA,kind=8) + Delf

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
subroutine Fission_levels(icomp)                
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine to fit matching energy given transition states above the 
!    fission barrier
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
   integer(kind=4), intent(in) :: icomp
   integer(kind=4) :: ia
   real(kind=8) :: spin, sg2cut, sum
   real(kind=8) :: chisq, chi_min
   real(kind=8) :: ematch
   real(kind=8) :: A
   integer(kind=4) :: nfit
   real(kind=8) :: ratio, sig2, sig2F
   real(kind=8), allocatable :: elv(:)
   real(kind=8), allocatable :: cum_rho(:),dcum_rho(:)
   real(kind=8), allocatable :: cum_fit(:)

   integer(kind=4) i,j, num
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   real(kind=8) :: sig2_param
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ia = nucleus(icomp)%A
   A = real(ia,kind=8)


!-------------------------------------------------------------------
   do i = 1, nucleus(icomp)%F_n_barr

!----   Find Temperature and E0 at start -  later fit discrete
      call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,           &
                     nucleus(icomp)%F_Barrier(i)%vib_enh,                  &
                     nucleus(icomp)%F_Barrier(i)%rot_enh)


!-----   Spin cut off parameter at Ecut
!-----   Taken from Reffo based on maximum liklihood arguments

      sg2cut = (0.83d0*A**0.26d0)**2
      num = nucleus(icomp)%F_Barrier(i)%num_discrete

      if(num > 0)then
         sg2cut=0.0d0
         sum = 0.0d0
         do j = 1, num
            spin = nucleus(icomp)%F_Barrier(i)%state_j(j)
            sg2cut = sg2cut + spin*(spin + 1.0d0)*(2.0d0*spin+1.0d0)
            sum = sum + (2.0d0*spin+1.0d0)
         end do
         sg2cut = sg2cut/(3.0d0*sum)

         nucleus(icomp)%F_Barrier(i)%level_param(9) = nucleus(icomp)%level_param(9)

         nucleus(icomp)%F_Barrier(i)%ecut = nucleus(icomp)%F_Barrier(i)%state_e(j) + 0.0001

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

         if(nfit > 5)then                                             !  There are so few discrete states, 
                                                                      !  no point in fitting to the cumulative 
                                                                      !  level density
            if(.not.allocated(elv))allocate(elv(nfit))                   !  allocate cumulative density array
            if(.not.allocated(cum_rho))allocate(cum_rho(nfit))           !  allocate cumulative density array
            if(.not.allocated(dcum_rho))allocate(dcum_rho(nfit))
            if(.not.allocated(cum_fit))allocate(cum_fit(nfit))           !  allocate cumulative density array
            do j = 1, nfit                                               !  Make cumulative rho array
               elv(j) = nucleus(icomp)%F_Barrier(i)%state_e(j)
               cum_rho(j)=real(j,kind=8)
               dcum_rho(j)=1.0d0/sqrt(cum_rho(j))                        !  assume error 1/sqrt
            end do


            ematch = 8.0
            nucleus(icomp)%F_Barrier(i)%level_param(6) = ematch
            call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,            &
                           nucleus(icomp)%F_Barrier(i)%vib_enh,                   &
                           nucleus(icomp)%F_Barrier(i)%rot_enh)

            call cumm_rho(nfit,elv,nucleus(icomp)%A,nucleus(icomp)%level_param,   &
                          nucleus(icomp)%vib_enh,                                 &
                          nucleus(icomp)%rot_enh,cum_fit)
            chisq = 0.0d0

            do j = 1, nfit
               chisq = chisq + (cum_rho(j)-cum_fit(j))**2/dcum_rho(j)**2
            end do
            ematch = nucleus(icomp)%F_Barrier(i)%level_param(6)
            chi_min = 1.0d10

 10         nucleus(icomp)%F_Barrier(i)%level_param(6) =                          &
                    nucleus(icomp)%F_Barrier(i)%level_param(6) - 0.05
            call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,            &
                           nucleus(icomp)%F_Barrier(i)%vib_enh,                   &
                           nucleus(icomp)%F_Barrier(i)%rot_enh)
            if(nucleus(icomp)%F_Barrier(i)%level_param(6) < 3.0)goto 11
            call cumm_rho(nfit,elv,nucleus(icomp)%A,nucleus(icomp)%level_param,   &
                          nucleus(icomp)%vib_enh,                                 &
                          nucleus(icomp)%rot_enh,cum_fit)
            chisq=0.0d0
            do j=1,nfit
               chisq = chisq + (cum_rho(j)-cum_fit(j))**2/dcum_rho(i)**2
            end do

            if(chisq < chi_min)then
               ematch = nucleus(icomp)%F_Barrier(i)%level_param(6)
               chi_min = chisq
            else
               goto 11
            end if
            goto 10
 11         nucleus(icomp)%F_Barrier(i)%level_param(6) = ematch
            call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,           &
                           nucleus(icomp)%F_Barrier(i)%vib_enh,                  &
                           nucleus(icomp)%F_Barrier(i)%rot_enh)
            if(allocated(elv))deallocate(elv)
            if(allocated(cum_rho))deallocate(cum_rho)
            if(allocated(dcum_rho))deallocate(dcum_rho)
            if(allocated(cum_fit))deallocate(cum_fit)
         else
!--------------------------------------------------------------------------------
!------   Fail safe, shouldn't get here, but need to do something to proceed
!------   too few levels to fit ematch, and ematch has not been specified
!------   punt and set to ematch of the compound nucleus
!--------------------------------------------------------------------------------
            ematch = nucleus(icomp)%F_Barrier(i)%level_param(6)
            if(ematch < 1.0d-6)then
               nucleus(icomp)%F_Barrier(i)%level_param(6) = nucleus(icomp)%level_param(6)
               ematch = nucleus(icomp)%F_Barrier(i)%level_param(6)
            end if
            call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,            &
                           nucleus(icomp)%F_Barrier(i)%vib_enh,                   &
                           nucleus(icomp)%F_Barrier(i)%rot_enh)
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
         Ematch = nucleus(icomp)%F_Barrier(i)%level_param(6)
         sig2 = sig2_param(Ematch,nucleus(icomp)%level_param,A)
         sig2F = sig2_param(Ematch,nucleus(icomp)%F_Barrier(i)%level_param,A)
         ratio = sig2F/sig2
         sg2cut = nucleus(icomp)%level_param(8)*ratio
         nucleus(icomp)%F_Barrier(i)%level_param(8) = sg2cut
         nucleus(icomp)%F_Barrier(i)%ecut = nucleus(icomp)%level_param(7)
         nucleus(icomp)%F_Barrier(i)%level_param(7) = nucleus(icomp)%F_Barrier(i)%ecut
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
   integer(kind=4) :: i, j, j_min
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
   integer(kind=4) :: ia
   real(kind=8) :: e_old
   logical converged
   real(kind=8) :: pmode, pe1, pbb
   save e_old

   real(kind=8) :: E0, T, E1
   real(kind=8) :: b, Max_J, F_Barrier, F_Barrier_hbw
   real(kind=8) :: aa, bb, cc
!--------   External Functions  ----------------------
   real(kind=8) :: HW_trans
   logical :: real8_equal
!-----------------------------------------------------

   par = 1.0d0
   if(ipar == 0)par = -1.0d0

   de = 0.005
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
!      if(E0 < 0.0d0)then
!         E1 = T*log(1.0d0-exp(E0/T))
!      else
!         E1 = -5.0d0
!      end if
      E1 = 0.0d0
      j_min = nint(E1/de)

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


      do j = 1,nucleus(icomp)%F_Barrier(i)%num_discrete
         if(real8_equal(nucleus(icomp)%F_barrier(i)%state_j(j),xji) .and.   &
            real8_equal(nucleus(icomp)%F_barrier(i)%state_pi(j),par))then
            Ef = nucleus(icomp)%F_barrier(i)%state_e(j)
            tt = HW_trans(Ex,Ef,F_Barrier,F_barrier_hbw)
            trans(i) = trans(i) + tt
         end if
      end do

      pmode = nucleus(icomp)%F_barrier(i)%level_param(16)
      pe1 = nucleus(icomp)%F_barrier(i)%level_param(17)
      pbb = nucleus(icomp)%F_barrier(i)%level_param(18)

      emin = E1
      if(nucleus(icomp)%F_Barrier(i)%num_discrete > 0)emin = ecut

      T_f = 0.0d0
      j = 0

      converged = .false.

      do while (.not. converged)
         Ef = real(j,kind=8)*de + emin 

         call rho_J_par_e(Ef, xji, ipar,                              &
                          nucleus(icomp)%F_barrier(i)%level_param,    &
                          nucleus(icomp)%F_barrier(i)%vib_enh,        &
                          nucleus(icomp)%F_barrier(i)%rot_enh,        &
                          ia, rho, apu, sig2, K_vib, K_rot)
         tt = HW_trans(Ex,Ef,F_Barrier,F_barrier_hbw)
         tt = tt*rho*de
         if(T_f > 0.0 .and. tt/T_f < 1.0d-6)converged = .true.
         if(Ef > Ex .and. tt < 1.0d-7)converged = .true.
         T_f = T_f + tt
         j = j + 1
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
         beta_2 = 3.7d0*nucleus(i)%beta(2)
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
   end do

   return
end subroutine init_barrier_data

