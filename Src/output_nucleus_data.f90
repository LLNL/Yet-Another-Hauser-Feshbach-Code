!
!*******************************************************************************
!
subroutine output_nucleus_data(j_max, itarget)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine prints out data for each nucleus in the decay network
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nodeinfo
!        options
!        print_control
!        useful_data
!        nuclei
!        Channel_info
!        particles_def
!        directory_structure
!        constants
!
!     Subroutines:
!
!        rho_J_par_E
!        rhoe
!        cumm_rho
!
!     External functions:
!
!        real(kind=8) :: parity_fac
!        real(kind=8) :: EL_f
!        real(kind=8) :: EL_f_component
!        real(kind=8) :: EL_trans
!        real(kind=8) :: ML_f
!        real(kind=8) :: ML_f_component
!        real(kind=8) :: ML_trans
!        real(kind=8) :: EL_absorption
!        real(kind=8) :: ML_absorption
!        integer(kind=4) :: find_ibin
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
   use print_control
   use useful_data
   use nuclei
   use Channel_info
   use particles_def
   use directory_structure
   use constants
   implicit none
!-------------------------------------------------------
   integer(kind=4), intent(in) :: j_max
   integer(kind=4), intent(in) :: itarget
!-------------------------------------------------------
   integer(kind=4) :: i, j, k, m, n, jj, ip, ii
   integer(kind=4) :: if1, ifi
   integer(kind=4) :: iA, iZ, iN
   integer(kind=4) :: l_radiation
   integer(kind=4) :: iprint,printZA(2,100)
   integer(kind=4) :: n_bin
   integer(kind=4) :: numcc
!---------------------------------------------------------------------
   integer(kind=4) :: nfit
   real(kind=8),allocatable :: cum_rho(:), cumm_fit(:), elv(:)
   real(kind=8) :: energy, prob, prob1(0:60), prob_jpi(0:60,0:1)
   real(kind=8) :: rho(0:60,0:1)
   real(kind=8) :: rho_Fm, rho_J_ip, sig2, apu
   real(kind=8) :: xj, pfac, dde
   real(kind=8) :: E0, T, E1
   real(kind=8) :: sum_rho
   real(kind=8) :: K_vib, K_rot
   real(kind=8) :: f_E, trans_E, f_M, trans_M
   real(kind=8), allocatable, dimension (:) :: gsf
   real(kind=8) :: EL_cs, ML_cs
   character(len=1) :: char
   character(len=1) :: char_pos,char_neg
   integer(kind=4) :: lev_below, char_start
   logical :: lprint
   real(kind=8) :: pmode, pe1, pbb
   real(kind=8) :: part_prob(0:7)
   integer(kind=4) :: num_points
   real(kind=8) :: Ex
   integer(kind=4) :: Ix
   character(len=132) :: temp_string
   character(len=200) :: fstring
   character(len=5) :: nuke_label
   integer(kind=4) :: inuke_end
   integer(kind=4) :: nnn
   character(len=7) :: f_units


!   real(kind=8) :: ematch, delta, Um, aparam
!   real(kind=8) :: sig2_em, sig2_min, sg2cut, sig
!   real(kind=8) :: deriv, ecut, shell, gamma
!   integer(kind=4) :: sig_model

!---------   External functions
   real(kind=8) :: parity_fac
   real(kind=8) :: EL_f
   real(kind=8) :: EL_f_component, EL_trans
   real(kind=8) :: ML_f, ML_f_component
   real(kind=8) :: ML_trans
   real(kind=8) :: EL_absorption
   real(kind=8) :: ML_absorption
   integer(kind=4) :: find_ibin
!   real(kind=8) :: sig2_param
!   real(kind=8) :: aparam_u

!--------------   Start subrotuine
   num_points = int(30.0d0/de,kind=4) + 1

   char_pos = '+'
   char_neg = '-'
   iprint = 0
   printZA(1:2,1:100) = 0
   do i = 1, num_comp
      lprint = .false.

      do j = 1, iprint           !   Check to see if this nucleus has been printed already
         if(nucleus(i)%Z == printZA(1,j) .and. nucleus(i)%A == printZA(2,j))then
            lprint=.true.            ! found in lest of previous printed nuclei
            exit
         end if
      end do
      if(lprint)cycle                 ! It has been printed before so exit
      iprint = iprint+1                ! It has been printed before so add to list
      printZA(1,iprint) = nucleus(i)%Z
      printZA(2,iprint) = nucleus(i)%A

      iA = nucleus(i)%A
      iZ = nucleus(i)%Z
      iN = iA - iZ
      write(13,*)'+++++++++++++++++++++++++++++++++++++++++++++++'
      write(13,*)
      write(13,*)'Information for compound nucleus #',i
      if(i == itarget)write(13,'(''Target Nucleus'')')
      char_start = 1
      if(nucleus(i)%atomic_symbol(1:1) == ' ')char_start = 2
      nuke_label(1:5) = '     '
      inuke_end = 5
      if(iA < 10)then
         write(nuke_label(1:1),'(i1)')iA
         inuke_end = 1
      elseif(iA < 100)then
         write(nuke_label(1:2),'(i2)')iA
         inuke_end = 2
      elseif(iA < 1000)then
         write(nuke_label(1:3),'(i3)')iA
         inuke_end = 3
      end if
      if(char_start == 1)then
         nuke_label(inuke_end+1:inuke_end+2) = nucleus(i)%atomic_symbol(char_start:2)
         inuke_end = inuke_end + 2
      else
         nuke_label(inuke_end+1:inuke_end+1) = nucleus(i)%atomic_symbol(char_start:2)
         inuke_end = inuke_end + 1
      end if

      write(13,'(a5)')nuke_label(1:inuke_end)
      write(13,'(''Mass excess = '',f10.3,'' MeV'')')nucleus(i)%ME
      write(13,'(''Binding energy = '',f10.3,'' MeV'')')nucleus(i)%BE
      write(13,'(''Pairing model used = '',i3)')nucleus(i)%pair_model
      if(lev_option > 0)then
         if(nucleus(i)%pair_model == 0)then
            write(13,'(''Pairing gap derived from Gilbert & Cameron pairing model'')')
            write(13,'(''Pairing Energy, Delta = '',f10.3,'' MeV'')')nucleus(i)%level_param(3)
         elseif(nucleus(i)%pair_model == 1)then
            write(13,'(''Pairing gap derived from systematics'')')
            if(iand(iZ,1) == 0 .and. iand(iN,1) == 0)write(13,'(''Even-even nucleus: Delta = 24.0/sqrt(A) MeV'')')
            if(iand(iZ,1) == 0 .and. iand(iN,1) == 1)write(13,'(''Even-odd nucleus: = 12.0/sqrt(A) MeV'')')
            if(iand(iZ,1) == 1 .and. iand(iN,1) == 0)write(13,'(''Odd-even nucleus: = 12.0/sqrt(A) MeV'')')
            if(iand(iZ,1) == 1 .and. iand(iN,1) == 1)write(13,'(''Odd-Odd nucleus: Delta = 0 MeV'')')
            write(13,'(''Pairing Energy, Delta = '',f10.3,'' MeV'')')nucleus(i)%level_param(3)
         elseif(nucleus(i)%pair_model == 2)then
            write(13,'(''Pairing gap derived from binding energies'')')
            if(iand(iZ,1) == 0 .and. iand(iN,1) == 0)write(13,'(''Even-even nucleus: Delta = (Delta_p + Delta_n)'')')
            if(iand(iZ,1) == 0 .and. iand(iN,1) == 1)write(13,'(''Even-odd nucleus: Delta = (Delta_p + Delta_n)/2'')')
            if(iand(iZ,1) == 1 .and. iand(iN,1) == 0)write(13,'(''Odd-even nucleus: Delta = (Delta_p + Delta_n)/2'')')
            if(iand(iZ,1) == 1 .and. iand(iN,1) == 1)write(13,'(''Odd-Odd nucleus: Delta = 0 MeV'')')
            write(13,'(''Proton Pairing Energy, Delta_p = '',f10.3,'' MeV'')')nucleus(i)%Delta_exp_p
            write(13,'(''Neutron Pairing Energy, Delta_n ='',f10.3,'' MeV'')')nucleus(i)%Delta_exp_n
            write(13,'(''Pairing Energy, Delta = '',f10.3,'' MeV'')')nucleus(i)%level_param(3)
         elseif(nucleus(i)%pair_model == 3)then
            write(13,'(''Pairing gap derived from user input'')')
            write(13,'(''Pairing Energy, Delta = '',f10.3,'' MeV'')')nucleus(i)%level_param(3)
         end if
      else
         write(13,'(''Pairing Energy, Delta = '',f10.3,'' MeV'')')nucleus(i)%level_param(3)
      end if
      write(13,'(''Deformation parameters'')')
      do k = 2, 6
          write(13,'('' beta('',i1,'')'','' = '',f10.5)') k,nucleus(i)%beta(k)
      end do
      write(13,*)'Particle separation energies '
      do k=1,6
         write(13,'(2x,a8,1x,''Separation energy = '',f10.3,'' MeV'')')              &
           particle(k)%name,nucleus(i)%sep_e(k)
      end do

      if(nucleus(i)%D0exp > 0.0d0)then
         write(13,'('' Experimental D0   ='',f10.3,'' +/- '',f10.3,'' eV'')')        &
                     nucleus(i)%D0exp,nucleus(i)%dD0exp
      else
         write(13,'('' Experimental D0   =  UNAVAILABLE'')')
      end if
      if(nucleus(i)%D0 > 0.0d0)then
          write(13,'('' Calculated D0     ='',f10.3,'' eV'')')                       &
                  nucleus(i)%D0
      else
         write(13,'('' D0 Not Calculated'')')
      end if

      if(nucleus(i)%D1exp > 0.0d0)then
         write(13,'('' Experimental D1 ='',f10.3,'' +/- '',f10.3,'' eV'')')          &
                     nucleus(i)%D1exp,nucleus(i)%dD1exp
      else
         write(13,'('' Experimental D1   =  UNAVAILABLE'')')
      end if
      if(nucleus(i)%D1 > 0.0d0)then
          write(13,'('' Calculated D1     ='',f10.3,'' eV'')')                       &
                  nucleus(i)%D1
      else
         write(13,'('' D1 Not Calculated'')')
      end if


      if(nucleus(i)%Gamma_g >= 0.0d0)write(13,'('' Calcuated Gamma_gamma(l=0)    = '',f12.3,'' meV'')')nucleus(i)%Gamma_g
      if(nucleus(i)%Gamma_g < 0.0d0)write(13,'('' Gamma_gamma(l=0) is NOT CALCULATED; Ex_max < Sn in this run'')')
      if(nucleus(i)%Gamma_g_exp > 0.0d0)then
         if(nucleus(i)%dGamma_g_exp > 0.0d0)then
            write(13,'('' Experimental Gamma_gamma(l=0) = '',f12.3,'' +/- '',f12.3,'' meV'')')  &
                               nucleus(i)%Gamma_g_exp,                                          &
                               nucleus(i)%dGamma_g_exp
         else
            write(13,'('' User set value Gamma_gamma(l=0) = '',f12.3,'' meV'')')                &
                            nucleus(i)%Gamma_g_exp
         end if
      else
         write(13,'('' Experimental Gamma_gamma(l=0) = UNAVAILABLE'')')
      end if
      if(nucleus(i)%Gamma_g_1 >= 0.0d0)write(13,'('' Calcuated Gamma_gamma(l=1)    = '',f12.3,'' meV'')')nucleus(i)%Gamma_g_1
      if(nucleus(i)%Gamma_g_1 < 0.0d0)write(13,'('' Gamma_gamma(l=1) is NOT CALCULATED; Ex_max < Sn in this run'')')
      if(nucleus(i)%Gamma_g_1_exp > 0.0d0)then
         write(13,'('' Experimental Gamma_gamma(l=1) = '',f12.3,'' +/- '',f12.3,'' meV'')')   &
                            nucleus(i)%Gamma_g_1_exp,                                         &
                            nucleus(i)%dGamma_g_1_exp
      else
         write(13,'('' Experimental Gamma_gamma(l=1) = UNAVAILABLE'')')
      end if
!------------------------------------------------------------------
      write(13,*)
      write(13,'(''******************************************************************'')')
      write(13,'(''  Data on discrete levels are based on the RIPL-3 data, with     *'')')
      write(13,'(''  modifcations by the author and users. RIPL-3 data is           *'')')
      write(13,'(''  available at https://www-nds.iaea.org/RIPL-3/. Discrete level  *'')')
      write(13,'(''  information is based on M. Verpelli and R. Capote Noy,         *'')')
      write(13,'(''  INDC(NDS)-0702, IAEA, 2015, and R. Capote, et al., Nuclear     *'')')
      write(13,'(''  Data Sheets 110 (2009), 3107-3214.                             *'')')
      write(13,'(''******************************************************************'')')
      write(13,*)
      write(13,*)'Discrete states used in calculation'
      write(13,'(''State'',6x,''Energy'',2x,''Spin'','//               &
                 '1x,''Par'',8x,''T1/2'',6x,'//                        &
                 '''Isomer'',2x,''Modified'')')
      lev_below = 0
      do k = 1,nucleus(i)%num_discrete
         if(nucleus(i)%state(k)%energy > nucleus(i)%level_param(7)     &
           .and. lev_below == 0)lev_below = 1
         if(lev_below /= 1)then
            if(i == itarget .and. k == target%istate)then
               write(13,'(2x,i3,2x,f10.5,1x,f5.1,1x,f3.0,1x,'//        &
                          '1x,e15.7,1x,l4,4x,l4,6x,''Target State'')') &
                 k,                                                    &
                 nucleus(i)%state(k)%energy,                           &
                 nucleus(i)%state(k)%spin,                             &
                 nucleus(i)%state(k)%parity,                           &
                 nucleus(i)%state(k)%t12,                              &
                 nucleus(i)%state(k)%isomer,                           &
                 nucleus(i)%state(k)%state_modified
            else
               write(13,'(2x,i3,2x,f10.5,1x,f5.1,1x,f3.0,1x,'//        &
                          '1x,e15.7,1x,l4,4x,l4)')                     &
                 k,                                                    &
                 nucleus(i)%state(k)%energy,                           &
                 nucleus(i)%state(k)%spin,                             &
                 nucleus(i)%state(k)%parity,                           &
                 nucleus(i)%state(k)%t12,                              &
                 nucleus(i)%state(k)%isomer,                           &
                 nucleus(i)%state(k)%state_modified
            end if
        elseif(lev_below == 1)then
            lev_below = 2
            write(13,'(''   --------------------   Ecut   ---------------------------------'')')
            write(13,'(2x,i3,2x,f10.5,1x,f5.1,1x,f3.0,1x,'//           &
                       '1x,e15.7,1x,l4,4x,l4)')                        &
              k,                                                       &
              nucleus(i)%state(k)%energy,                              &
              nucleus(i)%state(k)%spin,                                &
              nucleus(i)%state(k)%parity,                              &
              nucleus(i)%state(k)%t12,                                 &
              nucleus(i)%state(k)%isomer,                              &
              nucleus(i)%state(k)%state_modified
        end if
      end do

      write(13,*)'Decay properties of discrete states'
      do k = nucleus(i)%num_discrete, 2, -1
         write(13,*)'Decay of state',k
         write(13,*)'Number of transitions',nucleus(i)%state(k)%nbranch
         write(13,'(''       i --->   f'',12x,''branch'',6x,''prob_gamma'',9x,''prob_ic'',9x,''Modified'')')
         write(13,'(''     ---      ---'',3(5x,''-------------''),5x,''--------'')')

         prob = 0.0d0
         do m = 1, nucleus(i)%state(k)%nbranch
            prob = prob + nucleus(i)%state(k)%branch(m)
            write(13,'(4x,i4,'' --->'',i4,3(3x,e15.7),5x,l4)')         &
                k,nucleus(i)%state(k)%ibranch(m),                      &
                nucleus(i)%state(k)%branch(m),                         &
                nucleus(i)%state(k)%p_gamma(m),                        &
                nucleus(i)%state(k)%p_ic(m),                           &
                nucleus(i)%state(k)%branch_modified(m)

         end do
      end do

      if(i == itarget)then
         write(13,*)
         write(13,'(''Coupled-channels information for target nucleus'')')
         write(13,'(''Number of coupled-channels = '',i4)')OpticalCS%numcc
         numcc = 0
         do k = 1, OpticalCS%numcc
            if(OpticalCS%state(k)%state_type == 1)numcc = numcc + 1
         end do
         write(13,'(''States coupled to the ground-state band'')')
         write(13,*)
         write(13,'(''    N    J  PAR     Energy'')')
         write(13,'(''  ---  ---  ---     ------'')')
         do k = 1, numcc
            write(13,'(1x,i4,2(1x,f4.1),1x,f10.4)')k,OpticalCS%state(k)%spin,      &
                OpticalCS%state(k)%parity,OpticalCS%state(k)%energy
         end do
         write(13,*)
         write(13,'(''Additional DWBA states'')')
         write(13,*)
         write(13,'(''    N    J  PAR    K     Energy'')')
         write(13,'(''  ---  ---  ---  ---     ------'')')
         do k = numcc + 1, OpticalCS%numcc
            write(13,'(1x,i4,3(1x,f4.1),1x,f10.4)')k,OpticalCS%state(k)%spin,      &
                OpticalCS%state(k)%parity,OpticalCS%state(k)%K,OpticalCS%state(k)%energy
         end do
      end if


!------------------------------------------------------------------
      write(13,*)
      write(13,'(''Level Density information for '',a5)')nuke_label(1:inuke_end)
      if(nucleus(i)%fit_D0 .and. nucleus(i)%fit_aparam)then
         write(13,*)'Fitting to D0 by adjusting the a-parameter'
      elseif(nucleus(i)%fit_D0 .and. .not. nucleus(i)%fit_aparam)then
         write(13,*)'Fitting to D0 by adjusting the shell correction'
      elseif(.not. nucleus(i)%fit_D0)then
         write(13,*)'No fit to D0 is performed'
      end if
      if(nucleus(i)%D0exp > 0.0d0)then
         write(13,'('' Experimental D0 ='',f10.3,'' +/- '',f10.3,'' eV'')')    &
                     nucleus(i)%D0exp,nucleus(i)%dD0exp
      else
         write(13,'('' Experimental D0 =  UNAVAILABLE'')')
      end if
      if(nucleus(i)%D0 > 0.0d0)then
         write(13,'('' Calculated D0 ='',f10.3)')    &
                     nucleus(i)%D0
      else
         write(13,'('' D0 Not Calculated'')')
      end if
!-----------------------------------------------------------------------------
!------    Modeled or read in level densities
!------    First in the list is modeled
!-----------------------------------------------------------------------------
      if(.not. nucleus(i)%lev_den_read)then
         write(13,'(''E1 defined as the energy so that int(E1,0)rho(E)dE = 1'')')
         write(13,'(''Generally, E1 < 0 if E0 < 0'')')
         write(13,'(''Otherwise E1 undefined as int(-infty,0)rho(E)DE < 1'')')
         E0 = nucleus(i)%level_param(15)
         T = nucleus(i)%level_param(14)
         E1 = 0.0d0
         if(E0 < 0.0d0)E1 = T*log(1.0d0-exp(E0/T))
         write(13,*)'Level-density parameters'
         write(13,*)'Level-density model = ',nucleus(i)%level_model
         write(13,'('' a =        '',f12.7)')nucleus(i)%level_param(1)
         write(13,'('' lambda =   '',f12.7)')nucleus(i)%level_param(2)
         write(13,'('' Delta =    '',f12.7)')nucleus(i)%level_param(3)
         write(13,'('' Shell =    '',f12.7)')nucleus(i)%level_param(4)
         write(13,'('' gamma =    '',f12.7)')nucleus(i)%level_param(5)
         write(13,'('' E_match =  '',f12.7)')nucleus(i)%level_param(6)
         write(13,'('' E_cut =    '',f12.7)')nucleus(i)%level_param(7)
         write(13,'('' sigma_cut ='',f12.7)')nucleus(i)%level_param(8)
         write(13,'('' T =        '',f12.7)')nucleus(i)%level_param(14)
         write(13,'('' E0 =       '',f12.7)')nucleus(i)%level_param(15)
         write(13,'('' a(Sn) =    '',f12.7)')nucleus(i)%a_Sn
         write(13,'('' sig2(Sn) = '',f12.7)')nucleus(i)%sig2_Sn
         if(E1 < 0.0d0)then
             write(13,'('' E1 =       '',f12.7)')E1
         else
             write(13,'(''E1 is undefined as E0 > 0, E1 ~ -10.0'')')
         end if
         write(13,*)
         if(nint(nucleus(i)%level_param(10)) == 0)then
            write(13,*)'No rotational collective enhancement included'
            write(13,*)'K_rot = 1.0'
         else
            write(13,*)'Rotational Collective enhancement model = ',    &
                nucleus(i)%level_param(10)
            write(13,*)'Collective enhancement factors'
            write(13,*)'K_rot = Max[x(1)*(Factor-1)/(1+exp((E-x(2))/x(3))),0]+1'
!            if(nint(nucleus(i)%level_param(10)) == 1)then
!               write(13,*)'Factor = x(1)*sig2*(1+beta2/3)'
!            elseif(nint(nucleus(i)%level_param(10)) == 2)then
!               write(13,*)'Factor = 2.0*x(1)*sig2*(1+beta2/3)'
!            elseif(nint(nucleus(i)%level_param(10)) == 3)then
!               write(13,*)'Factor = x(1)*sig2**3/2*(1+beta2/3)*(1-2*beta2/3)'
!            elseif(nint(nucleus(i)%level_param(10)) == 4)then
!               write(13,*)'Factor = 2.0*x(1)*sig2**3/2*(1+beta2/3)*(1-2*beta2/3)'
!            end if

            if(nint(nucleus(i)%level_param(10)) == 1)then
                   if(nint(nucleus(i)%level_param(19)) == 2 .or.            &
                      nint(nucleus(i)%level_param(19)) == 3)then
       write(13,*)'Factor = x(1)*sig2*(1+beta2/3)'
       write(13,*)' --  Axially symmetric'
                   elseif(nint(nucleus(i)%level_param(19)) == 4 .or.            &
                      nint(nucleus(i)%level_param(19)) == 5)then
       write(13,*)'Factor = x(1)*sig2*(1+sqrt(5/4pi)*beta2)'
       write(13,*)' --  Axially symmetric'
                   end if
                elseif(nint(nucleus(i)%level_param(10)) == 2)then
                   if(nint(nucleus(i)%level_param(19)) == 2 .or.            &
                      nint(nucleus(i)%level_param(19)) == 3)then
       write(13,*)'Factor = 2.0*x(1)*sig2*(1+beta2/3)'
       write(13,*)' -- Left-right asymmetric'
                   elseif(nint(nucleus(i)%level_param(19)) == 4 .or.            &
                      nint(nucleus(i)%level_param(19)) == 5)then
       write(13,*)'Factor = 2.0*x(1)*sig2*(1+sqrt(5/4pi)*beta2)'
       write(13,*)' -- Left-right asymmetric'
                   end if
                elseif(nint(nucleus(i)%level_param(10)) == 3)then
                   if(nint(nucleus(i)%level_param(19)) == 2 .or.            &
                      nint(nucleus(i)%level_param(19)) == 3)then
       write(13,*)'Factor = x(1)*sqrt(pi/2)**sig2**(3/2)*(1+beta2/3)*(1-2*beta2/3)'
       write(13,*)' -- Triaxial and left-right asymmetric'
                   elseif(nint(nucleus(i)%level_param(19)) == 4 .or.            &
                      nint(nucleus(i)%level_param(19)) == 5)then
       write(13,*)'Factor = x(1)*sqrt(pi/2)**sig2**(3/2)*(1+sqrt(5/4pi)*beta2)*(1-sqrt(5/4pi)beta2/2)'
       write(13,*)' -- Triaxial and left-right asymmetric'
                   end if
                elseif(nint(nucleus(i)%level_param(10)) == 4)then
                   if(nint(nucleus(i)%level_param(19)) == 2 .or.            &
                      nint(nucleus(i)%level_param(19)) == 3)then
       write(13,*)'Factor = 2.0*x(1)*sqrt(pi/2)**sig2**(3/2)*(1+beta2/3)*(1-2*beta2/3)'
       write(13,*)' -- Triaxial and not left-right asymmetric'
                   elseif(nint(nucleus(i)%level_param(19)) == 4 .or.            &
                      nint(nucleus(i)%level_param(19)) == 5)then
       write(13,*)'Factor = 2.0*x(1)*sqrt(pi/2)**sig2**(3/2)*(1+sqrt(5/4pi)*beta2)*(1-sqrt(5/4pi)beta2/2)'
       write(13,*)' -- Triaxial and not left-right asymmetric'
                   end if
                end if



            write(13,*)'Rotational collective-enhancement parameters'

            write(13,'(''x(1) = '',f10.4)')nucleus(i)%rot_enh(1)
            write(13,'(''x(2) = '',f10.4)')nucleus(i)%rot_enh(2)
            write(13,'(''x(3) = '',f10.4)')nucleus(i)%rot_enh(3)
            write(13,*)
         end if
         if(nint(nucleus(i)%level_param(11)) == 0)then
            write(13,*)'No vibrational collective enhancement included'
            write(13,*)'K_vib = 1.0'
         elseif(nint(nucleus(i)%level_param(11)) == 1)then
            write(13,*)'Vibrational collective enhancement model 1 used'
            write(13,*)
            write(13,*)'K_vib = exp(0.05555*A**(2./3.)*T**(4./3.))'
            write(13,*)'K_vib = Max[x(1)*(K_vib-1)/(1+exp((E-x(2))/x(3))),0]+1'
            write(13,*)'T = sqrt(U/apu)'
            write(13,*)
            write(13,'(''x(1) = '',f10.4)')nucleus(i)%vib_enh(1)
            write(13,'(''x(2) = '',f10.4)')nucleus(i)%vib_enh(2)
            write(13,'(''x(3) = '',f10.4)')nucleus(i)%vib_enh(3)
            write(13,*)
         elseif(nint(nucleus(i)%level_param(11)) == 2)then
            write(13,*)'K_vib model 2 - approximated with Bose-gas relationship'
            write(13,*)
            write(13,*)'K_vib = exp(dS - dU/T)'
            write(13,*)'T = sqrt(U/apu)'
            write(13,*)
            write(13,*)'dS = Sum_i (2*L(i)+1)*[(1+n(i))*log(1+n(i)) - n(i)*log(n(i))]'
            write(13,*)'dU = Sum_i (2*L(i)+1)*w(i)*n(i)'
            write(13,*)'n(i) = exp(-g(i)/2*w(i))/[exp(w(i)/T) - 1]'
            write(13,*)'g(i) = C*(w(i)**2 + 4*pi**2T**2)'
            write(13,*)
            write(13,*)'C = 0.0075*A**(1/3)'
            write(13,*)'w(2) = 65*A**(-5/6)/(1*0.05*Shell)'
            write(13,*)'w(3) = 100*A**(-5/6)/(1*0.05*Shell)'
            write(13,*)
            write(13,*)'The sum extends over quadrupole and octupole vibrations, L=2,3'
         elseif(nint(nucleus(i)%level_param(11)) == 3)then
            write(13,*)'Vibrational collective enhancement model 1 used'
            write(13,*)
            write(13,*)'K_vib = exp(0.05555*A**(2./3.)*T**(4./3.))'
            write(13,*)'K_vib = Max[x(1)*(K_vib-1)/(1+exp((E-x(2))/x(3))),0]+1'
            write(13,*)'T = sqrt(U/a);  a = A/13'
            write(13,*)
            write(13,'(''x(1) = '',f10.4)')nucleus(i)%vib_enh(1)
            write(13,'(''x(2) = '',f10.4)')nucleus(i)%vib_enh(2)
            write(13,'(''x(3) = '',f10.4)')nucleus(i)%vib_enh(3)
            write(13,*)
         elseif(nint(nucleus(i)%level_param(11)) == 4)then
            write(13,*)'K_vib model 2 - approximated with Bose-gas relationship'
            write(13,*)
            write(13,*)'K_vib = exp(dS - dU/T)'
            write(13,*)'T = sqrt(U/a);  a = A/13'
            write(13,*)
            write(13,*)'dS = Sum_i (2*L(i)+1)*[(1+n(i))*log(1+n(i)) - n(i)*log(n(i))]'
            write(13,*)'dU = Sum_i (2*L(i)+1)*w(i)*n(i)'
            write(13,*)'n(i) = exp(-g(i)/2*w(i))/[exp(w(i)/T) - 1]'
            write(13,*)'g(i) = C*(w(i)**2 + 4*pi*T**2)'
            write(13,*)
            write(13,*)'C = 0.0075*A**(1/3)'
            write(13,*)'w(2) = 65*A**(-5/6)/(1*0.05*Shell)'
            write(13,*)'w(3) = 100*A**(-5/6)/(1*0.05*Shell)'
            write(13,*)
            write(13,*)'The sum extends over quadrupole and octupole vibrations, L=2,3'
         end if
         write(13,*)

         write(13,*)'Rho = K_rot*K_vib*Rho'

         write(13,*)
         if(nint(nucleus(i)%level_param(16)) == 0)then
             write(13,*)'Parity Factor = 0.5 for positive and negative parities'
         elseif(nint(nucleus(i)%level_param(16)) == 1)then
             write(13,*)'Parity Factor = 0.5*tanh(a*(E-E0)) for positive parity states'
             write(13,'(''E0 = '',f5.3)')nucleus(i)%level_param(17)
             write(13,'(''a = '',f5.3)')nucleus(i)%level_param(18)
         elseif(nint(nucleus(i)%level_param(16)) == 2)then
             write(13,*)'Parity Factor = 0.5*tanh(a*(E-E0)) for negative parity states'
             write(13,'(''E0 = '',f5.3)')nucleus(i)%level_param(17)
             write(13,'(''a = '',f5.3)')nucleus(i)%level_param(18)
         end if
         write(13,*)
 
         write(13,'(''Level density States/MeV for '',a5)')nuke_label(1:inuke_end)
         write(13,'(''Positive Parity '')')

         write(temp_string,*)min(j_max,60)+1
         fstring = "(6x,'E',9x,'a(U)',1x,'sqrt(sig2)',6x,'P_fac',"//                &
                   "6x,'K_vib',6x,'K_rot',8x,'Enh',9x,'Rho',5x,"                    &
                   //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

         write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
         fstring = "(1x,'--------',6(1x,'----------'),1x,'---------------',"        &
                   //trim(adjustl(temp_string))//"(1x,'---------------'))"
         write(13,fstring)
         pmode = nucleus(i)%level_param(16)
         pe1 = nucleus(i)%level_param(17)
         pbb = nucleus(i)%level_param(18)

         fstring = "(1x,f8.3,6(1x,f10.3),(1x,e15.7),"//trim(adjustl(temp_string))//"(1x,e15.7))"
         ip = 1
         do k = 1, nucleus(i)%nbin
            energy = nucleus(i)%e_grid(k)
             call rhoe(energy, nucleus(i)%level_param,                              &
                       nucleus(i)%vib_enh,                                          &
                       nucleus(i)%rot_enh,                                          &
                       ia, rho_Fm, apu, sig2, K_vib, K_rot)
             pfac = parity_fac(energy,xj,ip,pmode,pe1,pbb)
             write(13,fstring)                                                      &
                  energy,apu,sqrt(sig2),pfac,K_vib,K_rot,K_vib*K_rot,rho_FM*pfac,   &
                  (nucleus(i)%bins(jj,ip,k)%rho,jj = 0, min(j_max,60))
         end do


!----   Format to write j_max elements to file
         write(temp_string,*)min(j_max,60)+1
         write(13,'(''Negative Parity '')')
         fstring = "(6x,'E',9x,'a(U)',1x,'sqrt(sig2)',6x,'P_fac',"//                &
                   "6x,'K_vib',6x,'K_rot',8x,'Enh',9x,'Rho',5x,"                    &
                   //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"
         write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
         fstring = "(1x,'--------',6(1x,'----------'),1x,'---------------',"        &
                   //trim(adjustl(temp_string))//"(1x,'---------------'))"
         write(13,fstring)
         pmode = nucleus(i)%level_param(16)
         pe1 = nucleus(i)%level_param(17)
         pbb = nucleus(i)%level_param(18)

         fstring = "(1x,f8.3,6(1x,f10.3),(1x,e15.7),"//trim(adjustl(temp_string))//"(1x,e15.7))"

         ip = 0
         do k = 1, nucleus(i)%nbin
            energy = nucleus(i)%e_grid(k)
            sum_rho = 0.0d0
            call rhoe(energy,nucleus(i)%level_param,                                &
                      nucleus(i)%vib_enh,                                           &
                      nucleus(i)%rot_enh,                                           &
                      ia,rho_Fm,apu,sig2,K_vib,K_rot)
            pfac = parity_fac(energy,xj,ip,pmode,pe1,pbb)
            write(13,fstring)                                                       &
                 energy,apu,sqrt(sig2),pfac,K_vib,K_rot,K_vib*K_rot,rho_FM*pfac,    &
                 (nucleus(i)%bins(jj,ip,k)%rho,jj = 0, min(j_max,60))
         end do
!-----------------------------------------------------------------------------
!------   Next is the case where they are read in
!-----------------------------------------------------------------------------
      else
         write(13,*)
         write(13,*)'Level density read in from user input'
         write(13,'(''Level density States/MeV for '',a5)')nuke_label(1:inuke_end)
         write(13,'(''Positive Parity '')')

         write(temp_string,*)min(j_max,60)+1
         fstring = "(6x,'E',10x,'Rho',5x,"                                              &
                   //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"
         write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
         fstring = "(1x,'--------',1x,'---------------',"                           &
                   //trim(adjustl(temp_string))//"(1x,'---------------'))"
         write(13,fstring)

         fstring = "(1x,f8.3,(1x,e15.7),"//trim(adjustl(temp_string))//"(1x,e15.7))"
         ip = 1
         do k = 1, nucleus(i)%nbin
            energy = nucleus(i)%e_grid(k)
            sum_rho = 0.0d0
            do jj = 0, min(j_max,60)
               sum_rho = sum_rho + nucleus(i)%bins(jj,ip,k)%rho
            end do
            write(13,fstring)                                                       &
                  energy,sum_rho,                                    &
                  (nucleus(i)%bins(jj,ip,k)%rho,jj = 0, min(j_max,60))
         end do
         write(13,*)
         write(13,'(''Negative Parity '')')

         write(temp_string,*)min(j_max,60)+1
         fstring = "(6x,'E',10x,'Rho',5x,"                                              &
                   //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"
         write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
         fstring = "(1x,'--------',1x,'---------------',"                           &
                   //trim(adjustl(temp_string))//"(1x,'---------------'))"
         write(13,fstring)

         fstring = "(1x,f8.3,(1x,e15.7),"//trim(adjustl(temp_string))//"(1x,e15.7))"
         ip = 0
         do k = 1, nucleus(i)%nbin
            energy = nucleus(i)%e_grid(k)
            sum_rho = 0.0d0
            do jj = 0, min(j_max,60)
               sum_rho = sum_rho + nucleus(i)%bins(jj,ip,k)%rho
            end do
            write(13,fstring)                                                       &
                  energy,sum_rho,                                    &
                  (nucleus(i)%bins(jj,ip,k)%rho,jj = 0, min(j_max,60))
         end do

      end if

 !------------------------------------------------------------------
      if(.not. nucleus(i)%lev_den_read)then
         nfit = max(nucleus(i)%num_discrete,nucleus(i)%ncut)
         allocate(cum_rho(nfit))
         allocate(cumm_fit(nfit))
         allocate(elv(nfit))
         if(nucleus(i)%ncut >= 5) then
            if(nucleus(i)%fit_ematch)then
               write(13,*)'E_match was fit to the cummulative density for known discrete states'
               if(nucleus(i)%cum_rho_ratio > 1.25)then
                  write(13,*)
                  write(13,'(''******************************************************************'')')
                  write(13,'(''*****    CAUTION     CAUTION     CAUTION    CAUTION          *****'')')
                  write(13,'(''*****  The fitted cumulative level density for at            *****'')')
                  write(13,'(''*****  E_cut is more than 1.3 larger than the experimental   *****'')')
                  write(13,'(''*****  value. You should check this and possibly change      *****'')')
                  write(13,'(''*****  E_cut with the command "lev_ecut Z A  Value"          *****'')')
                  write(13,'(''*****  Problems of this nature usually means E_cut           *****'')')
                  write(13,'(''*****  is too high                                           *****'')')
                  write(13,'(''******************************************************************'')')
                  write(13,*)
               end if
               if(nucleus(i)%cum_rho_ratio < 0.75)then
                  write(13,*)
                  write(13,'(''******************************************************************'')')
                  write(13,'(''*****    CAUTION     CAUTION     CAUTION    CAUTION          *****'')')
                  write(13,'(''*****  The fitted cumulative level density for at            *****'')')
                  write(13,'(''*****  E_cut is more than 0.7 smaller than the experimental  *****'')')
                  write(13,'(''*****  value. You should check this and possibly change      *****'')')
                  write(13,'(''*****  E_cut with the command "lev_ecut Z A  Value"          *****'')')
                  write(13,'(''*****  Problems of this nature usually means E_cut           *****'')')
                  write(13,'(''*****  is too high                                           *****'')')
                  write(13,'(''******************************************************************'')')
                  write(13,*)
               end if
            else
               write(13,*)'E_match was not fit cumulative density for known discrete states'
               write(13,*)'User input over rides fitting or default was used'
            end if
            nfit = nucleus(i)%ncut

            write(13,*)
            write(13,*)'Cumulative level density up to E_cut'
            cum_rho(1) = 1.0d0
            elv(1) = nucleus(i)%state(1)%energy
            cum_rho(1) = 1.0d0
            do j = 2, nfit
               elv(j) = nucleus(i)%state(j)%energy
               cum_rho(j) = cum_rho(j-1) + 1.0d0
            end do
            do j = 1, nfit
               if(j > 1)cum_rho(j) = cum_rho(j-1) + 1.0d0
               write(13,'(f10.4,1x,f15.2)')elv(j), cum_rho(j)
               if(j < nfit)then
                  write(13,'(f10.4,1x,f15.2)')elv(j+1), cum_rho(j)
               end if
            end do
            if(nucleus(i)%ncut == nucleus(i)%num_discrete)then
               open(unit=23,file = nuke_label(1:inuke_end)//'-CLD-exp.dat',status='unknown')
               write(23,*)'#Cumulative level density up to E_cut'
               cum_rho(1) = 1.0d0
               elv(1) = nucleus(i)%state(1)%energy
               cum_rho(1) = 1.0d0
               do j = 2, nfit
                  elv(j) = nucleus(i)%state(j)%energy
                  cum_rho(j) = cum_rho(j-1) + 1.0d0
               end do
               do j = 1, nfit
                  if(j > 1)cum_rho(j) = cum_rho(j-1) + 1.0d0
                  write(23,'(f10.4,1x,f15.2)')elv(j), cum_rho(j)
                  if(j < nfit)then
                     write(23,'(f10.4,1x,f15.2)')elv(j+1), cum_rho(j)
                  end if
               end do
            end if

            if(nucleus(i)%ncut == nucleus(i)%num_discrete)then
               open(unit=23,file = nuke_label(1:inuke_end)//'-CLD-exp.dat',status='unknown')
                write(23,*)'#Modeled cumulative level density up to E_cut'
            end if
            write(13,*)
            write(13,*)'Modeled cumulative level density up to E_cut'
            cumm_fit(1:nfit) = 0.0d0
            call cumm_rho(nfit,elv,ia,nucleus(i)%level_param,             &
                          nucleus(i)%vib_enh,nucleus(i)%rot_enh,cumm_fit)
            do j = 1, nfit
               write(13,'(f10.4,1x,f15.2)')elv(j),cumm_fit(j)
               if(nucleus(i)%ncut == nucleus(i)%num_discrete)             &
                  write(23,'(f10.4,1x,f15.2)')elv(j),cumm_fit(j)
            end do

            if(nucleus(i)%ncut /= nucleus(i)%num_discrete)then
               nfit = nucleus(i)%num_discrete
!----   Added print to file - needed to edit .out file too often to get this info
               open(unit=23,file = nuke_label(1:inuke_end)//'-CLD-exp.dat',status='unknown')
               write(13,*)
               write(13,*)'Cumulative level density for all discrete states'
               write(23,*)'#Cumulative level density for all discrete states'
               write(23,*)'#E_cut = ',nucleus(i)%level_param(7)
               cum_rho(1:nfit) = 0.0d0
               cum_rho(1) = 1.0d0
               elv(1) = nucleus(i)%state(1)%energy
               cum_rho(1) = 1.0d0
               do j = 2, nfit
                  elv(j) = nucleus(i)%state(j)%energy
                  cum_rho(j) = cum_rho(j-1) + 1.0d0
               end do
               do j = 1, nfit
                  if(j > 1)cum_rho(j) = cum_rho(j-1) + 1.0d0
                  write(13,'(f10.4,1x,f15.2)')elv(j), cum_rho(j)
                  write(23,'(f10.4,1x,f15.2)')elv(j), cum_rho(j)
                  if(j < nfit)then
                     write(13,'(f10.4,1x,f15.2)')elv(j+1), cum_rho(j)
                     write(23,'(f10.4,1x,f15.2)')elv(j+1), cum_rho(j)
                  end if
               end do
               close(unit=23)
               open(unit=23,file = nuke_label(1:inuke_end)//'-CLD-modeled.dat',status='unknown')
               write(13,*)
               write(13,*)'Modeled cumulative level density up to maximum discrete state'
               write(23,*)'#Modeled cumulative level density up to maximum discrete state'
               write(23,*)'#E_cut = ',nucleus(i)%level_param(7)
               cumm_fit(1:nfit) = 0.0d0
               call cumm_rho(nfit,elv,ia,nucleus(i)%level_param,             &
                             nucleus(i)%vib_enh,nucleus(i)%rot_enh,cumm_fit)
               do j = 1, nfit
                  write(13,'(f10.4,1x,f15.2)')elv(j), cumm_fit(j)
                  write(23,'(f10.4,1x,f15.2)')elv(j), cumm_fit(j)
               end do
            end if
         end if
         deallocate(cum_rho)
         deallocate(cumm_fit)
         deallocate(elv)
         close(unit=23)
      end if

      if(.not. allocated(gsf))allocate(gsf(max_num_gsf))

      nucleus(i)%nbin_em = num_points

!------------------------------------------------------------------
!------------------------------------------------------------------
!-----------   Electric multipoles -
!------------------------------------------------------------------
      do l_radiation = 1, nucleus(i)%lmax_E
         write(13,*)
         if(nucleus(i)%EL_mode(l_radiation)%gsf_read)then
            write(13,*)'Gamma Strength Function obtained from manual input'
         else
            write(13,*)'Gamma Strength Function obtained with internal analytic models'
            if(l_radiation == 1)then
               if(nucleus(i)%e1_model == 1)then
                  write(13,'(''E1 model used = '',''Lorenztian'')')
               elseif(nucleus(i)%e1_model == 2)then
                  write(13,'(''E1 model used = '',''Kopecky-Uhl'')')
               end if
            end if
            write(13,'(''Electric-'',i1,'' resonance parameters for '',a5)')          &
            l_radiation,nuke_label(1:inuke_end)
            do k = 1, nucleus(i)%EL_mode(l_radiation)%num_gsf
               write(13,'(''Line ='',i2,'//                                           &
                         ''' Centroid = '',f16.6,'' MeV'','//                         &
                          ''' Width = '',f16.6,'' MeV'','//                           &
                          ''' Strength = '',e15.7,'' mb'')')                          &
                   k,                                                                 &
                   nucleus(i)%EL_mode(l_radiation)%gsf(k)%er,                         &
                   nucleus(i)%EL_mode(l_radiation)%gsf(k)%gr,                         &
                   nucleus(i)%EL_mode(l_radiation)%gsf(k)%sr
            end do
         end if

         write(temp_string,*)nucleus(i)%EL_mode(l_radiation)%num_gsf

         f_units(1:5) = 'MeV^-'
         if(2*l_radiation + 1 < 10)then
            write(f_units(6:6),'(i1)')2*l_radiation + 1
            f_units(7:7) = ' '
         else
            write(f_units(6:7),'(i2)')2*l_radiation + 1
         end if

         write(13,'(''E'',i1,'' Strength function, F, in units of '',a7)')l_radiation, f_units
         if(nucleus(i)%EL_mode(l_radiation)%num_gsf == 1)then
            write(13,'(''     Energy         F               T         absorption(b)'')')
            write(13,'(''  ---------   -------------   -------------   -------------'')')
         else
            fstring = "(5x,'Energy',"//trim(adjustl(temp_string))//                &
                      "(7x,'F(',i1,')',5x),9x,'F',15x,'T',9x,'absorption(b)')"
            write(13,fstring)(k,k=1,nucleus(i)%EL_mode(l_radiation)%num_gsf)
            fstring = "(2x,'---------',"//trim(adjustl(temp_string))//             &
                      "(3x,'-------------'),3(3x,'-------------'))"
            write(13,fstring)

         end if
         fstring = "(1x,f10.5,"//trim(adjustl(temp_string))//"(1x,e15.7),3(1x,e15.7))"
         do j = 0, nucleus(i)%nbin_em
            energy = real(j,kind=8)*de
            gsf(1:max_num_gsf) = 0.0d0
            do k = 1, nucleus(i)%EL_mode(l_radiation)%num_gsf
               gsf(k) = EL_f_component(i, l_radiation, k, energy, nucleus(i)%sep_e(1))
            end do
            f_E = EL_f(i, l_radiation, energy, nucleus(i)%sep_e(1))
            trans_E = EL_trans(i, l_radiation, energy, nucleus(i)%sep_e(1))
            EL_cs = EL_absorption(i, l_radiation, energy, nucleus(i)%sep_e(1))
            if(nucleus(i)%EL_mode(l_radiation)%num_gsf == 1)then
               write(13,'(1x,f10.5,(3(1x,e15.7)))')                                &
                    energy, f_E, trans_E, EL_cs
            else
               write(13,fstring)                                                   &
                    energy, (gsf(k),k=1,nucleus(i)%EL_mode(l_radiation)%num_gsf),  &
                    f_E, trans_E, EL_cs
            end if
         end do
      end do
!------------------------------------------------------------------
!------------------------------------------------------------------
!-----------   Magnetic multipoles -
!------------------------------------------------------------------
      do l_radiation = 1, nucleus(i)%lmax_M
         write(13,*)
         if(nucleus(i)%ML_mode(l_radiation)%gsf_read)then
            write(13,*)'Gamma Strength Function obtained from manual input'
         else
            write(13,*)'Gamma Strength Function obtained with internal analytic models'
            write(13,'(''Magnetic-'',i1,'' resonance parameters for '',a5)')          &
               l_radiation,nuke_label(1:inuke_end)
            do k = 1, nucleus(i)%ML_mode(l_radiation)%num_gsf
               write(13,'(''Line ='',i2,'//                                           &
                         ''' Centroid = '',f16.6,'' MeV'','//                         &
                          ''' Width = '',f16.6,'' MeV'','//                           &
                          ''' Strength = '',e15.7,'' mb'')')                          &
                   k,                                                                 &
                   nucleus(i)%ML_mode(l_radiation)%gsf(k)%er,                         &
                   nucleus(i)%ML_mode(l_radiation)%gsf(k)%gr,                         &
                   nucleus(i)%ML_mode(l_radiation)%gsf(k)%sr
            end do
         end if

         write(temp_string,*)nucleus(i)%ML_mode(l_radiation)%num_gsf

         f_units(1:5) = 'MeV^-'
         if(2*l_radiation + 1 < 10)then
            write(f_units(6:6),'(i1)')2*l_radiation + 1
            f_units(7:7) = ' '
         else
            write(f_units(6:7),'(i2)')2*l_radiation + 1
         end if

         write(13,'(''M'',i1,'' Strength function, F, in units of '',a7)')l_radiation, f_units
         if(nucleus(i)%ML_mode(l_radiation)%num_gsf == 1)then
            write(13,'(''     Energy         F               T         absorption(b)'')')
            write(13,'(''  ---------   -------------   -------------   -------------'')')
         else
            fstring = "(5x,'Energy',"//trim(adjustl(temp_string))//                &
                      "(7x,'F(',i1,')',5x),9x,'F',15x,'T',9x,'absorption(b)')"
            write(13,fstring)(k,k=1,nucleus(i)%ML_mode(1)%num_gsf)
            fstring = "(2x,'---------',"//trim(adjustl(temp_string))//             &
                      "(3x,'-------------'),3(3x,'-------------'))"
            write(13,fstring)

         end if
         fstring = "(1x,f10.5,"//trim(adjustl(temp_string))//"(1x,e15.7),3(1x,e15.7))"
         do j = 0, nucleus(i)%nbin_em
            energy = real(j,kind=8)*de
            gsf(1:max_num_gsf) = 0.0d0
            do k = 1, nucleus(i)%ML_mode(l_radiation)%num_gsf
               gsf(k) = ML_f_component(i, l_radiation, k, energy)
            end do
            f_M = ML_f(i, l_radiation, energy)
            trans_M = ML_trans(i, l_radiation, energy)
            ML_cs = ML_absorption(i, l_radiation, energy)
            if(nucleus(i)%ML_mode(l_radiation)%num_gsf == 1)then
               write(13,'(1x,f10.5,(3(1x,e15.7)))')                                &
                    energy, f_M, trans_M, ML_cs
            else
               write(13,fstring)                                                   &
                    energy,(gsf(k),k=1,nucleus(i)%ML_mode(l_radiation)%num_gsf),   &
                    f_M, trans_M, ML_cs
            end if
         end do
      end do
!------   Finished with EM strength functions
      if(allocated(gsf))deallocate(gsf)


      if(fission .and. .not. nucleus(i)%fission)then
         write(13,'(''+++++++++++++++++++++++++++++++++++++++++++++++++++++++'')')
         write(13,'(''This nucleus is not set up to undergo fission          '')')
         write(13,'(''This is most likely due to no initial fission barrier  '')')
         write(13,'(''data in file Fission-barrier.dat. To override this     '')')
         write(13,'(''use option f_num_barrier and manually set fission      '')')
         write(13,'(''barrier information.                                   '')')
      end if
      if(nucleus(i)%Fission)then
         write(13,*)
         write(13,'(''Fission parameters for '',a5)')nuke_label(1:inuke_end)
         do j = 1, nucleus(i)%F_n_barr
            E0 = nucleus(i)%F_barrier(j)%level_param(15)
            T = nucleus(i)%F_barrier(j)%level_param(14)
            E1 = -10.0d0
            if(E0 < 0.0d0)E1 = T*log(1.0d0-exp(E0/T))
            write(13,'(''Parameters for Fission Barrier #'',i3)')j
            write(13,'(''Height = '',f10.4)')nucleus(i)%F_Barrier(j)%barrier
            write(13,'(''Width  = '',f10.4)')nucleus(i)%F_Barrier(j)%hbw
            write(13,'(''Deformation parameter beta(2) = '',f10.6)')nucleus(i)%F_Barrier(j)%beta_2
            write(13,'(''Excitation energy dependent barrier, F = F*x(1)*exp(-((Ex-x(2))/x(3))**2)'')')
            write(13,'(''x(1) = '',f10.4,'' x(2) = '',f10.4,'' x(3) = '',e15.7)')            &
                  nucleus(i)%F_Barrier(j)%barrier_damp(1),                                   &
                  nucleus(i)%F_Barrier(j)%barrier_damp(2),                                   &
                  nucleus(i)%F_Barrier(j)%barrier_damp(3)
            write(13,*)
            write(13,'(''Level-density parameters above the barrier'')')
            write(13,*)'Level-density parameters'
            write(13,'('' a =        '',f12.7)')nucleus(i)%F_barrier(j)%level_param(1)
            write(13,'('' lambda =   '',f12.7)')nucleus(i)%F_barrier(j)%level_param(2)
            write(13,'('' Delta =    '',f12.7)')nucleus(i)%F_barrier(j)%level_param(3)
            write(13,'('' Shell =    '',f12.7)')nucleus(i)%F_barrier(j)%level_param(4)
            write(13,'('' gamma =    '',f12.7)')nucleus(i)%F_barrier(j)%level_param(5)
            write(13,'('' E_match =  '',f12.7)')nucleus(i)%F_barrier(j)%level_param(6)
            write(13,'('' E_cut =    '',f12.7)')nucleus(i)%F_barrier(j)%level_param(7)
            write(13,'('' T =        '',f12.7)')nucleus(i)%F_barrier(j)%level_param(14)
            write(13,'('' E0 =       '',f12.7)')nucleus(i)%F_barrier(j)%level_param(15)
            if(E1 < 0.0d0)then
               write(13,'('' E1 =       '',f12.7)')E1
            else
               write(13,'(''E1 is undefined as E0 > 0, E1 ~ -5.0'')')
            end if
            write(13,*)


            if(nint(nucleus(i)%F_Barrier(j)%level_param(10)) == 0)then
               write(13,*)'No rotational collective enhancement included'
               write(13,*)'K_rot = 1.0'
            else
               write(13,*)'Rotational Collective enhancement model = ',    &
                  nucleus(i)%F_Barrier(j)%level_param(10)
                write(13,*)'Collective enhancement factors'
                write(13,*)'K_rot = Max[x(1)*(Factor-1)/(1+exp((E-x(2))/x(3))),0]+1'
                if(nint(nucleus(i)%F_Barrier(j)%level_param(10)) == 1)then
                   if(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 2 .or.            &
                      nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 3)then
       write(13,*)'Factor = x(1)*sig2*(1+beta2/3)'
       write(13,*)' --  Axially symmetric'
                   elseif(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 4 .or.            &
                      nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 5)then
       write(13,*)'Factor = x(1)*sig2*(1+sqrt(5/4pi)*beta2)'
       write(13,*)' --  Axially symmetric'
                   end if
                elseif(nint(nucleus(i)%F_Barrier(j)%level_param(10)) == 2)then
                   if(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 2 .or.            &
                      nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 3)then
       write(13,*)'Factor = 2.0*x(1)*sig2*(1+beta2/3)'
       write(13,*)' -- Left-right asymmetric'
                   elseif(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 4 .or.            &
                      nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 5)then
       write(13,*)'Factor = 2.0*x(1)*sig2*(1+sqrt(5/4pi)*beta2)'
       write(13,*)' -- Left-right asymmetric'
                   end if
                elseif(nint(nucleus(i)%F_Barrier(j)%level_param(10)) == 3)then
                   if(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 2 .or.            &
                      nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 3)then
       write(13,*)'Factor = x(1)*sqrt(pi/2)**sig2**(3/2)*(1+beta2/3)*(1-2*beta2/3)'
       write(13,*)' -- Triaxial and left-right asymmetric'
                   elseif(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 4 .or.            &
                      nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 5)then
       write(13,*)'Factor = x(1)*sqrt(pi/2)**sig2**(3/2)*(1+sqrt(5/4pi)*beta2)*(1-sqrt(5/4pi)beta2/2)'
       write(13,*)' -- Triaxial and left-right asymmetric'
                   end if
                elseif(nint(nucleus(i)%F_Barrier(j)%level_param(10)) == 4)then
                   if(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 2 .or.            &
                      nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 3)then
       write(13,*)'Factor = 2.0*x(1)*sqrt(pi/2)**sig2**(3/2)*(1+beta2/3)*(1-2*beta2/3)'
       write(13,*)' -- Triaxial and not left-right asymmetric'
                   elseif(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 4 .or.            &
                      nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 5)then
       write(13,*)'Factor = 2.0*x(1)*sqrt(pi/2)**sig2**(3/2)*(1+sqrt(5/4pi)*beta2)*(1-sqrt(5/4pi)beta2/2)'
       write(13,*)' -- Triaxial and not left-right asymmetric'
                   end if
                end if

                if(nint(nucleus(i)%F_Barrier(j)%level_param(19)) == 3)then
                   write(13,*)'Deformation-dependent spin cutoff parameter is used to compute the spin distribution'
                   write(13,*)'sig2*(1+beta2/3)'
                end if


                write(13,'('' Rotational collective-enhancement parameters'')')

                write(13,'(''x(1) = '',f10.4)')nucleus(i)%F_Barrier(j)%rot_enh(1)
                write(13,'(''x(2) = '',f10.4)')nucleus(i)%F_Barrier(j)%rot_enh(2)
                write(13,'(''x(3) = '',f10.4)')nucleus(i)%F_Barrier(j)%rot_enh(3)
                write(13,*)
             end if
            if(nint(nucleus(i)%F_Barrier(j)%level_param(11)) == 0)then
               write(13,*)'No vibrational collective enhancement included'
               write(13,*)'K_vib = 1.0'
            elseif(nint(nucleus(i)%F_Barrier(j)%level_param(11)) == 1)then
               write(13,*)'Vibrational collective enhancement model 1 used'
               write(13,*)'K_vib = exp(0.05555*A**(2./3.)*T**(4./3.))'
               write(13,*)'K_vib = Max[x(1)*(K_vib-1)/(1+exp((E-x(2))/x(3))),0]+1'
               write(13,*)'T = sqrt(U/apu)'
               write(13,'(''x(1) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(1)
               write(13,'(''x(2) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(2)
               write(13,'(''x(3) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(3)
               write(13,*)
            elseif(nint(nucleus(i)%F_Barrier(j)%level_param(11)) == 2)then
               write(13,*)'K_vib model 2 - approximated with Bose-gas relationship'
               write(13,*)
               write(13,*)'K_vib = exp(dS - dU/T)'
               write(13,*)'T = sqrt(U/apu)'
               write(13,*)
               write(13,*)'dS = Sum_i (2*L(i)+1)*[(1+n(i))*log(1+n(i)) - n(i)*log(n(i))]'
               write(13,*)'dU = Sum_i (2*L(i)+1)*w(i)*n(i)'
               write(13,*)'n(i) = exp(-g(i)/2*w(i))/[exp(w(i)/T) - 1]'
               write(13,*)'g(i) = C*(w(i)**2 + 4*pi*T**2)'
               write(13,*)
               write(13,*)'C = 0.0075*A**(1/3)'
               write(13,*)'w(2) = 65*A**(-5/6)/(1*0.05*Shell)'
               write(13,*)'w(3) = 100*A**(-5/6)/(1*0.05*Shell)'
               write(13,*)
               write(13,*)'The sum extends over quadrupole and octupole vibrations, L=2,3'
            elseif(nint(nucleus(i)%F_Barrier(j)%level_param(11)) == 3)then
               write(13,*)'Vibrational collective enhancement model 1 used'
               write(13,*)'K_vib = exp(0.05555*A**(2./3.)*T**(4./3.))'
               write(13,*)'K_vib = Max[x(1)*(K_vib-1)/(1+exp((E-x(2))/x(3))),0]+1'
               write(13,*)'T = sqrt(U/a);  a = A/13'
               write(13,'(''x(1) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(1)
               write(13,'(''x(2) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(2)
               write(13,'(''x(3) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(3)
               write(13,*)
            elseif(nint(nucleus(i)%F_Barrier(j)%level_param(11)) == 4)then
               write(13,*)'K_vib model 2 - approximated with Bose-gas relationship'
               write(13,*)
               write(13,*)'K_vib = exp(dS - dU/T)'
               write(13,*)'T = sqrt(U/a);  a = A/13'
               write(13,*)
               write(13,*)'dS = Sum_i (2*L(i)+1)*[(1+n(i))*log(1+n(i)) - n(i)*log(n(i))]'
               write(13,*)'dU = Sum_i (2*L(i)+1)*w(i)*n(i)'
               write(13,*)'n(i) = exp(-g(i)/2*w(i))/[exp(w(i)/T) - 1]'
               write(13,*)'g(i) = C*(w(i)**2 + 4*pi*T**2)'
               write(13,*)
               write(13,*)'C = 0.0075*A**(1/3)'
               write(13,*)'w(2) = 65*A**(-5/6)/(1*0.05*Shell)'
               write(13,*)'w(3) = 100*A**(-5/6)/(1*0.05*Shell)'
               write(13,*)
               write(13,*)'The sum extends over quadrupole and octupole vibrations, L=2,3'
            end if
            write(13,*)
            write(13,*)'Rho = K_rot*K_vib*Rho'

            write(13,*)
            if(nint(nucleus(i)%F_Barrier(j)%level_param(16)) == 0)then
               write(13,*)'Parity Factor = 0.5 for positive and negative parities'
            elseif(nint(nucleus(i)%F_Barrier(j)%level_param(16)) == 1)then
               write(13,*)'Parity Factor = 0.5*tanh(a*(E-E0)) for positive parity states'
               write(13,'(''E0 = '',f5.3)')nucleus(i)%F_Barrier(j)%level_param(17)
               write(13,'(''a = '',f5.3)')nucleus(i)%F_Barrier(j)%level_param(18)
            elseif(nint(nucleus(i)%F_Barrier(j)%level_param(16)) == 2)then
               write(13,*)'Parity Factor = 0.5*tanh(a*(E-E0)) for negative parity states'
               write(13,'(''E0 = '',f5.3)')nucleus(i)%F_Barrier(j)%level_param(17)
               write(13,'(''a = '',f5.3)')nucleus(i)%F_Barrier(j)%level_param(18)
            end if
            write(13,*)

            write(13,*)
            write(13,'(''Number of discrete transition states = '', i5)')            &
                nucleus(i)%F_barrier(j)%num_discrete
            if(nucleus(i)%F_barrier(j)%num_discrete > 0)then
               write(13,'(''Transition states above barrier #'',i4)')j
               do k = 1, nucleus(i)%F_barrier(j)%num_discrete
                  char = '+'
                  if(nucleus(i)%F_barrier(j)%state_pi(k) < 0)char = '-'
                  write(13,'(i5,1x,f10.3,1x,f4.1,a1)')                               &
                     k,nucleus(i)%F_barrier(j)%state_e(k),                           &
                       nucleus(i)%F_barrier(j)%state_j(k),char
               end do
            end if
!-------   Fission level densities

            dde = 0.2d0
            if(nucleus(i)%nbin > 1)dde = nucleus(i)%e_grid(2) - nucleus(i)%e_grid(1)

            n_bin = int((20.0d0 - nucleus(i)%F_barrier(j)%ecut)/dde)



            pmode = nucleus(i)%F_barrier(j)%level_param(16)
            pe1 = nucleus(i)%F_barrier(j)%level_param(17)
            pbb = nucleus(i)%F_barrier(j)%level_param(18)

            ip = 0

            write(13,'(''Level density States/MeV'')')
            write(13,'(''Positive parity'')')
            write(temp_string,*)min(j_max,60)+1
            fstring = "(6x,'E',9x,'a(U)',1x,'sqrt(sig2)',6x,'P_fac',"//             &
                     "6x,'K_vib',6x,'K_rot',8x,'Enh',9x,'Rho',5x,"                  &
                     //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

            write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
            fstring = "(1x,'--------',6(1x,'----------'),1x,'---------------',"        &
                      //trim(adjustl(temp_string))//"(1x,'---------------'))"
            write(13,fstring)
            fstring = "(1x,f8.3,6(1x,f10.3),(1x,e15.7),"//trim(adjustl(temp_string))//"(1x,e15.7))"
            ip = 1
            do k = 1, n_bin
               energy = real(k,kind=8)*dde
               if(nucleus(i)%F_barrier(j)%num_discrete > 0 .and.         &
                  energy < nucleus(i)%F_barrier(j)%ecut) cycle
               do jj = 0, min(j_max,60)
                  xj = jj + nucleus(i)%jshift
                  call rho_J_par_e(energy, xj, ip,                                   &
                                   nucleus(i)%F_barrier(j)%level_param,              &
                                   nucleus(i)%F_barrier(j)%vib_enh,                  &
                                   nucleus(i)%F_barrier(j)%rot_enh,                  &
                                   ia,rho_Fm,apu,sig2,K_vib,K_rot)
                  rho(jj,ip) = rho_FM
               end do
               call rhoe(energy,nucleus(i)%F_barrier(j)%level_param,              &
                                nucleus(i)%F_barrier(j)%vib_enh,                  &
                                nucleus(i)%F_barrier(j)%rot_enh,                  &
                                ia,rho_Fm,apu,sig2,K_vib,K_rot)

!   write(6,*)energy,sig2_param(energy,nucleus(i)%F_barrier(j)%level_param,nucleus(i)%A)

!      sig = nucleus(i)%F_barrier(j)%level_param(12)
!      shell = nucleus(i)%F_barrier(j)%level_param(4)
!      gamma = nucleus(i)%F_barrier(j)%level_param(5)
!      sig_model = nint(nucleus(i)%F_barrier(j)%level_param(13))
!      ematch = nucleus(i)%F_barrier(j)%level_param(6)
!      delta = nucleus(i)%F_barrier(j)%level_param(3)
!      Um = Ematch - delta
!      aparam = nucleus(i)%F_barrier(j)%level_param(1)
!      sg2cut = nucleus(i)%F_barrier(j)%level_param(8)
!      ecut = nucleus(i)%F_barrier(j)%level_param(7)
!      apu = aparam_u(Um,aparam,shell,gamma)
!      sig2_em = sig*sqrt(max(0.2d0,Um*apu))/aparam
!      if(sig_model == 0)sig2_em = sig*sqrt(max(0.2d0,Um*apu))/aparam
!      if(sig_model == 1)sig2_em = sig*sqrt(max(0.2d0,Um)/apu)
!      sig2_em = max(sig2_em,sig2_min)
!      deriv = (sig2_em - sg2cut)/(ematch - ecut)
!      sig2 = sig2_em - deriv*(ematch - energy)
!      if(sig2 < sig2_min) sig2 = sig2_min

!   write(6,*)energy,Ematch, delta, Um, ecut
!   write(6,*)aparam, apu
!   write(6,*)sig2_em, sig2_min, sg2cut
!   write(6,*)deriv

!   write(6,*)'sig2 ',sig2

               pfac = parity_fac(energy,xj,ip,pmode,pe1,pbb)
               write(13,fstring)                                                     &
                   energy,apu,sqrt(sig2),pfac,K_vib,K_rot,K_vib*K_rot,rho_FM*pfac,   &
                   (rho(jj,ip),jj = 0, min(j_max,60))
            end do
            write(13,*)
!  stop
            write(13,'(''Level density States/MeV'')')
            write(13,'(''Negative parity'')')
            write(temp_string,*)min(j_max,60)+1
            fstring = "(6x,'E',9x,'a(U)',1x,'sqrt(sig2)',6x,'P_fac',"//             &
                     "6x,'K_vib',6x,'K_rot',8x,'Enh',9x,'Rho',5x,"                  &
                     //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

            write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
            fstring = "(1x,'--------',6(1x,'----------'),1x,'---------------',"        &
                      //trim(adjustl(temp_string))//"(1x,'---------------'))"
            write(13,fstring)
            fstring = "(1x,f8.3,6(1x,f10.3),(1x,e15.7),"//trim(adjustl(temp_string))//"(1x,e15.7))"
            ip = 0
            do k = 1, n_bin
               energy = real(k,kind=8)*dde
               if(nucleus(i)%F_barrier(j)%num_discrete > 0 .and.         &
                  energy < nucleus(i)%F_barrier(j)%ecut) cycle
               sum_rho = 0.0d0
               call rhoe(energy,nucleus(i)%F_barrier(j)%level_param,                 &
                         nucleus(i)%F_barrier(j)%vib_enh,                            &
                         nucleus(i)%F_barrier(j)%rot_enh,                            &
                         ia,rho_Fm,apu,sig2,K_vib,K_rot)
               pfac=parity_fac(energy,xj,ip,pmode,pe1,pbb)
               do jj = 0, min(j_max,60)
                  xj = jj + nucleus(i)%jshift
                  call rho_J_par_e(energy,xj, ip,                                    &
                                   nucleus(i)%F_barrier(j)%level_param,              &
                                   nucleus(i)%F_barrier(j)%vib_enh,                  &
                                   nucleus(i)%F_barrier(j)%rot_enh,                  &
                                   ia,rho_J_ip,apu,sig2,K_vib,K_rot)
                  rho(jj,ip) = rho_J_ip
               end do
               write(13,fstring)                                                     &
                   energy,apu,sqrt(sig2),pfac,K_vib,K_rot,K_vib*K_rot,rho_FM*pfac,   &
                   (rho(jj,ip),jj = 0, min(j_max,60))
            end do


         end do
         write(13,*)
      end if


!--------------   HF-denominators
      write(13,*)
      write(13,'(''Hauser-Feshbach denominators for '',a5)')nuke_label(1:inuke_end)
      write(13,*)'Positive Parity'
      write(temp_string,*)min(j_max,60)+1
      fstring = "('    E     ',"                      &
                //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

      write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
      fstring = "('----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
      write(13,fstring)
      fstring = '(f10.6,'//trim(adjustl(temp_string))//'(1x,e15.7))'
      do n = 1, nucleus(i)%nbin
         energy=nucleus(i)%e_grid(n)
         write(13,fstring)                                                           &
           energy,(nucleus(i)%bins(j,0,n)%HF_den,j=0,min(j_max,60))
      end do
      write(13,*)
      write(13,*)'Negative Parity'
      fstring = "('    E     ',"                      &
                //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

      write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
      fstring = "('----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
      write(13,fstring)
      fstring = '(f10.6,'//trim(adjustl(temp_string))//'(1x,e15.7))'
      do n = 1, nucleus(i)%nbin
         energy=nucleus(i)%e_grid(n)
         write(13,fstring)                                                           &
           energy,(nucleus(i)%bins(j,1,n)%HF_den,j=0,min(j_max,60))
      end do
      write(13,*)
      write(13,'(''Summed decay transmission coefficients for each channel for '',a5)')nuke_label(1:inuke_end)
      do if1 = 1, nucleus(i)%num_decay
         k = nucleus(i)%decay_particle(if1)
         write(13,*)
         write(13,*)'Particle type = ',k
         write(13,*)'Positive Parity'
         fstring = "('    E     ',"                      &
                   //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

         write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
         fstring = "('----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
         write(13,fstring)
         fstring = '(f10.6,'//trim(adjustl(temp_string))//'(1x,e15.7))'
         do n=1,nucleus(i)%nbin
            energy=nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob1(j) = nucleus(i)%bins(j,0,n)%HF_trans(if1)
            end do
            write(13,fstring)energy,(prob1(j),j=0,min(j_max,60))
         end do
         write(13,*)
         write(13,*)'Negative Parity'
         fstring = "('    E     ',"                      &
                   //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

         write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
         fstring = "('----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
         write(13,fstring)
         fstring = '(f10.6,'//trim(adjustl(temp_string))//'(1x,e15.7))'
         do n = 1, nucleus(i)%nbin
            energy=nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob1(j) = nucleus(i)%bins(j,1,n)%HF_trans(if1)
            end do
            write(13,fstring)energy,(prob1(j),j=0,min(j_max,60))
         end do
      end do

      if(nucleus(i)%Fission)then
         if1 = nucleus(i)%num_decay + 1
         write(13,*)
         write(13,*)'Positive Parity: Decay type = Fission'
         fstring = "('    E     ',"                      &
                   //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

         write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
         fstring = "('----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
         write(13,fstring)
         fstring = '(f10.6,'//trim(adjustl(temp_string))//'(1x,e15.7))'
         do n = 1, nucleus(i)%nbin
            energy = nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob1(j) = nucleus(i)%bins(j,0,n)%HF_trans(if1)
            end do
            write(13,fstring)energy,(prob1(j),j = 0, min(j_max,60))
         end do
         write(13,*)
         write(13,*)'Negative Parity: Decay type = Fission'
         fstring = "('    E     ',"                      &
                   //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

         write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
         fstring = "('----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
         write(13,fstring)
         fstring = '(f10.6,'//trim(adjustl(temp_string))//'(1x,e15.7))'
         do n = 1, nucleus(i)%nbin
            energy = nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob1(j) = nucleus(i)%bins(j,1,n)%HF_trans(if1)
            end do
            write(13,fstring)energy,(prob1(j),j = 0, min(j_max,60))
         end do

         do ii = 1, nucleus(i)%F_n_barr
            write(13,*)
            write(13,*)'Positive Parity: Trans-Coef for Fission Barrier #', ii
            fstring = "('    E     ',"                      &
                      //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

            write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
            fstring = "('----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
            write(13,fstring)
            fstring = '(f10.6,'//trim(adjustl(temp_string))//'(1x,e15.7))'
            do n = 1, nucleus(i)%nbin
               energy = nucleus(i)%e_grid(n)
               do j = 0, min(j_max,60)
                  prob1(j) = nucleus(i)%bins(j,0,n)%HF_trans(if1+ii)
               end do
               write(13,fstring)energy,(prob1(j),j = 0, min(j_max,60))
            end do
            write(13,*)
            write(13,*)'Negative Parity: Trans-Coef for Fission Barrier #', ii
            fstring = "('    E     ',"                      &
                      //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

            write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
            fstring = "('----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
            write(13,fstring)
            fstring = '(f10.6,'//trim(adjustl(temp_string))//'(1x,e15.7))'
            do n = 1, nucleus(i)%nbin
               energy = nucleus(i)%e_grid(n)
               do j = 0, min(j_max,60)
                  prob1(j) = nucleus(i)%bins(j,1,n)%HF_trans(if1+ii)
               end do
               write(13,fstring)energy,(prob1(j),j = 0, min(j_max,60))
            end do
         end do
      end if
      write(13,*)
      write(13,'(''Decay probabilities for each channel in decay from '',a5)')nuke_label(1:inuke_end)
      nnn = nucleus(i)%num_decay
      if(nucleus(i)%fission)nnn = nnn + 1
      do if1 = 1, nnn
         if(if1 <= nucleus(i)%num_decay)then
            k = nucleus(i)%decay_particle(if1)
         else            !   Fission
            k = 7
         end if
         write(13,*)
         if(k < 7)write(13,*)'Particle type = ',particle(k)%name
         if(k == 7)write(13,*)'Decay type = Fission'
         write(13,*)'Positive Parity'
         fstring = "('    E     ',"                                            &
                   //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

         write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
         fstring = "('----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
         write(13,fstring)
         fstring = '(f10.6,'//trim(adjustl(temp_string))//'(1x,e15.7))'
         do n = 1,nucleus(i)%nbin
            energy = nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob1(j) = 0.0d0
               do ifi = 1 , nucleus(i)%bins(j,0,n)%num_decay
                  if(nucleus(i)%bins(j,0,n)%decay_particle(ifi) == k)  &
                     prob1(j) = nucleus(i)%bins(j,0,n)%HF_prob(ifi)
               end do
            end do
            write(13,fstring)energy,(prob1(j),j = 0, min(j_max,60))
         end do
         write(13,*)'Negative Parity'
         fstring = "('    E     ',"                                             &
                   //trim(adjustl(temp_string))//"(5x,'J = ',f4.1,3x))"

         write(13,fstring)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
         fstring = "('----------',"//trim(adjustl(temp_string))//"(1x,'---------------'))"
         write(13,fstring)
         fstring = '(f10.6,'//trim(adjustl(temp_string))//'(1x,e15.7))'
         do n = 1,nucleus(i)%nbin
            energy = nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob1(j) = 0.0d0
               do ifi = 1 , nucleus(i)%bins(j,1,n)%num_decay
                  if(nucleus(i)%bins(j,1,n)%decay_particle(ifi) == k)  &
                     prob1(j) = nucleus(i)%bins(j,1,n)%HF_prob(ifi)
               end do
            end do
            write(13,fstring)energy,(prob1(j),j = 0, min(j_max,60))
         end do
      end do
   end do
   if(pop_calc)then
      write(13,*)'*****************************************************************'
      write(13,*)
      write(13,'(''Decay probabilities for each initial excitation Population for '')')
      write(13,'(3x,6x,''Ex'',8x,8(3x,8x,a8))')(particle(k)%name, k = 0, 7)
      write(13,'(''   ----------------'',8(''   ----------------''))')
!----   loop over initial population energies
      do i = 1, num_pop_e
         prob_jpi(0:60,0:1) = 0.0d0
         Ex = Pop_data(i)%Ex_pop
         if(.not. j_pop_calc)then
            do j = 1, Pop_data(i)%num_pop
               Ix = nint(Pop_data(i)%j_pop(j) - nucleus(1)%jshift)
               ip = nint(Pop_data(i)%par_pop(j) + 1.0d0)/2
               prob_jpi(Ix,ip) = Pop_data(i)%bin_pop(j)
            end do
         end if
         n = find_ibin(Ex,1)
         part_prob(0:7) = 0.0d0
         do if1 = 1, nucleus(1)%num_decay
            k = nucleus(1)%decay_particle(if1)
            do ip = 0, 1
               do j = 0, min(j_max,60)
                  do ifi = 1, nucleus(1)%bins(j,ip,n)%num_decay
                     if(nucleus(1)%bins(j,ip,n)%decay_particle(ifi) == k)then
                        part_prob(k) = part_prob(k) + nucleus(1)%bins(j,ip,n)%HF_prob(ifi)*  &
                                       prob_jpi(j,ip)
                     end if
                  end do
               end do
            end do
         end do
         write(13,'(9(3x,1pe16.7))')Ex,(part_prob(k),k=0,7)
      end do
      write(13,*)
      write(13,*)'*****************************************************************'
   end if
   flush(13)
   return
end subroutine output_nucleus_data

