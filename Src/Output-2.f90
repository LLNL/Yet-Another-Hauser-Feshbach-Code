!
!*******************************************************************************
!
subroutine start_IO(num_comp)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine prints out data defining the Hauser-Feshbach calculation,
!    basically printing all the otpions, etc. onto the main output file
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
   use nodeinfo
   use options
   use pre_equilibrium_no_1
   use print_control
   use useful_data
   use nuclei
   use Channel_info
   use particles_def
   use directory_structure
   use constants
   implicit none
!-------------------------------------------------------
   integer(kind=4), intent(in) :: num_comp
!-------------------------------------------------------
   integer(kind=4) :: i, ip
   integer(kind=4) char_start
   integer(kind=4) :: ilast

   ilast = index(version,' ')-1

   write(13,*)'**********************************************'
   write(13,*)'* Livermore Monte Carlo Hauser-Feshbach code'
   write(13,*)'* Version # ',version(1:ilast)
   write(13,*)'**********************************************'
   write(13,*)
   write(13,*)'Complete list of compound nuclei'
   write(13,*)'Complete list of compound nuclei'
   do i = 1, num_comp
      char_start = 1
      if(nucleus(i)%atomic_symbol(1:1) == ' ')char_start = 2
      if(char_start == 1 .and. print_me)write(6,'(''Nucleus # '',i3,2x,i3,a2)')     &
         i,nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
      if(char_start == 2 .and. print_me)write(6,'(''Nucleus # '',i3,2x,i3,a1)')     &
         i,nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
      if(print_me)write(6,*)'Maximum excitation energy = ',               &
         nucleus(i)%Ex_max
      if(char_start == 1)write(13,'(''Nucleus # '',i3,2x,i3,a2)')         &
         i,nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
      if(char_start == 2)write(13,'(''Nucleus # '',i3,2x,i3,a1)')         &
         i,nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
      write(13,*)'Maximum excitation energy = ',                          &
         nucleus(i)%Ex_max
      write(13,*)'Energy grid set up with de = ',de
   end do
   write(6,*)
   if(.not.pop_calc)then
      ip = projectile%particle_type
      if(print_me)then
         write(6,'(''Projectile = '',a10)')particle(ip)%name
         write(6,'(''Min and max incidient energies: Min ='',f16.6,'//       &
                   ''' MeV, Max = '',f16.6,'' MeV'')')                       &
                projectile%e_min,projectile%e_max
      end if
      write(13,'(''Projectile = '',a10)')particle(ip)%name
      write(13,'(''Min and max incidient energies: Min ='',f16.6,'//      &
                ''' MeV, Max = '',f16.6,'' MeV'')')                       &
                projectile%e_min,projectile%e_max
   end if
   i = target%icomp
   char_start = 1
   if(nucleus(i)%atomic_symbol(1:1) == ' ')char_start = 2
   if(.not.pop_calc)then
      if(print_me)then
         if(char_start == 1)write(6,'(''Target nucleus = '',i3,a2)')                    &
                 nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
         if(char_start == 2)write(6,'(''Target nucleus = '',i3,a1)')                    &
                 nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
      end if
      if(char_start == 1)write(13,'(''Target nucleus = '',i3,a2)')                      &
              nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
      if(char_start == 2)write(13,'(''Target nucleus = '',i3,a1)')                      &
              nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
   else
      if(print_me)then
         if(char_start == 1)write(6,'(''Initial nucleus = '',i3,a2)')                   &
                 nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
         if(char_start == 2)write(6,'(''Initial nucleus = '',i3,a1)')                      &
                 nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
      end if
      if(char_start == 1)write(13,'(''Initial nucleus = '',i3,a2)')                  &
                 nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
      if(char_start == 2)write(13,'(''Initial nucleus = '',i3,a1)')                     &
              nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
   end if
   write(6,'(''Tracking the following emitted particles'')')
   write(13,'(''Tracking the following emitted particles'')')
   do ip = 0, 6
      if(particle(ip)%in_decay)then
         if(print_me)write(6,'(a10)')particle(ip)%name
         write(13,'(a10)')particle(ip)%name
      end if
   end do
   write(13,*)
   write(13,*)'*******************************************************'
   write(13,*)
   write(13,*)'Calculation options set in this calculation:'
   write(13,*)
   write(13,*)'output_mode = ',output_mode
   write(13,*)'fission = ',fission
   write(13,*)'channels = ',channels
   write(13,*)'dump_events = ',dump_events
   write(13,*)'track_gammas = ',track_gammas
   write(13,*)'track_primary_gammas = ',track_primary_gammas
   write(13,*)'PREEQ_Model = ',PREEQ_Model
   write(13,*)'analytic_preeq = ',analytic_preeq
   write(13,*)'Preeq_V = ', Preeq_V
   write(13,*)'Preeq_V1:',Preeq_V1
   write(13,*)'Preeq_K:',Preeq_K
   write(13,*)'preeq_pair_model = ',preeq_pair_model
   write(13,*)'preeq_fwell = ',preeq_fwell
   write(13,*)'Preeq_g_a = ',Preeq_g_a
   write(13,*)'Preeq_g_div = ',preeq_g_div
   write(13,*)'M2_C1 = ',M2_C1
   write(13,*)'M2_C2 = ',M2_C2
   write(13,*)'M2_C3 = ',M2_C3
   write(13,*)'M2_Rpp = ',M2_Rpp
   write(13,*)'M2_Rnn = ',M2_Rnn
   write(13,*)'M2_Rpn = ',M2_Rpn
   write(13,*)'M2_Rnp = ',M2_Rnp
   write(13,*)'WF_model = ',WF_model
   write(13,*)'trans_p_cut = ',trans_p_cut
   write(13,*)'trans_e_cut = ',trans_e_cut
   write(13,*)'prob_cut = ',prob_cut
   write(13,*)'pop_check = ',pop_check
   write(13,*)'E1_model = ',E1_model
   write(13,*)'pop_calc = ',pop_calc
   write(13,*)'fit_Gamma_gamma = ',fit_Gamma_gamma
   write(13,*)'All_gammas = ',All_gammas
   write(13,*)'Out_gammas_vs_E = ',Out_gammas_vs_E
   write(13,*)'rho_cut = ',rho_cut
   write(13,*)'pair_model = ',pair_model
   write(13,*)'quasi_elastic ', quasi_elastic
   write(13,*)'quasi_e = ', quasi_e
   write(13,*)'Max_J_allowed = ',Max_J_allowed
   write(13,*)'num_mc_samp = ', num_mc_samp
   write(13,*)'biased_sampling = ', biased_sampling
   write(13,*)'trans_avg_l = ', trans_avg_l
   write(13,*)'explicit_channels = ', explicit_channels
   write(13,*)'num_theta_angles = ', num_theta_angles
   write(13,*)
   write(13,*)'*******************************************************'
   write(13,*)
   write(13,*)

   write(13,*)'***************************************************************'
   write(13,*)'* Built-in input parameter files, such as masses, level       *'
   write(13,*)'* density, discrete levels, fission barriers and gamma        *'
   write(13,*)'* strength functions based on the RIPL library [1]            *'   
   write(13,*)'*                                                             *'
   write(13,*)'* Levels used in this calculation are derived from the RIPL   *'
   write(13,*)'* evaluated libraries based on the ENSDF database.            *'
   write(13,*)'*                                                             *'
   write(13,*)'* [1] M. Verpelli and R. Capote Noy, INDC(NDS)-0702,          *'
   write(13,*)'*     IAEA, 2015. Available online at                         *'
   write(13,*)'* https://www-nds.iaea.org/publications/indc/indc-nds-0702/   *'
   write(13,*)'*                                                             *'
   write(13,*)'* [2] R.Capote, M.Herman, P.Oblozinsky, P.G.Young, S.Goriely, *'
   write(13,*)'*     T.Belgya, A.V.Ignatyuk, A.J.Koning, S.Hilaire,          *'
   write(13,*)'*     V.A.Plujko, M.Avrigeanu, Zhigang Ge, Yinlu Han,         *'
   write(13,*)'*     S.Kailas, J.Kopecky, V.M.Maslov, G.Reffo, M.Sin,        *'
   write(13,*)'*     E.Sh.Soukhovitskii and P. Talou, "RIPL - Reference      *'
   write(13,*)'*     Input Parameter Library for Calculation of Nuclear      *'
   write(13,*)'*     Reactions and Nuclear Data Evaluations", Nuclear Data   *'
   write(13,*)'*     Sheets 110 (2009) 3107-3214                             *'
   write(13,*)'*                                                             *'
   write(13,*)'**************************************************************'

   flush(13)
   return
end subroutine start_IO
!
!*******************************************************************************
!
subroutine output_trans_coef
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine prints out transmission coefficients on the calculation
!    grid
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
!-------------------------------------------------------------
   integer(kind=4) :: j, k, l
   integer(kind=4) :: i, isp
   real(kind=8) :: energy
   real(kind=8) :: xj, xi
!------------------------------------------------------------------
   write(13,'(''Transmission coefficients for '','//              &
              '''each emitted particle (COM energy)'')')
   do k=1,6
      if(.not.particle(k)%in_decay)cycle


      isp = nint(2.0d0*particle(k)%spin)
      xj = -particle(k)%spin
      do i = 0, isp
         xi = real(i,kind=8)
         write(13,*)
         write(13,'('' Transmission coefficients for '',a8)')       &
               particle(k)%name
         if(xj + xi < 0)then
            write(13,'(''For j = l -'',f4.1)')abs(xj+xi)
         else
            write(13,'(''For j = l +'',f4.1)')xj+xi
         end if
         write(13,'(''       T         l='',i3,5x,50(8x,i3,5x))')   &
            (l,l = 0, particle(k)%lmax)
         do j = 1, particle(k)%nbin
            energy = real(j,kind=8)*de
            write(13,'(1x,f10.5,(40(1x,e15.7)))')                   &
                  energy, (particle(k)%trans(i,l,j),l = 0, particle(k)%lmax)
      end do
      end do
   end do

   flush(13)

   return
end subroutine output_trans_coef
!
!*******************************************************************************
!
subroutine output_nucleus_data(num_comp, j_max, itarget)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine prints out data for each nucleus in the decay network
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
   integer(kind=4), intent(in) :: num_comp, j_max, itarget
!-------------------------------------------------------
   integer(kind=4) :: i, j, k, m, n, jj, ip, ii
   integer(kind=4) :: if1
   integer(kind=4) :: iA, iZ, iN
   integer(kind=4) :: l_radiation
   integer(kind=4) :: iprint,printZA(2,100)
   integer(kind=4) :: n_bin
   integer(kind=4) :: numcc
!---------------------------------------------------------------------
   integer(kind=4) :: nfit
   real(kind=8),allocatable :: cum_rho(:), cumm_fit(:), elv(:)
   real(kind=8) :: energy, prob, prob1(0:60), prob2(0:60)
   real(kind=8) :: rho(0:60,0:1)
   real(kind=8) :: rho_Fm, sig2, apu
   real(kind=8) :: xj, jfac, pfac, dde
   real(kind=8) :: E0, T, E1
   real(kind=8) :: sum_rho
   real(kind=8) :: K_vib, K_rot
   real(kind=8) :: f_E, str_E, f_M, str_M
   character(len=1) :: char
   character(len=1) :: char_pos,char_neg
   integer(kind=4) :: lev_below, char_start
   logical :: lprint
   real(kind=8) :: pmode, pe1, pbb
   integer(kind=4) :: num_points
!---------   External functions
   real(kind=8) :: spin_fac, parity_fac
   real(kind=8) :: EL_f, EL_trans
   real(kind=8) :: ML_f, ML_trans

!--------------   Start subrotuine
   num_points = int(30.0d0/de,kind=4) + 1

   char_pos='+'
   char_neg='-'
   iprint=0
   printZA(1:2,1:100)=0
   do i=1,num_comp
      lprint=.false.
      do j=1,iprint           !   Check to see if this nucleus has been printed already
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
      if(char_start == 1)write(13,'(i3,a2)')nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
      if(char_start == 2)write(13,'(i3,a1)')nucleus(i)%A,nucleus(i)%atomic_symbol(char_start:2)
      write(13,'(''Maximum excitation energy = '',f10.3,'' MeV'')')nucleus(i)%Ex_max
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
         write(13,'(2x,a8,1x,''Separation energy = '',f10.3,'' MeV'')')          &
           particle(k)%name,nucleus(i)%sep_e(k)
      end do

!      write(6,*)'Level density information for this compound nucleus - not target'
      if(nucleus(i)%D0exp > 0.0)then
         write(13,'('' Experimental D0   ='',f10.3,'' +/- '',f10.3,'' eV'')')    &
                     nucleus(i)%D0exp,nucleus(i)%dD0exp
      else
         write(13,'('' Experimental D0   =  UNAVAILABLE'')')
      end if
      write(13,'('' Calculated D0     ='',f10.3,'' eV'')')                       &
                  nucleus(i)%D0

      if(nucleus(i)%D1exp > 0.0)then
         write(13,'('' Experimental D1 ='',f10.3,'' +/- '',f10.3,'' eV'')')      &
                     nucleus(i)%D1exp,nucleus(i)%dD1exp
      else
         write(13,'('' Experimental D1   =  UNAVAILABLE'')')
      end if
      write(13,'('' Calculated D1     ='',f10.3,'' eV'')')                       &
                  nucleus(i)%D1


!      write(6,*)'Gamma_gamma data for this compound nucleus - not as a target'
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
      write(13,*)'Discrete states used in calculation'
      write(13,'(''State'',6x,''Energy'',2x,''Spin'','//               &
                 '1x,''Par'',8x,''T1/2'',6x,'//                        &
                 '''Isomer'')')
      lev_below = 0
      do k = 1,nucleus(i)%num_discrete
         if(nucleus(i)%state(k)%energy > nucleus(i)%level_param(7)     &
           .and. lev_below == 0)lev_below = 1 
         if(lev_below /= 1)then
            if(i == itarget .and. k == target%istate)then
               write(13,'(2x,i3,2x,f10.5,1x,f5.1,1x,f3.0,1x,'//        &
                          '1x,e15.7,1x,l4,4x,''Target State'')')       &
                 k,                                                    &
                 nucleus(i)%state(k)%energy,                           &
                 nucleus(i)%state(k)%spin,                             &
                 nucleus(i)%state(k)%parity,                           &
                 nucleus(i)%state(k)%t12,                              &
                 nucleus(i)%state(k)%isomer
            else
               write(13,'(2x,i3,2x,f10.5,1x,f5.1,1x,f3.0,1x,'//        &
                          '1x,e15.7,1x,l4)')                           &
                 k,                                                    &
                 nucleus(i)%state(k)%energy,                           &
                 nucleus(i)%state(k)%spin,                             &
                 nucleus(i)%state(k)%parity,                           &
                 nucleus(i)%state(k)%t12,                              &
                 nucleus(i)%state(k)%isomer
            end if
        elseif(lev_below == 1)then
            lev_below = 2
            write(13,'(''   ---------------   Ecut   ---------------------'')')
            write(13,'(2x,i3,2x,f10.5,1x,f5.1,1x,f3.0,1x,'//           &
                       '1x,e15.7,1x,l4)')                              &
              k,                                                       &
              nucleus(i)%state(k)%energy,                              &
              nucleus(i)%state(k)%spin,                                &
              nucleus(i)%state(k)%parity,                              &
              nucleus(i)%state(k)%t12,                                 &
              nucleus(i)%state(k)%isomer
        end if
      end do

      write(13,*)'Decay properties of discrete states'
      do k = nucleus(i)%num_discrete, 1, -1
         write(13,*)'Decay of state',k
         write(13,*)'Number of transitions',nucleus(i)%state(k)%nbranch
         write(13,'(''       i --->   f            branch      prob_gamma         prob_ic'')')
         write(13,'(''     ---      ---     -------------   -------------   -------------'')')

         prob = 0.0d0
         do m = 1, nucleus(i)%state(k)%nbranch
            prob = prob + nucleus(i)%state(k)%branch(m)
            write(13,'(4x,i4,'' --->'',i4,4(3x,e15.7))')k,nucleus(i)%state(k)%ibranch(m),    &
                                                          nucleus(i)%state(k)%branch(m),     &
                                                          nucleus(i)%state(k)%p_gamma(m),    &
                                                          nucleus(i)%state(k)%p_ic(m)
  
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
            write(13,'(1xi4,2(1x,f4.1),1x,f10.4)')k,OpticalCS%state(k)%spin,      &
                OpticalCS%state(k)%parity,OpticalCS%state(k)%energy
         end do
         write(13,*)
         write(13,'(''Additional DWBA states'')')
         write(13,*)
         write(13,'(''    N    J  PAR    K     Energy'')')
         write(13,'(''  ---  ---  ---  ---     ------'')')
         do k = numcc + 1, OpticalCS%numcc
            write(13,'(1xi4,3(1x,f4.1),1x,f10.4)')k,OpticalCS%state(k)%spin,      &
                OpticalCS%state(k)%parity,OpticalCS%state(k)%K,OpticalCS%state(k)%energy
         end do
      end if


!------------------------------------------------------------------
      E0 = nucleus(i)%level_param(15)
      T = nucleus(i)%level_param(14)
      E1 = -10.0d0
      if(E0 < 0.0d0)E1 = T*log(1.0d0-exp(E0/T))
      write(13,*)
      write(13,*)'Level Density information'
      write(13,*)'E1 defined as the energy so that int(E1,0)rho(E)dE = 1'
      write(13,*)'Generally, E1 < 0 if E0 < 0'
      write(13,*)'Otherwise E1 undefined as int(-infty,0)rho(E)DE < 1'
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
      write(13,'('' Calculated D0   ='',f10.3,'' eV'')')    &
                  nucleus(i)%D0
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
      if(nucleus(i)%level_param(10) == 0)then
         write(13,*)'No rotational collective enhancement included'
         write(13,*)'K_rot = 1.0'
      else
         write(13,*)'Rotational Collective enhancement model = ',    &
             nucleus(i)%level_param(10)
         write(13,*)'Collective enhancement factors'
         write(13,*)'K_rot = Max[x(1)*(Factor-1)/(1+exp((E-x(2))/x(3))),0]+1'
         if(nucleus(i)%level_param(10) == 1)then
            write(13,*)'Factor = x(1)*sig2*(1+beta2/3)'
         elseif(nucleus(i)%level_param(10) == 2)then
            write(13,*)'Factor = 2.0*x(1)*sig2*(1+beta2/3)'
         elseif(nucleus(i)%level_param(10) == 3)then
            write(13,*)'Factor = x(1)*sig2**3/2*(1+beta2/3)*(1-2*beta2/3)'
         elseif(nucleus(i)%level_param(10) == 4)then
            write(13,*)'Factor = 2.0*x(1)*sig2**3/2*(1+beta2/3)*(1-2*beta2/3)'
         end if

         write(13,*)'Rotational collective-enhancement parameters'

         write(13,'(''x(1) = '',f10.4)')nucleus(i)%rot_enh(1)
         write(13,'(''x(2) = '',f10.4)')nucleus(i)%rot_enh(2)
         write(13,'(''x(3) = '',f10.4)')nucleus(i)%rot_enh(3)
         write(13,*)
      end if
      if(nucleus(i)%level_param(11) == 0)then
         write(13,*)'No vibrational collective enhancement included'
         write(13,*)'K_vib = 1.0'
      elseif(nucleus(i)%level_param(11) == 1)then
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
      elseif(nucleus(i)%level_param(11) == 2)then
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
      elseif(nucleus(i)%level_param(11) == 3)then
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
      elseif(nucleus(i)%level_param(11) == 4)then
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

      write(13,'(''Level density States/MeV '')')
      write(13,'(''Positive Parity '')')
      write(13,1901)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
 1901 format('      E         a(U)    sqrt(sig2)   P_fac      K_vib      K_rot       Enh         Rho     ',     &
            40(5x,'J = ',f4.1,3x))
      pmode = nucleus(i)%level_param(16)
      pe1 = nucleus(i)%level_param(17)
      pbb = nucleus(i)%level_param(18)

      do k = 1, nucleus(i)%nbin
         energy = nucleus(i)%e_grid(k)
         sum_rho = 0.0d0
         call rhoe(energy,nucleus(i)%level_param,                      &
                   nucleus(i)%vib_enh,                                 &
                   nucleus(i)%rot_enh,                                 &
                   ia,rho_Fm,apu,sig2,K_vib,K_rot)
          ip = 1
          pfac=parity_fac(energy,xj,ip,pmode,pe1,pbb)
          do jj = 0, min(j_max,60)
             xj = jj + nucleus(i)%jshift
             jfac=spin_fac(xj,sig2)
             rho(jj,ip)=rho_FM*jfac*pfac
             sum_rho = sum_rho + nucleus(i)%bins(jj,ip,k)%rho
          end do
          write(13,'(1x,f8.3,6(1x,f10.3),60(1x,e15.7))')                         &
               energy,apu,sqrt(sig2),pfac,K_vib,K_rot,K_vib*K_rot,sum_rho,       &
               (nucleus(i)%bins(jj,ip,k)%rho,jj = 0, min(j_max,60))
      end do
      write(13,'(''Negative Parity '')')
      write(13,1902)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
 1902 format('      E         a(U)    sqrt(sig2)   P_fac      K_vib      K_rot       Enh         Rho     ',     &
            40(5x,'J = ',f4.1,3x))
      pmode = nucleus(i)%level_param(16)
      pe1 = nucleus(i)%level_param(17)
      pbb = nucleus(i)%level_param(18)

      do k = 1, nucleus(i)%nbin
         energy = nucleus(i)%e_grid(k)
         sum_rho = 0.0d0
         call rhoe(energy,nucleus(i)%level_param,                      &
                   nucleus(i)%vib_enh,                                 &
                   nucleus(i)%rot_enh,                                 &
                   ia,rho_Fm,apu,sig2,K_vib,K_rot)
          ip = 0
          pfac=parity_fac(energy,xj,ip,pmode,pe1,pbb)
          do jj = 0, min(j_max,60)
             xj = jj + nucleus(i)%jshift
             jfac=spin_fac(xj,sig2)
             rho(jj,ip)=rho_FM*jfac*pfac
             sum_rho = sum_rho + nucleus(i)%bins(jj,ip,k)%rho
          end do
          write(13,'(1x,f8.3,6(1x,f10.3),60(1x,e15.7))')                         &
               energy,apu,sqrt(sig2),pfac,K_vib,K_rot,K_vib*K_rot,sum_rho,       &
               (nucleus(i)%bins(jj,ip,k)%rho,jj = 0, min(j_max,60))
      end do
 !------------------------------------------------------------------
      nfit = max(nucleus(i)%num_discrete,nucleus(i)%ncut)
      allocate(cum_rho(nfit))
      allocate(cumm_fit(nfit))
      allocate(elv(nfit))
      if(nucleus(i)%ncut >= 5) then
         if(nucleus(i)%fit_ematch)then
            write(13,*)'E_match was fit cummulative density for known discrete states'
         else
            write(13,*)'E_match was not fit cummulative density for known discrete states'
            write(13,*)'User input over rides fitting or default was used'
         end if
         nfit = nucleus(i)%ncut
         write(13,*)
         write(13,*)'Cummulative level density up to E_cut'
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
         write(13,*)
         write(13,*)'Modeled cummulative level density up to E_cut'
         cumm_fit(1:nfit) = 0.0d0
         call cumm_rho(nfit,elv,ia,nucleus(i)%level_param,             &
                       nucleus(i)%vib_enh,nucleus(i)%rot_enh,cumm_fit)
         do j = 1, nfit
            write(13,'(f10.4,1x,f15.2)')elv(j),cumm_fit(j)
         end do

         if(nucleus(i)%ncut /= nucleus(i)%num_discrete)then
            nfit = nucleus(i)%num_discrete
            write(13,*)
            write(13,*)'Cummulative level density for all discrete states'
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
               if(j < nfit)then
                  write(13,'(f10.4,1x,f15.2)')elv(j+1), cum_rho(j)
               end if
            end do
            write(13,*)
            write(13,*)'Modeled cummulative level density up to maximum discrete state'
            cumm_fit(1:nfit) = 0.0d0
            call cumm_rho(nfit,elv,ia,nucleus(i)%level_param,             &
                          nucleus(i)%vib_enh,nucleus(i)%rot_enh,cumm_fit)
            do j = 1, nfit
               write(13,'(f10.4,1x,f15.2)')elv(j), cumm_fit(j)
            end do
         end if
      end if
      deallocate(cum_rho)
      deallocate(cumm_fit)
      deallocate(elv)

      nucleus(i)%nbin_em = num_points

      if(nucleus(i)%e1_model == 1)then
         write(13,'(''E1 model used = '',''Lorenztian'')')
      elseif(nucleus(i)%e1_model == 2)then
         write(13,'(''E1 model used = '',''Kopecky-Uhl'')')
      end if
!------------------------------------------------------------------
      write(13,*)
      write(13,*)'Electric-dipole resonance parameters'
      do k = 1, 3
         write(13,'(''Mode ='',i2,'//                             &
                    ''' Centroid = '',f16.6,'' MeV'','//          &
                    ''' Width = '',f16.6,'' MeV'','//             &
                    ''' Strength = '',f16.6,'' mb'')')            &
             k,                                                   &
             nucleus(i)%er_E1(k),                                 &
             nucleus(i)%gr_E1(k),                                 &
             nucleus(i)%sr_E1(k)
            end do
      write(13,'(''E1 strength function'')')
      write(13,'(''     Energy         F               T'')')
      write(13,'(''  ---------   -------------   -------------'')')
      l_radiation = 1
      do j = 0, nucleus(i)%nbin_em
         energy = dfloat(j)*de
         f_E = EL_f(i, l_radiation, energy, nucleus(i)%sep_e(1))
         str_E = EL_trans(i, l_radiation, energy, nucleus(i)%sep_e(1))
         write(13,'(1x,f10.5,(2(1x,e15.7)))')                     &
              energy, f_E, str_E
!              energy,nucleus(i)%f_E(j,1),nucleus(i)%str_E(j,1)
      end do
!------------------------------------------------------------------
!-----------   The other Electric multipoles - 
!-----------   different because E1 has up to three modes
!------------------------------------------------------------------
      do l_radiation = 2, nucleus(i)%lmax_E
         write(13,*)
         write(13,'(''Electric-'',i1,'' resonance parameters'')')  &
            l_radiation
            write(13,'('' Centroid = '',f16.6,'' MeV'','//         &
                       ''' Width = '',f16.6,'' MeV'','//           &
                       ''' Strength = '',e15.7,'' mb'')')          &
                nucleus(i)%er_E(l_radiation),                      &
                nucleus(i)%gr_E(l_radiation),                      &
                nucleus(i)%sr_E(l_radiation)
         write(13,'(''E'',i1,'' Strength function'')')l_radiation
         write(13,'(''     Energy         F               T'')')
         write(13,'(''  ---------   -------------   -------------'')')
         do j = 0, nucleus(i)%nbin_em
            energy = dfloat(j)*de
            f_E = EL_f(i, l_radiation, energy, nucleus(i)%sep_e(1))
            str_E = EL_trans(i, l_radiation, energy, nucleus(i)%sep_e(1))
            write(13,'(1x,f10.5,(2(1x,e15.7)))')                   &
                 energy, f_E, str_E
!            write(13,'(1x,f10.5,(20(1x,e15.7)))')     &
!                   energy,                            &
!                   nucleus(i)%f_E(j,l_radiation),     &
!                   nucleus(i)%str_E(j,l_radiation)
         end do
      end do
!------------------------------------------------------------------
      do l_radiation = 1, nucleus(i)%lmax_M
         write(13,*)
         write(13,'(''Magnetic-'',i1,'' resonance parameters'')')  &
            l_radiation
         write(13,'('' Centroid = '',f10.5,'' MeV'','//            &
                       ''' Width = '',f10.5,'' MeV'','//           &
                       ''' Strength = '',e12.4,'' mb'')')          &
                nucleus(i)%er_M(l_radiation),                      &
                nucleus(i)%gr_M(l_radiation),                      &
                nucleus(i)%sr_M(l_radiation)
         write(13,'(''M'',i1,'' Strength function'')')l_radiation
         write(13,'(''     Energy         F               T'')')
         write(13,'(''  ---------   -------------   -------------'')')
         do j = 0, nucleus(i)%nbin_em
            energy = dfloat(j)*de
            f_M = ML_f(i, l_radiation, energy)
            str_M = ML_trans(i, l_radiation, energy)
            write(13,'(1x,f10.5,(2(1x,e15.7)))')                   &
                 energy, f_M, str_M
!            write(13,'(1x,f10.5,(20(1x,e15.7)))')                  &
!                  energy,                                          &
!                  nucleus(i)%f_M(j,l_radiation),                   &
!                  nucleus(i)%str_M(j,l_radiation)
         end do
      end do
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
         write(13,*)'Fission parameters'
         do j = 1, nucleus(i)%F_n_barr
            E0 = nucleus(i)%F_barrier(j)%level_param(15)
            T = nucleus(i)%F_barrier(j)%level_param(14)
            E1 = -10.0d0
            if(E0 < 0.0d0)E1 = T*log(1.0d0-exp(E0/T))
            write(13,'(''Parameters for Fission Barrier #'',i3)')j
            write(13,'(''Height = '',f10.4)')nucleus(i)%F_Barrier(j)%barrier
            write(13,'(''Width  = '',f10.4)')nucleus(i)%F_Barrier(j)%hbw
            write(13,'(''Deformation parameter beta(2) = '',f10.6)')nucleus(i)%F_Barrier(j)%beta_2
            write(13,'(''Excitation energy dependent barrier, F = F*x(1)*exp(-x(3)**2*(Ex-x(2))**2)'')')
            write(13,'(''x(1) = '',f10.4,'' x(2) = '',f10.4,'' x(3) = '',f10.4)')            &
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


            if(nucleus(i)%F_Barrier(j)%level_param(10) == 0)then
               write(13,*)'No rotational collective enhancement included'
               write(13,*)'K_rot = 1.0'
            else
               write(13,*)'Rotational Collective enhancement model = ',    &
                  nucleus(i)%F_Barrier(j)%level_param(10)
                write(13,*)'Collective enhancement factors'
                write(13,*)'K_rot = Max[x(1)*(Factor-1)/(1+exp((E-x(2))/x(3))),0]+1'
                if(nucleus(i)%F_Barrier(j)%level_param(10) == 1)then
                   write(13,*)'Factor = x(1)*sig2*(1+beta2/3) --  Axially symmetric'
                elseif(nucleus(i)%F_Barrier(j)%level_param(10) == 2)then
                   write(13,*)'Factor = 2.0*x(1)*sig2*(1+beta2/3) -- Left-right asymmetric'
                elseif(nucleus(i)%F_Barrier(j)%level_param(10) == 3)then
                   write(13,*)'Factor = x(1)*sqrt(pi/2)*sig2**3/2*(1+beta2/3)*(1-2*beta2/3) -- Triaxial and left-right asymmetric'
                elseif(nucleus(i)%F_Barrier(j)%level_param(10) == 4)then
                   write(13,*)                     &
          'Factor = 2.0*x(1)*sqrt(pi/2)*sig2**3/2*(1+beta2/3)*(1-2*beta2/3) -- Triaxial and not left-right asymmetric'
                end if

                write(13,'('' Rotational collective-enhancement parameters'')')
  
                write(13,'(''x(1) = '',f10.4)')nucleus(i)%F_Barrier(j)%rot_enh(1)
                write(13,'(''x(2) = '',f10.4)')nucleus(i)%F_Barrier(j)%rot_enh(2)
                write(13,'(''x(3) = '',f10.4)')nucleus(i)%F_Barrier(j)%rot_enh(3)
                write(13,*)
             end if
            if(nucleus(i)%F_Barrier(j)%level_param(11) == 0)then
               write(13,*)'No vibrational collective enhancement included'
               write(13,*)'K_vib = 1.0'
            elseif(nucleus(i)%F_Barrier(j)%level_param(11) == 1)then
               write(13,*)'Vibrational collective enhancement model 1 used'
               write(13,*)'K_vib = exp(0.05555*A**(2./3.)*T**(4./3.))'
               write(13,*)'K_vib = Max[x(1)*(K_vib-1)/(1+exp((E-x(2))/x(3))),0]+1'
               write(13,*)'T = sqrt(U/apu)'
               write(13,'(''x(1) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(1)
               write(13,'(''x(2) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(2)
               write(13,'(''x(3) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(3)
               write(13,*)
            elseif(nucleus(i)%F_Barrier(j)%level_param(11) == 2)then
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
            elseif(nucleus(i)%F_Barrier(j)%level_param(11) == 3)then
               write(13,*)'Vibrational collective enhancement model 1 used'
               write(13,*)'K_vib = exp(0.05555*A**(2./3.)*T**(4./3.))'
               write(13,*)'K_vib = Max[x(1)*(K_vib-1)/(1+exp((E-x(2))/x(3))),0]+1'
               write(13,*)'T = sqrt(U/a);  a = A/13'
               write(13,'(''x(1) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(1)
               write(13,'(''x(2) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(2)
               write(13,'(''x(3) = '',f10.4)')nucleus(i)%F_Barrier(j)%vib_enh(3)
               write(13,*)
            elseif(nucleus(i)%F_Barrier(j)%level_param(11) == 4)then
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

            if(nucleus(i)%nbin > 1)then
               dde = nucleus(i)%e_grid(2) - nucleus(i)%e_grid(1)
            else
               dde = 0.2
            end if
            n_bin = int((20.0d0 - nucleus(i)%F_barrier(j)%ecut)/dde)
            write(13,'(''Level density States/MeV'')')
            write(13,'(''Positive parity only'')')
            write(13,1903)(real(jj) + nucleus(i)%jshift, jj = 0, min(j_max,60))
 1903       format('      E         a(U)    sqrt(sig2)   P_fac      K_vib      K_rot       Enh         Rho     ',     &
                  40(5x,'J = ',f4.1,3x))
            pmode = nucleus(i)%F_barrier(j)%level_param(16)
            pe1 = nucleus(i)%F_barrier(j)%level_param(17)
            pbb = nucleus(i)%F_barrier(j)%level_param(18)

            do k = 1, n_bin
               energy = real(k,kind=8)*dde
               if(energy < nucleus(i)%F_barrier(j)%ecut) cycle
               sum_rho = 0.0d0
               call rhoe(energy,nucleus(i)%F_barrier(j)%level_param,                 &
                         nucleus(i)%F_barrier(j)%vib_enh,                            &
                         nucleus(i)%F_barrier(j)%rot_enh,                            &
                         ia,rho_Fm,apu,sig2,K_vib,K_rot)
               ip=1
               pfac=parity_fac(energy,xj,ip,pmode,pe1,pbb)
               do jj = 0, min(j_max,60)
                  xj = jj + nucleus(i)%jshift
                  jfac=spin_fac(xj,sig2)
                  rho(jj,ip)=rho_FM*jfac*pfac
                  sum_rho = sum_rho + rho(jj,ip)
               end do
               write(13,'(1x,f8.3,6(1x,f10.3),60(1x,e15.7))')                        &
                   energy,apu,sqrt(sig2),pfac,K_vib,K_rot,K_vib*K_rot,rho_FM*pfac,   &
                   (rho(jj,ip),jj = 0, min(j_max,60))
            end do
         end do
         write(13,*)
      end if

!--------------   HF-denominators
      write(13,*)
      write(13,*)'Hauser-Feshbach denominators'
      write(13,*)'Positive Parity'
      do n = 1, nucleus(i)%nbin
         energy=nucleus(i)%e_grid(n)
         write(13,'(f10.6,60(1x,e15.7))')                                            &
           energy,(nucleus(i)%bins(j,0,n)%HF_den,j=0,j_max)
      end do
      write(13,*)
      write(13,*)'Negative Parity'
      do n = 1, nucleus(i)%nbin
         energy=nucleus(i)%e_grid(n)
         write(13,'(f10.6,60(1x,e15.7))')                                            &
           energy,(nucleus(i)%bins(j,1,n)%HF_den,j=0,j_max)
      end do
      write(13,*)
      write(13,*)'Summed decay transmission coefficients for each channel'
      do if1 = 1, nucleus(i)%num_decay
         k = nucleus(i)%decay_particle(if1)
         write(13,*)
         write(13,*)'Particle type = ',k
         write(13,*)'Positive Parity'
         do n=1,nucleus(i)%nbin
            energy=nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob1(j) = nucleus(i)%bins(j,0,n)%HF_trans(if1)
            end do
            write(13,'(f10.6,60(1x,e15.7))')                                         &
              energy,(prob1(j),j = 0, j_max)
         end do
         write(13,*)
         write(13,*)'Negative Parity'
         do n = 1, nucleus(i)%nbin
            energy=nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob2(j) = nucleus(i)%bins(j,1,n)%HF_trans(if1)
            end do
            write(13,'(f10.6,60(1x,e15.7))')                                         &
              energy,(prob2(j),j = 0, j_max)
         end do
      end do
      
      if(nucleus(i)%Fission)then
         if1 = nucleus(i)%num_decay + 1
         write(13,*)
         write(13,*)'Positive Parity: Decay type = Fission'
         do n = 1, nucleus(i)%nbin
            energy = nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob1(j) = nucleus(i)%bins(j,0,n)%HF_trans(if1)
            end do
            write(13,'(f10.6,60(1x,e15.7))')                                         &
              energy,(prob1(j),j = 0, min(j_max,60))
         end do
         write(13,*)
         write(13,*)'Negative Parity: Decay type = Fission'
         do n = 1, nucleus(i)%nbin
            energy = nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob2(j) = nucleus(i)%bins(j,1,n)%HF_trans(if1)
            end do
            write(13,'(f10.6,60(1x,e15.7))')                                         &
              energy,(prob2(j),j = 0, min(j_max,60))
         end do

         do ii = 1, nucleus(i)%F_n_barr
            write(13,*)
            write(13,*)'Positive Parity: Trans-Coef for Fission Barrier #', ii
            do n = 1, nucleus(i)%nbin
               energy = nucleus(i)%e_grid(n)
               do j = 0, min(j_max,60)
                  prob1(j) = nucleus(i)%bins(j,0,n)%HF_trans(if1+ii)
               end do
               write(13,'(f10.6,60(1x,e15.7))')                                      &
                 energy,(prob1(j),j = 0, min(j_max,60))
            end do
            write(13,*)
            write(13,*)'Negative Parity: Trans-Coef for Fission Barrier #', ii
            do n = 1, nucleus(i)%nbin
               energy = nucleus(i)%e_grid(n)
               do j = 0, min(j_max,60)
                  prob2(j) = nucleus(i)%bins(j,1,n)%HF_trans(if1+ii)
               end do
               write(13,'(f10.6,60(1x,e15.7))')                                      &
                 energy,(prob2(j),j = 0, min(j_max,60))
            end do
         end do
      end if
      write(13,*)
      write(13,*)'Decay probabilities for each channel'
      do if1 = 1, nucleus(i)%num_decay
         k = nucleus(i)%decay_particle(if1)
         write(13,*)
         if(k < 7)write(13,*)'Particle type = ',particle(k)%name
         if(k == 7)write(13,*)'Decay type = Fission'
         write(13,*)'Positive Parity'
         do n = 1,nucleus(i)%nbin
            energy = nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob1(j) = nucleus(i)%bins(j,0,n)%HF_prob(if1)
            end do
            write(13,'(f10.6,60(1x,e15.7))')                                         &
              energy,(prob1(j),j = 0, min(j_max,60))
         end do
         write(13,*)'Negative Parity'
         do n = 1,nucleus(i)%nbin
            energy = nucleus(i)%e_grid(n)
            do j = 0, min(j_max,60)
               prob2(j) = nucleus(i)%bins(j,1,n)%HF_prob(if1)
            end do
            write(13,'(f10.6,60(1x,e15.7))')                                         &
              energy,(prob2(j),j = 0, min(j_max,60))
         end do
      end do
   end do
   flush(13)
   return
end subroutine output_nucleus_data

