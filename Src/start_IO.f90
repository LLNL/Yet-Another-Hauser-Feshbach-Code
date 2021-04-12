!
!*******************************************************************************
!
subroutine start_IO
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
      write(13,*)'Energy grid set up with de = ',nucleus(i)%delta_e(1)
!      write(13,*)'Energy grid set up with de = ',de
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
!   write(13,*)'output_mode = ',output_mode
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
   write(13,*)'all_discrete_states = ',all_discrete_states
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
