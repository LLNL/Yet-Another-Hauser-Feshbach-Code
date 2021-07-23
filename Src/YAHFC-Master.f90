!
!*****************************************************************************
!
program YAHFC_MASTER
!
!*****************************************************************************
!
!  Discussion:
!
!    This program is the driver for the Hauser-Feshbach code system YAHFC
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        constants
!        nodeinfo
!        print_control
!        useful_data
!        log_factorial
!        nuclei
!        Channel_info
!        Scatter_info
!        particles_def
!        directory_structure
!        pre_equilibrium_no_1
!        Gauss_integration
!
!     Subroutines:
!
!        log_fact
!        gauss_quad
!        system_clock
!        getenv
!        particle_data
!        parse_string
!        parse_command
!        lower_case_word
!        set_up_decay_chain
!        get_spectrum
!        gdr_param
!        EM_str_param
!        get_lev_den
!        fit_lev_den
!        Fission_data
!        Fission_levels
!        finish_lev_den
!        Find_T_E0
!        Setup_Energy_bins
!        fix_incident_energies
!        fix_pop_energies
!        set_min_particle_energy
!        optical_setup
!        reaction_max_J
!        remove_states
!        modeled_level_density
!        read_level_density
!        fit_nuke_Gamma_gamma
!        HF_denominator
!        start_IO
!        output_trans_coef
!        output_nucleus_data
!        memory_used
!        nucleus_label
!        check_directories
!        compound_xs
!        gauss_quad
!        HF_primary_decay_setup
!        photo_xs
!        pre_equilibrium_1
!        PREEQ_sample
!        find_prob_point
!        MC_primary_decay
!        MC_decay_bin
!        force_decay
!        MC_decay_state
!        Boost_frame
!        print_preeq_spectra
!        print_direct_spectra
!        print_x_particle_spectra
!        print_primary_decay
!        Legendre_expand
!        print_channel_gammas
!        print_nuke_data
!        print_total_cs
!        print_reaction_cs
!        print_absorption_cs
!        print_preeq_cs
!        print_direct_cs
!        print_elastic
!        print_compound_elastic
!        print_inelastic
!        print_pop_mult_data
!        print_channels
!        print_fission_cs
!        read_saved_parameters
!        exit_YAHFC
!
!     External functions:
!
!        integer(kind=4) :: find_channel
!        integer(kind=4) :: find_ibin
!        real(kind=8) :: random_64
!        real(kind=8) :: random_32
!        real(kind=8) :: poly
!        real(kind=8) :: clebr
!        real(kind=8) :: interpolate_exp
!        real(kind=8) :: Gauss_var
!        integer(kind=4) :: rank_commands
!
!     MPI routines:
!
!        MPI_INIT
!        MPI_COMM_RANK
!        MPI_COMM_SIZE
!        MPI_Abort    -----   via exit_YAHFC
!        MPI_Bcast
!        MPI_Barrier
!        MPI_AllReduce
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
!*****************************************************************************
!
!----------------------------------------------------------------------
!-------    Use modules
!----------------------------------------------------------------------
      use variable_kinds
      use options
      use constants
      use nodeinfo
      use print_control
      use useful_data
      use log_factorial
      use nuclei
      use Channel_info
      use Scatter_info
      use particles_def
      use directory_structure
      use pre_equilibrium_no_1
      use Gauss_integration
!---------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------
#if(USE_MPI == 1)
      include 'mpif.h'
      integer(kind=4) :: num_data
#endif
!---------------------------------------------------------------------
      character(len=1) quote
      character(len=80) bigblank        ! 80 blank character
      character(len=132) blank132
      logical :: call_tco
      logical :: write_error
!-----------------------------------------------------------------------
      character(len=80) directory
      integer(kind=4) ifile, idir, core_file
      integer(kind=4) icount
      character(len=40) file_name
      character(len=132) out_buff
      integer(kind=4) :: icharr, il1
!      character(len=132) :: unix_command
!-----------------------------------------------------------------------
      integer(kind=4) :: num_e
      real(kind=8) :: ee_max

!-----------------------------------------------------------------------
!----------   Temporary storage for state data if we need to remove discrete states

      real(kind=8) :: E_in, e_rel
      real(kind=8) :: E_f, Ex_f
      real(kind=8) :: rel_factor
      integer(kind=4) :: iproj, itarget
      integer(kind=4) :: jproj

      real(kind=8) :: spin_target, spin_proj
      real(kind=8) :: sp1, sp2
      real(kind=8) :: cs_threshold
      integer(kind=4) :: isp
      integer(kind=4) :: isp_max
      integer(kind=4) :: istate, jstate
      real(kind=8) :: sum
!      real(kind=8) :: sum_check
      real(kind=8) :: ex_tot

      real(kind=8) :: temp_cs
      integer(kind=4) :: num_eff, num_elastic
      real(kind=8) :: xnum_eff, xnum_elastic
      real(kind=8) :: smear(0:6)
      real(kind=8) :: shift, shift_min, shift_max
      integer(kind=4) :: icc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Information needed to parse a 132-length character string into
!----   individual words (at most 66 possible words)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer(kind=4) :: numw
      integer(kind=4) :: startw(66), stopw(66)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Data required for Haser-Feshbach Monte Carlo probabilities
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer(kind=4) :: ii, nn
      integer(kind=4) :: l_i, is_i, l_f
      integer(kind=4) :: icomp_i, ip_i, nbin_i
      integer(kind=4) :: icomp_f, ip_f, nbin_f, idb
      real(kind=8) :: xip_f
      character(len=1) ch_par(0:1)
      integer(kind=4) :: i_f, isp_f_max


      integer(kind=4) :: num_part, n_dat, dim_part

      parameter (n_dat = 24)

      real(kind=8), allocatable :: part_data(:,:)
      real(kind=8), allocatable :: extra_angle_data(:,:)
      integer(kind=4) :: num_theta, nang, nang_tot
      integer(kind=4) :: Ang_L_max
      real(kind=8), allocatable :: part_Ang_data(:,:)
      real(kind=8), allocatable :: Leg_poly(:,:)

      integer(kind=4) ::  L_Ang
      real(kind=8) :: xL_Ang, xl_i

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----  Probabilities for incident channels
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real(kind=8) :: mass_target, mass_exit, Q
      real(kind=8) :: mass_proj
      integer(kind=4) :: num_channel
      real(kind=8), pointer, dimension(:) :: channel_prob
      integer(kind=4), pointer, dimension(:,:) :: ichannel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      Label for channels
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer(kind=4) :: ichann
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      Data to implement max_particle constraints during decays
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer(kind=4) :: num_part_type(0:6)
      real(kind=8) :: part_fact(0:7)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Random number data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real(kind=8) :: ran
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      integer(kind=4) :: l_max
      integer(kind=4) :: ifind, ilast
      integer(kind=4) :: num_s
      integer(kind=4) :: ictype
      real(kind=8) :: xI_f
      real(kind=8) ::  par_i, xnbin_i
      integer(kind=4) :: iX_i, Ix_f
      real(kind=8) :: preeq_prob

      real(kind=8) :: pop_sum, prev
      real(kind=8), allocatable, dimension (:) :: pop_prob
      integer(kind=4), allocatable, dimension (:) :: pop_j, pop_ip
      logical :: pop_input
      integer(kind=4) ijj

      real(kind=8), allocatable, dimension (:,:) :: decay_mult, decay_var

      integer(kind=4) :: num_energies

      real(kind=8) :: sum_e
      real(kind=8) :: e_diff, max_e_diff, avg_e_diff

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real(kind=8) :: KE
      real(kind=8) :: Ef
      real(kind=8) :: Ex
      integer(kind=4) :: i,j,k,l,m,is
      integer(kind=4) :: Z_i,A_i
      integer(kind=4) :: Z_f,A_f
      integer(kind=4) :: in,ip,ia
      integer(kind=4) :: icomp
      integer(kind=4) :: cum_fit
      character(len=1) overide
      real(kind=8) :: hang
      logical :: hung_up
      real(kind=8) :: sum_ie, sum_ce
      integer(kind=4) :: inuc
      integer(kind=4) :: max_num
!---------------------------------------------------------------------
      integer(kind=4) :: state_map(1000)
      integer(kind=4) :: ibranch, numd
      integer(kind=4) :: JJ
      logical :: check
!---------------------------------------------------------------------
      character(len=2) char_pos,char_neg
      integer(kind=4) :: max_part
!----------------------------------------------------------------------
      integer(kind=4) :: max_gammas
!----------------------------------------------------------------------
      integer(kind=4) :: nbin
      integer(kind=4) :: n, n_f, if1
      real(kind=8) :: xj_min, xj_max
      real(kind=8) :: xI_min, xI_max
      integer(kind=4) :: j_max
      real(kind=8) :: xj
      integer(kind=4) :: Ix, Ix_min, Ix_max
      real(kind=8) :: xnorm
!--------------------  Data to collect input information -------------
      integer(kind=4) :: read_err
      integer(kind=4) :: nchar
      character(len=132) :: command
      character(len=132) char_temp
      integer(kind=4) :: itemp
      integer(kind=4) :: num_command
      logical finish
!----------------------------------------------------------------
      real(kind=8), allocatable :: reaction_cs(:)
      real(kind=8), allocatable :: absorption_cs(:)
      real(kind=8), allocatable :: preeq_css(:,:)
      real(kind=8), allocatable :: preeq_spect(:,:)
      real(kind=8), allocatable :: preeq_spect_full(:,:)
      real(kind=8), allocatable :: fission_cs(:)
      character(len=5) target_label
      integer(kind=4) :: ntar
      integer(kind=4) :: ilab
      character(len=100) :: file_lab

!------------------------------------------------------------------
!--------   Data to reconstruct angular distributions
      integer(kind=4) :: ixx_max
      real(kind=8) :: delta_ix
      integer(kind=4) :: max_jx_10
      real(kind=8) :: delta_jx_10
      integer(kind=4) :: max_jx_20
      real(kind=8) :: delta_jx_20
      integer(kind=4) :: max_jx_50
      real(kind=8) :: delta_jx_50
      integer(kind=4) :: max_jx_100
      real(kind=8) :: delta_jx_100
      real(kind=8) :: avg, avg_diff, expected_diff
      integer(kind=4) :: jx
!------------------------------------------------------------------
!--------    Gauss-Legendre integration useful to construct "smooth"
!--------    angular distributions
      real(kind=8), allocatable :: xvalue(:)

!--------   Total Inelastic Scattering
      real(kind=8) :: tot_direct, sum_d
      logical direct, compound, preeq_decay, cc_decay, dwba_decay
      logical inelastic, compound_elastic, fission_decay
!-------    Sum of all processes in these different channels

!
!--------     Primary Gamma-Spectrum
      real(kind=8), allocatable :: Primary_Gamma_Spectrum(:)
!
!-------    Cross section and spectrum from x-particle emission
!
      real(kind=8), allocatable :: x_particle_cs(:,:)
      real(kind=8), allocatable :: x_particle_Spectrum(:,:)
      real(kind=8) :: de_spec2

      real(kind=8), allocatable :: direct_Spectrum(:)
      real(kind=8), allocatable :: dwba_Spectrum(:)
!---------   Statistics on Fissioning nuclei
      real(kind=8), allocatable :: Fiss_J_avg(:,:), Fiss_J_var(:,:)
      real(kind=8), allocatable :: Fiss_tally(:,:)

!--------     Spectrum of emitted particles in Lab frame
!--------     when pop_calc = .true.
      real(kind=8) :: angle_dist(0:90)

!--------   yrast  ---------------------------------
      real(kind=8) :: yrast(0:100,0:1)
      real(kind=8) :: yrast_actual(0:100,0:1)

      real(kind=8) :: x, x1, theta

!-------   Lorentz boost matrices to transform to Lab and COM frames
      real(kind=8) :: Boost_Lab(0:3,0:3), Boost_COM(0:3,0:3)
      real(kind=8) :: mass_1, mass_2
      real(kind=8) :: pp_res, p_res(0:3), v_res(0:3)
      real(kind=8) :: beta, gamma

      real(kind=8) :: theta_0, phi_0, T_1, phi, T_2, T_L, theta_L, phi_L

      integer(kind=4) :: n_min

      real(kind=8) :: check_sum, sum_inelastic

      real(kind=8) :: pnorm

      real(kind=8) :: tally_weight
      real(kind=8) :: tally_norm

!-------    logical variable to identify when an inelastic transition to continuous bins first decays 
!-------    by a gamma to a discrete state
      logical :: collect

      integer(kind=4),allocatable :: command_rank(:)
      character(len=132), dimension(:), allocatable :: command_buffer

      integer(kind=4) :: num_bad_samp_e, num_bad_spread_e

      integer :: cnt, cnt1, cnt_r, cnt_m

      real(kind=8) :: ratio

!-----------------------------------------------------
!---------    Directory structure of Library outputs
      character(len=132) :: lib_dir
      integer(kind=4) :: ilib_dir

!----   Interface for compound_xs, where pointers are allocated
   interface
      subroutine compound_xs(e_in,itarget,istate,iproj,sigma,   &
                             num_channel,channel_prob,ichannel)
         use variable_kinds
         use options
         use nuclei
         use particles_def
         use constants
         use nodeinfo
         implicit none
         real(kind=8) :: e_in
         integer(kind=4) :: itarget, istate, iproj
         real(kind=8) :: sigma
         integer(kind=4) :: num_channel
         real(kind=8), pointer, dimension(:) :: channel_prob
         integer(kind=4), pointer, dimension(:,:) :: ichannel
      end subroutine
   end interface

!----   Interface for photo_xs, where pointers are allocated

   interface
      subroutine photo_xs(e_in, itarget, istate, sigma,   &
                          num_channel, channel_prob, ichannel)
         use variable_kinds
         use options
         use nuclei
         use particles_def
         use constants
         use nodeinfo
         implicit none
         real(kind=8) :: e_in
         integer(kind=4) :: itarget, istate
         real(kind=8) :: sigma
         integer(kind=4) :: num_channel
         real(kind=8), pointer, dimension(:) :: channel_prob
         integer(kind=4), pointer, dimension(:,:) :: ichannel
      end subroutine
   end interface

   interface
      subroutine find_prob_point(num,prob_array,prob,ifind)
         implicit none
         integer(kind=4) :: num
         real(kind=8),pointer :: prob_array(:)
         real(kind=8) :: prob
         integer(kind=4) :: ifind
      end subroutine find_prob_point
   end interface

!--------------------External functions--------------------------------
   integer(kind=4) :: find_channel
   integer(kind=4) :: find_ibin
   real(kind=8) :: random_64
   real(kind=8) :: random_32
   real(kind=8) :: poly
   real(kind=8) :: clebr
   real(kind=8) :: interpolate_exp
   real(kind=8) :: Gauss_var
   integer(kind=4) :: rank_commands
!***********************************************************************
!***********************************************************************
!-----------------   Start Main body of calculation     ----------------
!***********************************************************************
!***********************************************************************

!----------------------------------------------------------------------+
!------   Setup MPI world                                              +
!----------------------------------------------------------------------+
#if(USE_MPI == 1)
   call MPI_INIT(ierr)

   icomm = MPI_COMM_WORLD
   call MPI_COMM_RANK(icomm, iproc, ierr)
   call MPI_COMM_SIZE(icomm, nproc, ierr)
#else
   nproc = 1
   iproc = 0
   icomm = 0
   ierr = 0
#endif


   ilast = index(version,' ') - 1
   if(iproc == 0)then
      write(6,'(''****************************************************************'')')
      write(6,'(''* Yet Another Hauser Monte Hauser-Feshbach code - Monte Carlo  *'')')
      write(6,'(''* Version # '',a8,43x,''*'')')version(1:ilast)
      write(6,'(''* Author W. E. Ormand, Lawrence Livermore National Laboratory  *'')')
!      write(6,'(''* This code system is distributed under the license            *'')')
      write(6,'(''*--------------------------------------------------------------*'')')
      write(6,'(''* This code makes use of RIPL-3 data, with modifcations by the *'')')
      write(6,'(''* author and users. RIPL-3 data is available at                *'')')
      write(6,'(''* https://www-nds.iaea.org/RIPL-3/. Discrete level information *'')')
      write(6,'(''* is based on M. Verpelli and R. Capote Noy, INDC(NDS)-0702,   *'')')
      write(6,'(''* IAEA, 2015, and R. Capote, et al., Nuclear Data Sheets 110   *'')')
      write(6,'(''* (2009), 3107-3214.                                           *'')')
      write(6,'(''****************************************************************'')')
      write(6,*)
   end if

      quote = "'"
      angle_dist(0:90) = 0.0d0

      out_file(1:132) = ' '

!----------  Default Options
      print_output = .true.
      print_spectra = .true.
      print_libraries = .true.
      call_tco = .true.
      ex_set = .false.
      track_gammas = .false.
      track_primary_gammas = .false.
      fit_gamma_gamma = .true.
      lev_fit_d0 = .true.
      fit_aparam = .true.
      fit_ematch = .true.
      Apply_Coulomb_Barrier = .true.
      all_discrete_states = .false.
      Preeq_g_a = .false.
      biased_sampling = .false.
      explicit_channels = .false.
      dump_events = .false.
      binary_event_file = .true.
      event_generator = .false.
      use_unequal_bins = .false.
      xs_only = .false.
      verbose = .true.
      if(nproc > 1) verbose = .false.
      primary_decay = .false.
      pop_calc_prob = .true.
      pop_calc = .false.
      file_energy_index = .false.
      read_saved_params = .false.
      refresh_library_directories = .true.
      do_dwba = .false.
      use_tran_states = .true.
      use_eval_levels = .true.

      lev_option = -1


      num_bad_samp_e = 0
      num_bad_spread_e = 0

      output_mode = 1
      t12_isomer = 1.0
      write_me = 0
      ex_pop_file(1:20)= ' '
      Max_J_allowed = 20
      Fiss_Max_J = 0.0d0
      rho_cut = 0.25
      T_norm = 1.0d0
      pair_model = 1               !   Default for pairing model
      preeq_delta = 0.0d0
      Preeq_V = 38.0d0
      Preeq_g_div = 15.d0
      optical = 'fresco'
      num_theta_angles = 10
      target%istate = 1
      max_num_gsf = 4
!      max_em_l = 3
      tally_norm = 0.0d0
      if(iproc == 0)print_me = .true.
      mold_cutoff = 50.0d0


!
!----   Start with no optical potentials being set
!----   later they may be set with a choice of an option
!----   but if not set with an option, they will be set to default
!
      spin_proj = 0.0d0
      tot_direct = 0.0d0

      cs_scale = 1.0d0
      cs_units = ' b'

      local_cc_file(1:132)=' '
      exist_cc_file = .false.
      cc_scale = 1.0d0

      trans_avg_l = .false.
      scale_elastic = .false.

      ch_par(0) = '-'
      ch_par(1) = '+'

      itarget = 0
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   useful constants stored in module constants
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      pop_check = 1.0d-6
      pi = 2.0*asin(1.0d0)
      two_pi = 2.0d0*pi
!      e_sq = 1.4399764d0                    !   MeV*fm
      hbar = 6.58211899d-22
!      hbar_c = 197.3269631d0                  !   MeV*fm
!----   Value used in FRESCO
      hbar_c = 197.32705d0                   !   MeV*fm
      fmsq_eq_barn = 0.01d0                  !  1 fm**2 = 0.01 b
      fmsq_eq_mbarn = 10.0d0                 !  1 fm**2 = 10 mb
      barn_eq_fmsq = 100.0d0                 !  1 b = 100 fm**2
      mass_u = 931.494095367d0
!      fine_structure = e_sq/hbar_c
!      fine_structure = 7.2973525376d-3
!----   Value used in FRESCO
      fine_structure = 1.0d0/137.03599d0
      e_sq = hbar_c*fine_structure                   !   MeV*fm

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------   relativistic metric
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      g_metric(0:3,0:3) = 0.0d0
      g_metric(0,0) = 1.0d0
      do i = 1, 3
         g_metric(i,i) = -1.0d0
      end do

!-----------   Cross sections below this threshold are skipped
      cs_threshold = 1.0d-9

      call log_fact(num_fac, fact)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------    Set up Gauss-Laguerre quadrature values
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      n_glag = 100
      allocate(x_glag(n_glag),w_glag(n_glag))
      alf = 0.0d0
      bet = 0.0d0
      call gauss_quad(n_glag, 2, alf, bet, x_glag, w_glag)
      allocate(gauss_laguerre(4))
      gauss_laguerre(1)%num = 100
      allocate(gauss_laguerre(1)%nodes(gauss_laguerre(1)%num))
      allocate(gauss_laguerre(1)%weights(gauss_laguerre(1)%num))
      call gauss_quad(gauss_laguerre(1)%num, 2, alf, bet,                &
                      gauss_laguerre(1)%nodes, gauss_laguerre(1)%weights)
      gauss_laguerre(2)%num = 50
      allocate(gauss_laguerre(2)%nodes(gauss_laguerre(2)%num))
      allocate(gauss_laguerre(2)%weights(gauss_laguerre(2)%num))
      call gauss_quad(gauss_laguerre(2)%num, 2, alf, bet,                &
                      gauss_laguerre(2)%nodes, gauss_laguerre(2)%weights)
      gauss_laguerre(3)%num = 25
      allocate(gauss_laguerre(3)%nodes(gauss_laguerre(3)%num))
      allocate(gauss_laguerre(3)%weights(gauss_laguerre(3)%num))
      call gauss_quad(gauss_laguerre(3)%num, 2, alf, bet,                &
                      gauss_laguerre(3)%nodes, gauss_laguerre(3)%weights)
      gauss_laguerre(4)%num = 10
      allocate(gauss_laguerre(4)%nodes(gauss_laguerre(4)%num))
      allocate(gauss_laguerre(4)%weights(gauss_laguerre(4)%num))
      call gauss_quad(gauss_laguerre(4)%num, 2, alf, bet,                &
                      gauss_laguerre(4)%nodes, gauss_laguerre(4)%weights)
      


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------    Set up Gauss-Legendre quadrature values
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Ang_L_max = 30
      n_gleg = Ang_L_max + 1
      allocate(x_gleg(n_gleg), w_gleg(n_gleg))
      alf = 0.0d0
      bet = 0.0d0
      call gauss_quad(n_gleg, 1, alf, bet, x_gleg, w_gleg)
      allocate(Gauss_Leg(0:Ang_L_max,n_gleg))
      do L = 0, Ang_L_max
         do ix = 1, n_gleg
            Gauss_Leg(L,ix) = poly(L,1,alf,beta,x_gleg(ix))
         end do
      end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------    ixx_max and jxx_max need to be even
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ixx_max = 50
      delta_ix = 2.0d0/real(ixx_max,kind=8)
!---   Fill Legendre polynomial array
      if(.not. allocated(Leg_poly))allocate(Leg_poly(0:Ang_L_max,0:ixx_max))
      alf = 0.0d0
      bet = 0.0d0
      do i = 0, ixx_max
         x = real(i,kind=8)*delta_ix - 1.0d0
         do L = 0, Ang_L_max
            Leg_poly(L,i) = poly(L,1,alf,bet,x)
         end do
      end do

      max_jx_10 = 10
      delta_jx_10 = 2.0d0/real(max_jx_10,kind=8)
      max_jx_20 = 20
      delta_jx_20 = 2.0d0/real(max_jx_20,kind=8)
      max_jx_50 = 50
      delta_jx_50 = 2.0d0/real(max_jx_50,kind=8)
      max_jx_100 = 100
      delta_jx_100 = 2.0d0/real(max_jx_100,kind=8)

      part_lmax = 20

      d_zero = 0.0d0
      d_half = 0.5d0
      d_one = 1.0d0
      d_two = 2.0d0
      d_three = 3.0d0
      d_four = 4.0d0
      factorial(0) = 1.0d0
      do i = 1, 50
         factorial(i) = factorial(i-1)*real(i,kind=8)
      end do
      lfactorial(0) = 0.0d0
      do i = 1, 200
         lfactorial(i) = lfactorial(i-1) + dlog(real(i,kind=8))
      end do

      ifresco_shape = 13


      trans_p_cut = 1.0d-7
      trans_e_cut = 1.0d-14
      prob_cut = 1.0d-7

      num_mc_samp = 1000000
!
!---   initialize iseed to int(pi*10^9)
!
      iseed_64 = 3141592654_int_64
      iseed_32 = 3141592_int_64
!---
!---  "Randomize" iseed_64 with cnt from system clock. Generally, this
!---  will be different for each runs.
!---
      call system_clock(COUNT = cnt, COUNT_RATE = cnt_r, COUNT_MAX = cnt_m)
      cnt = max(mod(cnt,100000),1)
      iseed_64 = iseed_64 + cnt
      iseed_64 = iseed_64 + iproc*31415_int_64
      if(iand(iseed_64,1_int_64) /= 1_int_64)iseed_64 = iseed_64 + 1_int_64
      iseed_64 = -iseed_64

      call system_clock(COUNT = cnt, COUNT_RATE = cnt_r, COUNT_MAX = cnt_m)
      cnt = max(mod(cnt,100000),1)
      iseed_32 = iseed_32 + cnt
      iseed_32 = iseed_32 + iproc*31415_int_32
      if(iand(iseed_32,1_int_32) /= 1_int_32)iseed_32 = iseed_32 + 1_int_32
      iseed_32 = -iseed_32

!----------------------------------------------------------------------
      do i=1,80
         bigblank(i:i)=' '
      end do
      do i=1,132
         blank132(i:i)=' '
      end do
      char_pos='+'
      char_neg='-'
      max_particle(1)=-1
      max_particle(2)=-1
      max_particle(3)=0
      max_particle(4)=0
      max_particle(5)=0
      max_particle(6)=0
      max_part=0
      do k=1,6
         if(max_particle(k) > max_part)max_part=max_particle(k)
      end do

!--------------------------------------------------------------------------------------
      call getenv('HOME',home_path)
      len_home=index(home_path,' ')
      home_path(len_home:len_home)='/'
      call getenv('YAHFC_DATA',data_path)
      len_path=index(data_path,' ')
      data_path(len_path:len_path)='/'
      if(len_path <= 0)then
         if(iproc == 0)then
            write(6,*)'Environment variable for data is not set correctly'
            write(6,*)'Program exiting'
         end if
         call exit_YAHFC(101)
      end if
      KE = 0.0d0

      call particle_data

      do k = 1, 6
         particle(k)%lmax = part_lmax
         particle(k)%do_dwba = .false.
      end do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------    maximum multipole for electric transitions
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      E1_model = 2
      e_l_max = 3
      m_l_max = 3
      primary_setup = .false.
      target%specified = .false.
      fission = .true.
      channels = .true.
      PREEQ_Model = 1
      preeq_pair_model = 1
      preeq_fwell = 0
      analytic_preeq = .false.
      WF_model = 1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   projectile data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      projectile%Z = -1
      projectile%A = -1
      projectile%particle_type = -1
      projectile%e_max = -10000.
      projectile%e_min = 10000.
      projectile%e_step = 0.20d0
      projectile%specified = .false.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------  default values for Pre-equilibrium sq Matrix element
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      M2_C1 = 1.0d0
      M2_C2 = 1.0d0
      M2_C3 = 1.0d0
      M2_Rpp = 1.0d0
      M2_Rnn = 1.5d0
      M2_Rpn = 1.0d0
      M2_Rnp = 1.0d0

      quasi_elastic =.false.
      quasi_e = 0.0d0

!-------------------------------------------------------------------
      i_bind = 0
      compound_setup = .false.
      de = 0.1d0
      mass_file = 'aw'
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!-------    set up buffer for commands to be parsed later
!-------    this offers the ability to enter comands in
!-------    any order whatsoever 
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    Get input commands.                                           +
!------    Store in temporary file to count, store in buffer, order      +
!------    their priority and then run                                   +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      energy_input = .false.
      num_command = 0
      finish = .false.
      if(iproc == 0)then
         open(unit=17,file='YAHFC-commands.txt',status='unknown')
         do while(.not. finish)
            command(1:132) = ' '
            read(5,'(a)',iostat = read_err)command
            if(read_err /= 0)then
               finish = .true.
               cycle
            end if
            if(command(1:1) == '#' .or. command(1:1) == '!' .or. command(1:1) == ' ')cycle
            call parse_string(command, numw, startw, stopw)
            if(numw < 1)cycle
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Convert string to lower case
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            nchar = stopw(1) - startw(1) + 1
            call lower_case_word(nchar,command(startw(1):stopw(1)))
!-----   check if energy input is in commands
            if(command(startw(1):stopw(1)) == 'proj_eminmax')energy_input = .true.
            if(command(startw(1):stopw(1)) == 'proj_e_file')energy_input = .true.
            num_command = num_command + 1
            write(17,'(a132)')command
            if(command(startw(1):stopw(1)) == 'end')then
                finish = .true.
                cycle
            end if
         end do
         close(unit=17)
      end if

!
!----------------------------------------------------------
!----   NEED TO MPI_BCAST num_command here
!-----------------------------------------------------------------------------------------
!
#if(USE_MPI == 1)
      if(nproc > 1)then
         call MPI_Barrier(icomm, ierr)
         call MPI_BCAST(num_command, 1, MPI_INTEGER, 0, icomm, ierr)
         call MPI_BCAST(energy_input, 1, MPI_LOGICAL, 0, icomm, ierr)
      end if
#endif
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Set up buffer for input commands
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(.not. allocated(command_buffer))allocate(command_buffer(num_command+1))
      if(.not. allocated(command_rank))allocate(command_rank(num_command+1))
      command_buffer(1:num_command) = blank132
      command_rank(1:num_command) = 10
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Read commands back in from temporary file
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      open(unit=17,file='YAHFC-commands.txt',status='unknown')
      do i = 1, num_command
         read(17,'(a)')command_buffer(i)
      end do
      close(unit = 17)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Set rank for command
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i = 1, num_command
         command_rank(i) = rank_commands(command_buffer(i))
      end do

!----------------------------------------------------------------------------------------+
!-------   order the commands by rank                                                    +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   ranks need to be adjusted.                                                    +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i = 1, num_command - 1
         do j = i + 1, num_command
            if(command_rank(j) < command_rank(i))then
               itemp = command_rank(i)
               char_temp = command_buffer(i)
               command_rank(i) = command_rank(j)
               command_buffer(i) = command_buffer(j)
               command_rank(j) = itemp
               command_buffer(j) = char_temp
            end if
         end do
      end do
!----------------------------------------------------------------------------------------+
!-------   Check if ran_seed is in the commands (would have highest rank)                +
!-------   If not add it and push other commands down, this is why buffer was allocated  +
!-------   to num_command + 1                                                            +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call parse_string(command_buffer(1), numw, startw, stopw)
      if(command_buffer(1)(startw(1):stopw(1)) /= 'ran_seed')then
         do i = num_command, 1, -1
            command_buffer(i+1) = command_buffer(i)
            command_rank(i+1) = command_rank(i)
         end do
         num_command = num_command + 1
         command_buffer(1)(1:9) = 'ran_seed '
!----------------------------------------------------------------------------------------+
!---  "Randomize" iseed_32 with cnt from system clock. Generally, this                   +
!---  will be different for each runs.                                                   +
!----------------------------------------------------------------------------------------+
         cnt = 0
         if(iproc == 0)then
             call system_clock(COUNT = cnt, COUNT_RATE = cnt_r, COUNT_MAX = cnt_m)
             cnt = max(mod(cnt,10000),1)
             call system_clock(COUNT = cnt1, COUNT_RATE = cnt_r, COUNT_MAX = cnt_m)
             cnt1 = max(mod(cnt,10000),1)
             cnt = cnt*cnt1
         end if
#if(USE_MPI == 1)
         call MPI_Barrier(icomm, ierr)
         call MPI_BCAST(cnt, 1, MPI_INTEGER, 0, icomm, ierr)
#endif
         iseed_32 = cnt
         if(iand(iseed_32,1) /= 1)iseed_32 = iseed_32 + 1
         iseed_start = iseed_32
         write(command_buffer(1)(10:19),'(i10)')iseed_32
         command_rank(1) = 0
      end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Write back out to unit=17 in the order of execution                           +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(iproc == 0)then
         open(unit=17,file='YAHFC-commands.txt',status='unknown')
         write(17,'(''#YAHFC commands in order of execution'')')
         write(17,'(''# Default Random iseed_32 = '',i16)')iseed_32
         write(17,'(''# Default Random iseed_64 = '',i16)')iseed_64
         do i = 1, num_command
            write(17,'(a132)')command_buffer(i)
         end do
         close(unit=17)
      end if

      num_comp=0
      finish=.false.
      iproj = -10
      do i = 1, num_command
         call parse_command(command_buffer(i), finish)
         if(target%specified)istate = target%istate
         if(.not. compound_setup .and. (target%specified .and.        &
               projectile%specified .and. ex_set))then
            iproj = projectile%particle_type

            Z_i = target%Z + projectile%Z
            A_i = target%A + projectile%A
            call set_up_decay_chain(projectile%Z,projectile%A,        &
                                    target%Z,target%A,num_comp)
!----------------------------------------------------------------------------------------+
!------   Now that decay chain is set up, also set up data for each                      +
!------   compound nucleus                                                               +
!----------------------------------------------------------------------------------------+
            overide = 'c'                             !   use discrete states marked with 'u' in RIPL-2 file
            do icomp = 1, num_comp
               Z_f = nucleus(icomp)%Z
               A_f = nucleus(icomp)%A
               if(Z_f == target%Z .and. A_f == target%A)target%icomp = icomp
               itarget = target%icomp
!----------------------------------------------------------------------------------------+
!--------   Read in discrete state data                                                  +
!----------------------------------------------------------------------------------------+
               call get_spectrum(data_path,len_path,                      &
                                 overide,symb(Z_f),Z_f,A_f,icomp)
!----------------------------------------------------------------------------------------+
!--------   level density setup                                                          +
!----------------------------------------------------------------------------------------+
               if(A_f > 20)then
                  lev_option = 1
                  if(A_f > 130) lev_option = 2
               else
                  lev_option = 0
               end if
               call get_lev_den(data_path,len_path,                       &
                                symb(Z_f),Z_f,A_f,icomp)
               if(nucleus(icomp)%fit_ematch)call fit_lev_den(icomp)
!----------------------------------------------------------------------------------------+
!--------   EM strength function parameters                                              +
!----------------------------------------------------------------------------------------+
               call gdr_param(data_path,len_path,icomp)
               call EM_str_param(icomp)
               nucleus(icomp)%fission = .false.
!----------------------------------------------------------------------------------------+
!---------    Fission barrier info                                                       +
!----------------------------------------------------------------------------------------+
               if(fission .and. nucleus(icomp)%Z >= 80)then
                  call Fission_data(data_path,len_path,icomp)
               end if
!----------------------------------------------------------------------------------------+
!--------   Read in saved parameters from data_path(1:len_path)//'Saved-Parameters.txt'  +
!----------------------------------------------------------------------------------------+
               if(read_saved_params)call read_saved_parameters(data_path,len_path,icomp)
            end do
            spin_target = nucleus(itarget)%state(target%istate)%spin
            spin_proj = 0.0
            if(.not. pop_calc)spin_proj = particle(iproj)%spin
            compound_setup=.true.
         end if
      end do
!----------------------------------------------------------------------------------------+
!------    Finished getting input commands                                               +
!----------------------------------------------------------------------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Check if we have sufficient data to proceed                       +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(.not.target%specified)then
         if(iproc == 0)write(6,*)'Target nucleus not specified, quitting'
         call exit_YAHFC(101)
      end if
      if(.not.projectile%specified)then
         if(iproc == 0)write(6,*)'Projectile nucleus not specified, quitting'
         call exit_YAHFC(101)
      end if
      mass_proj = particle(iproj)%mass
      mass_target = nucleus(itarget)%mass + nucleus(itarget)%state(target%istate)%energy
      rel_factor = mass_target/(mass_target + mass_proj)
!-------------------------------------------------------------------------------------
      num_theta = num_theta_angles
      if(xs_only)num_theta = 1
!----------------------------------------------------------------------------------------+
!------  After everything is set up, also check if preequilibrium model parameters       +
!------  makes sense. In particular, the finite well parameter, can't have               +
!------  Preq_V1 > Preeq_V                                                               +
!----------------------------------------------------------------------------------------+
      if(preeq_V <= Preeq_V1)then
         Preeq_V1 = max(Preeq_V - 5.0d0,0.0d0)
         if(iproc == 0)then
            write(6,*)'------------   WARNING  ------------------'
            write(6,*)'Problem with preequilibrium model - Preeq_V1 > Preeq_V'
            write(6,*)'Resetting to make it more physical'
            write(6,*)'Preeq_V1 = max(Preeq_V - 5.0,0.0)'
            write(6,*)'Preeq_V = ',Preeq_V
            write(6,*)'Preeq_V1 = ',Preeq_V1
         end if
      end if
!-------------------------------------------------------------------------------+
!------   First call to random number generator                                 +
!-------------------------------------------------------------------------------+
      ran = random_64(iseed_64)
      ran = random_32(iseed_32)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Make sure we have the correct number of strength function lines, 
!----   check str for the components
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i = 1, num_comp
         do L = 1, e_l_max
            do j = 1, max_num_gsf
               if(abs(nucleus(i)%EL_mode(L)%gsf(j)%sr) >= 1.0d-6)           &
                  nucleus(i)%EL_mode(L)%num_gsf = j
            end do
         end do
      end do
      do i = 1, num_comp
         do L = 1, m_l_max
            do j = 1, max_num_gsf
               if(abs(nucleus(i)%ML_mode(L)%gsf(j)%sr) >= 1.0d-6)           &
                  nucleus(i)%ML_mode(L)%num_gsf = j
            end do
         end do
      end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Check if we have any bad data set up that would                   +
!---------    mess things up.                                                   +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call fail_safe_check(de)

!----------------------------------------------------------------------------------------+
!------    Do clean up and finish setting up nuclei, level density, etc.                 +
!----------------------------------------------------------------------------------------+
      do icomp = 1, num_comp
         ia = nucleus(icomp)%A
!-------------------------------------------------------------------------------
!----    Set EM strength parameters again if they have been adjusted
!----    and some are still default
!-------------------------------------------------------------------------------
         call EM_str_param(icomp)
!----------------------------------------------------------------------------------------+
!-------------    Check if D0exp > 0 and requesting to fit it, if not
!-------------    make flag to fit it .false.
!----------------------------------------------------------------------------------------+
         if(nucleus(icomp)%D0exp <= 0.0d0 .and. nucleus(icomp)%fit_D0)   &
            nucleus(icomp)%fit_D0 = .false.
         call finish_lev_den(icomp)
!----------------------------------------------------------------------------------------+
!----   If fitted cumulative density differs from experiment by more than a
!----   factor of 1.25 or less than a factor of 0.75, warn users to check
!----------------------------------------------------------------------------------------+
        if(nucleus(icomp)%fit_ematch .and. nucleus(icomp)%cum_rho_ratio > 1.3d0)then
           if(iproc == 0)then
              write(6,*)
              write(6,'(''******************************************************************'')')
              write(6,'(''****** CAUTION ** CAUTION ** CAUTION ** CAUTION ** CAUTION  ******'')')
              write(6,'(''*****  The fitted cumulative level density for '',a5,'' at      *****'')')nucleus(icomp)%Label
              write(6,'(''*****  E_cut is more than 1.3 larger than the experimental   *****'')')
              write(6,'(''*****  value. You should check this. Possible causes are:    *****'')')
              write(6,'(''*****   1. Ecut is too high and the experimental spectrum    *****'')')
              write(6,'(''*****      is incomplete up to Ecut, so the modled spectrum  *****'')')
              write(6,'(''*****      is higher than experiment.                        *****'')')
              write(6,'(''*****   2. Level density parameters such as the pairing gap  *****'')')
              write(6,'(''*****      and/or the shell correction are such that the     *****'')')
              write(6,'(''*****      modeled level density is too large and cannot     *****'')')
              write(6,'(''*****      reproduce the experimental spectrum. This can     *****'')')
              write(6,'(''*****      can happen if D0 is not known, and in odd-odd     *****'')')
              write(6,'(''*****      nuclei where the default pairing gap is zero.     *****'')')
              write(6,'(''*****   3. There are too few discrete states to effectively  *****'')')
              write(6,'(''*****      fit the experimental cumulative density.          *****'')')
              write(6,'(''****** CAUTION ** CAUTION ** CAUTION ** CAUTION ** CAUTION  ******'')')
              write(6,'(''******************************************************************'')')
              write(6,*)
           end if
         end if
        if(nucleus(icomp)%fit_ematch .and. nucleus(icomp)%cum_rho_ratio < 0.7d0)then
           if(iproc == 0)then
              write(6,*)
              write(6,'(''****** CAUTION ** CAUTION ** CAUTION ** CAUTION ** CAUTION  ******'')')
              write(6,'(''*****    CAUTION   CAUTION   CAUTION  CAUTION   CAUTION      *****'')')
              write(6,'(''*****  The fitted cumulative level density for '',a5,'' at      *****'')')nucleus(icomp)%Label
              write(6,'(''*****  E_cut is more than 0.7 smaller than the experimental  *****'')')
              write(6,'(''*****  value. You should check this and possibly adjust      *****'')')
              write(6,'(''*****  level density parameters. In this case the            *****'')')
              write(6,'(''*****  modeled level density is too small and cannot be      *****'')')
              write(6,'(''*****  made to match the experimental spectrum.              *****'')')
              write(6,'(''*****  Possibly there are too few experimental states to     *****'')')
              write(6,'(''*****  effectively fit the experimental cumulative density   *****'')')
              write(6,'(''****** CAUTION ** CAUTION ** CAUTION ** CAUTION ** CAUTION  ******'')')
              write(6,'(''******************************************************************'')')
              write(6,*)
           end if
         end if

!----------------------------------------------------------------------------------------+
!----   Check if fission is actually possible, i.e., barriers are below maximum
!----   excitation energies
!----------------------------------------------------------------------------------------+
         if(nucleus(icomp)%fission)then
            call Fission_levels(data_path,len_path,icomp)
            do i = 1, nucleus(icomp)%F_n_barr
               if(nucleus(icomp)%F_Barrier(i)%level_param(6) >                        &
                  nucleus(icomp)%F_Barrier(i)%level_param(3))then
                  call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,          &
                                 nucleus(icomp)%F_Barrier(i)%vib_enh,                 &
                                 nucleus(icomp)%F_Barrier(i)%rot_enh)
               else
                  if(iproc == 0)then
                     write(6,*)'**************************************************'
                     write(6,*)'Attempt to set E_match < Delta'
                     write(6,*)'Nucleus icomp = ',icomp
                     write(6,*)'Z = ',nucleus(icomp)%Z,' A = ',nucleus(icomp)%A
                     write(6,*)'Fission Barrier # ',i
                     write(6,*)'Cannot continue with this calculation'
                     write(6,*)'**************************************************'
                  end if
                  call exit_YAHFC(101)
               end if
            end do
         end if
        if(.not. all_discrete_states)nucleus(icomp)%num_discrete = nucleus(icomp)%ncut
      end do

      do k = 1, 6
         if(particle(k)%in_decay .and. .not. particle(k)%opt_pot_set)then
            particle(k)%opt_pot = 1
            if(k == 1 .and. (nucleus(itarget)%Z >= 90 .and. nucleus(itarget)%Z <= 96)) &
                 particle(k)%opt_pot = 2
            if(k == 2 .and. (nucleus(itarget)%Z >= 90 .and. nucleus(itarget)%Z <= 96)) &
                 particle(k)%opt_pot = 2
            particle(k)%opt_pot_set = .true.
         end if
      end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------    Check to see if we need to overide fission                          +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(.not.fission)then
         do icomp = 1,num_comp
            if(nucleus(icomp)%fission)fission = .false.
         end do
      else
         fission = .false.
         do icomp = 1,num_comp
            if(nucleus(icomp)%fission)fission = .true.
         end do
      end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Set up excitation energy bins for each nucleus                  +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call Setup_Energy_bins(de)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------   Now that energy grids are set up, remap incident projectile   +
!--------   energies to fit on excitation energy grid                     +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(.not. pop_calc)then
         call fix_incident_energies(iproj, rel_factor)
      else
          call fix_pop_energies
      end if
      num_energies = projectile%num_e
!-------------------------------------------------------------------------------+
!----   Setup minimum energies for the optical model calculations               +
!-------------------------------------------------------------------------------+
      call set_min_particle_energy
!-------------------------------------------------------------------------------+
!----   Run FRESCOX to get transmission coefficients and optical model data     +
!-------------------------------------------------------------------------------+
      call optical_setup(data_path, len_path, iproj, itarget, istate, de, Ang_L_max)
!-------------------------------------------------------------------------------+
!----   Check for maximum angular momentum from the reaction dynamics           +
!-------------------------------------------------------------------------------+
      if(.not. pop_calc)call reaction_max_J
      j_max = Max_J_allowed
!-------------------------------------------------------------------------------+
!-----   Set maximum angular momentum for each nucleus                          +
!-------------------------------------------------------------------------------+
      do icomp = 1, num_comp
         nucleus(icomp)%j_max = max_J_allowed
      end do

      if(out_file(1:1) == ' ')out_file = 'YAHFC-Output'

      ifile = index(out_file,' ') - 1
      if(print_output)then
         if(iproc == 0)open(unit=13,file=out_file(1:ifile)//'.out',status='unknown')
      end if

      cum_fit=2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------   Calculate and print Q-values
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(print_me)write(6,*)'Q-values for each channel'
      do i = 1, num_channels
         inuc = Exit_Channel(i)%Final_nucleus
         mass_exit = 0.0d0
         if(explicit_channels)then
            do n = 1, Exit_channel(i)%num_particles
               k = Exit_channel(i)%decay_particles(n)
               mass_exit = mass_exit + particle(k)%mass
            end do
         else
            do k = 1, 6
               mass_exit = mass_exit + particle(k)%mass*Exit_channel(i)%num_part(k)
            end do
         end if
         Q = mass_target + particle(iproj)%mass - nucleus(inuc)%mass - mass_exit
         if(abs(Q) < 1.0d-6)Q = 0.0d0
         ilast = index(Exit_Channel(i)%Channel_Label,' ')
         if(print_me)write(6,'(''channel '',i4,1x,a20,f15.7)')i, Exit_Channel(i)%Channel_Label,Q
         Exit_Channel(i)%Q_value = Q
      end do

      pop_input = .false.
      if(pop_calc)then
         WF_model = 0
         PREEQ_Model = 0
         track_gammas = .true.
         pop_input = .true.
         num_pop = 2*(Max_J_allowed + 1)
         if(.not.allocated(pop_prob))allocate(pop_prob(num_pop))
         if(.not.allocated(pop_j))allocate(pop_j(num_pop))
         if(.not.allocated(pop_ip))allocate(pop_ip(num_pop))
         if(.not.allocated(decay_mult))allocate(decay_mult(0:6,num_energies))
         decay_mult(0:6,num_energies) = 0.0d0
         if(.not.allocated(decay_var))allocate(decay_var(0:6,num_energies))
         decay_var(0:6,num_energies) = 0.0d0
      else
         if(.not.allocated(pop_prob))allocate(pop_prob(1))
         if(.not.allocated(pop_j))allocate(pop_j(1))
         if(.not.allocated(pop_ip))allocate(pop_ip(1))
         if(.not.allocated(decay_mult))allocate(decay_mult(0:6,1))
         decay_mult(0:6,1) = 0.0d0
         if(.not.allocated(decay_var))allocate(decay_var(0:6,1))
         decay_var(0:6,1) = 0.0d0
      end if

!-------------------------------------------------------------------------+
!--------------    Allocate arrays for pre-equilibrium model              +
!-------------------------------------------------------------------------+
     if(.not. pop_calc .and. PREEQ_Model == 1)then
        nucleus(1)%PREEQ = .true.
        pp_max = 6
        pn_max = 6
        p0(2) = projectile%Z
        p0(1) = projectile%A - p0(2)
        Preeq_gam_fact = 1.0d0
!-------------------------------------------------------------------------+
!-------------    Allocate arrays for particle spectra                    +
!-------------------------------------------------------------------------+
        do icomp = 1, num_comp
           Ef = (4.0d0/real(nucleus(icomp)%A,kind=8))*Init_Kinetic_Energy
           Ex = (sqrt(nucleus(icomp)%Ex_max)+sqrt(Ef))**2
           nucleus(icomp)%nbin_part = int(Ex/de) + 1
           if(nucleus(icomp)%PREEQ .and. PREEQ_Model > 0)then
              allocate(nucleus(icomp)%PREEQ_part_spectrum(1:nucleus(icomp)%num_decay,0:nucleus(icomp)%nbin_part))
              nucleus(icomp)%PREEQ_part_spectrum(1:nucleus(icomp)%num_decay,0:nucleus(icomp)%nbin_part) = 0.0d0
              allocate(nucleus(icomp)%PREEQ_cs(1:num_energies))
              nucleus(icomp)%PREEQ_cs(1:num_energies) = 0.0d0
              allocate(nucleus(icomp)%PREEQ_part_cs(1:nucleus(icomp)%num_decay,1:num_energies))
              nucleus(icomp)%PREEQ_part_cs(1:nucleus(icomp)%num_decay,1:num_energies) = 0.0d0
           end if
        end do
        if(.not. allocated(dWk))allocate(dWk(0:pn_max+1,0:pp_max+1,0:nucleus(1)%nbin_part,0:6))
        if(.not. allocated(Wk))allocate(Wk(0:pn_max+1,0:pp_max+1,0:6))
        if(.not. allocated(W))allocate(W(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(tau))allocate(tau(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(taup))allocate(taup(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Pre_Prob))allocate(Pre_Prob(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Gam_p))allocate(Gam_p(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Gam_n))allocate(Gam_n(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Gam_pp))allocate(Gam_pp(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Gam_pn))allocate(Gam_pn(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Gam_0_pn))allocate(Gam_0_pn(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Gam_0_np))allocate(Gam_0_np(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Lamb_p_p))allocate(Lamb_p_p(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Lamb_p_n))allocate(Lamb_p_n(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Lamb_0_pn))allocate(Lamb_0_pn(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Lamb_0_np))allocate(Lamb_0_np(0:pn_max+1,0:pp_max+1))
        if(.not. allocated(Vwell))allocate(Vwell(0:(pn_max+pp_max)))
        if(.not. allocated(Vwell_g))allocate(Vwell_g(0:(pn_max+pp_max)))
     end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------    Allocate arrays for each nucleus to store calculated
!--------    cross sections
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(fission)then
         allocate(fission_cs(1:num_energies))
         fission_cs(1:num_energies) = 0.0d0
      else
         allocate(fission_cs(1))
         fission_cs(1) = 0.0d0
      end if
      if(.not. pop_calc)then
         allocate(direct_tot(1:num_energies))
         allocate(direct_cc(1:num_energies))
         allocate(direct_dwba(1:num_energies))
      end if
      allocate(reaction_cs(1:num_energies))
      allocate(absorption_cs(1:num_energies))
      do icomp = 1, num_comp
         nucleus(icomp)%num_cs = 1
         do i = 2,nucleus(icomp)%num_discrete            ! Loop over discrete states to check for isomers
            if(nucleus(icomp)%state(i)%isomer)nucleus(icomp)%num_cs = nucleus(icomp)%num_cs + 1
         end do
         allocate(nucleus(icomp)%hang_cs(1:num_energies))
         nucleus(icomp)%hang_cs(1:num_energies) = 0.0d0
         if(nucleus(icomp)%fission)then
            allocate(nucleus(icomp)%fission_cs(1:num_energies))
            nucleus(icomp)%fission_cs(1:num_energies) = 0.0d0
         end if
         allocate(nucleus(icomp)%state_lab(nucleus(icomp)%num_cs))
         nucleus(icomp)%state_lab(1) = 1                 !  index of the state being tracked
         icount = 1
         do i = 2, nucleus(icomp)%num_discrete
            if(nucleus(icomp)%state(i)%isomer)then
               icount = icount + 1
               nucleus(icomp)%state_lab(icount) = i
            end if
         end do
      end do
!--------------------------------------------------------------------------+
!------    Now allocate arrays to track fission                            +
!--------------------------------------------------------------------------+
      Ef = (4.0d0/real(nucleus(num_comp)%A,kind=8))*Init_Kinetic_Energy
      Ex = (sqrt(nucleus(1)%Ex_max)+sqrt(Ef))**2
!--------------------------------------------------------------------------+
!-----  Now check for any discrete states have J greater than max_J_allowed or any states
!-----  that decay to a state greater than max_J_allowed. If so, remove them the list of discrete states
!--------------------------------------------------------------------------+
      do i = 1, num_comp
         numd = 0
         do n = 1, nucleus(i)%num_discrete
            J = nint(nucleus(i)%state(n)%spin - nucleus(i)%jshift)
            state_map(n) = -1
            if(J <= max_J_allowed)then            !   J is below max_J_allowed
               check = .true.
               do m = 1, nucleus(i)%state(n)%nbranch
                  ibranch = nucleus(i)%state(n)%ibranch(m)
                  JJ = nint(nucleus(i)%state(ibranch)%spin - nucleus(i)%jshift)
                  if(JJ > max_J_allowed)check = .false.
!--------------------------------------------------------------------------+
!----   Check if the final state has been marked for removal
!----   if so, remove this one too
!--------------------------------------------------------------------------+
                  if(state_map(ibranch) == -1)check = .false.
               end do
               if(check)then
                  numd = numd + 1
                  state_map(n) = numd
               end if
            end if
         end do
         if(numd < nucleus(i)%num_discrete)call remove_states(i,nucleus(i)%num_discrete,state_map)
      end do
!----------------------------------------------------------------------------------------+
!-----  set up and initialze arrays for starting populations
!----------------------------------------------------------------------------------------+
      if(.not.allocated(target%pop_xjpi))allocate(target%pop_xjpi(0:j_max,0:1))

      sp1 = abs(spin_target-spin_proj)
      sp2 = spin_target+spin_proj
      isp = int(sp2-sp1)
!----------------------------------------------------------------------------------------+
!-------   Create level density arrays                                                   +
!----------------------------------------------------------------------------------------+
      do i = 1, num_comp
         nbin = nucleus(i)%nbin
         allocate(nucleus(i)%bins(0:j_max,0:1,1:nbin))
         nucleus(i)%bins(0:j_max,0:1,1:nbin)%rho = 0.0d0
         nucleus(i)%bins(0:j_max,0:1,1:nbin)%pop = 0.0d0
         nucleus(i)%state(1:nucleus(i)%num_discrete)%pop = 0.0d0
!----------------------------------------------------------------------------------------+
!----   Continuous level density with internal models                                    +
!----------------------------------------------------------------------------------------+
         if(.not. nucleus(i)%lev_den_read)call modeled_level_density(i, yrast, yrast_actual)
!----------------------------------------------------------------------------------------+
!----   Continuous level density read in from file                                       +
!----------------------------------------------------------------------------------------+
         if(nucleus(i)%lev_den_read)call read_level_density(i, yrast, yrast_actual)
      end do
!----------------------------------------------------------------------------------------+
!----  Fit to Gamma_gamma                                                                +
!----------------------------------------------------------------------------------------+
      call fit_nuke_Gamma_gamma

      if(.not.pop_calc)then
         e_rel = projectile%energy(num_energies)*rel_factor
         ee_max = e_rel + nucleus(1)%sep_e(iproj) +               &
                          nucleus(itarget)%state(target%istate)%energy
      else
         ee_max = 0.0
         do i = 1, num_pop_e
            if(Pop_data(i)%Ex_pop > ee_max) ee_max = Pop_data(i)%Ex_pop
         end do
      end if

!----------------------------------------------------------------------------------------+
!----   Setup the energy bin for the spectra: de_spec                                    +
!----------------------------------------------------------------------------------------+
      de_spec = de
      if(de <= 1.0d0)then
         ii = nint(1.0d0/de)
         de_spec = 1.0d0/real(ii,kind=8)
      else
         de_spec = real(nint(de),kind=8)
      end if

      de_spec2 = de_spec

      num_e = int(ee_max/de_spec) + 3


      if(print_me)write(6,'(''de = '',f16.9,''    de_spec = '',f16.9,''    num = '',i6)')de, de_spec, num_e

      if(.not. pop_calc)then
          if(PREEQ_Model > 0)then
             allocate(preeq_css(0:6,1:num_energies))
             preeq_css(0:6,1:num_energies) = 0.0d0
             allocate(preeq_spect(0:6,0:num_e))
             allocate(preeq_spect_full(0:6,0:num_e))
          else
             allocate(preeq_css(0:6,1))
             preeq_css(0:6,1:1) = 0.0d0
             allocate(preeq_spect(0:6,1))
             preeq_spect(0:6,1) = 0.0d0
             allocate(preeq_spect_full(0:6,1))
             preeq_spect_full(0:6,1) = 0.0d0
          end if
       else
          allocate(preeq_css(0:6,1))
          preeq_css(0:6,1) = 0.0d0
          allocate(preeq_spect(0:6,1))
          preeq_spect(0:6,1) = 0.0d0
          allocate(preeq_spect_full(0:6,1))
          preeq_spect_full(0:6,1) = 0.0d0
       end if
!-------------------------------------------------------------------------+
!------                                                                   +
!-------   Set up Hauser-Feshbach denominators                            +
!------                                                                   +
!-------------------------------------------------------------------------+
       if(iproc == 0)write(6,*)'Calculating HF probabilities'
       do i = 1, num_comp
          if(print_me)write(6,*)'HF-denominators for nucleus #',i
          call HF_denominator(i)
       end do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Now set up array tracking particles in the decay chain               +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      max_gammas = 0
      do icomp = 1, num_comp
         if(nucleus(icomp)%num_discrete > max_gammas) max_gammas = nucleus(icomp)%num_discrete
      end do
      dim_part = nucleus(1)%nbin + max_gammas*2

      nang_tot = num_theta
      if(.not. allocated(part_data))allocate(part_data(n_dat,dim_part))
      if(.not. allocated(part_Ang_data))allocate(part_Ang_data(0:Ang_L_max,dim_part))
      if(.not. allocated(extra_angle_data))allocate(extra_angle_data(3*num_theta,dim_part))

      if(.not.pop_calc .and. iproj > 0)then
         l_max = 0
         do k = 1, 6
            if(max_particle(k) > 0 .and. particle(k)%lmax > l_max)l_max = particle(k)%lmax
         end do
         if(.not.allocated(clb_l))allocate(clb_l(0:Ang_L_max,0:l_max))
         do L_Ang = 0, Ang_L_max
            xL_ang = real(L_Ang,kind=8)
            do l_i = 0, l_max
               xl_i = real(l_i,kind=8)
               clb_l(L_Ang,l_i) =  clebr(xl_i,0.0d0,xl_i,0.0d0,xL_ang,0.0d0)
            end do
         end do
      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Arrays tracking explicit Exit Channel Data                           +
!     Fill in remaining info not done in subroutine decay_chain            +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i = 1, num_channels
         inuc = Exit_Channel(i)%Final_nucleus
         num_s = 1
         do n = 2, nucleus(inuc)%num_discrete
            if(nucleus(inuc)%state(n)%isomer)num_s = num_s + 1
         end do
         Exit_Channel(i)%num_cs = num_s
         if(.not.allocated(Exit_Channel(i)%Channel_cs))                                            &
            allocate(Exit_Channel(i)%Channel_cs(-2:num_s,1:num_energies))
         Exit_Channel(i)%Channel_cs(-2:num_s,1:num_energies) = 0.0d0
         if(.not.allocated(Exit_Channel(i)%StateLabel))                                            &
            allocate(Exit_Channel(i)%StateLabel(num_s))
         Exit_Channel(i)%StateLabel(1) = 1
         num_s = 1
         do n = 2, nucleus(inuc)%num_discrete
            if(nucleus(inuc)%state(n)%isomer)then
               num_s = num_s + 1
               Exit_Channel(i)%StateLabel(num_s) = n
            end if
         end do

         if(track_gammas)then
            if(.not. allocated(Exit_channel(i)%state))                                             &
               allocate(Exit_channel(i)%state(nucleus(inuc)%num_discrete))
            do n = 1, nucleus(inuc)%num_discrete
               if(.not. allocated(Exit_channel(i)%state(n)%cs))                                    &
                  allocate(Exit_channel(i)%state(n)%cs(nucleus(inuc)%state(n)%nbranch,num_energies))
               Exit_channel(i)%state(n)%cs(1:nucleus(inuc)%state(n)%nbranch,1:num_energies) = 0.0d0
            end do
         end if

         if(.not.allocated(Exit_Channel(i)%part_mult))                                             &
                 allocate(Exit_Channel(i)%part_mult(0:6,-1:num_s))

         if(.not.allocated(Exit_Channel(i)%Spect))                                                 &
                 allocate(Exit_Channel(i)%Spect(0:6,-1:num_s))
         do j = -1, num_s
            do k = 0, 6
               if(.not.allocated(Exit_Channel(i)%Spect(k,j)%E_spec))                               &
                  allocate(Exit_Channel(i)%Spect(k,j)%E_spec(0:num_e))
!----
               if(.not.allocated(Exit_Channel(i)%Spect(k,j)%E_count))                              &
                  allocate(Exit_Channel(i)%Spect(k,j)%E_count(0:num_e))
!----
               if(.not.allocated(Exit_Channel(i)%Spect(k,j)%Ang_Dist))                             &
                  allocate(Exit_Channel(i)%Spect(k,j)%Ang_Dist(0:max_jx_10))
!----
               if(.not.allocated(Exit_Channel(i)%Spect(k,j)%Ang_L))                                &
                  allocate(Exit_Channel(i)%Spect(k,j)%Ang_L(0:Ang_L_max))
!----
               if(.not.allocated(Exit_Channel(i)%Spect(k,j)%E_Ang_Dist))                           &
                   allocate(Exit_Channel(i)%Spect(k,j)%E_Ang_Dist(0:max_jx_10,0:num_e))
!----
               if(.not.allocated(Exit_Channel(i)%Spect(k,j)%E_Ang_L))                              &
                   allocate(Exit_Channel(i)%Spect(k,j)%E_Ang_L(0:Ang_L_max,0:num_e))
!----
               if(.not.allocated(Exit_Channel(i)%Spect(k,j)%E_Ang_L_max))                          &
                   allocate(Exit_Channel(i)%Spect(k,j)%E_Ang_L_max(0:num_e))
            end do
         end do
      end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Arrays used to cache transmission coefficients used in               +
!     HF_primary_decay_setup and to compute the Molauer product            +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(.not. pop_calc)then
         do if1 = 1, nucleus(1)%num_decay                                      !  loop over nuclei in the decay chain
            i_f = nucleus(1)%decay_to(if1)
            k = nucleus(1)%decay_particle(if1)
            isp_f_max = nint(2.0d0*particle(k)%spin)
            if(k == 0)cycle
            if(.not. allocated(particle(k)%trans_bin))                                             &
               allocate(particle(k)%trans_bin(0:isp_f_max,0:particle(k)%lmax,nucleus(i_f)%nbin))
            if(.not. allocated(particle(k)%trans_discrete))                                        &
               allocate(particle(k)%trans_discrete(0:isp_f_max,0:particle(k)%lmax,nucleus(i_f)%num_discrete))
         end do
      end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Arrays tracking Inelastic and Compound Elastic Channels              +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(.not. pop_calc)then
         if(.not.allocated(Inelastic_cont_cs))                                                     &
             allocate(Inelastic_cont_cs(0:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_cont_cs(0:nucleus(itarget)%num_discrete,1:num_energies) = 0.0d0
         if(.not.allocated(Inelastic_cs))                                                          &
             allocate(Inelastic_cs(0:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_cs(0:nucleus(itarget)%num_discrete,1:num_energies) = 0.0d0
         if(.not.allocated(Inelastic_count))                                                       &
             allocate(Inelastic_count(0:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_count(0:nucleus(itarget)%num_discrete,1:num_energies) = 0
!-------------
         if(.not.allocated(Inelastic_L_max))                                                       &
             allocate(Inelastic_L_max(0:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_L_max(1:nucleus(itarget)%num_discrete,1:num_energies) = Ang_L_max
         if(.not.allocated(Inelastic_Ang_L))                                                       &
            allocate(Inelastic_Ang_L(0:Ang_L_max,0:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_Ang_L(0:Ang_L_max,0:nucleus(itarget)%num_discrete,1:num_energies) = 0.0d0
!-------------
         if(iproj > 0)then
            if(.not.allocated(direct_Spectrum))allocate(direct_Spectrum(0:num_e))
            if(.not.allocated(dwba_Spectrum))allocate(dwba_Spectrum(0:num_e))
            dwba_Spectrum(0:num_e) = 0.0d0
            if(quasi_elastic)then
               if(.not.allocated(QE_cs))allocate(QE_cs(1:num_energies))
               QE_cs(1:num_energies) = 0.0d0
               if(.not.allocated(QE_Spectrum))allocate(QE_Spectrum(0:num_e))
               if(.not.allocated(QE_Ang))allocate(QE_Ang(0:Ang_L_max,1:num_energies))
               QE_ang(0:Ang_L_max,1:num_energies)=0.0d0
            end if
            if(.not.allocated(Elastic_cs))allocate(Elastic_cs(1:num_energies))
            Elastic_cs(1:num_energies) = 0.0d0
            if(.not.allocated(Elastic_Ang))allocate(Elastic_ang(0:Ang_L_max,1:num_energies))
            Elastic_ang(0:Ang_L_max,1:num_energies)=0.0d0
!--------------------------------------------------------------------------------------------------
!------    Now arrays for excitations to each discrete state
!--------------------------------------------------------------------------------------------------
!------    Direct excitation
!--------------------------------------------------------------------------------------------------
            if(.not.allocated(direct_cs))allocate(direct_cs(1:OpticalCS%numcc))
            if(.not.allocated(direct_Ang))allocate(direct_Ang(0:Ang_L_max,1:OpticalCS%numcc))
            if(.not.allocated(direct_prob))allocate(direct_prob(0:ixx_max,1:OpticalCS%numcc))
            direct_cs(1:OpticalCS%numcc) = 0.0d0
            direct_Ang(0:Ang_L_max,OpticalCS%numcc) = 0.0d0
            direct_prob(0:ixx_max,1:OpticalCS%numcc) = 0.0d0

            if(.not.allocated(SE_cs))allocate(SE_cs(1:num_energies))
            SE_cs(1:num_energies) = 0.0d0
            if(.not.allocated(SE_Ang))allocate(SE_ang(0:Ang_L_max,1:num_energies))
            SE_ang(0:Ang_L_max,1:num_energies)=0.0d0
            if(.not.allocated(SE_prob))allocate(SE_prob(0:ixx_max))
         else            !    not used with photo-absorption, but prevents warnings
            if(.not.allocated(direct_Spectrum))allocate(direct_Spectrum(1))
            direct_Spectrum(1) = 0.0d0
            if(.not.allocated(dwba_Spectrum))allocate(dwba_Spectrum(1))
            dwba_Spectrum(1) = 0.0d0
            if(.not.allocated(SE_cs))allocate(SE_cs(1))
            SE_cs(1) = 0.0d0
            if(.not.allocated(SE_Ang))allocate(SE_ang(1,1))
            SE_ang(1,1)=0.0d0
            if(.not.allocated(SE_prob))allocate(SE_prob(1))
            SE_prob(1) = 0.0d0
            if(.not.allocated(direct_cs))allocate(direct_cs(1))
            if(.not.allocated(direct_Ang))allocate(direct_Ang(1,1))
            if(.not.allocated(direct_prob))allocate(direct_prob(1,1))
            direct_cs(1) = 0.0d0
            direct_Ang(1,1) = 0.0d0
            direct_prob(1,1) = 0.0d0
            if(.not.allocated(Elastic_cs))allocate(Elastic_cs(1))
            Elastic_cs(1) = 0.0d0
            if(.not.allocated(Elastic_Ang))allocate(Elastic_ang(1,1))
            Elastic_ang(1,1)=0.0d0
         end if

         if(.not.allocated(CE_cs))allocate(CE_cs(1:num_energies))
         CE_cs(1:num_energies) = 0.0d0
         if(.not.allocated(CE_Ang))allocate(CE_ang(0:Ang_L_max,1:num_energies))
         CE_ang(0:Ang_L_max,1:num_energies)=0.0d0
!--------------------------------------------------------------------------------------------------
!------   Arrays for storing results from inelasitic scattering
!--------------------------------------------------------------------------------------------------
!------    First is total of inelastic as a function of incident energy for each type,
!------    direct, compound to a discrete state, and compound through Continuous bins
!--------------------------------------------------------------------------------------------------
         if(.not.allocated(Inelastic_Total))allocate(Inelastic_Total(1:num_energies))
         inelastic_Total(1:num_energies) = 0.0d0
!--------------------------------------------------------------------------------------------------
!-----    Series of dummy allocations for variables only used for a reaction calculation
!-----    allocated because the gfortran compiler gives "erroneous" warnings of uninitialized
!-----    stride and offset for these variables later on. They are accessed within if blocks
!-----    specifying that it is NOT a population calculation (pop_calc = .false.), but this is
!-----    conditional check seems to be lost. Thus, allocating to size one, and initializing to zero.
!-----    This eliminates the warnings, which will allow for other checks if a variable is not
!-----    properly initialized.
!--------------------------------------------------------------------------------------------------
      else
         if(.not.allocated(Inelastic_cont_cs))allocate(Inelastic_cont_cs(1,1))
         Inelastic_cont_cs(1,1) = 0.0d0
         if(.not.allocated(Inelastic_cs))allocate(Inelastic_cs(1,1))
         Inelastic_cs(1,1) = 0.0d0
         if(.not.allocated(Inelastic_count))allocate(Inelastic_count(1,1))
         Inelastic_count(1,1) = 0
!-------------
         if(.not.allocated(Inelastic_L_max))allocate(Inelastic_L_max(1,1))
         Inelastic_L_max(1,1) = 0
         if(.not.allocated(Inelastic_Ang_L))allocate(Inelastic_Ang_L(1,1,1))
         Inelastic_Ang_L(1,1,1) = 0.0d0
!-------------
         if(.not.allocated(direct_Spectrum))allocate(direct_Spectrum(1))
         direct_Spectrum(1) = 0.0d0
         if(.not.allocated(dwba_Spectrum))allocate(dwba_Spectrum(1))
         dwba_Spectrum(1) = 0.0d0
         if(.not.allocated(QE_cs))allocate(QE_cs(1:num_energies))
         QE_cs(1:num_energies) = 0.0d0
         if(.not.allocated(QE_Spectrum))allocate(QE_Spectrum(1))
         if(.not.allocated(QE_Ang))allocate(QE_Ang(1,1))
         QE_ang(1,1)=0.0d0
         if(.not.allocated(Elastic_cs))allocate(Elastic_cs(1))
         Elastic_cs(1) = 0.0d0
         if(.not.allocated(Elastic_Ang))allocate(Elastic_ang(1,1))
         Elastic_ang(1,1)=0.0d0

         if(.not.allocated(CE_cs))allocate(CE_cs(1))
         CE_cs(1) = 0.0d0
         if(.not.allocated(CE_Ang))allocate(CE_ang(1,1))
         CE_ang(1,1)=0.0d0

         if(.not.allocated(SE_cs))allocate(SE_cs(1))
         SE_cs(1) = 0.0d0
         if(.not.allocated(SE_Ang))allocate(SE_ang(1,1))
         SE_ang(1,1)=0.0d0
         if(.not.allocated(SE_prob))allocate(SE_prob(1))
         SE_prob(1) = 0.0d0
!--------------------------------------------------------------------------------------------------
!------   Arrays for storing results from inelasitic scattering
!--------------------------------------------------------------------------------------------------
!------    First is total of inelastic as a function of incident energy for each type,
!------    direct, compound to a discrete state, and compound through Continuous bins
!--------------------------------------------------------------------------------------------------
         if(.not.allocated(Inelastic_Total))allocate(Inelastic_Total(1))
         inelastic_Total(1) = 0.0d0
!--------------------------------------------------------------------------------------------------
!------    Now arrays for excitations to each discrete state
!--------------------------------------------------------------------------------------------------
!------    Direct excitation
!--------------------------------------------------------------------------------------------------
         if(.not.allocated(direct_cs))allocate(direct_cs(1))
         if(.not.allocated(direct_Ang))allocate(direct_Ang(1,1))
         if(.not.allocated(direct_prob))allocate(direct_prob(1,1))
         direct_cs(1) = 0.0d0
         direct_Ang(1,1) = 0.0d0
         direct_prob(1,1) = 0.0d0
      end if


      if(.not.allocated(x_particle_cs))allocate(x_particle_cs(1:num_energies,0:6))
      if(.not.allocated(x_particle_Spectrum))allocate(x_particle_Spectrum(0:num_e,0:6))
      x_particle_cs(1:num_energies,0:6) = 0.0d0

      if(track_primary_gammas)then
         allocate(Primary_Gamma_Spectrum(0:num_e))
      else
         allocate(Primary_Gamma_Spectrum(0))
      end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Create core file name for output                               +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call nucleus_label(itarget, ntar, target_label)
      file_name(1:40) = ' '
      ifile = ntar
      if(.not. pop_calc)then
         k = projectile%particle_type
         file_name(1:ntar+3) = particle(k)%label//'+'//target_label(1:ntar)//'-'
         ifile = ntar + 3
      else
         file_name(1:ntar+1) = target_label(1:ntar)//'-'
         ifile = ntar + 1
      end if
      core_file = ifile

      if(iproc == 0)call check_directories(ntar, target_label, ilib_dir, lib_dir)

#if(USE_MPI == 1)
      if(nproc > 1)then
         call MPI_BCAST(ilib_dir, 1, MPI_INTEGER, 0, icomm, ierr)
         call MPI_BCAST(lib_dir, 132, MPI_CHARACTER, 0, icomm, ierr)
      end if
#endif
!----------------------------------------------------------------------------------------+
!-----     Print data to .out file                                                       +
!----------------------------------------------------------------------------------------+
      if(iproc == 0 .and. print_output)then
         call start_IO
         call output_trans_coef
!----------------------------------------------------------------------------------------+
!---------   Print out data used in the calculation for each                             +
!---------   compound nucleus to the .out file                                           +
!----------------------------------------------------------------------------------------+
         call output_nucleus_data(j_max, itarget)
      end if
!----------------------------------------------------------------------------------------+
!---------   Print out data related to decays for each nucleus                                                            +
!----------------------------------------------------------------------------------------+
     if(iproc == 0)call print_nuke_data(ilib_dir, lib_dir)
!
!**************************************************************************
!------    Write particle properties to file                           ---*
!**************************************************************************
!
      if(iproc == 0)then
         open(unit= 100, file = lib_dir(1:ilib_dir)//'/Particle-properties.dat', status = 'unknown')
         do k = -1, 6
            write(100,'(''#**************************'')')
            write(100,'(''# '',a)')particle(k)%name
            write(100,'(''Z = '',i3)')particle(k)%Z
            write(100,'(''A = '',i3)')particle(k)%A
            write(100,'(''J = '',f4.1)')particle(k)%spin
            write(100,'(''Parity = '', i2)')nint(particle(k)%par)
            write(100,'(''Mass = '',1pe16.7,'' amu'')')particle(k)%mass/mass_u
         end do
         close(unit = 100)
      end if

!-------------------------------------------------------------------------+
!---------                                                                +
!---------   Set up loop over incident energies                           +
!---------                                                                +
!-------------------------------------------------------------------------+
!----   removing fission check to make warning that stride might be initialized
      if(fission)then
         if(.not. allocated(Fiss_J_avg))allocate(Fiss_J_avg(num_comp,num_energies))
         Fiss_J_avg(1:num_comp,1:num_energies) = 0.0d0
         if(.not. allocated(Fiss_J_var))allocate(Fiss_J_var(num_comp,num_energies))
         Fiss_J_var(1:num_comp,1:num_energies) = 0.0d0
         if(.not. allocated(Fiss_tally))allocate(Fiss_tally(num_comp,num_energies))
         Fiss_tally(1:num_comp,1:num_energies) = 0.0d0
      else
         if(.not. allocated(Fiss_J_avg))allocate(Fiss_J_avg(1,1))
         Fiss_J_avg(1,1) = 0.0d0
         if(.not. allocated(Fiss_J_var))allocate(Fiss_J_var(1,1))
         Fiss_J_var(1,1) = 0.0d0
         if(.not. allocated(Fiss_tally))allocate(Fiss_tally(1,1))
         Fiss_tally(1,1) = 0.0d0
      end if

!----------------------------------------------------------------------------------------+
!---------   Print out memory used in this calculation                                   +
!----------------------------------------------------------------------------------------+
      call memory_used

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Start loop over initial Energies
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      if(iproc == 0)write(6,*)'Starting energy loop with ',num_energies,' values'

      do in = 1, num_energies
         hang = 0
         sum_ce = 0
         sum_ie = 0

         E_in = projectile%energy(in)

         do i = 1, num_channels
            do j = -1, Exit_Channel(i)%num_cs
               do k = 0, 6
                  Exit_Channel(i)%part_mult(k,j) = 0.0d0
                  Exit_Channel(i)%Spect(k,j)%E_spec(0:num_e) = 0.0d0
                  Exit_Channel(i)%Spect(k,j)%E_count(0:num_e) = 0
                  Exit_Channel(i)%Spect(k,j)%Ang_Dist(0:max_jx_10) = 0.0d0
                  Exit_Channel(i)%Spect(k,j)%Ang_L(0:Ang_L_max) = 0.0d0
                  Exit_Channel(i)%Spect(k,j)%Ang_L_Max = 0
                  Exit_Channel(i)%Spect(k,j)%E_Ang_Dist(0:max_jx_10,0:num_e) = 0.0d0
                  Exit_Channel(i)%Spect(k,j)%E_Ang_L(0:Ang_L_max,0:num_e) = 0.0d0
                  Exit_Channel(i)%Spect(k,j)%E_Ang_L_Max(0:num_e) = 0
               end do
            end do
         end do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Reset primary populations
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do i = 1, num_comp
            nbin = nucleus(i)%nbin
            j_max = nucleus(i)%j_max
            nucleus(i)%bins(0:j_max,0:1,1:nbin)%pop = 0.0d0
            do j = 1, nucleus(i)%num_discrete
               nucleus(i)%state(j)%pop = 0.0d0
            end do
         end do

         if(.not. pop_calc)then
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------    Can't do energy in center-of-momentum frame just yet as there is a small mismatch
!--------------    that leads to small neg energies -1.d-4. Use non-relativistic approximation to get
!--------------    initial excitation energy. OK to 1e-4 MeV which is smaller than binning.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            e_rel = E_in*mass_target/(mass_target + mass_proj)
            ex_tot = e_rel + nucleus(1)%sep_e(iproj) +                          &
                             nucleus(itarget)%state(target%istate)%energy

            if(ex_tot < nucleus(1)%e_grid(1))then
                if(iproc == 0)write(6,*)'This incident energy would populate discrete states only - skipping'
                cycle
            end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Set up process to boost emitted particles to the correct frame
!------   For a reaction with incident particle keep track of relative energy
!------   Hence lab frame and COM are not the same
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            KE = e_in - e_rel
            mass_2 = nucleus(1)%Mass + e_rel
            pp_res = sqrt(2.0d0*mass_2*KE)
            p_res(0) = mass_2 + KE
            pp_res = sqrt(p_res(0)**2 - mass_2**2)
            p_res(1) = pp_res
            p_res(2) = 0.0d0
            p_res(3) = 0.0d0
            do i = 1, 3
               v_res(i) = p_res(i)/mass_2
            end do
            beta = pp_res/mass_2
            gamma = 1.0d0/sqrt(1.0d0 - beta**2)
            Boost_Lab(0:3,0:3) = 0.0d0
            do i = 0, 3
               Boost_Lab(i,i) = 1.0d0
            end do
            Boost_Lab(0,0) = gamma
            Boost_Lab(0,1) = -gamma*v_res(1)
            Boost_Lab(1,0) = Boost_Lab(0,1)
            Boost_Lab(1,1) = 1.0d0 + (gamma - 1.0d0)*v_res(1)**2/beta**2
            if(print_me)then
               write(6,*)
               write(6,'(''-----------------------------------------------------------------------------'')')
               write(6,'(''*****************************************************************************'')')
               write(6,'(''-----------------------------------------------------------------------------'')')
               write(6,'(''Incident energy ='',1x,1pe15.7,'' MeV'')')E_in
               write(6,'(''Center-of-mass E ='',1x,1pe15.7,'' MeV'')')e_rel
               write(6,'(''Delta_E ='',1x,1pe15.7,'' MeV'')')de
            elseif(iproc == 0 .and. print_output)then
               write(13,*)
               write(13,'(''-----------------------------------------------------------------------------'')')
               write(13,'(''*****************************************************************************'')')
               write(13,'(''-----------------------------------------------------------------------------'')')
               write(13,'(''Incident energy ='',1x,1pe15.7,'' MeV'')')E_in
               write(13,'(''Center-of-mass E ='',1x,1pe15.7,'' MeV'')')e_rel
               write(13,'(''Delta_E ='',1x,1pe15.7,'' MeV'')')de
            end if
         else
            if(print_me)then
               write(6,*)
               write(6,'(''-----------------------------------------------------------------------------'')')
               write(6,'(''*****************************************************************************'')')
               write(6,'(''-----------------------------------------------------------------------------'')')
               write(6,'(''Initial energy ='',1x,1pe15.7,'' MeV'')')E_in
               write(6,'(''Delta_E ='',1x,1pe15.7,'' MeV'')')de
            elseif(iproc == 0 .and. print_output)then
               write(13,*)
               write(13,'(''-----------------------------------------------------------------------------'')')
               write(13,'(''*****************************************************************************'')')
               write(13,'(''-----------------------------------------------------------------------------'')')
               write(13,'(''Initial energy ='',1x,1pe15.7,'' MeV'')')E_in
               write(13,'(''Delta_E ='',1x,1pe15.7,'' MeV'')')de
            end if
         end if
!----------------------------------
!         nucleus(1:num_comp)%Fiss_cs = 0.0d0
         nbin = find_ibin(ex_tot,1)

         if(fission)fission_cs(in) = 0.0d0

         x_particle_Spectrum(0:num_e,0:6) = 0.0d0
         if(.not. pop_calc .and. iproj > 0)then
            if(PREEQ_Model > 0)then
               preeq_Spect(0:6,0:num_e) = 0.0d0
               preeq_Spect_full(0:6,0:num_e) = 0.0d0
            end if
            direct_Spectrum(0:num_e) = 0.0d0
            dwba_Spectrum(0:num_e) = 0.0d0
         end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------    Compound nuclear cross section
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(.not. pop_calc)then
            reaction_cs(in) = 0.0d0
            absorption_cs(in) = 0.0d0
            if(iproj > 0)then                  !*******   Particle projectile
               call compound_xs(e_in,itarget,istate,iproj,absorption_cs(in),  &
                                num_channel,channel_prob,ichannel)
               if(num_channel == 0)cycle

               if(print_me)write(6,'(''Reaction cross section = '',1x,1pe15.7,1x,a2)')absorption_cs(in)*cs_scale,cs_units

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------  Setup decay first nucleus in the chain. Apply width-fluctations    +
!------  if necessary. Thus, this nucleus is treated slightly differnet     +
!------  than others in the chain.                                          +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               if(WF_model == 1)then
!                  alf = 0.0d0
!                  bet = 0.0d0
!                  if(e_in < 0.5d0)then
!                     n_glag = 100
!                  elseif(e_in > 0.5d0 .and. e_in < 1.0d0)then
!                     n_glag = 50
!                  elseif(e_in >= 1.0d0 .and. e_in < 2.0d0)then
!                     n_glag = 25
!                  elseif(e_in >= 2.0d0)then
!                     n_glag = 10
!                  end if
!                  n_glag = 10
!                  call gauss_quad(n_glag, 2, alf, bet, x_glag, w_glag)
!               end if

               call HF_primary_decay_setup(e_in,iproj,itarget,1,istate,ex_tot,absorption_cs(in))
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    We have both the channel formation probability and all the decay
!------    probabilities. There are rare occasions with very low formation
!------    probability that also have no allowable decays. This can be caused by 
!------    a high spin and a set of discrete states where particle and gamma 
!------    emission are strongly hindered, or forbidden. Thus, there are no 
!------    primary decays for that channel. This creates a serious computational
!------    problem and is difficult to fix with the forced decay, etc.
!------    Here, the decays are checked and if there are no allowable decays,
!------    the channel is removed, by zeroing out its probability.
               pnorm = 0.0d0
               do i = num_channel, 1, -1
                  if(i > 1)channel_prob(i) = channel_prob(i) - channel_prob(i-1)
                  l_i = ichannel(1,i)
                  is_i = ichannel(2,i)
                  Ix_i = ichannel(3,i)
                   ip_i = ichannel(4,i)
                  if(channel(l_i,is_i,Ix_i)%num_decay == 0)channel_prob(i) = 0.0d0
                  pnorm = pnorm + channel_prob(i)
               end do
               channel_prob(1) = channel_prob(1)/pnorm
               do i = 2, num_channel
                  channel_prob(i) = channel_prob(i-1) + channel_prob(i)/pnorm
               end do          
 
            elseif(iproj == 0)then           !***********    Photo-absorption

               call photo_xs(e_in, itarget, istate, absorption_cs(in),           &
                             num_channel, channel_prob, ichannel)
            end if
         elseif(pop_calc)then              !*******    Population calculation
             reaction_cs(in) = Pop_data(in)%Norm_pop
             absorption_cs(in) = reaction_cs(in)
         end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Now set up arrays for pre-equilibrium emission                    +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         preeq_prob = 0.0d0
         if(.not. pop_calc .and. PREEQ_Model > 0 .and.                   &
            absorption_cs(in) >= cs_threshold .and. iproj > 0)then
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Delete pre-equilibrium model data so that we can start afresh with      +
!-------   new energy.                                                             +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            preeq_css(0:6,in) = 0.0d0
            preeq_spect(0:6,0:num_e) = 0.0d0
            preeq_spect_full(0:6,0:num_e) = 0.0d0
            preeq_prob = 0.0d0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Call pre-equilibrium model                                               +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            call pre_equilibrium_1(1,istate,in,e_in,ex_tot,de,absorption_cs(in))
            preeq_prob = nucleus(1)%PREEQ_cs(in)/absorption_cs(in)
         end if

         if(pop_calc)then
            if(.not. j_pop_calc)then
               pop_sum = 0.0d0
               do i = 1, Pop_data(in)%num_pop
                  pop_sum = pop_sum + Pop_data(in)%bin_pop(i)
               end do
               pop_prob(1:Pop_data(in)%num_pop) = 0.0d0
               pop_prob(1) = Pop_data(in)%bin_pop(1)/pop_sum
               target%pop_xjpi(0:j_max,0:1) = 0.0d0
               j =  nint(Pop_data(in)%j_pop(1) - nucleus(1)%jshift)
               ip = nint((Pop_data(in)%par_pop(1)+1.0d0)/2.0d0)
               pop_j(1) = j
               pop_ip(1) = ip
               target%pop_xjpi(j,ip) = Pop_data(in)%bin_pop(1)/pop_sum

               do i = 2, Pop_data(in)%num_pop
                  pop_prob(i) = pop_prob(i-1) + Pop_data(in)%bin_pop(i)/pop_sum
                  j =  nint(Pop_data(in)%j_pop(i) - nucleus(1)%jshift)
                  ip = nint((Pop_data(in)%par_pop(i)+1.0d0)/2.0d0)

                  pop_j(i) = j
                  pop_ip(i) = ip
                  target%pop_xjpi(j,ip) = Pop_data(in)%bin_pop(i)/pop_sum
               end do

               if(iproc == 0 .and. print_output)then
                  write(13,'(''*************************************************************'')')
                  write(13,'(''Initial populations for incident energy ='',1x,1pe15.7)')e_in
                  do ip = 0, 1
                     do j = 0, nucleus(1)%j_max
                        xj = j + nucleus(1)%jshift
                        write(13,'(2(1x,f4.1),2(1x,f10.6))')xj,2.0*real(ip)-1.0,target%pop_xjpi(j,ip)
                     end do
                  end do
               end if
            else
               pop_sum = 0.0d0
               do j = 0, nucleus(1)%j_max
                  do ip = 0, 1
                     pop_sum = pop_sum + nucleus(1)%bins(j,ip,nbin)%rho
                  end do
               end do
               ijj = 0
               pop_prob(1:num_pop) = 0.0d0
               prev = 0.0d0
               do j = 0, nucleus(1)%j_max
                  do ip = 0, 1
                     ijj = ijj + 1
                     pop_prob(ijj) = prev + nucleus(1)%bins(j,ip,nbin)%rho/pop_sum
                     prev = pop_prob(ijj)
                     pop_j(ijj) = j
                     pop_ip(ijj) = ip
                     target%pop_xjpi(j,ip) = nucleus(1)%bins(j,ip,nbin)%rho/pop_sum
                  end do
               end do

               if(iproc == 0 .and. print_output)then
                  write(13,'(''Initial populations for incident energy ='',1x,1pe15.7)')e_in
                  do ip =0, 1
                     do j = 0, nucleus(1)%j_max
                        xj = j + nucleus(1)%jshift
                        write(13,'(1x,f4.1,2(1x,f10.6))')xj,2.0*real(ip)-1.0,target%pop_xjpi(j,ip)
                     end do
                  end do
                  do ijj = 1, num_pop
                     write(13,'(3(1x,i5),f10.6)')ijj,pop_j(ijj), pop_ip(ijj), pop_prob(ijj)
                  end do
               end if
            end if
         end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Zero gamma cross section arrays
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do icomp = 1, num_comp
            do j = 1, nucleus(icomp)%num_discrete
               do k = 1, nucleus(icomp)%state(j)%nbranch
                  nucleus(icomp)%state(j)%cs(k) = 0.0d0
               end do
            end do
         end do

         if(track_primary_gammas)Primary_Gamma_Spectrum(0:num_e) = 0.0d0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Now, perform interpolation for Shape Elastic cross section
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(.not. pop_calc)then

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------  Allocate Direct component
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            reaction_cs(in) = absorption_cs(in)

            if(iproj > 0)then
               tot_direct = 0.0d0
               direct_tot(in) = 0.0d0
               direct_cc(in) = 0.0d0
               direct_dwba(in) = 0.0d0
               j = 1
               direct_cs(1:OpticalCS%numcc) = 0.0d0
               direct_Ang(0:Ang_L_max,1:OpticalCS%numcc) = 0.0d0
               direct_prob(0:ixx_max,1:OpticalCS%numcc) = 0.0d0

               do n = 1, OpticalCS%numcc
                  if(OpticalCS%state(n)%istate == target%istate .and. OpticalCS%state(n)%state_type == 1)then    !  This is shape elastic so skip this in decays
                     SE_cs(in) = interpolate_exp(0,E_in,OpticalCS%nume,OpticalCS%energy,OpticalCS%optical_cs(1,n))
                     do L = 0, Ang_L_max
                        SE_Ang(L,in) = interpolate_exp(0,E_in,OpticalCS%nume,OpticalCS%energy,OpticalCS%optical_leg(1,L,n))
                     end do
                     cycle
                  end if
                  if(e_rel < OpticalCS%state(n)%energy)cycle
                  direct_cs(n) = interpolate_exp(0,E_in,OpticalCS%nume,OpticalCS%energy,OpticalCS%optical_cs(1,n))
                  tot_direct = tot_direct + direct_cs(n)
                  reaction_cs(in) = reaction_cs(in) + direct_cs(n)
                  direct_tot(in) = direct_tot(in) + direct_cs(n)
                  if(OpticalCS%state(n)%state_type == 1)direct_cc(in) = direct_cc(in) + direct_cs(n)
                  if(OpticalCS%state(n)%state_type <= 0)direct_dwba(in) = direct_dwba(in) + direct_cs(n)
                  do L = 0, Ang_L_max
                     direct_Ang(L,n) = interpolate_exp(0,E_in,OpticalCS%nume,OpticalCS%energy,OpticalCS%optical_leg(1,L,n))
                  end do
                  direct_prob(0:ixx_max,n) = 0.0d0
                  do i = 1, ixx_max
                     x = real(i,kind=8)*delta_ix - 1.0d0
                     x1 = x - delta_ix
                     alf = 0.0d0
                     bet = 0.0d0
                     sum = 0.0d0
                     do L = 0, Ang_L_max
                        sum = sum + direct_Ang(L,n)*(poly(L,1,alf,bet,x) + poly(L,1,alf,bet,x1))
                     end do
                  direct_prob(i,n) = direct_prob(i-1,n) + 0.5d0*delta_ix*sum
                  end do
                  do i = 0, ixx_max
                  direct_prob(i,n) = direct_prob(i,n)/direct_prob(ixx_max,n)
                  end do
                  do L = 0, Ang_L_max
                     direct_Ang(L,n) = direct_Ang(L,n)/direct_prob(ixx_max,n)
                  end do
               end do
!---   Elastic scattering angle probability
               SE_prob(0:ixx_max) = 0.0d0
               do i = 1, ixx_max
                  x = real(i,kind=8)*delta_ix - 1.0d0
                  x1 = x - delta_ix
                  alf = 0.0d0
                  bet = 0.0d0
                  sum = 0.0d0
                  do L = 0, Ang_L_max
                     sum = sum + SE_Ang(L,in)*(poly(L,1,alf,bet,x) + poly(L,1,alf,bet,x1))
                  end do
                  SE_prob(i) = SE_prob(i-1) + 0.5d0*sum*delta_ix
               end do
               do i = 0, ixx_max
                  SE_prob(i) = SE_prob(i)/SE_prob(ixx_max)
               end do
            end if

            if(iproj == 1 .and. print_me)then
               write(6,'(''Total cross section ='',1x,1pe15.7,1x,a2)')(reaction_cs(in) + SE_cs(in))*cs_scale,cs_units
               write(6,'(''Total Reaction cross section ='',1x,1pe15.7,1x,a2)')reaction_cs(in)*cs_scale,cs_units
               write(6,'(''Compound cross section ='',1x,1pe15.7,1x,a2)')absorption_cs(in)*cs_scale,cs_units
               write(6,'(''Direct excitation cross section ='',1x,1pe15.7,1x,a2)')tot_direct*cs_scale,cs_units
               write(6,'(''Shape Elastic cross section ='',1x,1pe15.7,1x,a2)')SE_cs(in)*cs_scale,cs_units
            elseif(iproj == 0 .and. iproc == 0 .and. print_output)then
               write(13,'(''Total cross section ='',1x,1pe15.7,1x,a2)')reaction_cs(in)*cs_scale,cs_units
               write(13,'(''Total Reaction cross section ='',1x,1pe15.7,1x,a2)')reaction_cs(in)*cs_scale,cs_units
               write(13,'(''Compound cross section ='',1x,1pe15.7,1x,a2)')absorption_cs(in)*cs_scale,cs_units
            elseif(iproj > 1 .and. print_me)then
               write(6,'(''Total Reaction cross section ='',1x,1pe15.7,1x,a2)')reaction_cs(in)*cs_scale,cs_units
               write(6,'(''Compound cross section ='',1x,1pe15.7,1x,a2)')absorption_cs(in)*cs_scale,cs_units
               write(6,'(''Direct excitation cross section ='',1x,1pe15.7,1x,a2)')tot_direct*cs_scale,cs_units
            elseif(iproj > 1 .and. iproc == 0 .and. print_output)then
               write(13,'(''Total Reaction cross section ='',1x,1pe15.7,1x,a2)')reaction_cs(in)*cs_scale,cs_units
               write(13,'(''Compound cross section ='',1x,1pe15.7,1x,a2)')absorption_cs(in)*cs_scale,cs_units
               write(13,'(''Direct excitation cross section ='',1x,1pe15.7,1x,a2)')tot_direct*cs_scale,cs_units
            end if

!******************************************************************************
!------   If absorption_cs is too small skip all the decays as nothing can happen anyway
            if(iproj > 0 .and. absorption_cs(in) < cs_threshold)then
               if(print_me)then
                  write(6,'(''******************************************************************************'')')
                  write(6,'(''*  Absorption cross section < '',e12.5,'', skipping decay processes   *'')')cs_threshold
                  write(6,'(''******************************************************************************'')')
               elseif(iproc == 0 .and. print_output)then
                  write(13,'(''******************************************************************************'')')
                  write(13,'(''*  Absorption cross section < '',e12.5,'', skipping decay processes   *'')')cs_threshold
                  write(13,'(''******************************************************************************'')')
               end if
               goto 1901
            end if

!*****************************************************************************
!------
         else
            if(print_me)write(6,*)'Total population ',absorption_cs(in)
         end if

         if(.not. file_energy_index)then
            ifile = core_file
            file_name(ifile:ifile+4) = '_Ein_'
            ifile =ifile + 4
            if(e_in < 10.0)then
               file_name(ifile+1:ifile+2)='00'
               ifile=ifile+2
               write(file_name(ifile+1:ifile+6),'(f6.4)')e_in
               ifile=ifile+6
            elseif(e_in < 100.0)then
               file_name(ifile+1:ifile+2)='0'
               ifile=ifile+1
               write(file_name(ifile+1:ifile+7),'(f7.4)')e_in
               ifile=ifile+7
            elseif(e_in < 1000.0)then
               write(file_name(ifile+1:ifile+7),'(f8.4)')e_in
               ifile=ifile+8
            end if
         else
            ifile = core_file
            file_name(ifile:ifile+3) = '_iE_'
            ifile =ifile + 3
            
            file_name(ifile+1:ifile+3) = '000'
            if(in < 10)then
               write(file_name(ifile+3:ifile+3),'(i1)')in
            elseif(in < 100)then
               write(file_name(ifile+2:ifile+3),'(i2)')in
            else
               write(file_name(ifile+1:ifile+2),'(i3)')in
            end if
            ifile = ifile + 3
         end if

         directory(1:ilib_dir) = lib_dir(1:ilib_dir)
         idir = ilib_dir + 1
         directory(idir:idir) = '/'
         idir = idir + 1
         directory(idir:idir+10) = 'Event-files'
         idir = idir + 11
         directory(idir:idir) = '/'

!-----   If using MPI, then we need a different file name for each MPI process
#if(USE_MPI == 1)
         if(nproc > 1)then
            node_name(1:8) = '_node000'
            if(iproc < 10)then
               write(node_name(8:8),'(i1)')iproc
            else if(iproc >= 10 .and. iproc < 100)then
               write(node_name(7:8),'(i2)')iproc
            else if(iproc >= 100 .and. iproc < 1000)then
               write(node_name(6:8),'(i3)')iproc
            else if(iproc >= 1000)then
               write(6,*)'Error, many more MPI processes than expected'
               call exit_YAHFC(101)
            end if
         end if
#endif

         if(dump_events .and. binary_event_file)then
            if(nproc == 1)then
               open(unit=88,file=                                                                   &
                    directory(1:idir)//file_name(1:ifile)//'.Events.bin',                           &
                    form='unformatted',status='unknown')
               write(6,*)'Events written to binary file ', directory(1:idir)//file_name(1:ifile)//'.Events.bin'
            else
               open(unit=88,file=                                                                   &
                    directory(1:idir)//file_name(1:ifile)//node_name//'.Events.bin',                &
                    form='unformatted',status='unknown')
               write(6,*)'Events for MPI process ',iproc,' are written to binary file ',            &
                  directory(1:idir)//file_name(1:ifile)//node_name//'.Events.bin'
            end if
         end if
         if(dump_events .and. .not. binary_event_file)then
            if(nproc == 1)then
               open(unit=88,file=                                                                   &
                    directory(1:idir)//file_name(1:ifile)//'.Events.txt',                           &
                    status='unknown')
               write(6,*)'Events written to formatted file ', directory(1:idir)//file_name(1:ifile)//'.Events.txt'
            else
               open(unit=88,file=                                                                   &
                    directory(1:idir)//file_name(1:ifile)//node_name//'.Events.txt',                &
                    status='unknown')
               write(6,*)'Events for MPI process ',iproc,' are written to binary file ',            &
                    directory(1:idir)//file_name(1:ifile)//node_name//'.Events.txt'
            end if
         end if

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Start of Monte Carlo Sampling Loop     -----------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
         do i = 1, num_channels
            Exit_Channel(i)%num_event = 0
         end do

         max_e_diff = 0.0d0
         avg_e_diff = 0.0d0
         re_baseline = 0
         if(pop_calc)then
            decay_mult(0:6,in) = 0.0d0
            decay_var(0:6,in) = 0.0d0
         end if
!
!---------------------------------------------------------------------------------------
!--------    Monte Carlo sampling loop for each event                              -----
!---------------------------------------------------------------------------------------
!
         tally_norm = 0.0d0
         do nsamp = iproc, num_mc_samp - 1, nproc

            fission_decay = .false.

            num_part_type(0:6) = 0
            part_fact(0:7) = 1.0d0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------     Reset the boost matrix
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            Boost_COM(0:3,0:3) = 0.0d0
            do i = 0, 3
               Boost_COM(i,i) = 1.0d0
            end do
            hung_up = .false.

!------    It is a population calculation. Set initial kinetic energy for the decaying system
            if(pop_calc)then
               KE = -1.0
               do while (KE < 0.0)
                   KE = init_Kinetic_Energy + dInit_Kinetic_Energy*Gauss_var(iseed_32)
!                   KE = init_Kinetic_Energy + dInit_Kinetic_Energy*Gauss_var(iseed_64)
               end do
               nucleus(1)%Kinetic_Energy = KE
!------    Now get initial excitation energy. This is set in Pop_data(i)%Ex_pop
!------    with Gaussian spread in Pop_data(i)%dEx_pop
10001          E_in = 100000000.0d0
               do while (E_in > nucleus(1)%Ex_max .or. nbin <= 0)
                  E_in = Pop_data(in)%Ex_pop + Pop_data(in)%dEx_pop*Gauss_var(iseed_32)
!                  E_in = Pop_data(in)%Ex_pop + Pop_data(in)%dEx_pop*Gauss_var(iseed_64)
                  Ex_tot = E_in
                  nbin = find_ibin(Ex_tot, 1)
                  Ex_tot = nucleus(1)%e_grid(nbin)
               end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Setup boost matrix for this kinetic energy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               nucleus(1)%Kinetic_Energy = Kinetic_Energy
!               KE = nucleus(1)%Kinetic_Energy
!               mass_2 = nucleus(1)%Mass
               mass_2 = nucleus(1)%Mass + Ex_tot
               pp_res = sqrt(2.0d0*mass_2*KE)
               p_res(0) = mass_2 + KE
               pp_res = sqrt(abs(p_res(0)**2 - mass_2**2))
               p_res(1) = pp_res
               p_res(2) = 0.0d0
               p_res(3) = 0.0d0
               do i = 1, 3
                  v_res(i) = p_res(i)/mass_2
               end do
               beta = pp_res/mass_2
               gamma = 1.0d0/sqrt(1.0d0 - beta**2)
               Boost_Lab(0:3,0:3) = 0.0d0
               do i = 0, 3
                  Boost_Lab(i,i) = 1.0d0
               end do
               Boost_Lab(0,0) = gamma
               Boost_Lab(0,1) = -gamma*v_res(1)
               Boost_Lab(1,0) = Boost_Lab(0,1)
               if(beta > 1.0d-16)Boost_Lab(1,1) = 1.0d0 + (gamma - 1.0d0)*v_res(1)**2/beta**2
!------     Populations are based on level density, so recompute probabilities based on
!------     excitation energy with spread in the excitaiton energy
               if(j_pop_calc)then
                  pop_sum = 0.0
                  do j = 0, nucleus(1)%j_max
                     do ip = 0, 1
                        pop_sum = pop_sum + nucleus(1)%bins(j,ip,nbin)%rho
                     end do
                  end do
                  if(pop_sum <= 1.0d-7) goto 10001             !   No levels, get new excitation energy
                  ijj = 0
                  pop_prob(1:num_pop) = 0.0d0
                  prev = 0.0d0
                  do j = 0, nucleus(1)%j_max
                     do ip = 0, 1
                        ijj = ijj + 1
                        pop_prob(ijj) = prev + nucleus(1)%bins(j,ip,nbin)%rho/pop_sum
                        prev = pop_prob(ijj)
                        pop_j(ijj) = j
                        pop_ip(ijj) = ip
                        target%pop_xjpi(j,ip) = nucleus(1)%bins(j,ip,nbin)%rho/pop_sum
                     end do
                  end do
               end if
            end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Determine entrance Channel based on transmission coefficients     +
!-------   from the optical potential                                        +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
            num_part = 0
            part_data(1:n_dat,1:dim_part) = 0.0
            part_Ang_data(0:Ang_L_max,1:dim_part) = 0.0
            extra_angle_data(3*num_theta,1:dim_part) = 0.0d0
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Determine first decay, including possibility of pre-equilibrium   +
!-------   or direct reaction                                                +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

            compound = .false.
            direct = .false.
            preeq_decay = .false.
            fission_decay = .false.
            inelastic = .false.
            cc_decay = .false.
            dwba_decay = .false.
            if(.not.pop_calc .and. iproj > 0)then      !***********   particle projectile   *******
               ran = random_32(iseed_32)
               if(ran <= absorption_cs(in)/reaction_cs(in))then         !   Normal Compound formation
                  sp1 = abs(spin_target-spin_proj)
                  call find_prob_point(num_channel,channel_prob,ran,ifind)
                  icomp_i = 1
                  nbin_i = nbin
                  l_i = ichannel(1,ifind)
                  is_i = ichannel(2,ifind)
                  Ix_i = ichannel(3,ifind)
                  ip_i = ichannel(4,ifind)

                  xnbin_i = nbin_i
                  par_i = 2*ip_i - 1
                  ran = random_32(iseed_32)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------  Check if this part of reaction cross section is populated by pre-equilibirum
!-----------  or compound.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  if(ran <= preeq_prob .and. preeq_prob >= 1.0d-9 .and. E_in >= de/2.0)then     !  Pre-equilibrium decay
                     preeq_decay = .true.
                     compound = .false.

                     call PREEQ_sample(iproj, in, itarget, istate, e_in, ex_tot,         &
                                       l_i, is_i, Ix_i, ip_i, icomp_i, icomp_f,          &
                                       Ix_f, l_f, ip_f, nbin_f, idb,                     &
                                       n_dat, dim_part, num_part_type, part_fact,        &
                                       num_part, part_data,                              &
                                       Ang_L_max, part_Ang_data,                         &
                                       num_theta, extra_angle_data)

                  else                                                                          !  Normal compound nucleus decay
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------   Normal compound state. Decay first entrance state separately. Decay of
!-----------   of first energy bin is somewhat different as entrance formation
!-----------   matters because of width fluctuations. In addition, angular distribution
!-----------   is not isotropic for first decay.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                     compound = .true.
                     preeq_decay = .false.
                     num_part = 0
                     idb = 0
                     call MC_primary_decay(iproj, spin_target,                           &
                                           l_i, is_i, Ix_i, ex_tot, icomp_i,             &
                                           icomp_f, Ix_f, ip_f, nbin_f, idb,             &
                                           n_dat, dim_part, num_part_type, part_fact,    &
                                           num_part, part_data,                          &
                                           Ang_L_max, part_Ang_data,                     &
                                           ixx_max, delta_ix, Leg_poly,                  &
                                           num_theta, extra_angle_data)

                     if(nbin_f < 0)fission_decay = .true.
                     if(fission_decay)goto 101
                  end if
               else                                            !  Direct reaction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------    This is a direct component, not compound.
!--------    Find which discrete state it is excited to, and set decay to incident particle
!--------    for this specific reaction. Set idb == 1 so that we know it populates a
!--------    discrete state
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  preeq_decay = .false.
                  compound = .false.
                  direct = .true.                       !   Identify it as a direct excitation and not compound
                  cc_decay = .false.
                  dwba_decay = .false.
                  icomp_i = 1
                  ran = random_32(iseed_32)
                  sum_d = 0.0d0
                  j = 1
                  do n = 1, OpticalCS%numcc
                     j = OpticalCS%state(n)%istate
                     if(j == target%istate .and. OpticalCS%state(n)%state_type == 1)cycle
                     if(j == -1)cycle
                     sum_d = sum_d + direct_cs(n)/tot_direct
                     if(ran <= sum_d)then
                        num_part = 1
                        icomp_i = 1
                        icomp_f = itarget
                        nbin_f = j
                        xI_f = OpticalCS%state(n)%spin
                        xip_f = OpticalCS%state(n)%parity
                        Ex_f = 0.0d0
                        if(OpticalCS%state(n)%state_type == 1)then       !  Direct to discrete state
                           idb = 1
                           xI_f = OpticalCS%state(n)%spin
                           Ix_f = nint(xI_f - nucleus(icomp_f)%jshift)
                           xip_f = OpticalCS%state(n)%parity
                           ip_f = nint((xip_f + 1.0d0)/2.0d0)
                           E_f = ex_tot - nucleus(icomp_i)%sep_e(iproj) -           &
                                 nucleus(icomp_f)%state(nbin_f)%energy
                           Ex_f = nucleus(icomp_f)%state(nbin_f)%energy
                           cc_decay = .true.
                        elseif(OpticalCS%state(n)%state_type == 0)then    !  Direct to continuum bin with explicit J
                           idb = 0
                           xI_f = OpticalCS%state(n)%spin
                           Ix_f = nint(xI_f - nucleus(icomp_f)%jshift)
                           xip_f = OpticalCS%state(n)%parity
                           ip_f = nint((xip_f + 1.0d0)/2.0d0)
                           Ex = OpticalCS%state(n)%E_min +                          &
                                random_32(iseed_32)*OpticalCS%state(n)%Delta_E
                           nbin_f = find_ibin(Ex,itarget)
                           if(nbin_f < 1)then
                              if(iproc == 0)then
                                 write(6,*)
                                 write(6,*)'***************************************************************'
                                 write(6,*)'ERROR - Attempting to excite a continuous DWBA state below'
                                 write(6,*)'the first bin. To continue either rerun Frescox or manually'
                                 write(6,*)'set Ecut for the target nucleus below ',OpticalCS%state(n)%E_min
                                 write(6,*)'***************************************************************'
                                 write(6,*)
                              end if
                              call exit_YAHFC(902)
                           end if
                           Ex = nucleus(itarget)%e_grid(nbin_f)
                           E_f = ex_tot - nucleus(icomp_i)%sep_e(iproj) - Ex
                           Ex_f = nucleus(itarget)%e_grid(nbin_f)
                           dwba_decay = .true.
                        elseif(OpticalCS%state(n)%state_type == -1)then    !  Direct to continuum bin w/o explicit J
                           idb = 0
                           ran = random_32(iseed_32)
                           do Ix_f = OpticalCS%state(n)%Ix_min, OpticalCS%state(n)%Ix_max   !  find final spin for DWBA K state
                              if(ran <= OpticalCS%state(n)%spin_prob(Ix_f))exit
                           end do
                           xI_f = real(Ix_f,kind=8) + nucleus(icomp_f)%jshift
                           xip_f = OpticalCS%state(n)%parity
                           ip_f = nint((xip_f + 1.0d0)/2.0d0)
                           Ex = OpticalCS%state(n)%E_min +                          &
                                random_32(iseed_32)*OpticalCS%state(n)%Delta_E
                           nbin_f = find_ibin(Ex,itarget)
                           Ex = nucleus(itarget)%e_grid(nbin_f)
                           dwba_decay = .true.
                        end if
                        part_data(1,num_part) = real(icomp_f,kind=8)
                        part_data(2,num_part) = real(iproj,kind=8)
                        part_data(3,num_part) = xI_f
                        part_data(4,num_part) = xip_f
                        part_data(5,num_part) = real(nbin_f,kind=8)
                        part_data(6,num_part) = real(idb,kind=8)
                        part_data(7,num_part) = 0.0d0
                        part_data(8,num_part) = 1.0d0
                        part_data(9,num_part) = E_f
                        num_part_type(iproj) = num_part_type(iproj) + 1
                        if(num_part_type(iproj) >= max_particle(iproj))part_fact(iproj) = 0.0d0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------   get theta
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        theta_0 = 0.0d0
                        phi_0 = 0.0d0

                        if(.not. xs_only)then
                           do nang = 1, num_theta
                              ran = random_32(iseed_32)
                              do i = 1, ixx_max
                                 x = real(i,kind=8)*delta_ix - 1.0d0
                                 if(direct_prob(i,n) > ran)exit
                              end do
                              x = x - random_32(iseed_32)*delta_ix*0.999999d0
                              extra_angle_data(nang,num_part) = acos(x)
                           end do
                           theta_0 = extra_angle_data(1,num_part)
                           ran = random_32(iseed_32)
                           phi_0 = ran*2.0d0*pi
                        end if

                        part_data(10,num_part) = theta_0
                        part_data(11,num_part) = phi_0
                        part_data(19,num_part) = 1.0d0
                        part_data(20,num_part) = Ex_f
                        part_data(21,num_part) = icomp_i
                        part_data(22,num_part) = nucleus(itarget)%state(istate)%spin
                        part_data(23,num_part) = nucleus(itarget)%state(istate)%parity
                        part_data(24,num_part) = istate
                        do L = 0, Ang_L_max
                           part_Ang_data(L,num_part) = direct_Ang(L,n)
                        end do
                        exit
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Direct to bin, reset data to send to Monte carlo Decay
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                     end if
                  end do
               end if
            elseif(.not. pop_calc .and. iproj == 0)then  !*************  photo-absorption calculation   *******
               ran = random_32(iseed_32)
               call find_prob_point(num_channel,channel_prob,ran,ifind)
               icomp_i = 1
               nbin_i = nbin
               l_i = ichannel(1,ifind)
               is_i = ichannel(2,ifind)
               Ix_i = ichannel(3,ifind)
               ip_i = ichannel(4,ifind)
               call MC_decay_bin(icomp_i, Ix_i, ip_i, nbin_i,                 &
                                 icomp_f, Ix_f, ip_f, nbin_f, idb,            &
                                 n_dat, dim_part, num_part_type, part_fact,   &
                                 num_part, part_data,                         &
                                 Ang_L_max, part_Ang_data,                    &
                                 num_theta, extra_angle_data)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Check if decay is hung up in a bin with no place to go. If so, force
!----------   decay to a discrete state and collect statistics to rpeort later
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               if(nbin_f == 0)then      !   Hung up with no place to go

                  hung_up = .true.

                  call force_decay(icomp_i, nbin_i,                             &
                                   icomp_f, Ix_f, ip_f, nbin_f, idb,            &
                                   n_dat, dim_part, num_part, part_data,        &
                                   Ang_L_max, part_Ang_data,                    &
                                   num_theta, extra_angle_data)
               end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Check if decay is fission. If so, collect separately and exit decay loop
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               if(nbin_f < 0)fission_decay = .true.
               if(fission_decay)goto 101
            elseif(pop_calc)then         !*************   population calculation   *******
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   This is a population decay rather than a reaction
!---------   Population for J and parity were input for a given starting
!---------   exciation energy. Much like reaction calculation, conduct
!---------   first decay and then pass on to later routines to decay
!---------   to terminus at either the ground state of some nucleus or an
!---------   isomer where subsequent decay is forbidden.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

               compound = .true.
               ran = random_32(iseed_32)
               do i = 1, Pop_data(in)%num_pop
                  if(ran <= pop_prob(i))exit
               end do

               icomp_i = 1
               nbin_i = nbin
               Ix_i = pop_j(i)
               ip_i = pop_ip(i)


               call MC_decay_bin(icomp_i, Ix_i, ip_i, nbin_i,                 &
                                 icomp_f, Ix_f, ip_f, nbin_f, idb,            &
                                 n_dat, dim_part, num_part_type, part_fact,   &
                                 num_part, part_data,                         &
                                 Ang_L_max, part_Ang_data,                    &
                                 num_theta, extra_angle_data)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Check if decay is hung up in a bin with no place to go. If so, force
!----------   decay to a discrete state and collect statistics to rpeort later
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               if(nbin_f == 0)then      !   Hung up with no place to go

                  hung_up = .true.

                  call force_decay(icomp_i, nbin_i,                             &
                                   icomp_f, Ix_f, ip_f, nbin_f, idb,            &
                                   n_dat, dim_part, num_part, part_data,        &
                                   Ang_L_max, part_Ang_data,                    &
                                   num_theta, extra_angle_data)
               end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Check if decay is fission. If so, collect separately and exit decay loop
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               if(nbin_f < 0)fission_decay = .true.
               if(fission_decay)goto 101
            end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Decay entry states until they come to a discrete state. Finished when idb == 1
!---------   The decay can also finish with idb=0 if the contiuum bin hangs with no decay
!---------   path. This gets kicked out below with nbin_f == 0. Fission is also special
!---------   there are no further decays possible.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            icomp_i = icomp_f
            icomp_f = 0
            Ix_i = Ix_f
            Ix_f = 0
            ip_i = ip_f
            ip_f = 0
            nbin_i = nbin_f
            nbin_f =0

            if(Ix_i > nucleus(icomp_i)%j_max)then
               write(6,*)'Ix_i > j_max before decay loop, Ix_i =', Ix_i,' j_max = ',nucleus(icomp_i)%j_max
               write(6,*)'Processor # ',iproc
               write(6,*)'Check: pop_calc = ', pop_calc
               write(6,*)'Check: compound = ', compound
               write(6,*)'Check: preeq_decay = ', preeq_decay
               write(6,*)'Check: direct = ', direct
               write(6,*)'Check: cc_decay = ', cc_decay
               write(6,*)'Check: dwba_decay = ', dwba_decay
               write(6,*)'icomp = ', icomp_i
               write(6,*)'nbin_i = ', nbin_i
               write(6,*)'idb = ', idb
               write(6,'(2(1x,i10),3(1x,f10.5))')nsamp, num_part
               call exit_YAHFC(101)
            end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Option for doing only the primary decay
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Now decay continous energy bins until they reach a discrete state (occurs when idb=1)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            do while(idb /= 1)
               if(Ix_i > j_max)then
                  if(iproc == 0)write(6,*)'error Ix_i > j_max , Ix_i = ',Ix_i
                  flush(6)
                  call exit_YAHFC(101)
               end if
               call MC_decay_bin(icomp_i, Ix_i, ip_i, nbin_i,                   &
                                 icomp_f, Ix_f, ip_f, nbin_f, idb,              &
                                 n_dat, dim_part, num_part_type, part_fact,     &
                                 num_part, part_data,                           &
                                 Ang_L_max, part_Ang_data,                      &
                                 num_theta, extra_angle_data)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Check if decay is hung up in a bin with no place to go. If so, force       +
!----------   decay to a discrete state and collect statistics to rpeort later           +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              if(nbin_f == 0)then
                 hung_up = .true.
                 call force_decay(icomp_i, nbin_i,                             &
                                  icomp_f, Ix_f, ip_f, nbin_f, idb,            &
                                  n_dat, dim_part, num_part, part_data,        &
                                  Ang_L_max, part_Ang_data,                    &
                                  num_theta, extra_angle_data)
              end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Check if decay is fission. If so, collect separately and exit decay loop   +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

               if(nbin_f < 0)fission_decay = .true.
               if(fission_decay)goto 101

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------   reset initial and final state parameters
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               icomp_i = icomp_f
               icomp_f = 0
               Ix_i = Ix_f
               Ix_f = 0
               ip_i = ip_f
               ip_f = 0
               nbin_i = nbin_f
               nbin_f =0
            end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Finish decay discrete states all the way to the end     ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ichann = find_channel(n_dat, dim_part, num_part, part_data, num_part_type)
            call MC_decay_state(icomp_i, nbin_i,                             &
                                n_dat,dim_part,num_part,part_data,           &
                                Ang_L_max, part_Ang_data,                    &
                                num_theta, extra_angle_data, ichann, in)
!
101         continue
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------    Get channel number                                          +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------   ichann = label for this channel, now have only discrete gammas left
!----------------   in = index for incident energy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ichann = find_channel(n_dat, dim_part, num_part, part_data, num_part_type)

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------    Decay is finished     -------------------------------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------Channel_cs---    Determine the weight of this event     --------------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            tally_weight = 1.0d0
            sum_e = 0.0d0
            do nn = 1, num_part
               tally_weight = tally_weight*part_data(19,nn)
               if(.not. fission_decay)then
                   k = nint(part_data(2,nn))
                   icomp_i = nint(part_data(21,nn))
                   icomp_f = nint(part_data(1,nn))
                   idb = nint(part_data(6,nn))
                   nbin_f = nint(part_data(5,nn))
                   sum_e = sum_e + part_data(9,nn)
                   if(k > 0 .and. k < 7)sum_e = sum_e + nucleus(icomp_i)%sep_e(k)
                   if(nn == num_part .and. idb == 1)then
                      if(nucleus(icomp_f)%state(nbin_f)%isomer)sum_e = sum_e +              &
                         nucleus(icomp_f)%state(nbin_f)%energy
                   end if
                end if
            end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  The weight of this event is tally_weight      ------------------------------------+
!--------  Collect the total in tally_norm               ------------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            tally_norm = tally_norm + tally_weight
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  Collect data on primary decay               --------------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(primary_decay)then
               inuc = nint(part_data(1,1))
               iX_f = nint(part_data(3,1) - nucleus(inuc)%jshift)
               ip_f = nint(part_data(4,1) + 1.0d0)/2
               n_f = nint(part_data(5,1))
               idb = nint(part_data(6,1))
               if(idb == 0)then
                  nucleus(inuc)%bins(iX_f,ip_f,n_f)%pop =                             &
                      nucleus(inuc)%bins(iX_f,ip_f,n_f)%pop + tally_weight
               elseif(idb == 1)then
                  nucleus(inuc)%state(n_f)%pop =                                      &
                      nucleus(inuc)%state(n_f)%pop + tally_weight
               end if
            end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  Check if energy adds up as it should        --------------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(.not. fission_decay)then
               e_diff = abs(ex_tot - sum_e)
               if(e_diff > max_e_diff)max_e_diff = e_diff
               avg_e_diff = avg_e_diff + e_diff
               if(e_diff >= 1.0d-3)then
                  num_bad_samp_e = num_bad_samp_e + 1
                  write(6,*)'Processor #',iproc
                  write(6,*)'Energy not conserved to within 1 keV'
                  write(6,*)'E_in = ',E_in,' nsamp = ',nsamp,' iproc = ',iproc
                  write(6,*)'Ex_tot = ',ex_tot,' Sum_e ',sum_e
                  write(6,*)'Diff = ',e_diff
                  write(6,*)'Check: pop_cal = ',pop_calc
                  write(6,*)'Check: compound = ',compound
                  write(6,*)'Check: preeq_decay = ', preeq_decay
                  write(6,*)'Check: direct = ',direct
                  write(6,*)'Check: dwba = ',dwba_decay
!                  write(6,*)'Info written to ',out_file(1:ifile)//'.bad_samp_energy'
                  write(6,*)'Numer of particles = ', num_part
                  do nn = 1, num_part
                     write(6,'(1x,i5,21(2x,f12.5))')nn,(part_data(k,nn), k = 1, n_dat)
                  end do
!                  write(28,*)'Processor #',iproc
!                  write(28,*)'Energy not conserved to within 1 keV'
!                  write(28,*)'E_in = ',E_in,' nsamp = ',nsamp,' iproc = ',iproc
!                  write(28,*)'Ex_tot = ',ex_tot,' Sum_e ',sum_e
!                  write(28,*)'Diff = ',e_diff
!                  write(28,*)'Check: pop_cal = ',pop_calc
!                  write(28,*)'Check: compound = ',compound
!                  write(28,*)'Check: preeq_decay = ', preeq_decay
!                  write(28,*)'Check: direct = ',direct
!                  write(28,*)'Check: dwba = ',dwba_decay
!                  write(28,*)'Info written to ',out_file(1:ifile)//'.bad_samp_energy'
!                  do nn = 1, num_part
!                     write(28,'(1x,i5,21(2x,f12.5))')nn,(part_data(k,nn), k = 1, n_dat)
!                  end do
!                  flush(28)
               end if
            end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    Add "spread" to emitted energies based on size of the bins
!------    Helps to remove periodicity in the outgoing spectra as the transitions occur
!------    on a fixed grid.
!------    At each stage, we also need to transform to Lab and COM frames from the rest
!------   frame where each decay occurred
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(.not. xs_only)then
               sum_e = 0.0d0
               do nn = 1, num_part
                  k = nint(part_data(2,nn))
                  if(k == 7) cycle                     !   Fission event, essentially finished here
                  theta_0 = part_data(10,nn)
                  phi_0 = part_data(11,nn)
                  idb = nint(part_data(6,nn))
                  icomp_i = nint(part_data(21,nn))
                  nbin_i = nint(part_data(24,nn))
                  icomp_f = nint(part_data(1,nn))
                  nbin_f = nint(part_data(5,nn))
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Add spread +/- de/2 in emitted energies
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
                  if(nn < num_part-1 .and. idb == 0)then     !  Decay to a continuous energy bin
                     shift_max = 0.5d0*nucleus(icomp_f)%delta_e(nbin_f)
                     shift_min = -0.5d0*nucleus(icomp_f)%delta_e(nbin_f)
                     if(nint(part_data(6,nn+1)) == 1)shift_max = min(shift_max,part_data(9,nn+1))
!                     shift = random_64(iseed_64)*(shift_max-shift_min) + shift_min
                     shift = random_32(iseed_32)*(shift_max-shift_min) + shift_min
                     part_data(9,nn) = part_data(9,nn) + shift
                     part_data(9,nn+1) = part_data(9,nn+1) - shift
                  end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Transform to COM and Lab frame
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
                  e_rel = part_data(9,nn)
                  mass_1 = particle(k)%Mass
                  if(idb == 0)then
                     mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%e_grid(nbin_f)
                  else
                     mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%state(nbin_f)%energy
                  end if

                  call Boost_frame(e_rel, mass_1, mass_2, theta_0, phi_0,                 &
                                   Boost_Lab, Boost_COM, T_1, theta, phi,                 &
                                   T_2, T_L, theta_L, phi_L)

                  part_data(12,nn) = T_1
                  part_data(13,nn) = theta
                  part_data(14,nn) = phi
                  part_data(15,nn) = T_L
                  part_data(16,nn) = theta_L
                  part_data(17,nn) = phi_L
                  part_data(18,nn) = T_2
                  nucleus(icomp_f)%Kinetic_Energy = T_2

                  extra_angle_data(num_theta+1,nn) = theta               !   COM frame
                  extra_angle_data(2*num_theta+1,nn) = theta_L           !   Lab frame
!
!----    Sum energy from the Rest reference frames and check if "conserved"
!
                  sum_e = sum_e + part_data(9,nn)
                  if(k > 0 .and. k < 7)sum_e = sum_e + nucleus(icomp_i)%sep_e(k)
                  if(nn == num_part .and. idb == 1)then
                     if(nucleus(icomp_f)%state(nbin_f)%isomer)sum_e = sum_e +       &
                         nucleus(icomp_f)%state(nbin_f)%energy
                  end if
               end do

               do nang = 2, num_theta
                  do nn = 1, num_part
                     k = nint(part_data(2,nn))
                     if(k == 7) cycle                     !   Fission event, essentially finished here
                     theta_0 = extra_angle_data(nang,nn)
                     phi_0 = part_data(11,nn)
                     idb = nint(part_data(6,nn))
                     icomp_f = nint(part_data(1,nn))
                     nbin_f = nint(part_data(5,nn))
                     e_rel = part_data(9,nn)
                     mass_1 = particle(k)%Mass
                     if(idb == 0)then
                        mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%e_grid(nbin_f)
                     else
                        mass_2 = nucleus(icomp_f)%Mass + nucleus(icomp_f)%state(nbin_f)%energy
                     end if

                     call Boost_frame(e_rel, mass_1, mass_2, theta_0, phi_0,                 &
                                      Boost_Lab, Boost_COM, T_1, theta, phi,                 &
                                      T_2, T_L, theta_L, phi_L)
                     extra_angle_data(num_theta+nang,nn) = theta               !   COM frame
                     extra_angle_data(2*num_theta+nang,nn) = theta_L           !   Lab frame
                  end do
               end do
!
!----   if difference exceeds 1.0d-4 MeV write out to unit 400
!----   Count number of events and warn user at the end
!
               if(.not. fission .and. abs(ex_tot-sum_e) > 1.0d-4 )then
                   write(6,*)'Energy conservation issue detected with extra angles and smear in energy'
                   write(6,*)'iproc = ',iproc
                   write(6,*)'Info written to ',out_file(1:ifile)//'.bad_spread_energy'
                   num_bad_spread_e = num_bad_spread_e + 1
                   do nn = 1, num_part
                      k = nint(part_data(2,nn))
                      e_rel = part_data(9,nn)
                      theta = part_data(10,nn)
                      phi = part_data(11,nn)
                      idb = nint(part_data(6,nn))
                      icomp_i = nint(part_data(21,nn))
                      nbin_i = nint(part_data(24,nn))
                      icomp_f = nint(part_data(1,nn))
                      nbin_f = nint(part_data(5,nn))
                      write(6,*)i,k, idb, e_rel, theta, phi, icomp_i, nbin_i,icomp_f,nbin_f
                  end do
               end if
            end if


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------   For population decay, collect some useful statistics     -----------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             if(pop_calc)then
                do k = 0, 6
                   decay_mult(k,in) = decay_mult(k,in) + real(num_part_type(k),kind=8)*tally_weight
                   decay_var(k,in) = decay_var(k,in) + real(num_part_type(k),kind=8)**2*tally_weight
                end do
             end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------     Write out events to file     -----------------------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(dump_events .and. .not. binary_event_file)then                           !   Formatted output
               write(88,'(2(1x,i10),3(1x,f10.5))')nsamp, num_part, E_in,                  &
                    nucleus(1)%Kinetic_Energy, tally_weight
               do i = 1, num_part
                  write(88,'(1x,i5,30(2x,f12.5))')i,(part_data(k,i), k = 1, n_dat)
               end do
               flush(88)
            end if
            if(dump_events .and. binary_event_file)then                           !   Unformatted output
               write(88)nsamp, num_part, E_in, nucleus(1)%Kinetic_Energy,tally_weight,    &
                    ((part_data(k,i), k = 1, n_dat), i =1, num_part)
               flush(88)
            end if

            if(biased_sampling .and. dabs(tally_weight - 1.0d0) > 1.0d-6)then
               write(6,*)'Error in sampling weight for unbiased sampling with processor #',iproc
               write(6,*)nsamp,tally_weight
               call exit_YAHFC(101)
            end if

            if(fission_decay)then       !   Fission channel
               ictype = 0
               ichann = find_channel(n_dat, dim_part, num_part, part_data, num_part_type)
               Exit_Channel(ichann)%Channel_cs(ictype,in) =                               &
                      Exit_Channel(ichann)%Channel_cs(ictype,in) + tally_weight
               fission_cs(in) = fission_cs(in) + tally_weight
               do nn = 1, num_part
                  if(part_data(2,nn) < 0.0d0)then
                     icomp_i = nint(part_data(21,nn))
                     Fiss_J_avg(icomp_i,in) = Fiss_J_avg(icomp_i,in) + part_data(3,nn)
                     Fiss_J_var(icomp_i,in) = Fiss_J_var(icomp_i,in) + part_data(3,nn)**2
                     Fiss_tally(icomp_i,in) = Fiss_tally(icomp_i,in) + tally_weight
                  end if
               end do
            else
               nn = num_part
               ichann = find_channel(n_dat, dim_part, num_part, part_data, num_part_type)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    This is an error trap. If ichann < 1, then something went wrong and           +
!------    a valid channel has not been found
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               if(ichann < 1)then
                   write(6,*)'Processor # ',iproc
                  do k = 1, 6
                     write(6,*)k, num_part_type(k)
                  end do
                  write(6,*)num_part
                  write(6,*)'Particle types'
                  do k = 1, num_part
                     write(6,*)k, part_data(2,k)
                  end do
                  write(6,*)'ichan: pop_calc = ',pop_calc
                  write(6,*)'ichan: compound = ',compound
                  write(6,*)'ichan: preeq_decay = ', preeq_decay
                  write(6,*)'ichan: direct = ',direct
                  write(6,*)'ichan: cc_decay = ',cc_decay
                  write(6,*)'ichan: dwba_decay = ',dwba_decay
                  write(6,*)'idb = ',idb
               end if

               if(hung_up)Exit_Channel(ichann)%Channel_cs(-1,in) =                        &
                          Exit_Channel(ichann)%Channel_cs(-1,in) + tally_weight

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    This is an error trap. if num_part == 0 something went wrong
!------    No decay seems to have occurred
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

               if(num_part == 0)then
                  write(6,*)'num_part == 0'
                  write(6,*)'Check: pop_calc = ',pop_calc
                  write(6,*)'Check: compound = ',compound
                  write(6,*)'Check: preeq_decay = ', preeq_decay
                  write(6,*)'Check: direct = ',direct
                  write(6,*)'Check: cc_decay = ',cc_decay
                  write(6,*)'Check: dwba_decay = ',dwba_decay
                  write(6,*)'idb = ',idb
                  call exit_YAHFC(101)
               end if

               icomp_f = nint(part_data(1,num_part))
               n_f = nint(part_data(5,num_part))

               ictype = nucleus(icomp_f)%state(n_f)%exit_lab
               inelastic = .false.
               compound_elastic = .false.
               if(.not. pop_calc .and. icomp_f == itarget)then
                  nn = 1
                  jstate = nint(part_data(5,nn))
                  idb = nint(part_data(6,nn))
                  jproj = nint(part_data(2,nn))
                  if(jproj == iproj .and. idb == 1)then           !  Inelastic transition directly to discrete states
                     icc = int(part_data(12,nn)/de_spec) + 1
                     Inelastic_cs(jstate,in) = Inelastic_cs(jstate,in) + tally_weight
                     Inelastic_count(jstate,in) = Inelastic_count(jstate,in) + 1
                     Inelastic_total(in) = Inelastic_total(in) + tally_weight
                     if(.not. xs_only)then
                        do L = 0, Ang_L_max
                           Inelastic_Ang_L(L,jstate,in) = Inelastic_Ang_L(L,jstate,in) + part_Ang_data(L,nn)*tally_weight
                        end do
                     end if

!                     if(.not. xs_only)then
!                        do nang = 1, num_theta
!                           theta = extra_angle_data(num_theta+nang,nn)
!                           x = cos(theta)
!                           jx = nint((x+1.0d0)/delta_jx_10)
!                           if(jx < 0)jx = 0
!                           if(jx > max_jx_10)jx = max_jx_10
!                           Inelastic_Ang_Dist(jx,jstate,in) = Inelastic_Ang_Dist(jx,jstate,in) +      &
!                                                              tally_weight/real(num_theta,kind=8)
!                        end do
!                     end if
                  else                                          !  Inelastic in some other fashion
                     Inelastic_cs(0,in) = Inelastic_cs(0,in) + tally_weight
                     Inelastic_count(0,in) = Inelastic_count(0,in) + 1
                     Inelastic_total(in) = Inelastic_total(in) + tally_weight
                     Exit_channel(ichann)%Channel_cs(ictype,in) =                                  &
                         Exit_channel(ichann)%Channel_cs(ictype,in) + tally_weight
                     Exit_Channel(ichann)%num_event = Exit_Channel(ichann)%num_event + 1
                     collect = .true.
                     do nn = 1, num_part                           !  loop over particles and collect spectra
                        jstate = nint(part_data(5,nn))
                        idb = nint(part_data(6,nn))
                        k = nint(part_data(2,nn))
                        if(idb == 1 .and. collect)then
                           Inelastic_cont_cs(jstate,in) = Inelastic_cont_cs(jstate,in) + tally_weight
                           collect = .false.
                        end if
!------    If it is a discrete state, add a bit of smear to avoid accidental collision with boundary of spectrum
!------    bin
!                        e_shift = 0.0d0
!                        if(idb == 1)e_shift = (2.0d0*random_32(iseed_32) - 1.0d0)*de_spec*0.5d0
!                        icc = int((part_data(12,nn) + e_shift)/de_spec) + 1
                        icc = int(part_data(12,nn)/de_spec) + 1
                        if(k >= 0 .and. (icc >= 0 .and. icc <= num_e))then
                           Exit_Channel(ichann)%part_mult(k,ictype) =                           &
                                Exit_Channel(ichann)%part_mult(k,ictype) + tally_weight
                           if(.not. xs_only)then
                              Exit_Channel(ichann)%Spect(k,ictype)%E_spec(icc) =                &
                                  Exit_Channel(ichann)%Spect(k,ictype)%E_spec(icc) +            &
                                  tally_weight/de_spec
                              Exit_Channel(ichann)%Spect(k,ictype)%E_count(icc) =               &
                                  Exit_Channel(ichann)%Spect(k,ictype)%E_count(icc) + 1
                              do nang = 1, num_theta
                                 theta = extra_angle_data(num_theta+nang,nn)
                                 x = cos(theta)
                                 jx = nint((x+1.0d0)/delta_jx_10)
                                 if(jx < 0)jx = 0
                                 if(jx > max_jx_10)jx = max_jx_10
                                 Exit_Channel(ichann)%Spect(k,ictype)%E_Ang_Dist(jx,icc) =      &
                                     Exit_Channel(ichann)%Spect(k,ictype)%E_Ang_Dist(jx,icc) +  &
                                     tally_weight/(delta_jx_10*de_spec)/real(num_theta,kind=8)
                                 Exit_Channel(ichann)%Spect(k,ictype)%Ang_Dist(jx) =            &
                                     Exit_Channel(ichann)%Spect(k,ictype)%Ang_Dist(jx) +        &
                                     tally_weight/delta_jx_10/real(num_theta,kind=8)
                              end do
                           end if
                        end if
!*****************************************************************************************************************
!-----    This line will remove further discrete gammas from the spectra files --  if they are wanted comment out
!*****************************************************************************************************************
                        if(.not. collect)exit
                     end do
                  end if
               else

                  Exit_channel(ichann)%Channel_cs(ictype,in) =                                     &
                      Exit_channel(ichann)%Channel_cs(ictype,in) + tally_weight
                  Exit_Channel(ichann)%num_event = Exit_Channel(ichann)%num_event + 1
                  do nn = 1, num_part
                     k = nint(part_data(2,nn))
!------    If it is a discrete state, add a bit of smear to avoid accidental collision with boundary of spectrum
!------    bin
                     icc = int(part_data(12,nn)/de_spec) + 1
                     if(k >= 0 .and. (icc >= 0 .and. icc <= num_e))then
                        Exit_Channel(ichann)%part_mult(k,ictype) =                              &
                             Exit_Channel(ichann)%part_mult(k,ictype) + tally_weight
                        if(.not. xs_only)then
                           Exit_Channel(ichann)%Spect(k,ictype)%E_spec(icc) =                   &
                               Exit_Channel(ichann)%Spect(k,ictype)%E_spec(icc) +               &
                               tally_weight/de_spec
                           Exit_Channel(ichann)%Spect(k,ictype)%E_count(icc) =                  &
                               Exit_Channel(ichann)%Spect(k,ictype)%E_count(icc) + 1
                           do nang = 1, num_theta
                              theta = extra_angle_data(num_theta+nang,nn)
                              x = cos(theta)
                              jx = nint((x+1.0d0)/delta_jx_10)
                              if(jx < 0)jx = 0
                              if(jx > max_jx_10)jx = max_jx_10
                              Exit_Channel(ichann)%Spect(k,ictype)%E_Ang_Dist(jx,icc) =         &
                                  Exit_Channel(ichann)%Spect(k,ictype)%E_Ang_Dist(jx,icc) +     &
                                  tally_weight/(delta_jx_10*de_spec)/real(num_theta,kind=8)
                              Exit_Channel(ichann)%Spect(k,ictype)%Ang_Dist(jx) =               &
                                  Exit_Channel(ichann)%Spect(k,ictype)%Ang_Dist(jx) +           &
                                  tally_weight/delta_jx_10/real(num_theta,kind=8)
                           end do
                        end if
                     end if
                  end do
               end if
            end if

            if(.not. fission_decay)then
               do i = 1, num_part
                  k = nint(part_data(2,i))
                  if(k >= 0 .and. k <= 6)then                      !  k < 0 is internal conversion
                     icc = int(part_data(15,i)/de_spec2)
                     if(icc >= 0 .and. icc <= num_e)then
                        x_particle_spectrum(icc,k) = x_particle_spectrum(icc,k) + tally_weight
                     else
                        write(6,*)'iproc= ',iproc
                        write(6,*)'nsamp = ',nsamp
                        write(6,*)'num_part = ', num_part
                        write(6,*)'Check: pop_calc = ',pop_calc
                        write(6,*)'Check: compound = ',compound
                        write(6,*)'Check: preeq_decay = ', preeq_decay
                        write(6,*)'Check: direct = ',direct
                        write(6,*)'Check: cc_decay = ',cc_decay
                        write(6,*)'Check: dwba_decay = ',dwba_decay
                        write(6,*)'idb = ',idb
                        write(6,*)'paritcle type = ',int(part_data(2,i))
                        write(6,*)'energy = ',part_data(15,i)
                        write(6,*)'icc = ',icc,' num_e = ',num_e
                        write(6,*)'Issue with part_Lab_spectrum inde > num_spect_e'
                     end if
                  end if
               end do
            end if


            if(preeq_decay)then
               k = nint(part_data(2,1))
               if(k >= 0 .and. k <= 7)then
                  preeq_css(k,in) = preeq_css(k,in) + tally_weight
                  if(.not. xs_only)then
                     icc = int(part_data(15,1)/de_spec2)
                     if(icc >= 0 .and. icc <= num_e)then
                        preeq_spect_full(k,icc) = preeq_spect_full(k,icc) + tally_weight
                        if(.not. fission_decay)                                                          &
                           preeq_spect(k,icc) = preeq_spect(k,icc) + tally_weight
                     end if
                  end if
               end if
            end if
            if(direct .and. cc_decay .and. .not. fission_decay .and. .not. xs_only)then
               icc = int(part_data(15,1)/de_spec2)
               if(icc >= 0 .and. icc <= num_e)direct_spectrum(icc) = direct_spectrum(icc) + tally_weight
            end if
            if(direct .and. dwba_decay .and. .not. fission_decay .and. .not. xs_only)then
               icc = int(part_data(15,1)/de_spec2)
               if(icc >= 0 .and. icc <= num_e)dwba_spectrum(icc) = dwba_spectrum(icc) + tally_weight
            end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Finished Monte Carlo Sampling Loop     -----------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         end do                                   !   end do nsamp loop
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Finished Monte Carlo Sampling Loop     -----------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Now, if using MPI, sum over the MPI processes    -------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
#if(USE_MPI == 1)
      if(nproc > 1)then
         call MPI_Barrier(icomm, ierr)
         num_data = 1
         call MPI_Allreduce(MPI_IN_PLACE, tally_norm, num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
         num_data = (num_e+1)*7
         call MPI_Allreduce(MPI_IN_PLACE, x_particle_spectrum,                                     &
                            num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
!---  Data for primary_decay_only calculations, populations for each bin
         if(primary_decay)then
            do icomp = 1, num_comp
               do Ix_i = 1, nucleus(icomp)%j_max
                  do ip_i = 0, 1
                     do nn = 1, nucleus(icomp)%nbin
                        num_data = 1
                        call MPI_Allreduce(MPI_IN_PLACE, nucleus(icomp)%bins(Ix_f,ip_f,nn)%pop,    &
                                           num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
                     end do
                  end do
               end do
               do nn = 1, nucleus(icomp_f)%num_discrete
                  num_data = 1
                  call MPI_Allreduce(MPI_IN_PLACE, nucleus(icomp)%state(nn)%pop,                   &
                                     num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
               end do
            end do
         end if
!---  Fission information
         if(fission)then
            num_data = 1
            call MPI_Allreduce(MPI_IN_PLACE, fission_cs(in),                                       &
                                       num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
            num_data = num_comp
            call MPI_Allreduce(MPI_IN_PLACE, FIss_J_avg(1,in),                                     &
                               num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
            num_data = num_comp
            call MPI_Allreduce(MPI_IN_PLACE, FIss_J_var(1,in),                                     &
                               num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
            num_data = num_comp
            call MPI_Allreduce(MPI_IN_PLACE, FIss_Tally(1,in),                                     &
                               num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
         end if
         if(.not. pop_calc)then
!---  Inelastic scattering data
            num_data = 1
            call MPI_Allreduce(MPI_IN_PLACE, Inelastic_Total(in),                                  &
                               num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
            num_data = (nucleus(itarget)%num_discrete+1)
            call MPI_Allreduce(MPI_IN_PLACE, Inelastic_cs(0,in),                                   &
                               num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
            if(.not. xs_only)then
               do j = 0, nucleus(itarget)%num_discrete
                  num_data = Ang_L_max
                  call MPI_Allreduce(MPI_IN_PLACE, Inelastic_Ang_L(0,j,in),                           &
                                     num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
               end do
            end if

            num_data = (nucleus(itarget)%num_discrete+1)
            call MPI_Allreduce(MPI_IN_PLACE, Inelastic_count(0,in),                                &
                               num_data, MPI_INTEGER, MPI_SUM, icomm, ierr)
!---   Pre-equilibrium data
            num_data = 7
            call MPI_Allreduce(MPI_IN_PLACE, preeq_css(0,in),                                      &
                               num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
            if(.not. xs_only)then
               num_data = 7*(num_e+1)
               call MPI_Allreduce(MPI_IN_PLACE, preeq_spect,                                       &
                                  num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
               num_data = 7*(num_e+1)
               call MPI_Allreduce(MPI_IN_PLACE, preeq_spect_full,                                  &
                                  num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
!---   Direct and DWBA spectra
               num_data = (num_e+1)
               call MPI_Allreduce(MPI_IN_PLACE, direct_spectrum,                                   &
                                  num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
               num_data = (num_e+1)
               call MPI_Allreduce(MPI_IN_PLACE, dwba_spectrum,                                     &
                                  num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
            end if
         end if
!---   Exit_Channel data
         do i = 1, num_channels
            num_data = 1
            call MPI_Allreduce(MPI_IN_PLACE,Exit_Channel(i)%num_event,                             &
                               num_data, MPI_INTEGER, MPI_SUM, icomm, ierr)
            num_s = Exit_Channel(i)%num_cs
            num_data = num_s + 3                                !  allocated (-2:Exit_Channel(i)%num_cs)
            call MPI_Allreduce(MPI_IN_PLACE,Exit_Channel(i)%Channel_cs(-2,in),                     &
                               num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
            num_data = 7*(num_s+2)
            call MPI_Allreduce(MPI_IN_PLACE,Exit_Channel(i)%part_mult(0,-1),                    &
                               num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
            if(.not. xs_only)then
               do k = 0, 6
                  do ictype = -1, num_s
                     num_data = num_e + 1
                     call MPI_Allreduce(MPI_IN_PLACE,Exit_Channel(i)%Spect(k,ictype)%E_spec,    &
                                  num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
                     num_data = num_e + 1
                     call MPI_Allreduce(MPI_IN_PLACE,Exit_Channel(i)%Spect(k,ictype)%E_count,   &
                                  num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
                     num_data = (max_jx_10+1)*(num_s+2)
                     call MPI_Allreduce(MPI_IN_PLACE,Exit_Channel(i)%Spect(k,ictype)%E_Ang_Dist,&
                                  num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
                     num_data = (max_jx_10+1)
                     call MPI_Allreduce(MPI_IN_PLACE,Exit_Channel(i)%Spect(k,ictype)%Ang_Dist,  &
                                  num_data, MPI_REAL8, MPI_SUM, icomm, ierr)
                  end do
               end do
            end if
         end do
      end if
#endif
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Finished colelcting data from the MPI processes   -------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!---  Data for primary_decay_only calculations, populations for each bin
      if(primary_decay)then
         do icomp = 1, num_comp
            do jj = 0, nucleus(icomp)%j_max
               do ip = 0, 1
                  do nn = 1, nucleus(icomp)%nbin
                     nucleus(icomp)%bins(jj,ip,nn)%pop =                                       &
                        nucleus(icomp)%bins(jj,ip,nn)%pop/tally_norm
                  end do
               end do
            end do
            do nn = 1, nucleus(icomp)%num_discrete
               nucleus(icomp)%state(nn)%pop = nucleus(icomp)%state(nn)%pop/tally_norm
            end do
         end do
      end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------  Collect spectra                                    -------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         do k = 0, 6
            smear(k) = 0.2d0
         end do
!-----   Create core file name for spectra files

         ifile = core_file
         file_name(ifile:ifile+4) = '_Ein_'
         ifile =ifile + 4
         if(e_in < 10.0)then
            file_name(ifile+1:ifile+2)='00'
            ifile=ifile+2
            write(file_name(ifile+1:ifile+6),'(f6.4)')e_in
            ifile=ifile+6
         elseif(e_in < 100.0)then
            file_name(ifile+1:ifile+2)='0'
            ifile=ifile+1
            write(file_name(ifile+1:ifile+7),'(f7.4)')e_in
            ifile=ifile+7
         elseif(e_in < 1000.0)then
            write(file_name(ifile+1:ifile+7),'(f8.4)')e_in
            ifile=ifile+8
         end if

         max_num = nint(ee_max/de_spec2) + 1
         if(max_num > num_e) max_num = num_e

         temp_cs = 0.0d0
         if(pop_calc)then
            temp_cs = reaction_cs(in)
         else
            if(iproj == 1)then
               temp_cs = reaction_cs(in) + SE_cs(in)
            elseif(iproj > 1)then
               temp_cs = reaction_cs(in)
            end if
         end if

         xnum_eff = tally_norm*temp_cs/reaction_cs(in)                 !    number of effective events with this cross section
         xnum_elastic = xnum_eff - tally_norm
         num_eff = int(xnum_eff)
         num_elastic = int(xnum_elastic)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Here, collect the spectra for particles emitted and do not distinguish by channel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(.not. xs_only)then
            if(iproj == 1 .and. .not. pop_calc)then
!-------   Throw num_elastic events to simulate elastic scattering
               do i = 1, num_elastic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------     Reset the boost matrix
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  Boost_COM(0:3,0:3) = 0.0d0
                  do j = 0, 3
                     Boost_COM(j,j) = 1.0d0
                  end do
                  nsamp = nsamp + 1
                  num_part = 1
                  E_f = E_in*mass_target/(mass_target + mass_proj)
                  part_data(1,num_part) = real(icomp_f,kind=8)
                  part_data(2,num_part) = real(iproj,kind=8)
                  part_data(3,num_part) = nucleus(itarget)%state(istate)%spin
                  part_data(4,num_part) = nucleus(itarget)%state(istate)%parity
                  part_data(5,num_part) = real(0,kind=8)
                  part_data(6,num_part) = real(1,kind=8)
                  part_data(7,num_part) = 0.0d0
                  part_data(8,num_part) = 1.0d0
                  part_data(9,num_part) = E_f
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------   get theta
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               ran = random_64(iseed_64)
                  ran = random_32(iseed_32)
                  do j = 1, ixx_max
                     x = real(j,kind=8)*delta_ix - 1.0d0
                     if(SE_prob(j) > ran)exit
                  end do
!                  x = x - random_64(iseed_64)*delta_ix*0.9999999d0
                  x = x - random_32(iseed_32)*delta_ix*0.9999999d0
                  theta_0 = acos(x)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------   Now phi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                  ran = random_64(iseed_64)
                  ran = random_32(iseed_32)
                  phi_0 = ran*2.0d0*pi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------   Boost to COM and Lab frame
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  mass_1 = particle(iproj)%Mass
                  call Boost_frame(e_f, mass_1, mass_2, theta_0, phi_0,                            &
                                   Boost_Lab, Boost_COM, T_1, theta, phi,                          &
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
                  part_data(19,num_part) = 1.0d0
                  part_data(20,num_part) = nucleus(itarget)%state(istate)%energy
                  part_data(21,num_part) = icomp_i
                  part_data(21,num_part) = icomp_i
                  part_data(22,num_part) = nucleus(itarget)%state(istate)%spin
                  part_data(23,num_part) = nucleus(itarget)%state(istate)%parity
                  part_data(24,num_part) = istate
                  nucleus(icomp_f)%Kinetic_Energy = T_2

                  icc = int(T_L/de_spec2)
                  if(icc <= num_e)x_particle_Spectrum(icc,iproj) = x_particle_Spectrum(icc,iproj) + 1.0d0
                  x_particle_cs(in,iproj) = x_particle_cs(in,iproj) + 1.0d0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------     Write out events to file     -----------------------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  if(dump_events .and. .not. binary_event_file .and. iproc == 0)then                           !   Formatted output
                     write(88,'(2(1x,i10),3(1x,f10.5))')nsamp, num_part, E_in,                     &
                         nucleus(1)%Kinetic_Energy, tally_weight
                     do j = 1, num_part
                        write(88,'(1x,i5,21(2x,f12.5))')j,(part_data(m,j), m = 1, n_dat)
                     end do
                     flush(88)
                  end if
                  if(dump_events .and. binary_event_file .and. iproc == 0)then                           !   Unformatted output
                     write(88)nsamp, num_part, E_in, nucleus(1)%Kinetic_Energy,tally_weight,       &
                             ((part_data(m,j), m = 1, n_dat), j =1, num_part)
                     flush(88)
                  end if
               end do
            end if
            do k = 0, 6
               x_particle_cs(in,k) = x_particle_cs(in,k)*temp_cs/xnum_eff
               xnorm = 0.0d0
               do i = 0, num_e
                  x_particle_Spectrum(i,k) = x_particle_Spectrum(i,k)*temp_cs/xnum_eff/de_spec2
                  xnorm = xnorm + x_particle_Spectrum(i,k)*de_spec2
               end do
            end do
         end if

         if(dump_events)close(unit=88)


!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------  Collect and print particle spectra ssociated with various reaction types,  ----+
!-------  pre-equilibirum, direct, dwba, etc.                                        ----+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(.not. pop_calc .and. iproj > 0 .and. .not. xs_only)then
            do i = 0, num_e
               direct_spectrum(i) = direct_spectrum(i)*reaction_cs(in)/tally_norm/de_spec2
               dwba_spectrum(i) = dwba_spectrum(i)*reaction_cs(in)/tally_norm/de_spec2
               if(PREEQ_Model > 0)then
                  do k = 0, 6
                     preeq_spect(k,i) = preeq_spect(k,i)*reaction_cs(in)/tally_norm/de_spec2
                     preeq_spect_full(k,i) = preeq_spect_full(k,i)*reaction_cs(in)/tally_norm/de_spec2
                  end do
               end if
            end do
            if(iproc == 0)then
               if(iproj > 0)then

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Pre-equilibrium decays                                                         +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  if(PREEQ_Model > 0)call print_preeq_spectra(ilib_dir, lib_dir, ifile, file_name, &
                                                              e_in, de_spec2, num_e,               &
                                                              Preeq_spect, Preeq_spect_full, smear)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Direct and DWBA states decays                                                  +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  if(OpticalCS%numcc > 1)call print_direct_spectra(ilib_dir, lib_dir,              &
                                                                   ifile, file_name,               &
                                                                   e_in, de_spec2, max_num, num_e, &
                                                                   direct_Spectrum, dwba_Spectrum, &
                                                                   smear)
               end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Particle spectra (iproj,Xk)                                                    +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               call print_x_particle_spectra(ilib_dir, lib_dir, ifile, file_name,                  &
                                             e_in, de_spec2, max_num, num_e,                       &
                                             x_particle_spectrum, smear)
            end if
         end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   If primary_decay_only, Print populations following primary decay               +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(primary_decay)then
            call print_primary_decay(ilib_dir, lib_dir, ifile, file_name, e_in, reaction_cs(in))
         end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ---- finished with printing spectra                                                    +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------------   Distribute cs hung in discrete states to final states            +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(biased_sampling)then
            if(print_me)then
               write(6,*)'Statistics on continuous bins with hung decays'
               write(6,'(''                                         #hung              #total      fraction'')')
               write(6,'(''                              ----------------    ----------------    ----------'')')
            end if
            if(iproc == 0 .and. print_output)then
               write(13,*)'Statistics on continuous bins with hung decays'
               write(13,'(''                                         #hung              #total      fraction'')')
               write(13,'(''                              ----------------    ----------------    ----------'')')
            end if
         else
            if(print_me)then
               write(6,*)'Statistics on continuous bins with hung decays'
               write(6,'(''                                         #hung              #total      fraction       Num events'')')
               write(6,'(''                              ----------------    ----------------    ----------     ------------'')')
            end if
            if(iproc == 0 .and. print_output)then
               write(13,*)'Statistics on continuous bins with hung decays'
               write(13,'(''                                         #hung              #total      fraction       Num events'')')
               write(13,'(''                              ----------------    ----------------    ----------     ------------'')')
            end if
         end if

         do ichann = 1, num_channels
            hang = Exit_Channel(ichann)%Channel_cs(-1,in)
            sum = 0.0d0
            do j = 1, Exit_Channel(ichann)%num_cs
               sum = sum + Exit_Channel(ichann)%Channel_cs(j,in)
            end do
            ratio = 0.0d0
            ilast = index(Exit_Channel(ichann)%Channel_Label,' ') - 1
            if(sum > 0.0d0)ratio = Exit_Channel(ichann)%Channel_cs(-1,in)/sum
            if(biased_sampling)then
               if(print_me)write(6,'(''Channel: '',i4,1x,a12,2(4x,f16.3),4x,f10.4,''%'')')          &
                   ichann, Exit_Channel(ichann)%Channel_Label(1:ilast), hang, sum, ratio*100.
               if(iproc == 0 .and. print_output)write(13,'(''Channel: '',i4,1x,a12,2(4x,f16.3),4x,f10.4,''%'')')         &
                   ichann, Exit_Channel(ichann)%Channel_Label(1:ilast), hang, sum, ratio*100.
            else
               if(print_me)write(6,'(''Channel: '',i4,1x,a12,2(4x,f16.3),4x,f10.4,''%'',4x,i12)')   &
                   ichann, Exit_Channel(ichann)%Channel_Label(1:ilast),                             &
                   hang, sum, ratio*100., Exit_Channel(ichann)%num_event
               if(iproc == 0 .and. print_output)write(13,'(''Channel: '',i4,1x,a12,2(4x,f16.3),4x,f10.4,''%'',4x,i12)')  &
                   ichann, Exit_Channel(ichann)%Channel_Label(1:ilast),                             &
                   hang, sum, ratio*100., Exit_Channel(ichann)%num_event
            end if
         end do
         if(print_me)write(6,*)'Hung bins are forced to decay to discrete states randomly'
         if(iproc == 0 .and. print_output)write(13,*)'Hung bins are forced to decay to discrete states randomly'

         icharr = 32
         out_buff(1:icharr) = 'For more details consult file : '
         il1 = index(out_file,' ')-1
         out_buff(icharr+1:icharr+il1)=out_file(1:il1)
         icharr = icharr + il1
         out_buff(icharr+1:icharr+15) = '.sample_trouble'
         icharr = icharr + 15
         if(re_baseline > 0 .and. print_me)then
           write(6,*)
           write(6,'(''*********************************************************************************'')')
           write(6,'(i10,'' Decays out of '',i10,'' had to be re-baselined due to energy'')')re_baseline, num_mc_samp
           write(6,'(''conservation requirements'')')
           write(6,'(a)')out_buff(1:icharr)
           write(6,'(''This is not an indication of failure, but is a consequence of MC decay in'')')
           write(6,'(''the bin structure so that a previously energetically allowed decay becomes '')')
           write(6,'(''forbidden. This is more common when discrete states are allowed above E_cut'')')
           write(6,'(''*********************************************************************************'')')
         end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    Print particle spectra for each Channel                                        +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        if(fission)fission_cs(in) = fission_cs(in)/tally_norm
        if(PREEQ_Model > 0)then
           do k = 0, 6
              preeq_css(k,in) = preeq_css(k,in)/tally_norm
           end do
        end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Compute cross sections for inelastic channels and normalize angular distributions
!---------    Inelastic cross sections to each state
!---------    j = 0 is inelastic to continuous energy bins
!---------    j = istate is compound elastic
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        sum_inelastic = 0.0d0
        if(.not. pop_calc)then
           do j = 1, nucleus(itarget)%num_discrete
              sum_inelastic = sum_inelastic + Inelastic_cs(j,in)/tally_norm
              if(.not. xs_only .and. Inelastic_cs(j,in) > 1.0d-8)then
                  do L = 0, Ang_L_max
                     Inelastic_Ang_L(L,j,in) = Inelastic_Ang_L(L,j,in)/Inelastic_cs(j,in)
                  end do
              end if
              Inelastic_cs(j,in) = Inelastic_cs(j,in)*reaction_cs(in)/tally_norm
              Inelastic_cont_cs(j,in) = Inelastic_cont_cs(j,in)*reaction_cs(in)/tally_norm
   
           end do
           Inelastic_cs(0,in) = Inelastic_cs(0,in)*reaction_cs(in)/tally_norm
        end if



!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------    Now the other Channels    -------------------------------------------------
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        check_sum = 0.0d0
!
        if(.not. event_generator .and. .not. xs_only)then
           if(.not.allocated(xvalue))allocate(xvalue(0:max_jx_10))
           do jx = 0, max_jx_10
              xvalue(jx) = real(jx,kind=8)*delta_jx_10 - 1.0d0
           end do
!---  Loop over channels
           do i = 1, num_channels
              inuc = Exit_Channel(i)%Final_nucleus
              n_min = 1
              if(fission)n_min = 0
              do n = n_min, Exit_Channel(i)%num_cs
                 do k = 0, 6                                  !   Loop over particle types
!-----   New approach to Spectrum and Angular Distributions
                    do icc = 0, num_e
                       if(Exit_Channel(i)%Spect(k,n)%E_count(icc) == 0)cycle
                       Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(0,icc) =                            &
                          Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(0,icc)*2.0d0
                       Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(max_jx_10,icc) =                    &
                          Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(max_jx_10,icc)*2.0d0
                       sum = 0.0d0
!------    Normalize for each energy
                       do jx = 1, max_jx_10
                          sum = sum + (Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(jx,icc) +           &
                                       Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(jx-1,icc))*         &
                                      delta_jx_10*0.5d0
                       end do
                       if(sum > 1.0d-8)then
                          do jx = 0, max_jx_10
                             Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(jx,icc) =                     &
                                  Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(jx,icc)/sum
                          end do
                       end if

                       sum = Exit_Channel(i)%Spect(k,n)%E_count(icc)
!-----   Set max L based on number of counts in energy bin
                       Exit_Channel(i)%Spect(k,n)%E_Ang_L_max(icc) = 8
                       if(sum < 9000)Exit_Channel(i)%Spect(k,n)%E_Ang_L_max(icc) = 6
                       if(sum < 6000)Exit_Channel(i)%Spect(k,n)%E_Ang_L_max(icc) = 4
                       if(sum < 3000)Exit_Channel(i)%Spect(k,n)%E_Ang_L_max(icc) = 2
                       if(sum < 1000)Exit_Channel(i)%Spect(k,n)%E_Ang_L_max(icc) = 0
!-----   Calculate average, avgerage difference from average, and expected difference
!-----   Check if data have enough statistics to warrant more Legendre coefficients
                       if(Exit_Channel(i)%Spect(k,n)%E_Ang_L_max(icc) > 0)then
                          avg = 0.0d0
                          do jx = 0, max_jx_10
                             avg = avg + Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(jx,icc)
                          end do
                          avg = avg/real(max_jx_10 + 1,kind=8)
                          avg_diff = 0.0d0
                          do jx = 0, max_jx_10
                             avg_diff = avg_diff + abs(Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(jx,icc) - avg)
                          end do
                          avg_diff = avg_diff/real(max_jx_10+1,kind=8)
                          avg_diff = avg_diff/avg
                          expected_diff = real(Exit_Channel(i)%Spect(k,n)%E_count(icc),kind=8)/  &
                                          real(max_jx_10 + 1,kind=8)
                          expected_diff = sqrt(expected_diff)/expected_diff/sqrt(real(num_theta,kind=8))
                          if(avg_diff < 2.5d0*expected_diff)Exit_Channel(i)%Spect(k,n)%E_Ang_L_max(icc) = 4
                          if(avg_diff < 1.75d0*expected_diff)Exit_Channel(i)%Spect(k,n)%E_Ang_L_max(icc) = 2
                          if(avg_diff < expected_diff)Exit_Channel(i)%Spect(k,n)%E_Ang_L_max(icc) = 0
                       end if

                       call Legendre_expand(max_jx_10+1,xvalue(0),                                  &
                                            Exit_Channel(i)%Spect(k,n)%E_Ang_Dist(0,icc),           &
                                            Exit_Channel(i)%Spect(k,n)%E_Ang_L_max(icc),            &
                                            Exit_Channel(i)%Spect(k,n)%E_Ang_L(0,icc))

!----------------------------------------------------------------------------------------------------
                    end do                       !  End of Loop over icc (outgoing energies
                 end do                          !  End of Loop over k, particles

              end do
           end do
!------    Deallocate temporary array with x-values
           if(allocated(xvalue))deallocate(xvalue)
        end if

        if(.not. event_generator)then
           do i = 1, num_channels
              inuc = Exit_Channel(i)%Final_nucleus
              n_min = 1
              if(fission)n_min = 0
              do n = n_min, Exit_Channel(i)%num_cs
                 do k = 0, 6                                  !   Loop over particle types
                    if(Exit_Channel(i)%Channel_cs(n,in) >= 1.0d-10)then
                       Exit_Channel(i)%part_mult(k,n) =                                            &
                          Exit_Channel(i)%part_mult(k,n)/Exit_Channel(i)%Channel_cs(n,in)
                    else
                       Exit_Channel(i)%part_mult(k,n) = 0.0d0
                    end if
                 end do
              end do
              if(track_gammas)then
                 do n = 1, nucleus(inuc)%num_discrete
                    do j = 1, nucleus(inuc)%state(n)%nbranch
                       Exit_channel(i)%state(n)%cs(j,in) =                                         &
                              Exit_channel(i)%state(n)%cs(j,in)/tally_norm
                   end do
                 end do
              end if
           end do
        end if

!-----   Normalize the Channel cross sections and collect check_sum to make sure all
!-----   events are accounted for
        do i = 1, num_channels
           n_min = 1
           if(fission)n_min = 0
           do n = n_min, Exit_Channel(i)%num_cs
              Exit_Channel(i)%Channel_cs(n,in) =                                                   &
                   Exit_Channel(i)%Channel_cs(n,in)/tally_norm
              if(n > 0)check_sum = check_sum +                                                     &
                           Exit_Channel(i)%Channel_cs(n,in)
           end do
        end do


        check_sum = check_sum + sum_inelastic
        if(fission)check_sum = check_sum + fission_cs(in)
        if(print_me)then
           write(6,*)
           write(6,*)'Check that all events are accounted for'
           write(6,'('' Check_sum = '',f12.8)') check_sum
           if(dabs(check_sum - 1.0d0) <= 1.0d-8)write(6,*)'All events are accounted for'
        end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------  if we are tracking gammas, print them out now, rather than with other channels
!-------  after we finish all the incident energies. Print to a separate file for each
!-------  incident energy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(track_gammas .and. .not. event_generator .and. iproc == 0)                          &
             call print_channel_gammas(itarget, ilab, file_lab, ilib_dir, lib_dir,              &
                                       ch_par, in, e_in, reaction_cs(in), write_error)

 1901    continue

         if(.not.pop_calc)then
            if(associated(channel_prob))nullify(channel_prob)
            if(associated(ichannel))nullify(ichannel)
            if(iproj > 0)then
               spin_target=nucleus(itarget)%state(istate)%spin
               spin_proj=particle(iproj)%spin
               sp1=abs(spin_target-spin_proj)
               sp2=spin_target+spin_proj
               isp=nint(sp2-sp1)
               isp_max = nint(2.0d0*spin_proj)
               Ix = 0
               do l = 0, particle(iproj)%lmax
                  xj_min = abs(real(l,kind=8) - spin_proj)
                  xj_max = abs(real(l,kind=8) + spin_proj)
                  isp = nint(xj_max - xj_min)
                  xj = real(l,kind=8) - spin_proj
                  do is = 0, isp_max
                     xj = xj + real(is,kind=8)
                     if(xj < 0.0d0)cycle
                     xI_min = abs(xj - spin_target)
                     xI_max = xj + spin_target
                     Ix_min = max(nint(xI_min-nucleus(1)%jshift),0)
                     Ix_max = min(nint(xI_max-nucleus(1)%jshift),nucleus(1)%j_max)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------      Loop over possible entrance orbital angular momenta      +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                     do Ix = Ix_min, Ix_max
                        do if1 = 1, nucleus(1)%num_decay + 1
                           if(allocated(Channel(l,is,Ix)%Channel_decay(if1)%decay_prob))              &
                              deallocate(Channel(l,is,Ix)%Channel_decay(if1)%decay_prob)
                           if(allocated(Channel(l,is,Ix)%Channel_decay(if1)%decay_list))              &
                              deallocate(Channel(l,is,Ix)%Channel_decay(if1)%decay_list)
                        end do
                        if(allocated(Channel(l,is,Ix)%Channel_decay))                                 &
                              deallocate(Channel(l,is,Ix)%Channel_decay)
                        if(allocated(Channel(l,is,Ix)%Channel_prob))                                  &
                              deallocate(Channel(l,is,Ix)%Channel_prob)
                     end do
                  end do
               end do
               if(allocated(Channel))deallocate(Channel)
            end if
         end if

         if(pop_calc)then
            do k = 0, 6
               decay_mult(k,in) = decay_mult(k,in)/tally_norm
               decay_var(k,in) = decay_var(k,in)/tally_norm
               decay_var(k,in) = sqrt(decay_var(k,in) - decay_mult(k,in)**2)
            end do
         end if

         ilab = ntar
         file_lab(1:ntar) = target_label(1:ntar)
         ilab = ilab + 1
         file_lab(ilab:ilab) = '('
         ilab = ilab + 1
         file_lab(ilab:ilab) = particle(iproj)%label
         ilab = ilab + 1
         file_lab(ilab:ilab) = ','

         directory(1:ilib_dir) = lib_dir(1:ilib_dir)
         idir = ilib_dir + 1
         directory(idir:idir) = '/'

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Channel cross section data      ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         write_error = .false.
         if(iproc == 0)call print_channels(in, itarget, istate, ilab, file_lab,                 &
                                           ilib_dir, lib_dir, ch_par,                           &
                                           num_energies, num_e, max_jx_20, delta_jx_20,         &
                                           de_spec, cs_threshold, reaction_cs, write_error)
#if(USE_MPI == 1)
         if(nproc > 1)then
            call MPI_Barrier(icomm, ierr)
            call MPI_BCAST(write_error, 1, MPI_LOGICAL, 0, icomm, ierr)
         end if
#endif
         if(write_error)call exit_YAHFC(51)

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Finished incident Energy Loop     ----------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      end do               !   End of incident energy loop
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Finished incident Energy Loop     ----------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

      if(iproc == 0)write(6,*)'Finished energy loop'


!-------   Notify user about events that failed energy conservation limits
!     close(unit=28)     !   Close units where energy nonconservation is written
!     close(unit=400)

     if(num_bad_samp_e > 0 .or. num_bad_spread_e > 0)then
        ifile = index(out_file,' ') - 1
        if(num_bad_samp_e > 0)then
           if(print_me)write(6,*)'There were ',num_bad_samp_e,' events violated energy conservation during primary sampling'
!           if(print_me)write(6,*)'Events written to ',out_file(1:ifile)//'.bad_samp_energy'
        end if
        if(num_bad_spread_e > 0)then
           if(print_me)write(6,*)'There were ',num_bad_samp_e,' events violated energy conservation applying extra angles'
!           if(print_me)write(6,*)'Events written to ',out_file(1:ifile)//'.bad_spread_energy'
        end if
     else
        if(print_me)write(6,*)'There were no events that violated energy conservation at 0.1 keV level'
!----   Delete the bad-energy files as they are empty
!        unix_command(1:132) = ' '
!        unix_command(1:3) = ' rm '
!        ilast = 3
!        unix_command(ilast+1:ilast+ifile) = out_file(1:ifile)
!        ilast = ilast + ifile
!        unix_command(ilast+1:ilast+16) = '.bad_samp_energy'
!        ilast = ilast + 16
!        call system(unix_command(1:ilast)
!        ilast = ifile + 3
!        unix_command(ilast+1:ilast+18) = '.bad_spread_energy'
!        ilast = ilast + 18
!        call system(unix_command(1:ilast)
     end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Deallocate arrays needed for pre-equilibrium model
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     if(.not. pop_calc .and. PREEQ_Model == 1)then
        if(allocated(dWk))deallocate(dWk)
        if(allocated(Wk))deallocate(Wk)
        if(allocated(W))deallocate(W)
        if(allocated(tau))deallocate(tau)
        if(allocated(taup))deallocate(taup)
        if(allocated(Pre_Prob))deallocate(Pre_Prob)
        if(allocated(Gam_p))deallocate(Gam_p)
        if(allocated(Gam_n))deallocate(Gam_n)
        if(allocated(Gam_pp))deallocate(Gam_pp)
        if(allocated(Gam_pn))deallocate(Gam_pn)
        if(allocated(Gam_0_pn))deallocate(Gam_0_pn)
        if(allocated(Gam_0_np))deallocate(Gam_0_np)
        if(allocated(Lamb_p_p))deallocate(Lamb_p_p)
        if(allocated(Lamb_p_n))deallocate(Lamb_p_n)
        if(allocated(Lamb_0_pn))deallocate(Lamb_0_pn)
        if(allocated(Lamb_0_np))deallocate(Lamb_0_np)
        if(allocated(Vwell))deallocate(Vwell)
     end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  Output data to Library Files     -----------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(event_generator)then
         write(6,*)'****************************************************'
         write(6,*)'----   Finished simulating and writing events.   ---'
         write(6,*)'----   All finished.                             ---'
         write(6,*)'****************************************************'
#if(USE_MPI == 1)
         call MPI_Finalize(ierr)
         call exit
#else
         call exit_YAHFC(901)
#endif
      end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Reaction Label       -----------------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ilab = ntar
      file_lab(1:ntar) = target_label(1:ntar)
      ilab = ilab + 1
      file_lab(ilab:ilab) = '('
      ilab = ilab + 1
      file_lab(ilab:ilab) = particle(iproj)%label
      ilab = ilab + 1
      file_lab(ilab:ilab) = ','

      directory(1:ilib_dir) = lib_dir(1:ilib_dir)
      idir = ilib_dir + 1
      directory(idir:idir) = '/'

      if(print_libraries)then
         if(iproc == 0)write(6,*)'Writing data libraries'
         if(.not. pop_calc)then
            if(iproj > 0)then        !   don't print for photons
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Total cross section
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               write_error = .false.
               if(iproc == 0)call print_total_cs(itarget, istate, ilab, file_lab,                &
                                                 ilib_dir,lib_dir, ch_par, num_energies,         &
                                                 reaction_cs, SE_cs, write_error)
#if(USE_MPI == 1)
               if(nproc > 1)then
                  call MPI_Barrier(icomm, ierr)
                  call MPI_BCAST(write_error, 1, MPI_LOGICAL, 0, icomm, ierr)
               end if
#endif
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Reaction cross section
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               if(write_error)call exit_YAHFC(51)
               write_error = .false.
               if(iproc == 0)call print_reaction_cs(itarget, istate, ilab, file_lab,               &
                                                   ilib_dir,lib_dir, ch_par, num_energies,         &
                                                   reaction_cs, write_error)
#if(USE_MPI == 1)
               if(nproc > 1)then
                  call MPI_Barrier(icomm, ierr)
                  call MPI_BCAST(write_error, 1, MPI_LOGICAL, 0, icomm, ierr)
               end if
#endif
               if(write_error)call exit_YAHFC(51)
            end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Absorption cross section
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            write_error = .false.
            if(iproc == 0)call print_absorption_cs(itarget, istate, ilab, file_lab,                &
                                                   ilib_dir,lib_dir, ch_par, num_energies,         &
                                                   absorption_cs, write_error)
#if(USE_MPI == 1)
            if(nproc > 1)then
               call MPI_Barrier(icomm, ierr)
               call MPI_BCAST(write_error, 1, MPI_LOGICAL, 0, icomm, ierr)
            end if
#endif
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Preequilibrium cross section       ----------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            write_error = .false.
            if(PREEQ_Model > 0 .and. iproc == 0)                                                &
                               call print_preeq_cs(itarget, istate, ilab, file_lab,             &
                                                   ilib_dir, lib_dir, ch_par,                   &
                                                   num_energies, reaction_cs, preeq_css,        &
                                                   write_error)
#if(USE_MPI == 1)
            if(nproc > 1)then
               call MPI_Barrier(icomm, ierr)
               call MPI_BCAST(write_error, 1, MPI_LOGICAL, 0, icomm, ierr)
            end if
#endif
            if(write_error)call exit_YAHFC(51)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Direct cross sections           ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            write_error = .false.
            if(iproc == 0)call print_direct_cs(itarget, istate, ilab, file_lab,                 &
                                               ilib_dir, lib_dir, ch_par,                       &
                                               num_energies, direct_cc, direct_dwba,            &
                                               direct_tot, write_error)
#if(USE_MPI == 1)
            if(nproc > 1)then
               call MPI_Barrier(icomm, ierr)
               call MPI_BCAST(write_error, 1, MPI_LOGICAL, 0, icomm, ierr)
            end if
#endif
            if(write_error)call exit_YAHFC(51)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------  Elastic channels               ---------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            write_error = .false.
            if(iproc == 0)then
               if(iproj > 0)then
                   call print_elastic(itarget, istate, ilab, file_lab,                   &
                                      ilib_dir, lib_dir, ch_par,                         &
                                      num_energies, Ang_L_max, max_jx_100, delta_Jx_100, &
                                      SE_cs, SE_Ang, Elastic_cs, Elastic_Ang,            &
                                      nucleus(itarget)%num_discrete, Inelastic_cs,       &
                                      Inelastic_Ang_L, Inelastic_L_max, write_error)
                end if
                call print_compound_elastic(itarget, istate, ilab, file_lab,                   &
                                            ilib_dir, lib_dir, ch_par,                         &
                                            num_energies, Ang_L_max, max_jx_100, delta_Jx_100, &
                                            cs_threshold, nucleus(itarget)%num_discrete,       &
                                            Inelastic_cs, Inelastic_Ang_L, Inelastic_L_max,    &
                                            write_error)
            end if
#if(USE_MPI == 1)
            if(nproc > 1)then
               call MPI_Barrier(icomm, ierr)
               call MPI_BCAST(write_error, 1, MPI_LOGICAL, 0, icomm, ierr)
            end if
#endif
            if(write_error)call exit_YAHFC(51)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------  Inelastic channels                 -----------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            write_error = .false.
            if(iproc == 0)call print_inelastic(itarget, istate, ilab, file_lab,                 &
                                               ilib_dir, lib_dir, ch_par,                       &
                                               num_energies, Ang_L_max, max_jx_50, delta_jx_50, &
                                               cs_threshold, nucleus(itarget)%num_discrete,     &
                                               absorption_cs, Inelastic_cs, Inelastic_Ang_L,    &
                                               Inelastic_L_max, Inelastic_cont_cs,              &
                                               reaction_cs, write_error)
#if(USE_MPI == 1)
            if(nproc > 1)then
               call MPI_Barrier(icomm, ierr)
               call MPI_BCAST(write_error, 1, MPI_LOGICAL, 0, icomm, ierr)
            end if
#endif
            if(write_error)call exit_YAHFC(51)
         end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------  Population calculation statistics      -------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(pop_calc)then
            write_error = .false.
            if(iproc == 0)call print_pop_mult_data(itarget, istate, ilab, file_lab,                &
                                                   ilib_dir, lib_dir, ch_par,                      &
                                                   num_energies, decay_mult, decay_var,            &
                                                   write_error)
#if(USE_MPI == 1)
            if(nproc > 1)then
               call MPI_Barrier(icomm, ierr)
               call MPI_BCAST(write_error, 1, MPI_LOGICAL, 0, icomm, ierr)
            end if
#endif
            if(write_error)call exit_YAHFC(51)
         end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Fission cross section           ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         write_error = .false.
         if(iproc == 0 .and. fission)call print_fission_cs(itarget, istate, ilab, file_lab,        &
                                                           ilib_dir, lib_dir, ch_par,              &
                                                           num_energies, reaction_cs, fission_cs,  &
                                                           num_comp,Fiss_J_avg, Fiss_J_var,        &
                                                           Fiss_tally, write_error)

#if(USE_MPI == 1)
            if(nproc > 1)then
               call MPI_Barrier(icomm, ierr)
               call MPI_BCAST(write_error, 1, MPI_LOGICAL, 0, icomm, ierr)
            end if
#endif
         if(write_error)call exit_YAHFC(51)
      end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Discrete gammas                 ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     if(iproc == 0 .and. track_gammas)                                                             &
           call print_channel_gammas_all_ein(itarget, istate, ilab, file_lab, ilib_dir, lib_dir,   &
                                             ch_par, num_energies, reaction_cs)


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------    Finished Library output    -----------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------------------------------------------------------------------+
     if(iproc == 0)write(6,*)'That'//quote//'s all Folks'
  
#if(USE_MPI == 1)
     call MPI_Finalize(ierr)
#endif

end program YAHFC_MASTER
!
!*****************************************************************************80
!
subroutine nucleus_label(icomp,length,label)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This routine puts nucleus symbol and A into character string
!
!   Dependencies:
!
!     Modules:
!
!        nuclei
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
!*****************************************************************************80
!
   use nuclei
   implicit none
!--------------------------------------------------------------
   integer(kind=4) icomp,length
   character(len=5) label
   if(nucleus(icomp)%A > 99)then
      write(label(1:3),'(i3)')nucleus(icomp)%A
      length=3
   elseif(nucleus(icomp)%A > 9)then
      write(label(1:2),'(i2)')nucleus(icomp)%A
      length=2
   else
      write(label(1:1),'(i1)')nucleus(icomp)%A
      length=1
   end if
   if(nucleus(icomp)%atomic_symbol(1:1) == ' ')then
      label(length+1:length+1)=nucleus(icomp)%atomic_symbol(2:2)
      length=length+1
   else
      label(length+1:length+2)=nucleus(icomp)%atomic_symbol(1:2)
      length=length+2
   end if
end subroutine nucleus_label
!
!*****************************************************************************80
!
real(kind=8) function jhat(xj)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This function computes 2*J+1
!
!   Dependencies:
!
!     Modules:
!
!        None
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
!*****************************************************************************80
!
   implicit none
   real(kind=8) xj
!--------------------------------------------------------------
   jhat=2.0d0*xj+1.0d0
   return
end function jhat
!
!*****************************************************************************80
!
real(kind=8) function interpolate_exp(itype, x_in, num, grid, vec)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This function interpolates between two points on an exponential grid
!    here the x-grid was constructed by x(i) = x(1)*factor**(i-1)
!    routine finds a quadratic "fit" along the grid, namely y = a + b*x + c*x**2
!    for three points along the grid.
!    The function will do an interpolation that is either linear in y (itype == 0)
!    or log in y (itype == 1).
!
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: det3_3
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
!*****************************************************************************80
!
   use nodeinfo
   implicit none
   integer(kind=4) :: itype          !   if itype == 0, linear in y, itype == 1, log in y
   real(kind=8) :: x_in
   integer(kind=4) :: num
   real(kind=8) :: grid(num), vec(num)
!------------------------------------------------------
   real(kind=8) :: mat(3,3), mat0(3,3), vecy(3), a(0:2)
   integer(kind=4) :: i,j,k
   integer(kind=4) im1, i0, ip1
   real(kind=8) :: lxm1, lx0, lxp1, lx
   real(kind=8) :: ym1, y0, yp1, y
   real(kind=8) :: denom,numer
   real(kind=8) :: diff1, diff2
!---------   External Functions        ------------------------
   real(kind=8) :: det3_3
!--------------------------------------------------------------

  i0 = 0
  if(x_in < grid(1))then
     if(iproc == 0)write(6,*) 'E_in smaller than input grid'
     call exit_YAHFC(101)
  end if
  if(x_in > grid(num))then
     if(iproc == 0)then
        write(6,*)num
        write(6,*)x_in,grid(num)
        write(6,*)'E_in larger than input grid'
     end if
     call exit_YAHFC(101)
  end if
  do i = 1, num - 1                                  !  start from the bottom of the grid
     if(abs(x_in - grid(i)) <= 1.0d-6)then                  !  it is exactly on a grid point
        interpolate_exp = vec(i)
        return
     elseif(x_in > grid(i) .and. x_in < grid(i+1))then          ! grid(i) <= e <= grid(i+1)
!-----    Sandwiched between two points. Now which one is it closest too?
!-----    use the closest point on the grid for the center of the quadratic fit
        diff1 = abs(log(grid(i)) - log(x_in))
        diff2 = abs(log(grid(i+1)) - log(x_in))
        if(diff1 < diff2)then
           i0 = i
        else
           i0 = i + 1
        end if
        exit
     end if
  end do
  if(i == 1)then                   ! at the start, so we need to use i0 = 2
     im1 = 1
     i0 = 2
     ip1 = 3
  elseif(i == num)then             !  at the end, so we use i0 = num - 1
     ip1 = num
     i0 = num - 1
     im1 = num - 2
  else                             !  in the middle
     im1 = i0 - 1
     ip1 = i0 + 1
  end if
  if(ip1 >= num)then
     ip1 = num
     i0 = ip1 - 1
     im1 = i0 - 1
  end if
  lxm1 = log(grid(im1))
  lx0 = log(grid(i0))
  lxp1 = log(grid(ip1))
  ym1 = vec(im1)
  y0 = vec(i0)
  yp1 = vec(ip1)
  if(itype == 1)then
     ym1 = log(vec(im1))
     y0 = log(vec(i0))
     yp1 = log(vec(ip1))
  end if
  mat0(1,1) = 1.0d0
  mat0(2,1) = 1.0d0
  mat0(3,1) = 1.0d0
  mat0(1,2) = lxm1
  mat0(2,2) = lx0
  mat0(3,2) = lxp1
  mat0(1,3) = lxm1**2
  mat0(2,3) = lx0**2
  mat0(3,3) = lxp1**2
  vecy(1) = ym1
  vecy(2) = y0
  vecy(3) = yp1

  denom = det3_3(mat0)

  lx = log(x_in)
  y = 0.0d0
  do k = 1, 3
     do i = 1, 3
        do j = 1, 3
           mat(i,j) = mat0(i,j)
        end do
     end do
     do i = 1, 3
        mat(i,k) = vecy(i)
     end do
     numer = det3_3(mat)
     a(k-1) = numer/denom
     y = y + a(k-1)*lx**(k-1)
  end do
  interpolate_exp = 0.0d0
  if(itype == 0)then
     interpolate_exp = y
  elseif(itype == 1)then
     interpolate_exp = exp(y)
  end if
end function interpolate_exp
!
!*****************************************************************************80
!
real(kind=8) function det3_3(mat)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This functions computes the detemrinant of a 3x3 matrix
!
!   Dependencies:
!
!     Modules:
!
!        None
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
!*****************************************************************************80
!
   implicit none
   real(kind=8) :: mat(3,3)
   det3_3 = mat(1,1)*(mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) -     &
            mat(1,2)*(mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) +     &
            mat(1,3)*(mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))
   return
end function det3_3
!
!*****************************************************************************80
!
real(kind=8) function Gauss_var(iseed_32)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This function returns Gaussian distributed variables
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
!        real(kind=8) :: random_32
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
!*****************************************************************************80
!
    use variable_kinds
    use constants
    implicit none
    integer(kind=int_32), intent(inout) :: iseed_32
    real(kind=8) :: u, v
!--------   External functions
    real(kind=8) :: random_32
!------------------------------------------------------------------------------
    u = random_32(iseed_32)
    v = random_32(iseed_32)
    if(random_32(iseed_32) < 0.5) then
       Gauss_var = sqrt(-2.0*log(u))*cos(two_pi*v)
    else
       Gauss_var = sqrt(-2.0*log(u))*sin(two_pi*v)
    end if
    return
end function Gauss_var
!
!*****************************************************************************80
!
integer(kind=4) function find_channel(n_dat, dim_part, num_part, part_data, num_part_type)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This function returns the channel number given a set of decay particles
!
!   Dependencies:
!
!     Modules:
!
!        options
!        Channel_info
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
!*****************************************************************************80
!
   use options
   use Channel_info
   implicit none
   integer(kind=4), intent(in) :: n_dat, dim_part, num_part
   real(kind=8), intent(in) :: part_data(n_dat, dim_part)
   integer(kind=4), intent(in) :: num_part_type(0:6)
!----------------------------------------------------------------------------------
   integer(kind=8) :: channel_code
   integer(kind=4) :: iparticle
   integer(kind=4) :: i, n, k
   integer(kind=4) :: ichannel

! -----    Set up code to determine channel - based on emitted particles, order does matter

   if(explicit_channels)then
      channel_code = 0
      n = 0
      do i = 1, num_part
         iparticle = nint(part_data(2,i))
         if(iparticle > 0 .and. iparticle <= 6)then
            n = n + 1
            channel_code = ior(channel_code,ishft(int(iparticle,kind=8),(n-1)*3))
         end if
      end do
   else
      channel_code = 0
      do k = 1, 6
         channel_code = ior(channel_code,ishft(int(num_part_type(k),kind=8),(k-1)*5))
      end do
   end if
!
! ------   Find it in the list of channels
!
  ichannel = 0
   do i = 1, num_channels
      if(channel_code == Exit_channel(i)%Channel_code)then
         ichannel = i
         exit
      end if
   end do

   find_channel = ichannel

  return
end function find_channel
!
!*****************************************************************************80
!
integer(kind=4) function find_ibin(Ex, inuc)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This function returns the index for the bin that brackets energy Ex in
!    nucleus inuc
!
!   Dependencies:
!
!     Modules:
!
!        options
!        nuclei
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
!*****************************************************************************80
!
   use options
   use nuclei
   implicit none
   real(kind=8), intent(in) :: Ex
   integer(kind=4), intent(in) :: inuc
!-------------------------------------------------------------------------------
   integer(kind=4) :: i
!-------------------------------------------------------------------------------
   find_ibin = 0
   do i = 1, nucleus(inuc)%nbin
      if(ex >= nucleus(inuc)%e_grid(i) - 0.5d0*nucleus(inuc)%delta_e(i) .and.   &
         ex < nucleus(inuc)%e_grid(i) + 0.5d0*nucleus(inuc)%delta_e(i))then
         find_ibin = i
         exit
      end if
   end do
   return
end function find_ibin
!
!*****************************************************************************80
!
real(kind=8) function delta_e_value(Ex, ex_base, de)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This function returns the value of the bin width delta_e based on the
!    current excitation energy. Increasing the value with increasing
!    excitation energy
!
!   Dependencies:
!
!     Modules:
!
!        None
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
!*****************************************************************************80
!
   implicit none
   real(kind=8) :: Ex, ex_base, de
!---------------------------------------------------
   delta_e_value = de
   if(Ex > ex_base .and. Ex <= 2.0d0*ex_base)then
      delta_e_value = 2.0d0*de
   elseif(Ex > 2.0d0*ex_base .and. Ex <= 4.0d0*ex_base)then
      delta_e_value = 4.0d0*de
   elseif(Ex > 4.0d0*ex_base)then
      delta_e_value = 8.0d0*de
   end if
   return
end function delta_e_value
!
!*****************************************************************************80
!
subroutine reaction_max_J
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine calculates the maximum angular momentum 
!    in the calculation based on the cross section
!
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!        options
!        useful_data
!        Channel_info
!        particles_def
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: compound_cs
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
!*****************************************************************************80
!
   use nodeinfo
   use options
   use useful_data
   use Channel_info
   use particles_def
   use nuclei
   implicit none
   integer (kind=4) :: j, n
   integer(kind=4) :: ipar
   integer(kind=4) :: max_J
   integer(kind=4) :: iproj
   integer(kind=4) :: itarget, istate
   real(kind=8) :: e_in
   real(kind=8) :: xJ_max
   real(kind=8) :: xji, xI_f
   real(kind=8) :: test_sigma
   real(kind=8) :: sum
!------------------------------------------------------
   real(kind=8) :: compound_cs
!------------------------------------------------------
   iproj = projectile%particle_type
   itarget = target%icomp
   istate = target%istate
   if(projectile%particle_type > 0)then
      e_in = projectile%energy(projectile%num_e)

      xJ_max = real(particle(iproj)%lmax,kind=8) + particle(iproj)%spin + nucleus(1)%state(istate)%spin
      max_J = nint(xJ_max - nucleus(1)%jshift)

      test_sigma = 0.0d0
      do j = 0, max_J                                          !   loop over J values
         do ipar = 0, 1                                              !   parity of compound nucleus
            xji = real(j,kind=8) + nucleus(1)%jshift
            test_sigma = test_sigma + compound_cs(e_in,ipar,xji,itarget,istate,iproj)
         end do
      end do

      sum = 0.0d0
      do J = 0, max_J                                          !   loop over J values
         do ipar = 0, 1                                              !   parity of compound nucleus
            xji = real(J,kind=8) + nucleus(1)%jshift
            sum = sum + compound_cs(e_in,ipar,xji,itarget,istate,iproj)
         end do
         if(abs(sum - test_sigma) < 1.d-7)exit
      end do
      max_J = J
      do n = 1, OpticalCS%numcc
         xI_f = OpticalCS%state(n)%spin
         J = nint(xI_f - nucleus(itarget)%jshift)
         if(J > max_J)max_J = J
      end do
      max_J_allowed = min(max_J,20)
   else
      max_J_allowed = 20
   end if
   return
end subroutine reaction_max_J
!
!*****************************************************************************80
!
subroutine set_min_particle_energy
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine sets up the minimum energy for the incident energy
!    grid to set up the optical model calculations
!
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!        options
!        particles_def
!        nuclei
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
!*****************************************************************************80
!
   use nodeinfo
   use options
   use particles_def
   use nuclei
   implicit none
   integer(kind=4) :: k, in, iproj
   integer(kind=4) :: num_energies
   real(kind=8) :: min_e
!-----------------------------------------------------------------
   iproj = projectile%particle_type
   num_energies = projectile%num_e
   if(.not.pop_calc)then
      min_e = 1000.
      do in = 1, num_energies
         if(projectile%energy(in) < min_e) min_e = projectile%energy(in)
      end do
      do k = 1, 6
         if(particle(k)%in_decay)then
            if(k == 1)particle(k)%min_e = 1.0d-5
            if(k == 2)particle(k)%min_e = 1.0d-3
            if(k == 3)particle(k)%min_e = 1.0d-3
            if(k == 4)particle(k)%min_e = 1.0d-3
            if(k == 5)particle(k)%min_e = 1.0d-3
            if(k == 6)particle(k)%min_e = 1.0d-3
         end if
      end do
      if(iproj == 1 .and. min_e < particle(iproj)%min_e)then
         if(iproc == 0)then
            write(6,'(''*********************************************'')')
            write(6,'(''*  ERROR!!!  ERROR!!!  ERROR!!!  ERROR!!!   *'')')
            write(6,'(''*  Minimum energy projectile is below a     *'')')
            write(6,'(''*  safe value. Restart with                 *'')')
            write(6,'(''*  e_min >= '',1pe15.7,''                 *'')')particle(iproj)%min_e
            write(6,'(''*********************************************'')')
         end if
         call exit_YAHFC(101)
      end if
   else
      do k = 1, 6
         if(particle(k)%in_decay)then
            if(k == 1)particle(k)%min_e = 1.0d-4
            if(k == 2)particle(k)%min_e = 1.0d-2
            if(k == 3)particle(k)%min_e = 1.0d-2
            if(k == 4)particle(k)%min_e = 1.0d-2
            if(k == 5)particle(k)%min_e = 1.0d-2
            if(k == 6)particle(k)%min_e = 1.0d-2
         end if
      end do
   end if
   return
end subroutine set_min_particle_energy
!
!*****************************************************************************80
!
subroutine modeled_level_density(inuc, yrast, yrast_actual)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine sets up continuous energy level density using internal models
!
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!        options
!        nuclei
!
!     Subroutines:
!
!        rho_J_par_e
!
!     External functions:
!
!        integer(kind=4) :: find_ibin
!        real(kind=8) :: parity_fac
!        real(kind=8) :: spin_fac
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
!*****************************************************************************80
!
   use nodeinfo
   use options
   use nuclei
   implicit none
   integer(kind=4) :: inuc
   real(kind=8) :: yrast(0:100,0:1), yrast_actual(0:100,0:1)
!-----------------------------------------------------------------------
   integer(kind=4) :: num_states, nbin
   integer(kind=4) :: j, k, ii, ipar
   integer(kind=4) :: kstart
   real(kind=8) :: pmode, bb, pe1
   integer(kind=4) :: j_max
   real(kind=8) :: xj, xnum, dde, e_lev_min, energy
   real(kind=8) :: rho
   real(kind=8) :: apu, sig, K_vib, K_rot
!---------    External Functions   -------------------------------------------
   integer(kind=4) :: find_ibin
!-----------------------------------------------------------------------------
   yrast(0:100,0:1) = -1.0d0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Find yrast energies for each spin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   num_states = nucleus(inuc)%ncut
   nbin = nucleus(inuc)%nbin
   if(all_discrete_states)num_states = nucleus(inuc)%num_discrete
   do ii = 1, num_states
       j = nint(nucleus(inuc)%state(ii)%spin - nucleus(inuc)%jshift)
       ipar = nint((nucleus(inuc)%state(ii)%parity + 1.0d0)/2.0d0)
       if(yrast(j,ipar) < 0.0d0)then
          yrast(j,ipar) = nucleus(inuc)%state(ii)%energy
          yrast_actual(j,ipar) = nucleus(inuc)%state(ii)%energy
       end if
   end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---- Calculate yrast energy for each spin with modeled level density
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(nbin > 0)then
      j_max = nucleus(inuc)%j_max
      pmode = nucleus(inuc)%level_param(16)
      pe1 = nucleus(inuc)%level_param(17)
      bb = nucleus(inuc)%level_param(18)
      do j = 0, j_max
         xj = real(j,kind=8) + nucleus(inuc)%jshift
         do ipar = 0, 1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Start level density arrays when cummlative level density has reached xnum >= rho_cut
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            xnum = 0.0d0
            dde = 0.0025d0
            e_lev_min = 0.0d0
            do k = 1, 100000
               energy = real(k,kind=8)*dde
!               call rhoe(energy,nucleus(inuc)%level_param,                          &
!                                nucleus(inuc)%vib_enh,                              &
!                                nucleus(inuc)%rot_enh,                              &
!                                nucleus(inuc)%A,rho,apu,sig,K_vib,K_rot)
!               xnum = xnum + rho*spin_fac(xj,sig)*                                  &
!                             parity_fac(energy,xj,ipar,pmode,pe1,bb)*dde
               call rho_J_par_e(energy,xj,ipar,nucleus(inuc)%level_param,           &
                                nucleus(inuc)%vib_enh,                              &
                                nucleus(inuc)%rot_enh,                              &
                                nucleus(inuc)%A,rho,apu,sig,K_vib,K_rot)
               xnum = xnum + rho*dde
               if(xnum >= rho_cut)then
                  e_lev_min = energy
                  exit
               end if
            end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Reset. If yrast is determined by discrete level, use this.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(yrast(j,ipar) > -1.0d0)then
               e_lev_min = yrast(j,ipar)
            else
               if(j <= j_max .and. yrast(j,ipar) < 0.0)then
                  yrast(j,ipar) = e_lev_min
                  yrast_actual(j,ipar) = max(e_lev_min,nucleus(inuc)%e_grid(1))
               end if
            end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----     First bin to start filling the level density at this spin and parity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            kstart = max(find_ibin(e_lev_min,inuc),1)
            pmode = nucleus(inuc)%level_param(16)
            pe1 = nucleus(inuc)%level_param(17)
            bb = nucleus(inuc)%level_param(18)
            do k = 1, nbin
               energy = nucleus(inuc)%e_grid(k)
               if(energy >= e_lev_min)then         !  check if above computed yrast line
                  call rho_J_par_e(energy, xj, ipar, nucleus(inuc)%level_param,         &
                                   nucleus(inuc)%vib_enh,                               &
                                   nucleus(inuc)%rot_enh,                               &
                                   nucleus(inuc)%A, rho, apu, sig, K_vib, K_rot)
                  nucleus(inuc)%bins(j,ipar,k)%rho = max(rho,0.0d0)
               end if
            end do
         end do
      end do
   end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   With discrete states above E_cut, reduce continuous level density
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(all_discrete_states)then
      do ii = 1, nucleus(inuc)%num_discrete
         j = nint(nucleus(inuc)%state(ii)%spin - nucleus(inuc)%jshift)
         ipar = nint((nucleus(inuc)%state(ii)%parity + 1.0d0)/2.0d0)
         if(j <= Max_J_allowed .and. ii > nucleus(inuc)%ncut)then
            k = find_ibin(nucleus(inuc)%state(ii)%energy,inuc)
            if(k <= nbin .and. k > 0)nucleus(inuc)%bins(j,ipar,k)%rho =                       &
                                     max(nucleus(inuc)%bins(j,ipar,k)%rho -                   &
                                     1.0d0/nucleus(inuc)%delta_e(k), 0.0d0)
         end if
      end do
   end if
   return
end subroutine modeled_level_density
!
!*****************************************************************************80
!
subroutine read_level_density(inuc, yrast, yrast_actual)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine sets up continuous energy level density using internal models
!
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!        options
!        nuclei
!
!     Subroutines:
!
!        none
!
!     External functions:
!
!        integer(kind=4) :: count_words
!        integer(kind=4) :: find_ibin
!        real(kind=8) :: parity_fac
!        real(kind=8) :: spin_fac
!        real(kind=8) :: sig2_param
!        real(kind=8) :: interpolate_exp
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
!*****************************************************************************80
!
   use nodeinfo
   use options
   use nuclei
   implicit none
   integer(kind=4) :: inuc
   real(kind=8) :: yrast(0:100,0:1), yrast_actual(0:100,0:1)
!-----------------------------------------------------------------------
   integer(kind=4) :: num_states, nbin
   integer(kind=4) :: i, ii, j, iJ, jj, k, l, ll
   integer(kind=4) :: ip
   integer(kind=4) :: ipar, ipar_max, n
   real(kind=8) :: xI, xJ
   real(kind=8) :: pmode, bb, pe1
   integer(kind=4) :: num_read, j_max
   real(kind=8) :: energy
   real(kind=8) :: rho
   real(kind=8) :: sig
   real(kind=8), allocatable, dimension (:,:) :: rho_temp
   real(kind=8), allocatable, dimension (:,:) :: rho_sep
   real(kind=8), allocatable, dimension (:) :: e_grid
   character(len=2000) :: line
   logical :: exist_den_file
   integer(kind=4) :: io_error
   integer(kind=4) :: ilast
!----------------------------------------------------------------------------------------+
!-------------    External functions                                                     +
!----------------------------------------------------------------------------------------+
   integer(kind=4) :: count_words
   integer(kind=4) :: find_ibin
   real(kind=8) :: sig2_param
   real(kind=8) :: parity_fac
   real(kind=8) :: spin_fac
   real(kind=8) :: interpolate_exp
!----------------------------------------------------------------------------------------+
!-----     Start exectution of subroutine                                                +
!----------------------------------------------------------------------------------------+
   if(.not. allocated(rho_sep))allocate(rho_sep(0:nucleus(inuc)%j_max,0:1))
   j_max = -1
   nbin = nucleus(inuc)%nbin
   ipar_max = 1
   pmode = nucleus(inuc)%level_param(16)
   pe1 = nucleus(inuc)%level_param(17)
   bb = nucleus(inuc)%level_param(18)
   if(nucleus(inuc)%lev_den_both)ipar_max = 0
!----------------------------------------------------------------------------------------+
!-----     Loop over parities                                                            +
!----------------------------------------------------------------------------------------+
   do ipar = 0, ipar_max
      exist_den_file = .false.
      ilast = index(nucleus(inuc)%lev_den_file(ipar),' ') - 1
      inquire(file = nucleus(inuc)%lev_den_file(ipar)(1:ilast), exist = exist_den_file)
      if(.not. exist_den_file)then
         if(iproc == 0)write(6,*)'Error in subroutine read_level_density'
         if(iproc == 0)write(6,*)'Level density file does not exist. Looking for ',  &
             nucleus(inuc)%lev_den_file(ipar)(1:ilast)
         call exit_YAHFC(101)
      end if
      open(unit=50, file = nucleus(inuc)%lev_den_file(ipar)(1:ilast), status= 'unknown')
!----------------------------------------------------------------------------------------+
!-----    Find out how many data elements to read num_read, and number of J j_max        +
!----------------------------------------------------------------------------------------+
      num_read = 0
      io_error = 0
      do while(io_error == 0)
         line(1:2000) = ' '
         read(50,'(a)',iostat = io_error)line
         if(line(1:1) == '#' .or. line(1:1) == '!' .or. io_error /= 0)cycle
         num_read = num_read + 1
         if(j_max < 0)j_max = count_words(line) - 2
      end do
      if(j_max > 0 .and. j_max < nucleus(inuc)%j_max)then
         if(iproc == 0)write(6,*)'Error in subroutine read_level_density'
         if(iproc == 0)write(6,*)'Not enough j-values. Maximum J >= ',     &
            real(nucleus(inuc)%j_max,kind=8) + nucleus(inuc)%jshift
         if(iproc == 0)write(6,*)'Number of entries musy be >= ',nucleus(inuc)%j_max + 1
         call exit_YAHFC(101)
      end if
      if(.not. allocated(e_grid))allocate(e_grid(num_read))
      if(.not. allocated(rho_temp))allocate(rho_temp(num_read,0:j_max))
!----------------------------------------------------------------------------------------+
!-----    Read level density data in                                                     +
!----------------------------------------------------------------------------------------+
      rewind(50)
      io_error = 0
      n = 0
      do while(io_error == 0)
         line(1:2000) = ' '
         read(50,'(a)', iostat = io_error)line
            if(line(1:1) == '#' .or. line(1:1) == '!'.or. io_error /= 0)cycle
            n = n + 1
            read(line,*)e_grid(n),(rho_temp(n,j),j=0,j_max)
!----------------------------------------------------------------------------------------+
!----   If rho_temp < 1.0d-6, set to 1.0d-6, otherwise, interpolation will get crazy     +
!----------------------------------------------------------------------------------------+
         do j = 0, j_max
            if(rho_temp(n,j) < 1.0d-9)rho_temp(n,j) = 1.0d-9
         end do
      end do
!----------------------------------------------------------------------------------------+
!----   Some error checks                                                                +
!----------------------------------------------------------------------------------------+
      if(e_grid(1) > nucleus(inuc)%e_grid(1))then
         if(iproc == 0)write(6,*)'Error in subroutine read_level_density'
         if(iproc == 0)write(6,*)'Min energy in read is greater than in energy in level-density grid'
         if(iproc == 0)write(6,*)'Minimum energy must be <= ',nucleus(inuc)%e_grid(1)
         call exit_YAHFC(101)
      end if
      if(e_grid(num_read) < nucleus(inuc)%e_grid(nbin))then
         if(iproc == 0)write(6,*)'Error in subroutine read_level_density'
         if(iproc == 0)write(6,*)'Max energy in read is less than max energy in level-density grid'
         if(iproc == 0)write(6,*)'Maximum energy must be >= ',nucleus(inuc)%e_grid(nbin)
         call exit_YAHFC(101)
      end if
!----------------------------------------------------------------------------------------+
!----   Fill level density array nucleus(inuc)%bins(j,ipar,n)%rho                        +
!----   Start with normal, that is all j values were read in.                            +
!----   Follow with case when j_max == 0, and use spin cutoff parameter to calculate     +
!----   spin distribution                                                                +
!----------------------------------------------------------------------------------------+
      if(j_max > 0)then
         do n = 1, nucleus(inuc)%nbin
            do j = 0, nucleus(inuc)%j_max
               xj = real(j,kind=8) + nucleus(inuc)%jshift
               rho = interpolate_exp(1, nucleus(inuc)%e_grid(n), num_read, e_grid, rho_temp(1,j))
               if(ipar_max == 1)then
                  nucleus(inuc)%bins(j,ipar,n)%rho = rho
               else              !----   Use internal model for parity distribution
                  nucleus(inuc)%bins(j,ipar,n)%rho = rho*parity_fac(energy,xj,ipar,pmode,pe1,bb)
                  nucleus(inuc)%bins(j,ipar+1,n)%rho = rho*parity_fac(energy,xj,ipar+1,pmode,pe1,bb)
               end if
            end do
!----------------------------------------------------------------------------------------+
!----   Store level density at the separation energy                                     +
!----------------------------------------------------------------------------------------+
            energy = nucleus(inuc)%sep_e(1)
            do j = 0, nucleus(inuc)%j_max
               xj = real(j,kind=8) + nucleus(inuc)%jshift
               rho = interpolate_exp(1, energy, num_read, e_grid, rho_temp(1,j))
               if(ipar_max == 1)then
                  rho_sep(j,ipar) = rho
               else              !----   Use internal model for parity distribution
                  rho_sep(j,ipar) = rho*parity_fac(energy, xj, ipar, pmode, pe1, bb)
                  rho_sep(j,ipar+1) = rho*parity_fac(energy, xj, ipar+1, pmode, pe1, bb)
               end if
            end do
         end do
!----------------------------------------------------------------------------------------+
!----   Use internal model for spin distribution                                         +
!----------------------------------------------------------------------------------------+
      elseif(j_max == 0)then
         do n = 1, nucleus(inuc)%nbin
            energy = nucleus(inuc)%e_grid(n)
            rho = interpolate_exp(1, nucleus(inuc)%e_grid(n), num_read, e_grid, rho_temp(1,0))

            sig = sig2_param(energy, nucleus(inuc)%level_param,nucleus(inuc)%A)

            do j = 0, nucleus(inuc)%j_max
               xj = real(j,kind=8) + nucleus(inuc)%jshift
               nucleus(inuc)%bins(j,ipar,n)%rho = rho*spin_fac(xj,sig)
               if(ipar_max == 1)then
                  nucleus(inuc)%bins(j,ipar,n)%rho = rho*spin_fac(xj,sig)
               else              !----   Use internal model for parity distribution
                  nucleus(inuc)%bins(j,ipar,n)%rho = rho*spin_fac(xj,sig)*                 &
                      parity_fac(energy,xj,ipar,pmode,pe1,bb)
                  nucleus(inuc)%bins(j,ipar+1,n)%rho = rho*spin_fac(xj,sig)*               &
                      parity_fac(energy,xj,ipar+1,pmode,pe1,bb)
               end if
            end do         
         end do
!----------------------------------------------------------------------------------------+
!----   Store level density at the separation energy                                     +
!----------------------------------------------------------------------------------------+
         energy = nucleus(inuc)%sep_e(1)
         sig = sig2_param(energy, nucleus(inuc)%level_param,nucleus(inuc)%A)
         do j = 0, nucleus(inuc)%j_max
            xj = real(j,kind=8) + nucleus(inuc)%jshift
            rho = interpolate_exp(1, energy, num_read, e_grid, rho_temp(1,j))
            if(ipar_max == 1)then
               rho_sep(j,ipar) = rho
            else              !----   Use internal model for parity distribution
               rho_sep(j,ipar) = rho*parity_fac(energy,xj,ipar,pmode,pe1,bb)
               rho_sep(j,ipar+1) = rho*parity_fac(energy,xj,ipar+1,pmode,pe1,bb)
            end if
         end do
      end if
      if(allocated(rho_temp))deallocate(e_grid)
      if(allocated(rho_temp))deallocate(rho_temp)
   end do
   yrast(0:100,0:1) = -1.0
!----------------------------------------------------------------------------------------+
!------   Find yrast energies for each spin                                              +
!----------------------------------------------------------------------------------------+
   num_states = nucleus(inuc)%ncut
   nbin = nucleus(inuc)%nbin
   if(all_discrete_states)num_states = nucleus(inuc)%num_discrete
   do ii = 1, num_states
      j = nint(nucleus(inuc)%state(ii)%spin - nucleus(inuc)%jshift)
      ipar = nint((nucleus(inuc)%state(ii)%parity + 1.0d0)/2.0d0)
      if(yrast(j,ipar) < 0.0)then
         yrast(j,ipar) = nucleus(inuc)%state(ii)%energy
         yrast_actual(j,ipar) = nucleus(inuc)%state(ii)%energy
      end if
   end do
!----------------------------------------------------------------------------------------+
!------   Impose yrast energies on the level density based on discete states, even if    +
!------   the level density was read in                                                  +
!----------------------------------------------------------------------------------------+
   do j = 0, nucleus(inuc)%j_max
      do ipar = 0, 1
!----------------------------------------------------------------------------------------+
!----  Yrast line for a discrete state was found, so make sure level density below this  +
!----  energy is minimal                                                                 +
!----------------------------------------------------------------------------------------+
         if(yrast(j,ipar) > -1.0d0)then 
            do n = 1, nbin
               if(nucleus(inuc)%e_grid(n) < yrast(j,ipar))then
                  nucleus(inuc)%bins(j,ipar,n)%rho = 1.0d-9
               else
                  exit
               end if
            end do
         end if
      end do
   end do
!----------------------------------------------------------------------------------------+
!------   Set yrast line for j,ipar not set with discrete states.                        +
!------   Might be nucleus(inuc)%e_grid(1)                                               +
!----------------------------------------------------------------------------------------+
   do j = 0, nucleus(inuc)%j_max
      do ipar = 0, 1
         if(yrast(j,ipar) < 0.0d0)then
            yrast(j,ipar) = nucleus(inuc)%e_grid(1)
            do n = 1, nbin
               if(nucleus(inuc)%bins(j,ipar,n)%rho <= 1.0d-6)then
                  yrast(j,ipar) = nucleus(inuc)%e_grid(n)
                  yrast_actual(j,ipar) = nucleus(inuc)%e_grid(n)
               else
                  exit
               end if
            end do
         end if
      end do
   end do
!----------------------------------------------------------------------------------------+
!------   With discrete states above E_cut, reduce continuous level density              +
!----------------------------------------------------------------------------------------+
   if(all_discrete_states)then
      do ii = 1, nucleus(inuc)%num_discrete
         j = nint(nucleus(inuc)%state(ii)%spin - nucleus(inuc)%jshift)
         ipar = nint((nucleus(inuc)%state(ii)%parity + 1.0d0)/2.0d0)
         if(j <= Max_J_allowed .and. ii > nucleus(inuc)%ncut)then
            k = find_ibin(nucleus(inuc)%state(ii)%energy,inuc)
            if(k <= nbin .and. k > 0)nucleus(inuc)%bins(j,ipar,k)%rho =                       &
                                     max(nucleus(inuc)%bins(j,ipar,k)%rho -                   &
                                     1.0d0/nucleus(inuc)%delta_e(k), 0.0d0)
         end if
      end do
   end if
!----------------------------------------------------------------------------------------+
!------   Calculate D0 and D1                                                            +
!----------------------------------------------------------------------------------------+
   xI = nucleus(inuc)%target_spin
   ipar = nucleus(inuc)%target_ipar
   iJ = nint(2.0d0*xI)
!----   D0
   l = 0
   ip = ((2*ipar-1)*(-1)**l+1)/2
   ll = 2*l
   rho = 0.0d0
   do i = iabs(iJ-ll), (iJ + ll), 2
      do j = -1, 1, 2
         k = i + j
         if(k < 0)cycle
         xJ = real(k,kind=8)/2.0d0
         jj = nint(xj - nucleus(inuc)%jshift)
         rho = rho + rho_sep(jj,ip)
      end do
   end do
   nucleus(inuc)%D0 = 1.0d6/rho
!----   D1
   l = 1
   ip = ((2*ipar-1)*(-1)**l+1)/2
   ll = 2*l
   rho = 0.0d0
   do i = iabs(iJ-ll), (iJ + ll), 2
      do j = -1, 1, 2
         k = i + j
         if(k < 0)cycle
         xJ = real(k,kind=8)/2.0d0
         jj = nint(xj - nucleus(inuc)%jshift)
         rho = rho + rho_sep(jj,ip)
      end do
   end do
   nucleus(inuc)%D1 = 1.0d6/rho

   if(allocated(rho_sep))deallocate(rho_sep)

   return
end subroutine read_level_density
!
!*******************************************************************************
!
integer function count_words(input)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function counts the number distinct word elements there are in the
!    character string input (len = 2000). 
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
   implicit none
   character(len=2000), intent(in) :: input
!----------------------------------------------------------------------
   integer(kind=4) :: numw
   integer(kind=4) :: k
   logical :: word
!----------------------------------------------------------------------
   numw = 0
   word = .false.
   k = 1
 1 if(k > 1 .and. (input(k:k) == '!' .or. input(k:k) == '#'))then
      goto 4
   elseif(input(k:k) /= ' ')then
      if(.not.word)then
         numw = numw + 1
         word = .true.
         if(k == 2000)goto 4
      end if
   elseif(input(k:k) == ' ')then
      if(word)word = .false.
      if(k == 2000) goto 4
   end if
   k = k + 1
   if(k > 2000)goto 4
   goto 1
 4 continue
   count_words = numw
   return
end function count_words
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine positions an open file to the end-of-file so that
!    subsequent writing will be appended
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        none
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
!*****************************************************************************80
!
subroutine get_to_eof(iunit)
   implicit none
!-------------------------------------------------------------------------
   integer(kind=4), intent(in) :: iunit
!-------------------------------------------------------------------------
   integer(kind=4) :: read_error
   read_error = 0
   do while(read_error == 0)
      read(iunit,*,iostat = read_error)
   end do
   backspace(iunit)
   return
end subroutine get_to_eof
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine terminates execution with either call exit or
!    MPI_Abort if MPI is being used. Exit routine is controlled via the
!    pre-processor variable USE_MPI.
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        none
!
!     External functions:
!
!        None
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
!*****************************************************************************80
!
subroutine exit_YAHFC(error_code)
   use nodeinfo
   implicit none
!-------------------------------------------------------------------------
   integer(kind=4), intent(in) :: error_code
!-------------------------------------------------------------------------
#if(USE_MPI == 1)
   call MPI_Abort(icomm,error_code,ierr)
#else
   call exit(error_code)
#endif
   return
end subroutine exit_YAHFC




