program YAHFC_MASTER
!
!*****************************************************************************
!
!  Discussion:
!
!    This program is the driver for the Hauser-Feshbach code system YAHFC
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
      use nuclei
      use Channel_info
      use particles_def
      use directory_structure
      use pre_equilibrium_no_1
      use Gauss_integration
!---------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------
!      integer(kind=int_64) :: one_int_64
      character(len=1) quote
      character(len=80) bigblank        ! 80 blank character
      character(len=132) blank132
      logical call_tco
      logical e_in_problem
!-----------------------------------------------------------------------
      character(len=80) outfile         !   moved to module useful_data
      character(len=80) directory
      integer(kind=4) ifile, idir, core_file
      integer(kind=4) icount
      character(len=40) file_name
      character(len=132) out_buff
      integer(kind=4) :: icharr, il1
!-----------------------------------------------------------------------
      integer(kind=4) :: num_e
      real(kind=8) :: ee_max

      logical converged

!-----------------------------------------------------------------------
!----------   Temporary storage for state data if we need to remove discrete states

      real(kind=8) :: test_sigma, xji
      real(kind=8) :: E_in, e_rel, e_x, e_in2, E_out
      real(kind=8) :: E_f, Ex_f, e_shift
      real(kind=8) :: rel_factor
      real(kind=8) :: xmu
      integer(kind=4) :: iproj, itarget
      integer(kind=4) :: jproj

      real(kind=8) :: spin_target, spin_proj
      real(kind=8) :: sp1, sp2
      real(kind=8) :: cs_threshold
      real(kind=8) :: cs
!      real(kind=8) :: D0, D1
      integer(kind=4) :: isp
      integer(kind=4) :: isp_max
      integer(kind=4) :: istate, jstate
      real(kind=8) :: sum
      real(kind=8) :: ex_tot

      integer(kind=4) :: max_tries

      real(kind=8) :: temp_cs
      integer(kind=4) :: num_eff, num_elastic
      real(kind=8) :: xnum_eff, xnum_elastic
      real(kind=8) :: smear(0:6)
      real(kind=8) :: e1, e2
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
      integer(kind=4) :: ipi, ipf, nf
      character(len=1) ch_par(0:1)


      integer(kind=4) :: num_part, n_dat, dim_part

      parameter (n_dat = 24)

      real(kind=8), allocatable :: part_data(:,:)
      real(kind=8), allocatable :: extra_angle_data(:,:)
      integer(kind=4) :: num_theta, nang, nang_tot
      integer(kind=4) :: Ang_L_max
      integer(kind=4) :: L_max
      real(kind=8), allocatable :: part_Ang_data(:,:)

      integer(kind=4) ::  L_Ang
      real(kind=8) :: xL_Ang, xl_i
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Data required for Haser-Feshbach Monte Carlo probabilities 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real(kind=8) :: Gamma_g, Gamma_g_exp, dGamma_g_exp, Gamma_g_old, diff
      real(kind=8) :: step
      integer(kind=4) :: g_error
      integer(kind=4) :: num_res, num_res_old

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----  Probabilities for incident channels
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real(kind=8) :: mass_target, mass_exit, Q
      real(kind=8) :: mass_proj
      integer(kind=4) :: num_channel
      real(kind=8), pointer, dimension(:) :: channel_prob
      integer(kind=4), pointer, dimension(:,:) :: ichannel
      integer(kind=4) :: ie_nuc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      Label for channels
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer(kind=4) :: ichann
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!      Data to implement max_particle constraints during decays
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer(kind=4) :: num_part_type(6)
      real(kind=8) :: part_fact(0:7)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Random number data 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real(kind=8) :: ran
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer(kind=4) :: ifind, ilast
      integer(kind=4) :: num_s
      integer(kind=4) :: ictype
      real(kind=8) :: xI_f
      real(kind=8) ::  par_i, xnbin_i
      integer(kind=4) :: iX_i, Ix_f
      real(kind=8) :: preeq_prob
      real(kind=8) :: xnum, dde, e_lev_min
      integer(kind=4) :: kstart

      real(kind=8) :: pop_sum, prev
      real(kind=8), allocatable, dimension (:) :: pop_prob
      integer(kind=4), allocatable, dimension (:) :: pop_j, pop_ip
      logical :: pop_input
      integer(kind=4) ijj

      integer(kind=4) :: num_energies

      real(kind=8) :: tolerance

      real(kind=8) :: sum_e
      real(kind=8) :: e_diff, max_e_diff, avg_e_diff
      integer(kind=4) :: num_no_fiss

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      real(kind=8) :: KE
      real(kind=8) :: energy
      real(kind=8) :: Ef
      real(kind=8) :: Ex
      integer(kind=4) :: i,j,k,l,m,is
      integer(kind=4) :: Z_i,A_i
      integer(kind=4) :: Z_f,A_f
      integer(kind=4) :: in,ip,ia
      integer(kind=4) :: in2,in3
      integer(kind=4) :: icomp
      integer(kind=4) :: cum_fit
      real(kind=8) :: min_e
      character(len=1) overide
      real(kind=8) :: hang
      logical :: hung_up
      real(kind=8) :: sum_ie, sum_ce
      real(kind=8), allocatable, dimension(:) :: check_cs
      integer(kind=4) :: inuc
      integer(kind=4) :: max_num
!---------------------------------------------------------------------
      integer(kind=4) :: ich
      character(len=2) char_pos,char_neg
      real(kind=8) :: t12_isomer
      real(kind=8) :: f,f1
      integer(kind=4) :: l_radiation
      real(kind=8) :: apu
      real(kind=8) :: rho, sig
      real(kind=8) :: xa,xz
!      real(kind=8) :: ratio
      integer(kind=4) :: max_part
      real(kind=8) :: K_vib, K_rot
!----------------------------------------------------------------------
      integer(kind=4) :: num_comp
      integer(kind=4) :: max_gammas
!----------------------------------------------------------------------
      integer(kind=4) :: ibin
      integer(kind=4) :: nbin
      integer(kind=4) :: n, n_f, if1, i_f
      real(kind=8) :: xj_min, xj_max
      real(kind=8) :: xI_min, xI_max
      integer(kind=4) :: max_J
      integer(kind=4) :: j_max
      real(kind=8) :: xj, xI
      integer(kind=4) :: Ix, Ix_min, Ix_max
      integer(kind=4) :: ipar
      real(kind=8) :: xnorm
      real(kind=8) :: xnorm1, xnorm2, xnorm3, xnorm4
!--------------------  Data to collect input information -------------
      integer(kind=4) :: read_err
      integer(kind=4) :: nchar
      integer(kind=4) :: icommand
      character(len=132) :: command
      character(len=132) char_temp
      integer(kind=4) :: itemp
      integer(kind=4) :: num_command
      logical finish
!----------------------------------------------------------------------
!---------------------   Level information ----------------------------
      real(kind=8) :: ex_max, S_part
      real(kind=8) :: dem, dep
      real(kind=8), allocatable :: comp_avg(:,:,:),comp_var(:,:,:)
      real(kind=8), allocatable :: test_e(:)
      integer(kind=4) :: ishift
!----------------------------------------------------------------
      real(kind=8), allocatable :: reaction_cs(:)
      real(kind=8), allocatable :: absorption_cs(:)
      real(kind=8), allocatable :: preeq_css(:,:)
      real(kind=8), allocatable :: preeq_spect(:,:)
      real(kind=8), allocatable :: preeq_spect_full(:,:)
      real(kind=8), allocatable :: fission_cs(:)
      character(len=5) target_label
      integer(kind=4) :: ntar
      integer(kind=4) :: ilab, ilab2
      character(len=100) :: file_lab
      character(len=20) file_lab2


      real(kind=8) :: population
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
      real(kind=8), allocatable :: Ang_Dist(:)
      real(kind=8), allocatable, dimension(:,:,:) :: Temp_Ang_dist
      real(kind=8) :: Temp
      real(kind=8) :: avg, avg_diff, expected_diff
      integer(kind=4) :: jx
!--------   Used to compute Rutherford Crossection

      real(kind=8) :: Sig_C, Coul, P_L, shape, comp

      real(kind=8) :: egrid, delta_e


      character(len=5) resid_label
      integer(kind=4) :: nres

!------------------------------------------------------------------
!--------    Gauss-Legendre integration useful to construct "smooth"
!--------    angular distributions
      integer(kind=4) :: LL_max
      real(kind=8), allocatable :: xvalue(:)

!--------   Total Inelastic Scattering
      logical direct, compound, preeq_decay,cc_decay, dwba_decay
      logical inelastic,compound_elastic, fission_decay
!-------    Sum of all processes in these different channels
      real(kind=8), allocatable :: Inelastic_Total(:)           !  Total inelastic cross section as a function of energy
      real(kind=8), allocatable :: Inelastic_Continuous(:)       !  Total inelastic from compound going to contiuum bins
      real(kind=8), allocatable :: Inelastic_C_Discrete(:)      !  Total inelastic from compound going to discrete states 
      real(kind=8), allocatable :: Inelastic_P_Discrete(:)      !  Total inelastic from pre-equilibirum going to discrete states
      real(kind=8), allocatable :: Inelastic_Direct(:)          !  Total direct excitation
      real(kind=8), allocatable :: Inelastic_Preeq(:)           !  Total direct excitation

!-------    Inelastic state by state starts here 
!-------    First sum of compound and compound and direct
      real(kind=8), allocatable :: Inelastic_cs(:,:)
      integer(kind=4), allocatable :: Inelastic_count(:,:)
      integer(kind=4), allocatable :: Inelastic_L_max(:,:)
      real(kind=8), allocatable :: Inelastic_Ang_L_avg(:,:,:)
      real(kind=8), allocatable :: Inelastic_Ang_L(:,:,:)
      real(kind=8), allocatable :: Inelastic_Ang_dist(:,:,:)
      real(kind=8), allocatable :: Inelastic_g_mult(:)
      real(kind=8), allocatable :: Inelastic_g_mult_var(:)
      real(kind=8) :: event_g_mult
!-------   Direct Inelastic Scattering
      real(kind=8) :: tot_direct, sum_d
      real(kind=8), allocatable :: direct_cs(:)
      real(kind=8), allocatable :: direct_Ang(:,:)
      real(kind=8), allocatable :: direct_prob(:,:)
      real(kind=8), allocatable :: direct_tot(:)
      real(kind=8), allocatable :: direct_cc(:)
      real(kind=8), allocatable :: direct_dwba(:)
      real(kind=8), allocatable :: SE_prob(:)
      real(kind=8), allocatable :: Elastic_cs(:)
      real(kind=8), allocatable :: Elastic_Ang(:,:)
!-------    Shape Elastic Scattering
      real(kind=8), allocatable :: SE_cs(:)
      real(kind=8), allocatable :: SE_Ang(:,:)
!-------    Compound Elastic Scattering
      real(kind=8), allocatable :: CE_cs(:)
      real(kind=8), allocatable :: CE_Ang(:,:)
!
!--------     Quasi-elastic Scattering
      real(kind=8), allocatable :: QE_cs(:)
      real(kind=8), allocatable :: QE_Spectrum(:)
      real(kind=8), allocatable :: QE_Ang(:,:)
!
!--------     Primary Gamma-Spectrum
      real(kind=8), allocatable :: Primary_Gamma_Spectrum(:)
!
!-------    Cross section and spectrum from x-particle emission
!
      real(kind=8), allocatable :: x_particle_cs(:,:)
      real(kind=8), allocatable :: x_particle_Spectrum(:,:)
      real(kind=8), allocatable :: x_particle_Spectrum_dist(:,:)
      real(kind=8), allocatable :: x_particle_Spectrum_res(:,:)
      real(kind=8), allocatable :: x_particle_Spectrum_res_dist(:,:)
      real(kind=8) :: de_spec2

      real(kind=8) :: preeq_test_cs
      real(kind=8), allocatable :: preeq_Spectrum(:)
      real(kind=8), allocatable :: direct_Spectrum(:)
      real(kind=8), allocatable :: dwba_Spectrum(:)
      real(kind=8), allocatable :: res1(:), res2(:), resp(:,:), resp2(:,:)
 !     real(kind=8) :: max_n_e
!---------   Statistics on Fissioning nuclei
      real(kind=8), allocatable :: Fiss_J_avg(:,:), Fiss_J_var(:,:)
      real(kind=8), allocatable :: Fiss_tally(:,:)

!
!--------     Avg particles when using pop_calc = .true.
      real(kind=8) :: avg_part(0:6), avg_part_e(0:6)
      real(kind=8) :: var_part(0:6), var_part_e(0:6)
!
!--------     Spectrum of emitted particles in Lab frame 
!--------     when pop_calc = .true.
!      integer(kind=4) :: num_spect_e, ind_e, ile
!      real(kind=8), allocatable :: part_Lab_spectrum(:,:)
      real(kind=8) :: angle_dist(0:90)

      integer(kind=4) :: num_compound, num_pre_equilibrium, num_fission
      real(kind=8) :: test_count

!--------   yrast  ---------------------------------
      real(kind=8) :: yrast(0:100,0:1)
      real(kind=8) :: yrast_actual(0:40,0:1)
      character(len=9) yrast_file

      real(kind=8) :: x, x1, theta

      real(kind=8) :: pmode, pe1, bb

      integer(kind=4) :: num_events, num_inelastic

      real(kind=8) :: Boost_Lab(0:3,0:3), Boost_COM(0:3,0:3)
      real(kind=8) :: mass_1, mass_2
      real(kind=8) :: pp_res, p_res(0:3), v_res(0:3)
      real(kind=8) :: beta, gamma

      real(kind=8) :: theta_0, phi_0, T_1, phi, T_2, T_L, theta_L, phi_L
      integer(kind=4) :: num_bad_e

      integer(kind=4) :: n_min

      real(kind=8) :: check_sum, sum_inelastic

      real(kind=8) :: tally_weight
      real(kind=8) :: tally_norm

      integer(kind=4),allocatable :: command_rank(:)
      character(len=132), dimension(:), allocatable :: command_buffer

      integer :: cnt, cnt_r, cnt_m

      real(kind=8) :: ratio
      real(kind=8) :: thresh
      logical :: print_cs
!-----------------------------------------------------
!---------    Directory structure of Library outputs
      logical f_exist
      character(len=132) unix_cmd
      integer(kind=4) :: icmd
      character(len=132) file_check
      integer(kind=4) :: if_check
      character(len=132) :: lib_dir
      integer(kind=4) :: ilib_dir
!--------   Data needed to output data libraries
      character(len=14) iso_label
      real(kind=8), allocatable :: Ang_dis(:,:)
      integer(kind=4) :: icc_max

!-------------------- Memory variables
    real(kind=8) :: mem,mem_icomp,mem_tot, mem_bins

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
      real(kind=8) :: spin_fac
      real(kind=8) :: parity_fac
      real(kind=8) :: random_64
      real(kind=8) :: poly
      real(kind=8) :: clebr
      real(kind=8) :: interpolate
      real(kind=8) :: Gauss_var
      real(kind=8) :: EL_f
      real(kind=8) :: EL_trans
      real(kind=8) :: ML_f
      real(kind=8) :: ML_trans
      integer(kind=4) :: rank_commands
      real(kind=8) :: compound_cs

!-----------------   Start Main body of calculation     ----------------

   ilast = index(version,' ')-1
   write(6,*)'**********************************************'
   write(6,*)'* Livermore Monte Carlo Hauser-Feshbach code'
   write(6,*)'* Version # ',version(1:ilast)
   write(6,*)'**********************************************'
   write(6,*)

      quote = "'"
!------   dummy set up with MPI
      angle_dist(0:90)=0.0d0

      iproc=0
      nproc=1
!----------  Default Options
      print_output=.true.
      print_spectra=.true.
      output_mode = 1
      call_tco=.true.
      t12_isomer=1.0
      write_me = 0
      ex_pop_file(1:20)=' '
      ex_set = .false.
      Max_J_allowed = 20
      Fiss_Max_J = 0.0d0
      pop_calc = .false.
      track_gammas = .false.
      Out_gammas_vs_E = .false.
      track_primary_gammas = .false.
      rho_cut = 0.25
      fit_Gamma_gamma = .false.
      Apply_Coulomb_Barrier = .true.
      All_gammas = .false.
      T_norm = 1.0d0
      pair_model = 1               !   Default for pairing model
      preeq_delta = 0.0d0
      Preeq_V = 38.0d0
      Preeq_g_div = 15.d0
      Preeq_g_a = .false.
      biased_sampling = .false.
      optical = 'fresco'
      explicit_channels = .false.
      num_theta_angles = 10
      target%istate = 1
      dump_events = .false.
      binary_event_file = .true.
      event_generator = .false.
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
      hbar_c = 197.32705d0                  !   MeV*fm
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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------    ixx_max and jxx_max need to be even
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ixx_max = 50
      delta_ix = 2.0d0/real(ixx_max,kind=8)

      max_jx_10 = 10
      delta_jx_10 = 2.0d0/real(max_jx_10,kind=8)
      max_jx_20 = 20
      delta_jx_20 = 2.0d0/real(max_jx_20,kind=8)
      max_jx_50 = 50
      delta_jx_50 = 2.0d0/real(max_jx_50,kind=8)
      max_jx_100 = 100
      delta_jx_100 = 2.0d0/real(max_jx_100,kind=8)

      pp_max = 6
      pn_max = 6

      part_lmax = 20

      d_zero = 0.0d0
      d_half = 0.5d0
      d_one = 1.0d0
      d_two = 2.0d0
      d_three = 3.0d0
      d_four = 4.0d0
      factorial(0) = 1.0d0
      do i = 1, 50
         factorial(i) = factorial(i-1)*dfloat(i)
      end do
      lfactorial(0) = 0.0d0
      do i = 1, 101
         lfactorial(i) = lfactorial(i-1) + dlog(dfloat(i))
      end do

      ifresco_shape = 13

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------    Set up Gauss-Laguerre quadrature values
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      n_glag = 50
      allocate(x_glag(n_glag),w_glag(n_glag))
      alf = 0.0d0
      bet = 0.0d0
!      call cdgqf(n_glag, 5, alf, 0.0d0, x_glag, w_glag)
      call gauss_quad(n_glag, 2, alf, bet, x_glag, w_glag)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------    Set up Gauss-Legendre quadrature values
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      Ang_L_max = 30
      n_gleg = Ang_L_max + 1
      allocate(x_gleg(n_gleg), w_gleg(n_gleg))
      alf = 0.0d0
      bet = 0.0d0
      call gauss_quad(n_gleg, 1, alf, bet, x_gleg, w_gleg)


      trans_p_cut = 1.0d-6
      trans_e_cut = 1.0d-15
      prob_cut = 1.0d-6

      num_mc_samp = 1000000
!
!---   initialize iseed to int(pi*10^9)
!
      iseed = 3141592654_int_64
!---
!---  "Randomize" iseed with cnt from system clock. Generally, this 
!---  will be different for each runs.
!---
      call system_clock(COUNT = cnt, COUNT_RATE = cnt_r, COUNT_MAX = cnt_m)
      cnt = max(mod(cnt,100000),1)
      iseed = iseed + cnt
      if(iand(iseed,1_int_64) /= 1_int_64)iseed = iseed + 1_int_64

      iseed = -iseed


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
         write(6,*)'Environment variable for data is not set correctly'
         write(6,*)'Program exiting'
         stop
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
      m_l_max = 2
      out_file = 'YAHFC-Calc'
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
      M2_Rnn = 1.0d0
      M2_Rpn = 1.5d0
      M2_Rnp = 1.5d0

      quasi_elastic =.false.
      quasi_e = 0.0d0

!-------------------------------------------------------------------  
      i_bind = 0
      compound_setup = .false.
      de = 0.1d0
      mass_file = 'aw'
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!-------    set up buffer for commands to parsed later
!-------    this offers the ability to enter comands in
!-------    any order whatsoever - maybe not a great idea
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    Get input commands.                                           +
!------    Store in temporary file to count, store in buffer, order      +
!------    their priority and then run                                   +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      num_command = 0
      finish = .false.
      open(unit=17,file='YAHFC-commands.txt',status='unknown')
      do while(.not. finish)
         read(5,'(a)',iostat = read_err)command
         if(read_err /= 0)then
            finish = .true.
            cycle
         end if
         if(command(1:1) == '#' .or. command(1:1) == '!')cycle
         call parse_string(command,numw,startw,stopw)
         if(numw < 1)cycle
         nchar = stopw(numw)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Convert string to lower case
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         nchar = stopw(1) - startw(1) + 1
         call lower_case_word(nchar,command(startw(1):stopw(1)))
         if(command(startw(1):stopw(1)) == 'end')then
             finish = .true.
             cycle
         end if
         num_command = num_command + 1
         write(17,'(a132)')command
      end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Set up buffer for input commands
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(.not. allocated(command_buffer))allocate(command_buffer(num_command))
      if(.not. allocated(command_rank))allocate(command_rank(num_command))
      command_buffer(1:num_command) = blank132
      command_rank(1:num_command) = 10
      rewind(17)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Read commands back in from temporary file
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
!-------   order the commands by rank
!-------   Lowest number is highest rank! This is give extra room later if
!-------   ranks need to be adjusted.
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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Write back out to unit=17 in the order of execution
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      open(unit=17,file='YAHFC-commands.txt',status='unknown')
      write(17,'(''YAHFC commands in order of execution'')')
      do i = 1, num_command
         write(17,'(a132)')command_buffer(i)
      end do
      close(unit=17)

      num_comp=0
      icommand=1
      finish=.false.

      iproj = -10
      do i = 1, num_command
         call parse_command(num_comp,icommand,command_buffer(i),finish)
         if(target%specified)istate = target%istate
         if(.not. compound_setup .and. (target%specified .and.        &
               projectile%specified .and. ex_set))then
            iproj = projectile%particle_type

            Z_i = target%Z + projectile%Z
            A_i = target%A + projectile%A
            call set_up_decay_chain(projectile%Z,projectile%A,  &
                                    target%Z,target%A,num_comp)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   now that decay chain is set up, also set up data for each
!------   compound nucleus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            overide = 'c'                             !   use discrete states marked with 'u' in RIPL-2 file
            do icomp = 1, num_comp
               Z_f = nucleus(icomp)%Z
               A_f = nucleus(icomp)%A
               if(Z_f == target%Z .and. A_f == target%A)target%icomp = icomp
               itarget = target%icomp
               call get_spectrum(data_path,len_path,               &
                                 overide,t12_isomer,               &
                                 symb(Z_f),Z_f,A_f,icomp)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Find E1 strength function parameters
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               nucleus(icomp)%E1_model = E1_model
               nucleus(icomp)%E1_default = .true.
               nucleus(icomp)%er_E1(1:4) = 0.0d0
               nucleus(icomp)%gr_E1(1:4) = 0.0d0
               nucleus(icomp)%sr_E1(1:4) = 0.0d0
               nucleus(icomp)%Gamma_g = 0.0d0
               nucleus(icomp)%Gamma_g_exp = 0.0d0
               nucleus(icomp)%dGamma_g_exp = 0.0d0

               call gdr_param(data_path,len_path,                  &
                              nucleus(icomp)%Z,nucleus(icomp)%A,   &
                              nucleus(icomp)%er_E1,                &
                              nucleus(icomp)%gr_E1,                &
                              nucleus(icomp)%sr_E1)     

               if(A_f > 20)then
                  lev_option = 1
                  if(A_f > 130) lev_option = 2
               else
                  lev_option = 0
               end if
               call get_lev_den(data_path,len_path,                &
                                symb(Z_f),Z_f,A_f,icomp)
               if(nucleus(icomp)%fit_ematch)call fit_lev_den(icomp)
               nucleus(icomp)%fission=.false.
               if(fission .and. nucleus(icomp)%Z >= 80)then
                  call Fission_data(data_path,len_path,icomp)
               end if
            end do
            spin_target=nucleus(itarget)%state(target%istate)%spin
            spin_proj = 0.0
            if(.not.pop_calc)spin_proj=particle(iproj)%spin
            compound_setup=.true.
         end if
      end do

      num_theta = num_theta_angles

!------  After everything is set up, also check if preequilibrium model parameters
!------  makes sense. In particular, the finite well parameter, can't have
!------  Preq_V1 > Preeq_V
      if(preeq_V <= Preeq_V1)then
         write(6,*)'------------   WARNING  ------------------'
         write(6,*)'Problem with preequilibrium model - Preeq_V1 > Preeq_V'
         write(6,*)'Resetting to make it more physical'
         Preeq_V1 = max(Preeq_V - 5.0d0,0.0d0)
         write(6,*)'Preeq_V1 = max(Preeq_V - 5.0,0.0)'
         write(6,*)'Preeq_V = ',Preeq_V
         write(6,*)'Preeq_V1 = ',Preeq_V1
      end if
!
!------   First call to random number generator
!
      ran = random_64(iseed)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    Finished getting input commands      -------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      do icomp = 1, num_comp

         ia = nucleus(icomp)%A
!-------------    Check if D0exp > 0 and requesting to fit it, if not
!-------------    make flag to fit it .false.
         if(nucleus(icomp)%D0exp <= 0.0d0 .and. nucleus(icomp)%fit_D0)   &
            nucleus(icomp)%fit_D0 = .false.

         call finish_lev_den(icomp)

         if(nucleus(icomp)%fission)then
            do i = 1, nucleus(icomp)%F_n_barr

               if(nucleus(icomp)%F_Barrier(i)%level_param(6) >                        &
                  nucleus(icomp)%F_Barrier(i)%level_param(3))then
                  call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,          &
                                 nucleus(icomp)%F_Barrier(i)%vib_enh,                 &
                                 nucleus(icomp)%F_Barrier(i)%rot_enh)
               else
                  write(6,*)'**************************************************'
                  write(6,*)'Attempt to set E_match < Delta'
                  write(6,*)'Nucleus icomp = ',icomp
                  write(6,*)'Z = ',nucleus(icomp)%Z,' A = ',nucleus(icomp)%A
                  write(6,*)'Fission Barrier # ',i
                  write(6,*)'Cannot continue with this calculation'
                  write(6,*)'**************************************************'
                  stop
               end if
            end do
            call Fission_levels(icomp)
         end if

        if(.not. all_gammas)nucleus(icomp)%num_discrete = nucleus(icomp)%ncut

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


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Check if we have sufficient data to proceed
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(.not.target%specified)then
         write(6,*)'Target nucleus not specified, quitting'
         stop
      end if
      if(.not.projectile%specified)then
         write(6,*)'Projectile nucleus not specified, quitting'
         stop
      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------    Check to see if we need to overide fission
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

      nucleus(1:num_comp)%PREEQ=.false.
      if(.not. pop_calc .and. PREEQ_Model == 1)then
         p0(2) = projectile%Z
         p0(1) = projectile%A - p0(2)
         nucleus(1)%PREEQ = .true.
      end if



      ifile=index(out_file,' ') - 1
      if(print_output)then
         open(unit=13,file=out_file(1:ifile)//'.out',status='unknown')
         open(unit=400,file=out_file(1:ifile)//'.bad_energy',status='unknown')
         open(unit=15,file=out_file(1:ifile)//'_e_avg.dat',status='unknown')
         open(unit=16,file=out_file(1:ifile)//'_j_avg.dat',status='unknown')
      end if

      cum_fit=2

      mass_proj = particle(iproj)%mass
      mass_target = nucleus(itarget)%mass + nucleus(itarget)%state(target%istate)%energy
      write(6,*)'Q-values for each channel'
      do i = 1, num_channels
         inuc = Exit_Channel(i)%Final_nucleus
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------   Calculate Q-value
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

!    write(6,*)mass_target
!    write(6,*)particle(iproj)%mass
!    write(6,*)inuc, nucleus(inuc)%mass
!    write(6,*)mass_exit

         Q = mass_target + particle(iproj)%mass - nucleus(inuc)%mass - mass_exit
         if(abs(Q) < 1.0d-6)Q = 0.0d0

         ilast = index(Exit_Channel(i)%Channel_Label,' ')
         write(6,'(''channel '',i4,1x,a20,f15.7)')i, Exit_Channel(i)%Channel_Label,Q
 
         Exit_Channel(i)%Q_value = Q

      end do


      if(pop_calc)then
         WF_model = 0
         PREEQ_Model = 0
         track_gammas = .true.
      end if
      num_energies = projectile%num_e


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!--------    Allocate arrays for each nucleus to store calculated 
!--------    cross sections
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(fission)then
         allocate(fission_cs(1:num_energies))
         fission_cs(1:num_energies)=0.0d0
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
            if(nucleus(icomp)%state(i)%isomer)nucleus(icomp)%num_cs=nucleus(icomp)%num_cs+1
         end do
         allocate(nucleus(icomp)%channel_cs(1:num_energies,1:nucleus(icomp)%num_cs))
         nucleus(icomp)%channel_cs(1:num_energies,1:nucleus(icomp)%num_cs) = 0.0d0
         allocate(nucleus(icomp)%hang_cs(1:num_energies))
         nucleus(icomp)%hang_cs(1:num_energies) = 0.0d0
         if(nucleus(icomp)%PREEQ .and. PREEQ_Model > 0)then
            allocate(nucleus(icomp)%PREEQ_cs(1:num_energies))
            nucleus(icomp)%PREEQ_cs(1:num_energies) = 0.0d0
            allocate(nucleus(icomp)%PREEQ_part_cs(1:nucleus(icomp)%num_decay,1:num_energies))
            nucleus(icomp)%PREEQ_part_cs(1:nucleus(icomp)%num_decay,1:num_energies) = 0.0d0
         end if
         if(nucleus(icomp)%fission)then
            allocate(nucleus(icomp)%fission_cs(1:num_energies))
            nucleus(icomp)%fission_cs(1:num_energies)=0.0d0
         end if
         allocate(nucleus(icomp)%state_lab(nucleus(icomp)%num_cs))
         nucleus(icomp)%state_lab(1)=1                 !  index of the state being tracked
         icount=1
         do i=2,nucleus(icomp)%num_discrete
            if(nucleus(icomp)%state(i)%isomer)then
               icount=icount+1
               nucleus(icomp)%state_lab(icount)=i
            end if
         end do
      end do

      pop_input = .false.
      if(pop_calc)then
         pop_input = .true.
         num_pop = 2*(Max_J_allowed + 1)
         if(.not.allocated(pop_prob))allocate(pop_prob(num_pop))
         if(.not.allocated(pop_j))allocate(pop_j(num_pop))
         if(.not.allocated(pop_ip))allocate(pop_ip(num_pop))
      end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------                                                                +
!---------   Count and catalog the compound nuclei in each chain          +
!---------   Which compound nuclei does nucleus(i) decay to and by        +
!---------   emitting what particle                                       +
!---------                                                                +
!-------------------------------------------------------------------------+
!-------------------------------------------------------------------------+
!---------  Set up energy grid and de  -  The goal is to have the         +
!---------  excitation energy of the first compound nucleus correspond    +
!---------  to the separation energy of the projectile                    +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Ex_max = nucleus(1)%Ex_max
      S_part = 0.0
      if(.not.pop_calc)then
         S_part = nucleus(1)%sep_e(iproj)
         if(Ex_max - S_part > de/2.0d0)then
            nbin = max(int((Ex_max-S_part)/de+0.5d0),1)
            dep = (Ex_max-S_part)/(dfloat(nbin)+0.5d0)
            dem = dep
            if(Ex_max - S_part > de/2.0d0)then
               dem = (Ex_max - S_part)/(real(nbin)-0.5d0)
            end if
            de = min(dep,dem)
            nbin = max(int((Ex_max - S_part - de/2.0d0)/de + 1.0d0),1)
         else
            nbin = 1
         end if
      else
         nbin = max(int((Ex_max - S_part)/de+0.5d0),1)
         dep = (Ex_max - S_part)/(dfloat(nbin)+0.5d0)
         dem = dep
         if(Ex_max - S_part > de/2.0d0)then
            dem = (Ex_max - S_part)/(real(nbin)-0.5d0)
         end if
         de = min(dep,dem)
         nbin = max(int((Ex_max - S_part - de/2.0d0)/de + 1.0d0),1)
      end if     
!-------------------------------------------------------------------------+
!---------                                                                +
!---------   Check energy levels and set up energy grid across            +
!---------   all nuclei involved in the decay chain                       +
!---------                                                                +
!-------------------------------------------------------------------------+
!----------------------    Allocate arrays for particle spectra
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do icomp = 1, num_comp
         Ef = (4.0d0/real(nucleus(icomp)%A,kind=8))*Init_Kinetic_Energy
         Ex = (sqrt(nucleus(icomp)%Ex_max)+sqrt(Ef))**2
         nucleus(icomp)%nbin_part = int(Ex/de) + 1
         if(nucleus(icomp)%PREEQ .and. PREEQ_Model > 0)then
            allocate(nucleus(icomp)%PREEQ_part_spectrum(1:nucleus(icomp)%num_decay,0:nucleus(icomp)%nbin_part))
            nucleus(icomp)%PREEQ_part_spectrum(1:nucleus(icomp)%num_decay,0:nucleus(icomp)%nbin_part)=0.0d0
         end if
      end do

      do i = 1, num_comp                  !   initial number of bins
         if(nucleus(i)%Ex_max > nucleus(i)%level_param(7) + de/2.)then   
            nbin = int((nucleus(i)%Ex_max -                                    &
                        nucleus(i)%level_param(7))/de+0.5)     ! energy difference level_param(7)=E_cut
            nbin = max(nbin,1)
         else
            nbin = 1                      !  Ex_max is less than the highest discrete state - no bins needed
         end if
         nucleus(i)%j_max = Max_J_allowed
         nucleus(i)%nbin = nbin
         if(nbin > 0)then
            allocate(nucleus(i)%e_grid(nbin))
            allocate(nucleus(i)%delta_e(nbin))
         end if
      end do
!
!-------------------   Set up energy grids for each nucleus
!
      do i = 1, num_comp
         nbin = nucleus(i)%nbin
         if(nbin == 0)cycle
         do j = 1, nbin
            nucleus(i)%e_grid(j) = nucleus(i)%Ex_max - (nbin-j)*de
            nucleus(i)%delta_e(j) = de
         end do
      end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------                                                                   +
!--------   Now that energy grids are set up, remap incident projectile   +
!--------   energies to fit on excitation energy grid                     +
!------                                                                   +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      num_bad_e = 0

      if(.not.allocated(test_e))allocate(test_e(num_energies))
  
      rel_factor = mass_target/(mass_target + mass_proj)
      e_in_problem = .false.
      xmu = dfloat(target%A)/dfloat(target%A + projectile%A)
      itarget = target%icomp
      do in = 1, num_energies
         e_in = projectile%energy(in)

         test_e(in) = projectile%energy(in)
         if(.not.pop_calc)then
            if(e_in < de/2.0)cycle                  !   keep it if energy is less than de/2
            e_rel = e_in*rel_factor
            e_x = e_rel + nucleus(1)%sep_e(projectile%particle_type)
!----   Find bin associated with this energy, not quite as simple as
!----   it used to be as we are considering the possiblity of unequal bins
            j = 0
            do i = 1, nucleus(1)%nbin
               if(e_x >= nucleus(1)%e_grid(i) - nucleus(1)%delta_e(i) .and.   &
                  e_x < nucleus(1)%e_grid(i) + nucleus(1)%delta_e(i))then
                  j = i
                  exit
               end if
            end do
            if(j > 0)then
               e_x = nucleus(1)%e_grid(j)
               e_rel = e_x - nucleus(1)%sep_e(projectile%particle_type)
               e_in = e_rel/rel_factor
               projectile%energy(in) = e_in
            else
               e_x = nucleus(1)%e_grid(1)
               e_rel = e_x - nucleus(1)%sep_e(projectile%particle_type)
               e_in = e_rel/rel_factor
               projectile%energy(in) = e_in
               e_in_problem = .true.
            end if
         else
            e_x = e_in
            j = nucleus(1)%nbin - nint((nucleus(1)%Ex_max - e_x)/de)
            e_x = nucleus(1)%e_grid(j)
            projectile%energy(in) = e_in
         end if
      end do
      if(e_in_problem)then
         write(6,*)'**************************************************************'
         write(6,*)'* Warning, compound excitation extends to low energy with    *'
         write(6,*)'* with low level density. Calculaitons will only be          *'
         write(6,*)'* for incident energies that populate the continuous energy  *'
         write(6,*)'* bins.                                                      *'
         write(6,*)'* Continue at your own risk!                                 *'
         write(6,*)'**************************************************************'
      end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Now look to see if we have to remove any energies from the list  +
!------   due to potential overlap                                         +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do in = 1, num_energies-1
         e_in = projectile%energy(in)
         if(e_in < de/2.0)cycle                  !   keep it if energy is less than de/2
 99      continue      
         do in2 = in + 1, num_energies
            e_in2 = projectile%energy(in2)
            if(abs(e_in - e_in2) < 1.0d-4)then          !  next energy is the same
               do in3 = in2 + 1, num_energies         !  push energies in list down
                  ishift = in3 - in2 - 1
                  projectile%energy(in2+ishift) = projectile%energy(in3)
               end do
               num_energies=num_energies-1            !  reduce number of energies in list
               goto 99
            end if
         end do
      end do
      projectile%num_e = num_energies
         
!--------------------------------------------------------------------------+
!------    Now allocate arrays to track average excitation energy and spin +
!------    for each nucleus before it decays- useful book-keeping          +
!--------------------------------------------------------------------------+
      allocate(comp_avg(2,num_comp,num_energies),comp_var(2,num_comp,num_energies))
!--------------------------------------------------------------------------+
!------    Now allocate arrays to track fission                            +
!--------------------------------------------------------------------------+

      Ef = (4.0d0/real(nucleus(num_comp)%A,kind=8))*Init_Kinetic_Energy
      Ex = (sqrt(nucleus(1)%Ex_max)+sqrt(Ef))**2

      do k=0,6
         particle(k)%nbin_spectrum = int(particle(k)%max_e/de) + 1
         if(PREEQ_Model > 0)then
            allocate(particle(k)%PREEQ_spectrum(0:particle(k)%nbin_spectrum))
            particle(k)%PREEQ_spectrum(0:particle(k)%nbin_spectrum)=0.0d0
         end if
      end do

      if(.not.pop_calc)then
         min_e = 1000.
         do in = 1, num_energies
            if(projectile%energy(in) < min_e) min_e = projectile%energy(in)
         end do
         do k = 1, 6
            if(particle(k)%in_decay)then
               if(k == 1)particle(k)%min_e = 1.0d-5
               if(k == 2)particle(k)%min_e = 1.0d-2
               if(k == 3)particle(k)%min_e = 1.0d-2
               if(k == 4)particle(k)%min_e = 1.0d-2
               if(k == 5)particle(k)%min_e = 1.0d-1
               if(k == 6)particle(k)%min_e = 1.0d-1
            end if
         end do
         if(min_e < particle(iproj)%min_e)then
            write(6,'(''*********************************************'')')
            write(6,'(''*  WARNING!WARNING!WARNING!WARNING!WARNING  *'')')
            write(6,'(''*  Minimum energy projectile is below what  *'')')
            write(6,'(''*  safe value. Restart with                 *'')')
            write(6,'(''*  e_min >= '',1pe15.7,''               *'')')particle(iproj)%min_e
            write(6,'(''*********************************************'')')
            stop
         end if
      else
         do k = 1, 6
            if(particle(k)%in_decay)then
               if(k == 1)particle(k)%min_e = 1.0d-4
               if(k == 2)particle(k)%min_e = 1.0d-2
               if(k == 3)particle(k)%min_e = 1.0d-2
               if(k == 4)particle(k)%min_e = 1.0d-2
               if(k == 5)particle(k)%min_e = 1.0d-1
               if(k == 6)particle(k)%min_e = 1.0d-1
            end if
         end do
      end if


!-------------------------------------------------------------------------+
!------                                                                   +
!------    Set up transmission coefficients for particle emission         +
!------                                                                   +
!-------------------------------------------------------------------------+

      call optical_setup(data_path, len_path, iproj, itarget, istate, de, num_comp, Ang_L_max)

!
!----   Check for maximum angular momentum from the reaction dynamics
!
      if(.not. pop_calc)then
         e_in = projectile%energy(num_energies)
  
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
         do j = 0, max_J                                          !   loop over J values
            do ipar = 0, 1                                              !   parity of compound nucleus
               xji = real(j,kind=8) + nucleus(1)%jshift
               sum = sum + compound_cs(e_in,ipar,xji,itarget,istate,iproj)
            end do
            if(abs(sum - test_sigma) < 1.d-5)exit
!            if(sum/test_sigma > 0.9999d0)exit
         end do
         max_J = J
!----
         max_J_allowed = min(max_J,40)
      end if        

!
!-----  Now check if any discrete are higher, and set to maximum spin of discrete states
!
      do i = 1, num_comp
         do n = 1, nucleus(i)%num_discrete
            J = nint(nucleus(i)%state(n)%spin - nucleus(i)%jshift)
            if(J > max_J_allowed)max_J_allowed = J
         end do
      end do


      write(6,*)'Maximum J-value set to ',max_J_allowed
      do i = 1, num_comp
         nucleus(i)%j_max = max_J_allowed
      end do


!-------------------------------------------------------------------------+
!------                                                                   +
!------   We have the transmission coefficients with maximum l-value      +
!------   so set up largest spin we need to track in the calculation      +
!------                                                                   +
!-------------------------------------------------------------------------+

      j_max = Max_J_allowed

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----  set up and initialze arrays for starting populations
      if(.not.allocated(target%pop_xjpi))allocate(target%pop_xjpi(0:j_max,0:1))

      sp1 = abs(spin_target-spin_proj)
      sp2 = spin_target+spin_proj
      isp = int(sp2-sp1)

      call start_IO(num_comp)

      call output_trans_coef


!-------------------------------------------------------------------------+
!------                                                                   +
!-------   Create level density arrays                                    +
!------                                                                   +
!-------------------------------------------------------------------------+


      do i = 1, num_comp
         nbin = nucleus(i)%nbin
         allocate(nucleus(i)%HF_den(0:j_max,0:1,1:nbin))
         allocate(nucleus(i)%bins(0:j_max,0:1,1:nbin))
         nucleus(i)%HF_den(0:j_max,0:1,1:nbin)=0.0d0
         nucleus(i)%bins(0:j_max,0:1,1:nbin)%rho = 0.0d0

         yrast_file(1:1)='Z'
         if(nucleus(i)%z < 10)then
            yrast_file(2:3)='00'
            write(yrast_file(4:4),'(i1)')nucleus(i)%z
         elseif(nucleus(i)%z < 100)then
            yrast_file(2:2)='0'
            write(yrast_file(3:4),'(i2)')nucleus(i)%z
         else
            write(yrast_file(2:4),'(i3)')nucleus(i)%z
         end if
         yrast_file(5:6)='-A'
         if(nucleus(i)%A < 10)then
            yrast_file(7:8)='00'
            write(yrast_file(9:9),'(i1)')nucleus(i)%A
         elseif(nucleus(i)%A < 100)then
            yrast_file(7:7)='0'
            write(yrast_file(8:9),'(i2)')nucleus(i)%A
         else
            write(yrast_file(7:9),'(i3)')nucleus(i)%A
         end if
         open(unit=99,file=yrast_file//'.yrast.dat',status='unknown')

         yrast(0:100,0:1) = -1.0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Find yrast energies for each spin
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do ii = 1, nucleus(i)%ncut
             j = nint(nucleus(i)%state(ii)%spin - nucleus(i)%jshift)
             ipar = nint((nucleus(i)%state(ii)%parity + 1.0d0)/2.0d0)
             if(yrast(j,ipar) < 0.0)then
                yrast(j,ipar) = nucleus(i)%state(ii)%energy
                yrast_actual(j,ipar) = nucleus(i)%state(ii)%energy
             end if
         end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   With discrete states above E_cut, reduce continuous level density
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do ii = 1, nucleus(i)%num_discrete
             j = nint(nucleus(i)%state(ii)%spin - nucleus(i)%jshift)
             ipar = nint((nucleus(i)%state(ii)%parity + 1.0d0)/2.0d0)
             if(j <= 20 .and. ii > nucleus(i)%ncut)then
                k = max(nint((nucleus(i)%state(ii)%energy - nucleus(i)%e_grid(1))/de)+1,1)
                if(k <= nbin)nucleus(i)%bins(j,ipar,k)%rho = nucleus(i)%bins(j,ipar,k)%rho - 1.0d0/de
             end if
         end do

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---- Calculate yrast energy for each spin with modeled level density
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(nbin > 0)then
            j_max = nucleus(i)%j_max
            pmode = nucleus(i)%level_param(16)
            pe1 = nucleus(i)%level_param(17)
            bb = nucleus(i)%level_param(18)
            do j = 0, j_max
               xj = dfloat(j) + nucleus(i)%jshift
               do ipar = 0, 1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Start level density arrays when cummlative level density has reached xnum >= rho_cut
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  xnum = 0.0d0
                  dde = 0.0025d0
                  e_lev_min = 0.0d0
                  do k = 1, 100000
                     energy = real(k,kind=8)*dde
                     call rhoe(energy,nucleus(i)%level_param,                   &
                                      nucleus(i)%vib_enh,                       &
                                      nucleus(i)%rot_enh,                       &
                                      nucleus(i)%A,rho,apu,sig,K_vib,K_rot)
                     xnum = xnum + rho*spin_fac(xj,sig)*                        &
                                   parity_fac(energy,xj,ipar,pmode,pe1,bb)*dde
                     if(xnum >= rho_cut)then
                        e_lev_min = energy
                        exit
                     end if
                  end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Reset. If yrast is determined by discrete level, use this.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  if(yrast(j,ipar) > 0.0)then
                      e_lev_min = yrast(j,ipar)
                  else
                     if(j <= 20 .and. yrast(j,ipar) < 0.0)then
                          yrast(j,ipar) = e_lev_min
                          yrast_actual(j,ipar) = max(e_lev_min,nucleus(i)%e_grid(1))
                     end if 
                  end if
                  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----     First bin to start filling the level density at this spin and parity
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  kstart = max(nint((e_lev_min - nucleus(i)%e_grid(1))/de),1)
                  pmode = nucleus(i)%level_param(16)
                  pe1 = nucleus(i)%level_param(17)
                  bb = nucleus(i)%level_param(18)
                  do k = 1, nbin
                     energy = nucleus(i)%e_grid(k)
                     if(k >= kstart)then         !  check if above computed yrast line
                        call rhoe(energy,nucleus(i)%level_param,                   &
                                         nucleus(i)%vib_enh,                       &
                                         nucleus(i)%rot_enh,                       &
                                         nucleus(i)%A,rho,apu,sig,K_vib,K_rot)
                        nucleus(i)%bins(j,ipar,k)%rho = max(rho*spin_fac(xj,sig)*  &
                                parity_fac(energy,xj,ipar,pmode,pe1,bb),0.0d0)
                     end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----      Still check if rho is < 0 from embedded discrete states, if so, set to zero.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                     if(nucleus(i)%bins(j,ipar,k)%rho < 0.0d0) nucleus(i)%bins(j,ipar,k)%rho = 1.0d-49
                  end do
               end do
            end do
         end if

         write(99,'(''#               Negative Parity          Positive Parity'')')
         write(99,'(''#  J           Yrast  Yrast-lev          Yrast  Yrast-lev'')')

         do j = 0, 20
            xj=dfloat(j)+nucleus(i)%jshift
            write(99,'(1x,f4.1,2(4x,2(1x,f10.4)))')xj,          &
               yrast_actual(j,0),yrast(j,0),yrast_actual(j,1),yrast(j,1)
         end do
         close(unit=99)
      end do


      do i = 1, num_comp
         nucleus(i)%lmax_E=e_l_max                     !   maximum Electro-magnetic multipole 
!----------------------------------------------------------------------
         xa=dfloat(nucleus(i)%A)
         xz=dfloat(nucleus(i)%Z)
         nucleus(i)%nbin_em = int(nucleus(i)%Ex_max/de) + 1
         allocate(nucleus(i)%str_E(0:nucleus(i)%nbin_em,nucleus(i)%lmax_E))
         allocate(nucleus(i)%f_E(0:nucleus(i)%nbin_em,nucleus(i)%lmax_E))
         nucleus(i)%str_E(0:nucleus(i)%nbin_em,1:nucleus(i)%lmax_E) = 0.0d0
         nucleus(i)%f_E(0:nucleus(i)%nbin_em,1:nucleus(i)%lmax_E) = 0.0d0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------   Electric dipole strength function
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         num_res = 0
         do k = 1,3
            if(nucleus(i)%sr_E1(k) > 1.0d-6)num_res = num_res +1
         end do
         nucleus(i)%num_res = num_res

         l_radiation = 1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Magnetic dipole
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(m_l_max >= 1)then
            nucleus(i)%lmax_M = m_l_max
            allocate(nucleus(i)%er_M(1:nucleus(i)%lmax_M))
            allocate(nucleus(i)%gr_M(1:nucleus(i)%lmax_M))
            allocate(nucleus(i)%sr_M(1:nucleus(i)%lmax_M))
            allocate(nucleus(i)%str_M(0:nucleus(i)%nbin_em,nucleus(i)%lmax_M))
            nucleus(i)%str_M(0:nucleus(i)%nbin_em,1:nucleus(i)%lmax_M) = 0.0d0
            allocate(nucleus(i)%f_M(0:nucleus(i)%nbin_em,nucleus(i)%lmax_M))
            nucleus(i)%f_M(0:nucleus(i)%nbin_em,1:nucleus(i)%lmax_M) = 0.0d0
            l_radiation = 1
            nucleus(i)%er_M(1) = 41.0d0/xa**(1.0d0/3.0d0)
            nucleus(i)%gr_M(1) = 4.0d0

!-----------   Get M1 value for f at 7 MeV

            f1 = EL_f(i, 1, 7.0d0, nucleus(i)%sep_e(1) )

!----------    Now M1 value for f at 7 MeV
 
            nucleus(i)%sr_M(1) = 1.0d0
            f = ML_f(i, l_radiation, 7.0d0)

            ratio=f1/f

            nucleus(i)%sr_M(1)=ratio/(0.0588*xa**(0.878))

         end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Higher electric multipoles
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(e_l_max > 1)then
            allocate(nucleus(i)%er_E(2:e_l_max))
            allocate(nucleus(i)%gr_E(2:e_l_max))
            allocate(nucleus(i)%sr_E(2:e_l_max))
         end if
         if(e_l_max >= 2)then           !  Include electric quadrupole
            nucleus(i)%er_E(2) = 63.0d0/xa**(1.0d0/3.0d0)
            nucleus(i)%gr_E(2) = 6.11-0.012*xa
            nucleus(i)%sr_E(2) = 0.00014d0*xz**2*nucleus(i)%er_E(2)/        &
                                 (xa**(1.0d0/3.0d0)*nucleus(i)%gr_E(2))

         end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------    Electric multipoles above E2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do l_radiation=3,e_l_max
            nucleus(i)%er_E(l_radiation) = nucleus(i)%er_E(l_radiation-1)
            nucleus(i)%gr_E(l_radiation) = nucleus(i)%gr_E(l_radiation-1)
            nucleus(i)%sr_E(l_radiation) = 8.0d-4*nucleus(i)%sr_E(l_radiation-1)
         end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------    Magnetic multipoles above M1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do l_radiation=2,m_l_max
            nucleus(i)%er_M(l_radiation) = nucleus(i)%er_M(l_radiation-1)
            nucleus(i)%gr_M(l_radiation) = nucleus(i)%gr_M(l_radiation-1)
            nucleus(i)%sr_M(l_radiation) = 8.0d-4*nucleus(i)%sr_M(l_radiation-1)
         end do
      end do

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Now that we have level densities and strength functions.       +
!-------   We examine Gamma_gamma and fit to it if necessary              +
!-------   Fit by adding an additional E1 strength function and adjust    +
!-------   the strength as necessary to reproduce Gamma_gamma for l=0.    +
!-------   The Additional E1 strength function is assumed to have         +
!-------   a centroid at 5.0 MeV with a width of 5 MeV. Note that this    +
!-------   additional mode can have a negative strength. Most likely      +
!-------   will not have a negative total strength funciton.              +
!-------   If this is the case, most likely the E1 strength is too        +
!-------   large to begin with.                                           +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      icomp = 1
      do icomp = 1, num_comp
         converged = .false.
         Gamma_g_exp = nucleus(icomp)%Gamma_g_exp
         call Gamma_gamma(icomp,0,de,Gamma_g,g_error)
         if(g_error > 0) converged = .true.                 !  Don't try to fit as Gamma_g can't be computed (Ex_max < Sn)
!----   Don't fit if Gamma_g_exp < 0.0, no data
         if(Gamma_g_exp < 0.0d0)converged = .true.

         if(nucleus(icomp)%dGamma_g_exp > 0.0d0)then
            dGamma_g_exp = nucleus(icomp)%dGamma_g_exp
         else
            dGamma_g_exp = Gamma_g_exp*0.05d0              !  just in case there is data, but error bar is messed up
         end if
         tolerance = min(dGamma_g_exp,0.1d0)             !  tolerance on fit. Make it small so that we can actually get Exp
         if(tolerance < 1.0d-3)tolerance = 0.1d0
         diff = 0.0d0

         if(Gamma_g_exp > 0.0d0)diff = abs(Gamma_g - Gamma_g_exp)                 !  current difference

         num_res_old = nucleus(icomp)%num_res

         if(diff > tolerance .and. fit_Gamma_gamma)then
            nucleus(icomp)%num_res = nucleus(icomp)%num_res + 1
            num_res = nucleus(icomp)%num_res
            nucleus(icomp)%er_E1(num_res) = 5.0d0
            nucleus(icomp)%gr_E1(num_res) = 5.0d0
            nucleus(icomp)%sr_E1(num_res) = 0.0d0
            step = nucleus(icomp)%sr_E1(1)/5000.0
            if(Gamma_g > Gamma_g_exp)step = - step       !   Calculated Gamma_g is too big, make strength negative
            nucleus(icomp)%sr_E1(num_res) = nucleus(icomp)%sr_E1(num_res) + step
            Gamma_g_old = Gamma_g
            call Gamma_gamma(icomp,0,de,Gamma_g,g_error)
            diff = abs(Gamma_g - Gamma_g_exp)
            if(diff < tolerance)converged = .true.
            if(step > 0.0 .and. Gamma_g > Gamma_g_exp) converged = .true.       !   Shouldn't really get here, but exit a trap
            if(step < 0.0 .and. Gamma_g < Gamma_g_exp) converged = .true.
            do while(.not. converged)
               nucleus(icomp)%sr_E1(num_res) = nucleus(icomp)%sr_E1(num_res) + step
               call Gamma_gamma(icomp,0,de,Gamma_g,g_error)
               diff = abs(Gamma_g - Gamma_g_exp)
               if(diff < tolerance)then
                  converged = .true.
               else
                  if(step < 0.0d0 .and. (Gamma_g_old > Gamma_g_exp .and. Gamma_g < Gamma_g_exp))then
                     step = -step/10.0
                  elseif(step > 0.0d0 .and. (Gamma_g_old < Gamma_g_exp .and. Gamma_g > Gamma_g_exp))then
                     step = -step/10.0
                  end if
               end if
               Gamma_g_old = Gamma_g
            end do
         end if

         nucleus(icomp)%Gamma_g = Gamma_g
         call Gamma_gamma(icomp,1,de,Gamma_g,g_error)
         nucleus(icomp)%Gamma_g_1 = Gamma_g

         if(nucleus(icomp)%Gamma_g > -1.0)then
            iA = nucleus(icomp)%A
            num_res = nucleus(icomp)%num_res
            write(6,*)
            ich = 1
            if(nucleus(icomp)%atomic_symbol(1:1) == ' ')then
               ich = 2
               write(6,'(5x,i3,a1)')iA,nucleus(icomp)%atomic_symbol(ich:2)
            else
               write(6,'(5x,i3,a2)')iA,nucleus(icomp)%atomic_symbol(ich:2)
            end if
	 
            if(Gamma_g > 0.0d0)then
               if(nucleus(icomp)%Gamma_g_exp > 0.0d0)then
                  write(6,'(''Gamma_gamma (l=0)'','' Calc = '',f8.3,'', Exp = '',f8.3,'' +/- '',f8.3)')    &
                     nucleus(icomp)%Gamma_g,nucleus(icomp)%Gamma_g_exp, nucleus(icomp)%dGamma_g_exp
               else
                  write(6,'(''Gamma_gamma (l=0)'','' Calc = '',f8.3,'', Exp = UNAVAILABLE'')')    &
                     nucleus(icomp)%Gamma_g
               end if
               if(nucleus(icomp)%Gamma_g_1_exp > 0.0d0)then
                  write(6,'(''Gamma_gamma (l=1)'','' Calc = '',f8.3,'', Exp = '',f8.3,'' +/- '',f8.3)')    &
                     nucleus(icomp)%Gamma_g_1,nucleus(icomp)%Gamma_g_1_exp, nucleus(icomp)%dGamma_g_1_exp
               else
                  write(6,'(''Gamma_gamma (l=1)'','' Calc = '',f8.3,'', Exp = UNAVAILABLE'')')    &
                     nucleus(icomp)%Gamma_g_1
               end if
               if(num_res > num_res_old)then
                  write(6,*)'In order to reproduce Gamma_gamma(l=0), and addition E1 resonance with paramters'
                  write(6,*)'Centroid = ',nucleus(icomp)%er_E1(num_res)
                  write(6,*)'EWidth   = ',nucleus(icomp)%gr_E1(num_res)
                  write(6,*)'Strength = ',nucleus(icomp)%sr_E1(num_res)
               end if
            elseif(Gamma_g < 0.0d0)then
               if(nucleus(icomp)%Gamma_g_exp > 0.0d0)then
                  write(6,'(''Gamma_gamma (l=0)'','' Calc = Not Calculated, Exp = '',f8.3,'' +/- '',f8.3)')    &
                     nucleus(icomp)%Gamma_g,nucleus(icomp)%Gamma_g_exp, nucleus(icomp)%dGamma_g_exp
               else
                  write(6,'(''Gamma_gamma (l=0)'','' Calc = Not Calculated, Exp = UNAVAILABLE'')')    &
                     nucleus(icomp)%Gamma_g
               end if
               if(nucleus(icomp)%Gamma_g_1_exp > 0.0d0)then
                  write(6,'(''Gamma_gamma (l=1)'','' Calc = Not Calculated, Exp = '',f8.3,'' +/- '',f8.3)')    &
                     nucleus(icomp)%Gamma_g_1,nucleus(icomp)%Gamma_g_1_exp, nucleus(icomp)%dGamma_g_1_exp
               else
                  write(6,'(''Gamma_gamma (l=1)'','' Calc = Not Calculated, Exp = UNAVAILABLE'')')    &
                     nucleus(icomp)%Gamma_g_1
               end if

            end if

            write(6,*)
         end if
      end do


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------                                                                   +
!-------   Now that we fit Gamma_gamma, set up arrays                     +
!-------   for strength functions                           +
!------                                                                   +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do i=1,num_comp
         nucleus(i)%lmax_E=e_l_max                     !   maximum Electro-magnetic multipole 
!----------------------------------------------------------------------
         xa=dfloat(nucleus(i)%A)
         xz=dfloat(nucleus(i)%Z)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------   Electric multipole strength functions
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do l_radiation = 1, e_l_max
            do j=0,nucleus(i)%nbin_em
               energy=dfloat(j)*de
               nucleus(i)%f_E(j,l_radiation) = EL_f(i,l_radiation, energy, nucleus(i)%sep_e(1))
               nucleus(i)%str_E(j,l_radiation) = EL_trans(i,l_radiation, energy, nucleus(i)%sep_e(1))
            end do
        end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Magnetic multipole strength functions
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do l_radiation = 1, m_l_max
            do j=0,nucleus(i)%nbin_em
               energy=dfloat(j)*de
               nucleus(i)%f_M(j,l_radiation) = ML_f(i, l_radiation, energy)
               nucleus(i)%str_M(j,l_radiation) = ML_trans(i,l_radiation, energy)
            end do
        end do
      end do



      if(.not.pop_calc)then
         e_rel = projectile%energy(num_energies)*rel_factor
         ee_max = e_rel + nucleus(itarget)%sep_e(iproj) +               &
                          nucleus(itarget)%state(target%istate)%energy
      else
         ee_max = 0.0
         do i = 1, num_pop_e
            if(Pop_data(i)%Ex_pop > ee_max) ee_max = Pop_data(i)%Ex_pop
         end do
      end if

      de_spec = de
      if(de <= 1.0d0)then
         ii = nint(1.0d0/de)
         de_spec = 1.0d0/real(ii,kind=8)
      else
         de_spec = real(nint(de),kind=8)    
      end if

      de_spec2 = de_spec
       
      num_e = int(ee_max/de_spec) + 3

      write(6,*)'de = ',de, 'de_spec = ',de_spec,'num_e = ', num_e

      if(.not. pop_calc)then
          if(PREEQ_Model > 0)then
             allocate(preeq_css(0:6,1:num_energies))
             preeq_css(0:6,1:num_energies) = 0.0d0
             allocate(preeq_spect(0:6,0:num_e))
             allocate(preeq_spect_full(0:6,0:num_e))
          end if
       end if



!-------------------------------------------------------------------------+
!------                                                                   +
!-------   Set up Hauser-Feshbach denominators                            +
!------                                                                   +
!-------------------------------------------------------------------------+

       do i=1,num_comp
          write(6,*)'HF-denominators for nucleus #',i
          call HF_denominator(i)
       end do

!-------------------------------------------------------------------------+
!---------                                                                +
!---------   Print out data used in the calculation for each              +
!---------   compound nucleus                                             +
!---------                                                                +
!-------------------------------------------------------------------------+

      call output_nucleus_data(num_comp, j_max, itarget)

!-------------------------------------------------------------------------+
!---------                                                                +
!---------      Cacluate memory being used for each compound nucleus      +
!---------                                                                +
!-------------------------------------------------------------------------+
      mem_tot=0.0d0
      do icomp = 1,num_comp
         mem_icomp = 0.0d0
         if(allocated(nucleus(icomp)%pop))then
             mem = dfloat(size(nucleus(icomp)%pop)*8)
             mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%PREEQ_pop))then
            mem = dfloat(size(nucleus(icomp)%PREEQ_pop)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%rho))then
            mem = dfloat(size(nucleus(icomp)%rho)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%HF_den))then
            mem = dfloat(size(nucleus(icomp)%HF_den)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%PREEQ_cs))then
            mem = dfloat(size(nucleus(icomp)%PREEQ_cs)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%PREEQ_part_cs))then
            mem = dfloat(size(nucleus(icomp)%PREEQ_part_cs)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%channel_cs))then
            mem = dfloat(size(nucleus(icomp)%channel_cs)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%er_M))then
            mem = dfloat(size(nucleus(icomp)%er_M)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%gr_M))then
            mem = dfloat(size(nucleus(icomp)%gr_M)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%sr_M))then
            mem = dfloat(size(nucleus(icomp)%sr_M)*8)
            mem_icomp = mem_icomp + mem
         end if
         mem_tot = mem_tot + mem_icomp
         if(allocated(nucleus(icomp)%er_E))then
            mem = dfloat(size(nucleus(icomp)%er_E)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%gr_E))then
            mem = dfloat(size(nucleus(icomp)%gr_E)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%sr_E))then
            mem = dfloat(size(nucleus(icomp)%sr_E)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%f_E))then
            mem = dfloat(size(nucleus(icomp)%f_E)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%f_M))then
            mem = dfloat(size(nucleus(icomp)%f_M)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%str_M))then
            mem = dfloat(size(nucleus(icomp)%str_M)*8)
            mem_icomp = mem_icomp + mem
         end if
         if(allocated(nucleus(icomp)%str_E))then
            mem = dfloat(size(nucleus(icomp)%str_E)*8)
            mem_icomp = mem_icomp + mem
         end if

!-------   Size of derived type bins

         mem_bins = 0.0d0
         if(allocated(nucleus(icomp)%bins))then
            nbin = nucleus(icomp)%nbin
            do n = 1, nbin
               do ip = 0, 1
                  do j = 0, j_max
                     do nn = 1, nucleus(icomp)%num_decay + 1      !  HF_prob element
                          mem_bins = mem_bins + 8.0d0
                     end do
                     do nn = 1, nucleus(icomp)%num_decay          !  nuke_decay element
                        mem_bins = mem_bins +                   &
                          real(nucleus(icomp)%bins(j,ip,n)%nuke_decay(nn)%num_decay,kind=8)*(8.0d0+4.0d0)
                     end do
                  end do
               end do
            end do
         end if
         mem_icomp = mem_icomp + mem_bins

         mem_tot = mem_tot + mem_icomp
         write(6,'(''Memory used for compound nucleus#'',i3,'' = '',f14.2,'' kbytes'')')icomp,mem_icomp/1.0d3 
      end do
      write(6,'(''Memory used = '',f14.3,'' kbytes'')')mem_tot/1.0d3


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
      


      if(.not.pop_calc)then
         if(.not.allocated(clb_l))allocate(clb_l(0:Ang_L_max,0:particle(iproj)%lmax))
         do L_Ang = 0, Ang_L_max
            xL_ang = real(L_Ang,kind=8)
            do l_i = 0, particle(iproj)%lmax
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
            allocate(Exit_Channel(i)%Channel_cs(1:num_energies,-2:num_s))
         Exit_Channel(i)%Channel_cs(1:num_energies,-2:num_s) = 0.0d0
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
                 allocate(Exit_Channel(i)%part_mult(0:6,-1:num_s,1:num_energies))
         Exit_Channel(i)%part_mult(0:6,-1:num_s,1:num_energies) = 0.0d0
         if(.not.allocated(Exit_Channel(i)%Ang_L))                                                 &
                 allocate(Exit_Channel(i)%Ang_L(0:Ang_L_max,0:num_e,                               &
                        0:6,-1:num_s,1:num_energies))
         Exit_Channel(i)%Ang_L(0:Ang_L_max,0:num_e,0:6,-1:num_s,1:num_energies) = 0.0d0
         if(.not.allocated(Exit_Channel(i)%Max_L))                                                 &
                 allocate(Exit_Channel(i)%Max_L(0:num_e,0:6,-1:num_s,1:num_energies))
         Exit_Channel(i)%Max_L(0:num_e,0:6,-1:num_s,1:num_energies) = 0.0d0

         if(.not.allocated(Exit_Channel(i)%Spect))                                                 &
                 allocate(Exit_Channel(i)%Spect(0:6,-1:num_s,1:num_energies))
         do in = 1, num_energies
            do j = -1, num_s
               do k = 0, 6
                  if(.not.allocated(Exit_Channel(i)%Spect(k,j,in)%E_spec))                         &
                     allocate(Exit_Channel(i)%Spect(k,j,in)%E_spec(0:num_e))
                  Exit_Channel(i)%Spect(k,j,in)%E_spec(0:num_e) = 0.0d0
!----
                  if(.not.allocated(Exit_Channel(i)%Spect(k,j,in)%E_count))                        &
                     allocate(Exit_Channel(i)%Spect(k,j,in)%E_count(0:num_e))
                  Exit_Channel(i)%Spect(k,j,in)%E_count(0:num_e) = 0
!----
                  if(.not.allocated(Exit_Channel(i)%Spect(k,j,in)%Ang_Dist))                       &
                     allocate(Exit_Channel(i)%Spect(k,j,in)%Ang_Dist(0:max_jx_10))
                  Exit_Channel(i)%Spect(k,j,in)%Ang_Dist(0:max_jx_10) = 0.0d0
!----
                  if(.not.allocated(Exit_Channel(i)%Spect(k,j,in)%Ang_L))                          &
                     allocate(Exit_Channel(i)%Spect(k,j,in)%Ang_L(0:Ang_L_max))
                  Exit_Channel(i)%Spect(k,j,in)%Ang_L(0:Ang_L_max) = 0.0d0
                  Exit_Channel(i)%Spect(k,j,in)%Ang_L_Max = 0
!----
                  if(.not.allocated(Exit_Channel(i)%Spect(k,j,in)%E_Ang_Dist))                     &
                      allocate(Exit_Channel(i)%Spect(k,j,in)%E_Ang_Dist(0:max_jx_10,0:num_e))
                  Exit_Channel(i)%Spect(k,j,in)%E_Ang_Dist(0:max_jx_10,0:num_e) = 0.0d0
!----
                  if(.not.allocated(Exit_Channel(i)%Spect(k,j,in)%E_Ang_L))                        &
                      allocate(Exit_Channel(i)%Spect(k,j,in)%E_Ang_L(0:Ang_L_max,0:num_e))
                  Exit_Channel(i)%Spect(k,j,in)%E_Ang_L(0:Ang_L_max,0:num_e) = 0.0d0
!----
                  if(.not.allocated(Exit_Channel(i)%Spect(k,j,in)%E_Ang_L_max))                    &
                      allocate(Exit_Channel(i)%Spect(k,j,in)%E_Ang_L_max(0:num_e))
                  Exit_Channel(i)%Spect(k,j,in)%E_Ang_L_Max(0:num_e) = 0
               end do
            end do
         end do


      end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Arrays tracking Inelastic and Compound Elastic Channels              + 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(.not. pop_calc)then
         if(.not.allocated(Inelastic_cs))                                                          &
             allocate(Inelastic_cs(0:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_cs(0:nucleus(itarget)%num_discrete,1:num_energies) = 0.0d0
         if(.not.allocated(Inelastic_count))                                                       &
             allocate(Inelastic_count(0:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_count(0:nucleus(itarget)%num_discrete,1:num_energies) = 0
         if(.not.allocated(Inelastic_g_mult))                                                      &
             allocate(Inelastic_g_mult(1:num_energies))
         Inelastic_g_mult(1:num_energies) = 0.0d0
         if(.not.allocated(Inelastic_g_mult_var))                                                  &
             allocate(Inelastic_g_mult_var(1:num_energies))
         Inelastic_g_mult_var(1:num_energies) = 0.0d0
!-------------
         if(.not.allocated(Inelastic_L_max))                                                       &
             allocate(Inelastic_L_max(0:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_L_max(1:nucleus(itarget)%num_discrete,1:num_energies) = 0
         if(.not.allocated(Inelastic_Ang_L))                                                       &
            allocate(Inelastic_Ang_L(0:Ang_L_max,0:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_Ang_L(0:Ang_L_max,0:nucleus(itarget)%num_discrete,1:num_energies) = 0.0d0
         if(.not.allocated(Inelastic_Ang_L_avg))                                                   &
            allocate(Inelastic_Ang_L_avg(0:Ang_L_max,0:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_Ang_L_avg(0:Ang_L_max,0:nucleus(itarget)%num_discrete,1:num_energies) = 0.0d0
         if(.not.allocated(Inelastic_Ang_dist))                                                    &
             allocate(Inelastic_Ang_dist(0:max_jx_10,1:nucleus(itarget)%num_discrete,1:num_energies))
         Inelastic_Ang_dist(0:max_jx_10,1:nucleus(itarget)%num_discrete,1:num_energies) = 0.0d0
!-------------
         if(.not.allocated(direct_Spectrum))allocate(direct_Spectrum(0:num_e))
         if(.not.allocated(dwba_Spectrum))allocate(dwba_Spectrum(0:num_e))
         dwba_Spectrum(0:num_e) = 0.0d0
         if(.not.allocated(res1))allocate(res1(0:num_e))
         res1(0:num_e) = 0.0d0
         if(.not.allocated(res2))allocate(res2(0:num_e))
         res2(0:num_e) = 0.0d0
         if(.not.allocated(resp))allocate(resp(0:6,0:num_e))
         resp(0:6,0:num_e) = 0.0d0
         if(.not.allocated(resp2))allocate(resp2(0:6,0:num_e))
         resp2(0:6,0:num_e) = 0.0d0
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

         if(.not.allocated(CE_cs))allocate(CE_cs(1:num_energies))
         CE_cs(1:num_energies) = 0.0d0
         if(.not.allocated(CE_Ang))allocate(CE_ang(0:Ang_L_max,1:num_energies))
         CE_ang(0:Ang_L_max,1:num_energies)=0.0d0

         if(.not.allocated(SE_cs))allocate(SE_cs(1:num_energies))
         SE_cs(1:num_energies) = 0.0d0
         if(.not.allocated(SE_Ang))allocate(SE_ang(0:Ang_L_max,1:num_energies))
         SE_ang(0:Ang_L_max,1:num_energies)=0.0d0
         if(.not.allocated(SE_prob))allocate(SE_prob(0:ixx_max))
!--------------------------------------------------------------------------------------------------
!------   Arrays for storing results from inelasitic scattering   
!--------------------------------------------------------------------------------------------------
!------    First is total of inelastic as a function of incident energy for each type,
!------    direct, compound to a discrete state, and compound through Continuous bins
!--------------------------------------------------------------------------------------------------
         if(.not.allocated(Inelastic_Total))allocate(Inelastic_Total(1:num_energies)) 
         if(.not.allocated(Inelastic_Continuous))allocate(Inelastic_Continuous(1:num_energies)) 
         if(.not.allocated(Inelastic_C_Discrete))allocate(Inelastic_C_Discrete(1:num_energies))
         if(.not.allocated(Inelastic_P_Discrete))allocate(Inelastic_P_Discrete(1:num_energies))
         if(.not.allocated(Inelastic_Direct))allocate(Inelastic_Direct(1:num_energies))
         inelastic_Total(1:num_energies) = 0.0d0
         inelastic_Continuous(1:num_energies) = 0.0d0
         inelastic_C_Discrete(1:num_energies) = 0.0d0
         inelastic_P_Discrete(1:num_energies) = 0.0d0
         inelastic_Direct(1:num_energies) = 0.0d0
         if(PREEQ_Model > 0)then
            if(.not.allocated(preeq_Spectrum))allocate(preeq_Spectrum(0:num_e))
            preeq_Spectrum(0:num_e) = 0.0d0
            if(.not. allocated(Inelastic_Preeq))allocate(Inelastic_Preeq(1:num_energies))
            inelastic_Preeq(1:num_energies) = 0.0d0
         end if
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
!
!-----    Series of dummy allocations for variables only used for a reaction calculation
!-----    allocated because the gfortran compiler gives "erroneous" warnings of uninitialized
!-----    stride and offset for these variables later on. They are accessed within if blocks 
!-----    specifying that it is NOT a population calculation (pop_calc = .false.), but this is 
!-----    conditional check seems to be lost. Thus, allocating to size one, and initializing to zero.
!-----    This eliminates the warnings, which will allow for other checks if a variable is not
!-----    properly initialized.
!
      else
         if(.not.allocated(Inelastic_cs))                                                          &
             allocate(Inelastic_cs(1,1))
         Inelastic_cs(1,1) = 0.0d0
         if(.not.allocated(Inelastic_count))                                                       &
             allocate(Inelastic_count(1,1))
         Inelastic_count(1,1) = 0
         if(.not.allocated(Inelastic_g_mult))                                                      &
             allocate(Inelastic_g_mult(1))
         Inelastic_g_mult(1) = 0.0d0
         if(.not.allocated(Inelastic_g_mult_var))                                                  &
             allocate(Inelastic_g_mult_var(1))
         Inelastic_g_mult_var(1) = 0.0d0
!-------------
         if(.not.allocated(Inelastic_L_max))                                                       &
             allocate(Inelastic_L_max(1,1))
         Inelastic_L_max(1,1) = 0
         if(.not.allocated(Inelastic_Ang_L))                                                       &
            allocate(Inelastic_Ang_L(1,1,1))
         Inelastic_Ang_L(1,1,1) = 0.0d0
         if(.not.allocated(Inelastic_Ang_L_avg))                                                   &
            allocate(Inelastic_Ang_L_avg(1,1,1))
         Inelastic_Ang_L_avg(1,1,1) = 0.0d0
         if(.not.allocated(Inelastic_Ang_dist))                                                    &
             allocate(Inelastic_Ang_dist(1,1,1))
         Inelastic_Ang_dist(1,1,1) = 0.0d0
!-------------
         if(.not.allocated(preeq_Spectrum))allocate(preeq_Spectrum(1))
         preeq_Spectrum(1) = 0.0d0
         if(.not.allocated(direct_Spectrum))allocate(direct_Spectrum(1))
         direct_Spectrum(1) = 0.0d0
         if(.not.allocated(dwba_Spectrum))allocate(dwba_Spectrum(1))
         dwba_Spectrum(1) = 0.0d0
         if(.not.allocated(res1))allocate(res1(1))
         res1(1) = 0.0d0
         if(.not.allocated(res2))allocate(res2(1))
         res2(1) = 0.0d0
         if(.not.allocated(resp))allocate(resp(1,1))
         resp(1,1) = 0.0d0
         if(.not.allocated(resp2))allocate(resp2(1,1))
         resp2(1,1) = 0.0d0
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
         if(.not.allocated(Inelastic_Continuous))allocate(Inelastic_Continuous(1)) 
         if(.not.allocated(Inelastic_C_Discrete))allocate(Inelastic_C_Discrete(1))
         if(.not.allocated(Inelastic_P_Discrete))allocate(Inelastic_P_Discrete(1))
         if(.not.allocated(Inelastic_Direct))allocate(Inelastic_Direct(1))
         inelastic_Total(1) = 0.0d0
         inelastic_Continuous(1) = 0.0d0
         inelastic_C_Discrete(1) = 0.0d0
         inelastic_P_Discrete(1) = 0.0d0
         inelastic_Direct(1) = 0.0d0
         if(PREEQ_Model > 0)then
            if(.not. allocated(Inelastic_Preeq))allocate(Inelastic_Preeq(1))
            inelastic_Preeq(1) = 0.0d0
         end if
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
      end if


      if(.not.allocated(x_particle_cs))allocate(x_particle_cs(1:num_energies,0:6))
      if(.not.allocated(x_particle_Spectrum))allocate(x_particle_Spectrum(0:num_e,0:6))
      if(.not.allocated(x_particle_Spectrum_dist))allocate(x_particle_Spectrum_dist(0:num_e,0:6))
      if(.not.allocated(x_particle_Spectrum_res))allocate(x_particle_Spectrum_res(0:num_e,0:6))
      if(.not.allocated(x_particle_Spectrum_res_dist))allocate(x_particle_Spectrum_res_dist(0:num_e,0:6))
      x_particle_cs(1:num_energies,0:6) = 0.0d0


      if(track_primary_gammas)allocate(Primary_Gamma_Spectrum(0:num_e))


      do i = 1, nucleus(1)%num_decay
         i_f = nucleus(1)%decay_to(i)
         if(i_f == itarget)ie_nuc = i_f
      end do


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Create core file name for output                               +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      call nucleus_label(itarget,ntar,target_label)
      ifile = ntar
      if(.not.pop_calc)then
         k=projectile%particle_type
         file_name(1:ntar+3)=particle(k)%label//'+'//target_label(1:ntar)//'-'
         ifile=ntar+3
      else
         file_name(1:ntar+1)=target_label(1:ntar)//'-'
         ifile = ntar + 1
      end if
      core_file=ifile
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------------------   Set up directories for output libraries
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(6,*)
      write(6,*)'Checking on directory structure for Data Libraries'
      f_exist = .false.
      lib_dir(1:ntar) = target_label(1:ntar)
      ilib_dir = ntar
      inquire(file = lib_dir(1:ilib_dir), exist = f_exist)
      if(.not. f_exist)then
         icmd = 6
         unix_cmd(1:icmd) = 'mkdir '
         icmd = icmd + 1
         unix_cmd(icmd:icmd + ilib_dir) = lib_dir(1:ilib_dir)
         icmd = icmd + ilib_dir
         write(6,*)unix_cmd(1:icmd)
         call system(unix_cmd(1:icmd))
      end if
      ilib_dir = ilib_dir + 1
      lib_dir(ilib_dir:ilib_dir) = '/'
      if(.not. pop_calc)then
         ilib_dir = ilib_dir + 1
         lib_dir(ilib_dir:ilib_dir) = particle(projectile%particle_type)%label
      else
         ilib_dir = ilib_dir + 1
         lib_dir(ilib_dir:ilib_dir+2) = 'Pop'
         ilib_dir = ilib_dir + 2
      end if
      inquire(file = lib_dir(1:ilib_dir), exist = f_exist)
      if(.not. f_exist)then
         icmd = 6
         unix_cmd(1:icmd) = 'mkdir '
         icmd = icmd + 1
         unix_cmd(icmd:icmd + ilib_dir) = lib_dir(1:ilib_dir)
         icmd = icmd + ilib_dir
         write(6,*)unix_cmd(1:icmd)
         call system(unix_cmd(1:icmd))
      end if
      if(.not. event_generator)then
         if(fission)then
            if_check = ilib_dir
            file_check(1:ilib_dir) = lib_dir
            if_check = if_check + 1
            file_check(if_check:if_check) = '/'
            if_check = if_check + 1
            file_check(if_check:if_check) = 'f'
            inquire(file = file_check(1:if_check), exist = f_exist)
            if(.not. f_exist)then
               icmd = 6
               unix_cmd(1:icmd) = 'mkdir '
               icmd = icmd + 1
               unix_cmd(icmd:icmd + if_check) = file_check(1:if_check)
               icmd = icmd + if_check
               write(6,*)unix_cmd(1:icmd)
               call system(unix_cmd(1:icmd))
            end if
         end if
         do i = 1, num_channels
            if_check = ilib_dir
            file_check(1:ilib_dir) = lib_dir
            if_check = if_check + 1
            file_check(if_check:if_check) = '/'
            if_check = if_check + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_check(if_check:if_check+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            if_check = if_check + ilast
            inquire(file = file_check(1:if_check), exist = f_exist)
            if(.not. f_exist)then
               icmd = 6
               unix_cmd(1:icmd) = 'mkdir '
               icmd = icmd + 1
               unix_cmd(icmd:icmd + if_check) = file_check(1:if_check)
               icmd = icmd + if_check
               write(6,*)unix_cmd(1:icmd)
                call system(unix_cmd(1:icmd))
            end if
         end do
      end if
      if_check = ilib_dir
      file_check(1:ilib_dir) = lib_dir
      if_check = if_check + 1
      file_check(if_check:if_check) = '/'
      if_check = if_check + 1
      ilast = 6
      file_check(if_check:if_check+ilast) = 'Spectra'
      if_check = if_check + ilast
      inquire(file = file_check(1:if_check), exist = f_exist)
      if(.not. f_exist)then
         icmd = 6
         unix_cmd(1:icmd) = 'mkdir '
         icmd = icmd + 1
         unix_cmd(icmd:icmd + if_check) = file_check(1:if_check)
         icmd = icmd + if_check
         write(6,*)unix_cmd(1:icmd)
         call system(unix_cmd(1:icmd))
      end if
!-----   Direct printing of event to their own directory


      if(dump_events)then
         if_check = ilib_dir
         file_check(1:ilib_dir) = lib_dir
         if_check = if_check + 1
         file_check(if_check:if_check) = '/'
         if_check = if_check + 1
         ilast = 10
         file_check(if_check:if_check+ilast) = 'Event-files'
         if_check = if_check + ilast
         inquire(file = file_check(1:if_check), exist = f_exist)
         if(.not. f_exist)then
            icmd = 6
            unix_cmd(1:icmd) = 'mkdir '
            icmd = icmd + 1
            unix_cmd(icmd:icmd + if_check) = file_check(1:if_check)
            icmd = icmd + if_check
            write(6,*)unix_cmd(1:icmd)
            call system(unix_cmd(1:icmd))
         end if
      end if

!**************************************************************************
!------    Write particle properties to file                           ---*
!**************************************************************************

      open(unit= 100, file = lib_dir(1:ilib_dir)//'/Particle-properties.dat', status = 'unknown')
      do k = -1, 6
         write(100,'(''#**************************'')')
!         ilast = index(particle(k)%name,' ')-1
!         if(ilast <= 0)ilast = 8
         write(100,'(''# '',a)')particle(k)%name
         write(100,'(''Z = '',i3)')particle(k)%Z
         write(100,'(''A = '',i3)')particle(k)%A
         write(100,'(''J = '',f4.1)')particle(k)%spin
         write(100,'(''Parity = '', i2)')nint(particle(k)%par)
         write(100,'(''Mass = '',1pe16.7,'' amu'')')particle(k)%mass/mass_u
      end do
      close(unit = 100)

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


      allocate(check_cs(1:num_comp))

      do k = 1, 6
         isp = nint(2*particle(k)%spin)
         write(100,*)particle(k)%name
         write(100,*)'lmax = ',particle(k)%lmax
         do i = 0, isp
            write(100,*)'isp = ', i
            do m = 1, particle(k)%nume
               write(100,*)particle(k)%e_grid(m),(particle(k)%trans_read(m,i,l),l = 0, particle(k)%lmax)
            end do
         end do
      end do








      do in = 1, num_energies

         check_cs(1:num_comp) = 0.0d0
         hang = 0
         sum_ce = 0
         sum_ie = 0


         E_in = projectile%energy(in)

         ibin = int(e_in/de)
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
                write(6,*)'This incident energy would populate discrete states only - skipping'
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

         end if
!-------------------------------------------------
            write(6,*)
         if(.not.pop_calc)then
            write(6,'(''-----------------------------------------------------------------------------'')')
            write(6,'(''*****************************************************************************'')')
            write(6,'(''-----------------------------------------------------------------------------'')')
            write(6,'(''Incident energy ='',1x,1pe15.7,'' MeV'')')E_in
            write(6,'(''Center-of-mass E ='',1x,1pe15.7,'' MeV'')')e_rel
            write(6,'(''Delta_E ='',1x,1pe15.7,'' MeV'')')de         
         else
            write(6,'(''-----------------------------------------------------------------------------'')')
            write(6,'(''*****************************************************************************'')')
            write(6,'(''-----------------------------------------------------------------------------'')')
            write(6,'(''Initial energy ='',1x,1pe15.7,'' MeV'')')E_in
            write(6,'(''Delta_E ='',1x,1pe15.7,'' MeV'')')de
         end if
!----------------------------------
         nucleus(1:num_comp)%Fiss_cs=0.0d0
         nbin = (nint((ex_tot-nucleus(1)%e_grid(1))/de) + 1)
         do i = 1, nucleus(1)%nbin
            if(ex_tot > nucleus(1)%e_grid(i) + nucleus(1)%delta_e(i) .and.   &
               nucleus(1)%e_grid(i) + nucleus(1)%delta_e(i) > ex_tot)then
               nbin = i
               exit
            end if
         end do

         if(fission)fission_cs(in) = 0.0d0

         x_particle_Spectrum(0:num_e,0:6) = 0.0d0
         if(.not. pop_calc .and. PREEQ_Model > 0)then
            preeq_Spect(0:6,0:num_e) = 0.0d0
            preeq_Spect_full(0:6,0:num_e) = 0.0d0
            preeq_Spectrum(0:num_e) = 0.0d0
            preeq_test_cs = 0.0d0
         end if
         if(.not. pop_calc)then
            direct_Spectrum(0:num_e) = 0.0d0
            dwba_Spectrum(0:num_e) = 0.0d0
         end if

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------    Compound nuclear cross section
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(.not. pop_calc)then
            reaction_cs(in) = 0.0d0
            absorption_cs(in) = 0.0d0
            call compound_xs(e_in,itarget,istate,iproj,absorption_cs(in),  &
                             num_channel,channel_prob,ichannel)

            write(6,'(''Reaction cross section = '',1x,1pe15.7,1x,a2)')absorption_cs(in)*cs_scale,cs_units
            write(13,'(''Initial populations for incident energy = '',1pe15.7,''MeV'')')e_in
            write(13,'(''   l     J  par    Pop'')')
            write(13,'(''  --- ----  ---  -------------'')')
            do nn = 1, num_channel
               ipar = 2*ichannel(4,nn)-1
               population = 0.0
               if(nn == 1)then
                  population = channel_prob(nn)
               else
                  population = channel_prob(nn) - channel_prob(nn-1)
               end if
               xI = ichannel(3,nn) + nucleus(1)%jshift
               write(13,'(1x,i3,1x,f5.1,1x,i4,2(1x,e15.6))')ichannel(1,nn),xI,ipar,channel_prob(nn),population
            end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------  Setup decay first nucleus in the chain. Apply width-fluctations    +
!------  if necessary. Thus, this nucleus is treated slightly differnet     +
!------  than others in the chain.                                          +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(WF_model == 1)then
               if(e_in > 0.5d0 .and. e_in < 1.0d0 .and. n_glag .ne. 40)then
                  n_glag = 40
                  alf = 0.0d0
                  bet = 0.0d0
                  call gauss_quad(n_glag, 2, alf, bet, x_glag, w_glag)
               elseif(e_in >= 1.0d0 .and. e_in < 1.75d0 .and. n_glag .ne. 25)then
                  n_glag = 25
                  alf = 0.0d0
                  bet = 0.0d0
                  call gauss_quad(n_glag, 2, alf, bet, x_glag, w_glag)
               elseif(e_in >= 1.75d0 .and. n_glag .ne. 20)then
                  n_glag = 20
                  alf = 0.0d0
                  bet = 0.0d0
                  call gauss_quad(n_glag, 2, alf, bet, x_glag, w_glag)
               end if
            end if 

            call HF_primary_decay_setup(e_in,iproj,itarget,1,istate,        &
                                        ex_tot,nbin,de)
         else

            reaction_cs(in) = 1.0d0
            absorption_cs(in) = 1.0d0

         end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Now set up arrays for pre-equilibrium emission                    +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         preeq_prob = 0.0d0
         if(.not. pop_calc .and. PREEQ_Model > 0 .and. absorption_cs(in) >= cs_threshold)then
            preeq_css(0:6,in) = 0.0d0
            preeq_spect(0:6,0:num_e) = 0.0d0
            preeq_spect_full(0:6,0:num_e) = 0.0d0
            preeq_prob = 0.0d0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----    Call pre-equilibrium model                                               +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            call pre_equilibrium_1(1,istate,in,e_in,ex_tot,de,absorption_cs(in))


            preeq_prob = nucleus(1)%PREEQ_cs(in)/absorption_cs(in)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Delete pre-equilibrium model data so that we can start afresh with      +
!-------   new energy.                                                             +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         end if


         if(pop_calc)then
!             Kinetic_Energy = nucleus(1)%Kinetic_Energy
            if(.not. j_pop_calc)then
               pop_sum = 0.0d0
               do i = 1, Pop_data(in)%num_pop
                  pop_sum = pop_sum + Pop_data(in)%bin_pop(i)
               end do
               absorption_cs(in) = pop_sum
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

               write(13,'(''Initial populations for incident energy ='',1x,1pe15.7)')e_in
               do ip = 0, 1
                  do j = 0, nucleus(1)%j_max
                     xj = j + nucleus(1)%jshift
                     write(13,'(2(1x,f4.1),2(1x,f10.6))')xj,2.0*real(ip)-1.0,target%pop_xjpi(j,ip)
                  end do
               end do
               do ijj = 1, num_pop
                  write(13,'(3(1x,i5),f10.6)')ijj,pop_j(ijj), pop_ip(ijj), pop_prob(ijj)
               end do
            else
               pop_sum = 0.0d0
               do j = 0, nucleus(1)%j_max
                  do ip = 0, 1
                     pop_sum = pop_sum + nucleus(1)%bins(j,ip,nbin)%rho
                  end do
               end do
               absorption_cs(in) = 1.0d0
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
          
         if(.not.pop_calc)then

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------  Allocate Direct component
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            reaction_cs(in) = absorption_cs(in)
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
                  SE_cs(in) = interpolate(0,E_in,OpticalCS%nume,OpticalCS%energy,OpticalCS%optical_cs(1,n))
                  do L = 0, Ang_L_max
                     SE_Ang(L,in) = interpolate(0,E_in,OpticalCS%nume,OpticalCS%energy,OpticalCS%optical_leg(1,L,n))
                  end do
                  cycle
               end if
               if(e_rel < OpticalCS%state(n)%energy)cycle 
               direct_cs(n) = interpolate(0,E_in,OpticalCS%nume,OpticalCS%energy,OpticalCS%optical_cs(1,n))
               tot_direct = tot_direct + direct_cs(n)
               reaction_cs(in) = reaction_cs(in) + direct_cs(n)
               direct_tot(in) = direct_tot(in) + direct_cs(n)
               if(OpticalCS%state(n)%state_type == 1)direct_cc(in) = direct_cc(in) + direct_cs(n)
               if(OpticalCS%state(n)%state_type <= 0)direct_dwba(in) = direct_dwba(in) + direct_cs(n)
               do L = 0, Ang_L_max
                  direct_Ang(L,n) = interpolate(0,E_in,OpticalCS%nume,OpticalCS%energy,OpticalCS%optical_leg(1,L,n))
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


            if(iproj == 1)then
               write(6,'(''Total cross section ='',1x,1pe15.7,1x,a2)')(reaction_cs(in) + SE_cs(in))*cs_scale,cs_units
               write(6,'(''Total Reaction cross section ='',1x,1pe15.7,1x,a2)')reaction_cs(in)*cs_scale,cs_units
               write(6,'(''Compound cross section ='',1x,1pe15.7,1x,a2)')absorption_cs(in)*cs_scale,cs_units
               write(6,'(''Direct excitation cross section ='',1x,1pe15.7,1x,a2)')tot_direct*cs_scale,cs_units
               write(6,'(''Shape Elastic cross section ='',1x,1pe15.7,1x,a2)')SE_cs(in)*cs_scale,cs_units
            elseif(iproj > 1)then
               write(6,'(''Total Reaction cross section ='',1x,1pe15.7,1x,a2)')reaction_cs(in)*cs_scale,cs_units
               write(6,'(''Compound cross section ='',1x,1pe15.7,1x,a2)')absorption_cs(in)*cs_scale,cs_units
               write(6,'(''Direct excitation cross section ='',1x,1pe15.7,1x,a2)')tot_direct*cs_scale,cs_units
            end if

!******************************************************************************
!------   If absorption_cs is too small skip all the decays as nothing can happen anyway
            if(absorption_cs(in) < cs_threshold)then
               write(6,'(''******************************************************************************'')')
               write(6,'(''*  Absorption cross section < '',e12.5,'', skipping decay processes   *'')')cs_threshold
               write(6,'(''******************************************************************************'')')
               goto 1901
            end if

!*****************************************************************************
!------
         else
            write(6,*)'Total population',absorption_cs(in)
         end if


         if(pop_calc)then
            avg_part(0:6) = 0.0d0
            avg_part_e(0:6) = 0.0d0
            var_part(0:6) = 0.0d0
            var_part_e(0:6) = 0.0d0
         end if




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


         directory(1:ilib_dir) = lib_dir(1:ilib_dir)
         idir = ilib_dir + 1
         directory(idir:idir) = '/'
         idir = idir + 1
         directory(idir:idir+10) = 'Event-files'
         idir = idir + 11
         directory(idir:idir) = '/'


         if(dump_events .and. binary_event_file)then
            open(unit=88,file=                                                                     &
                 directory(1:idir)//file_name(1:ifile)//'.Events.bin',                        &
                 form='unformatted',status='unknown')
            write(6,*)'Events written to binary file ', directory(1:idir)//file_name(1:ifile)//'.Events.bin'
         end if
         if(dump_events .and. .not. binary_event_file)then
            open(unit=88,file=                                                                     &
                 directory(1:idir)//file_name(1:ifile)//'.Events.txt',                        &
                 status='unknown')
            write(6,*)'Events written to formatted file ', directory(1:idir)//file_name(1:ifile)//'.Events.txt'
         end if           
        
         num_events = 0
         num_inelastic = 0

         num_compound = 0
         num_pre_equilibrium = 0
         num_fission = 0
         test_count = 0.0d0
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Start of Monte Carlo Sampling Loop     -----------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
         do i = 1, num_channels
            Exit_Channel(i)%num_event = 0
         end do

         passes = 0
         total_tries = 0
         max_tries = 0
         num_no_fiss = 0
         num_force_decay = 0
         tally_norm = 0.0d0
         max_e_diff = 0.0d0
         avg_e_diff = 0.0d0
         re_baseline = 0
!
!---------------------------------------------------------------------------------------
!--------    Monte Carlo sampling loop for each event                              -----
!---------------------------------------------------------------------------------------
!
         do nsamp = 1, num_mc_samp

            fission_decay = .false.

            num_part_type(1:6) = 0
            part_fact(0:7) = 1.0d0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------     Reset the boost matrix 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            Boost_COM(0:3,0:3) = 0.0d0
            do i = 0, 3
               Boost_COM(i,i) = 1.0d0
            end do
            hung_up = .false.

            event_g_mult = 0.0d0

!------    It is a population calculation. Set initial kinetic energy for the decaying system
            if(pop_calc)then
               KE = -1.0
               do while (KE < 0.0)
                   KE = init_Kinetic_Energy + dInit_Kinetic_Energy*Gauss_var(iseed)
               end do
               nucleus(1)%Kinetic_Energy = KE
!------    Now get initial excitation energy. This is set in Pop_data(i)%Ex_pop
!------    with Gaussian spread in Pop_data(i)%dEx_pop
10001          E_in = 100000000.
               do while (E_in > nucleus(1)%Ex_max .or. nbin <= 0)
                  E_in = Pop_data(in)%Ex_pop + Pop_data(in)%dEx_pop*Gauss_var(iseed)
                  Ex_tot = E_in
!                  ibin = int(E_in/de)
!                  nbin = (nint((Ex_tot-nucleus(1)%e_grid(1))/de) + 1)
                  nbin = 0
                  do i = 1, nucleus(1)%nbin
                     egrid = nucleus(icomp_f)%e_grid(i)
                     delta_e = nucleus(icomp_f)%delta_e(i)
                     if((Ex_tot >= egrid - 0.5d0*delta_e) .and.                    &
                         (Ex_tot < egrid + 0.5d0*delta_e))then
                          nbin = i
                          exit
                     end if
                  end do
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
                  absorption_cs(in) = 1.0
                  ijj = 0
                  pop_prob(1:num_pop) = 0.0
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
!            if(.not.pop_calc)part_Ang_data(0:Ang_L_max,1:dim_part) = 0.0
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
            if(.not.pop_calc)then
               ran = random_64(iseed)
               if(ran <= absorption_cs(in)/reaction_cs(in))then         !   Normal Compound formation
                  sp1=abs(spin_target-spin_proj)
                  call find_prob_point(num_channel,channel_prob,ran,ifind) 
                  icomp_i = 1
                  nbin_i = nbin
                  l_i = ichannel(1,ifind)
                  is_i = ichannel(2,ifind)
                  Ix_i = ichannel(3,ifind)
                  ip_i = ichannel(4,ifind)
                  
                  xnbin_i = nbin_i
                  par_i = 2*ip_i - 1
                  ran = random_64(iseed)
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

                     num_pre_equilibrium = num_pre_equilibrium + 1
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
                     num_compound = num_compound + 1
                     call MC_primary_decay(iproj, spin_target,                           &
                                           l_i, is_i, Ix_i, ex_tot, icomp_i,             &
                                           icomp_f, Ix_f, ip_f, nbin_f, idb,             &
                                           n_dat, dim_part, num_part_type, part_fact,    &
                                           num_part, part_data,                          &
                                           Ang_L_max, part_Ang_data,                     &
                                           ixx_max, delta_ix,                            &
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
                  ran = random_64(iseed)
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
!---- Check this xip_f
                           xip_f = OpticalCS%state(n)%parity
                           ip_f = nint((xip_f + 1.0d0)/2.0d0)
!
                           Ex = OpticalCS%state(n)%E_min +                          &
                                random_64(iseed)*OpticalCS%state(n)%Delta_E

                           nbin_f = int((Ex - nucleus(itarget)%e_grid(1) - de/2.0d0)/de) + 1
                           if(nbin_f < 0)nbin_f = 1
                           Ex = nucleus(itarget)%e_grid(nbin_f)
 
                           E_f = ex_tot - nucleus(icomp_i)%sep_e(iproj) - Ex
                           Ex_f = nucleus(itarget)%e_grid(nbin_f)
                           dwba_decay = .true.
                        elseif(OpticalCS%state(n)%state_type == -1)then    !  Direct to continuum bin w/o explicit J
                           idb = 0
                           ran = random_64(iseed)
                           do Ix_f = OpticalCS%state(n)%Ix_min, OpticalCS%state(n)%Ix_max   !  find final spin for DWBA K state
                              if(ran <= OpticalCS%state(n)%spin_prob(Ix_f))exit
                           end do
                           xI_f = real(Ix_f,kind=8) + nucleus(icomp_f)%jshift
                           xip_f = OpticalCS%state(n)%parity
                           ip_f = nint((xip_f + 1.0d0)/2.0d0)
                           Ex = OpticalCS%state(n)%E_min +                          &
                                random_64(iseed)*OpticalCS%state(n)%Delta_E

                           nbin_f = int((Ex - nucleus(itarget)%e_grid(1) - de/2.0d0)/de) + 1
                           if(nbin_f < 0)nbin_f = 1
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
                        do nang = 1, num_theta
                           ran = random_64(iseed)
                           do i = 1, ixx_max
                              x = real(i,kind=8)*delta_ix - 1.0d0 -                                   &
                                  random_64(iseed)*delta_ix*0.999999d0
                              if(direct_prob(i,n) > ran)exit
                           end do
                           extra_angle_data(nang,num_part) = acos(x)
                        end do
                        
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------   Now phi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        ran = random_64(iseed)
                        phi_0 = ran*2.0d0*pi

                        part_data(10,num_part) = extra_angle_data(1,num_part)
                        part_data(11,num_part) = phi_0
!                        part_data(12,num_part) = T_1
!                        part_data(13,num_part) = theta
!                        part_data(14,num_part) = phi
!                        part_data(15,num_part) = T_L
!                        part_data(16,num_part) = theta_L
!                        part_data(17,num_part) = phi_L
!                        part_data(18,num_part) = T_2
                        part_data(19,num_part) = 1.0d0
                        part_data(20,num_part) = Ex_f
                        part_data(21,num_part) = icomp_i
                        part_data(22,num_part) = nucleus(itarget)%state(istate)%spin
                        part_data(23,num_part) = nucleus(itarget)%state(istate)%parity
                        part_data(24,num_part) = istate
!                        nucleus(icomp_f)%Kinetic_Energy = T_2

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
            else                                   !   population calculation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   This is a population decay rather than a reaction 
!---------   Population for J and parity were input for a given starting
!---------   exciation energy. Much like reaction calculation, conduct 
!---------   first decay and then pass on to later routines to decay 
!---------   to terminus at either the ground state of some nucleus or an
!---------   isomer where subsequent decay is forbidden.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

               compound = .true.
               ran = random_64(iseed)
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
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   Check if decay is hung up in a bin with no place to go. If so, force
!----------   decay to a discrete state and collect statistics to rpeort later   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
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
               write(6,*)'Ix_i > j_max before decay loop, Ix_i =', Ix_i
               write(6,*)'Check: pop_calc = ',pop_calc
               write(6,*)'Check: compound = ',compound
               write(6,*)'Check: preeq_decay = ', preeq_decay
               write(6,*)'Check: direct = ',direct
               write(6,*)'Check: cc_decay = ',cc_decay
               write(6,*)'Check: dwba_decay = ',dwba_decay
               write(6,*)'idb = ',idb
               write(89,'(2(1x,i10),3(1x,f10.5))')nsamp, num_part, E_in,                  &
                    nucleus(1)%Kinetic_Energy, tally_weight
               do i = 1, num_part
                  write(89,'(1x,i5,21(2x,f12.5))')i,(part_data(k,i), k = 1, n_dat)
               end do
               flush(89)
               stop
            end if


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Now decay continous energy bins until they reach a discrete state (occurs when idb=1)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            do while(idb /= 1)
               if(Ix_i > j_max)then
                  write(6,*)'error Ix_i > j_max , Ix_i = ',Ix_i
                  flush(6)
                  stop
               end if
               call MC_decay_bin(icomp_i, Ix_i, ip_i, nbin_i,                   &
                                 icomp_f, Ix_f, ip_f, nbin_f, idb,              &
                                 n_dat, dim_part, num_part_type, part_fact,     &
                                 num_part, part_data,                           &
                                 Ang_L_max, part_Ang_data,                      &
                                 num_theta, extra_angle_data)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Following is information relating to the number of times the decay "failed" due to
!------   energy conservation. namely due to random noise put into emission energy, the excitation
!------   energy of the decaying state fell below normally allowed decay. If this happens, the 
!------   subroutine will first try several times to bypass it, but if the decay is more probably than the 
!------   rest, for example neutron emission relative to gamma, then the subroutine will attempt to 
!------   re-baseline the decays, removing the "forbidden" decay and notify user. This is something somewhat 
!------   more prominant near thresholds and when discrete states are allowed above e_cut.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              total_tries = total_tries + tries
              passes = passes + 1
              if(tries > max_tries)max_tries = tries
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


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------    Get channel number                                          +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------    Now decay all discrete states to the bottom, but only if
!-------    we are not collecting inelastic to explicit discrete states
!-------    They do not get decayed here.
!-------    Decay if this is a population calculation.
!-------    
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------   ichann = label for this channel, now have only discrete gammas left
!----------------   in = index for incident energy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
            ichann = find_channel(n_dat, dim_part, num_part, part_data, num_part_type)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Finish decay all the way to the end     ----------------------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
            call MC_decay_state(icomp_i, nbin_i,                             &
                                n_dat,dim_part,num_part,part_data,           & 
                                Ang_L_max, part_Ang_data,                    &
                                num_theta, extra_angle_data, ichann, in)

101         continue

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
                      if(nucleus(icomp_f)%state(nbin_f)%isomer)sum_e = sum_e +       &
                         nucleus(icomp_f)%state(nbin_f)%energy
                   end if
                end if               
            end do          
            tally_norm = tally_norm + tally_weight

            if(.not. fission_decay)then
               num_no_fiss = num_no_fiss + 1
               e_diff = abs(ex_tot - sum_e)
               if(e_diff > max_e_diff)max_e_diff = e_diff
               avg_e_diff = avg_e_diff + e_diff
               if(e_diff >= 1.0d-2)then
                  write(28,*)'Energy not conserved to within 10 keV'
                  write(28,*)'E_in = ',E_in,' nsamp = ',nsamp
                  write(28,*)'Ex_tot = ',ex_tot,' Sum_e ',sum_e
                  write(28,*)'Diff = ',e_diff
                  write(28,*)'Ix_I > j_max before decay loop, Ix_i =', Ix_i
                  write(28,*)'Check: pop_cal = ',pop_calc
                  write(28,*)'Check: compound = ',compound
                  write(28,*)'Check: preeq_decay = ', preeq_decay
                  write(28,*)'Check: direct = ',direct
                  write(28,*)'Check: dwba = ',dwba_decay
                  do nn = 1, num_part
                     write(28,'(1x,i5,21(2x,f12.5))')nn,(part_data(k,nn), k = 1, n_dat)
                  end do
                  flush(28)
               end if
            end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    Add "spread" to emitted energies based on size of the bins
!------    Helps to remove periodicity in the outgoing spectra as the transitions occur
!------    on a fixed grid.
!------    At each stage, we also need to transform to Lab and COM frames from the rest
!------   frame where each decay occurred
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
                  shift = random_64(iseed)*(shift_max-shift_min) + shift_min
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
                  extra_angle_data(num_theta+1,nn) = theta               !   COM frame
                  extra_angle_data(2*num_theta+1,nn) = theta_L           !   Lab frame
               end do 
            end do

!
!----   if difference exceeds 1.0d-4 MeV write out to unit 400
!----   Count number of events and warn user at the end
!
            if(.not. fission .and. abs(ex_tot-sum_e) > 1.0d-4)then
                num_bad_e = num_bad_e + 1
                write(400,*)'Energy conservation issue',ex_tot,sum_e,abs(ex_tot-sum_e)
                write(400,*)'Incident energy = ', e_in, 'Sample # ',nsamp
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
                   write(400,*)i,k, idb, e_rel, theta, phi, icomp_i, nbin_i,icomp_f,nbin_f
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
                write(6,*)'Error in sampling weight for unbiased sampling'
                write(6,*)nsamp,tally_weight
                stop
            end if


            if(fission_decay)then       !   Fission channel
               ictype = 0
               ichann = find_channel(n_dat, dim_part, num_part, part_data, num_part_type)
               Exit_Channel(ichann)%Channel_cs(in,ictype) =                               &
                      Exit_Channel(ichann)%Channel_cs(in,ictype) + tally_weight
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


    if(ichann < 1)then
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

               if(hung_up)Exit_Channel(ichann)%Channel_cs(in,-1) =                        &
                          Exit_Channel(ichann)%Channel_cs(in,-1) + tally_weight

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    This is an error trap. if num_part == 0 something went wrong
!------    No decay seems to have occurred
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

               if(num_part == 0)then
                  write(6,*)'num_part == 0'
                  write(6,*)'Ix_I > j_max before decay loop, Ix_i =', Ix_i
                  write(6,*)'Check: pop_calc = ',pop_calc
                  write(6,*)'Check: compound = ',compound
                  write(6,*)'Check: preeq_decay = ', preeq_decay
                  write(6,*)'Check: direct = ',direct
                  write(6,*)'Check: cc_decay = ',cc_decay
                  write(6,*)'Check: dwba_decay = ',dwba_decay
                  write(6,*)'idb = ',idb
                  stop
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
                     do nang = 1, num_theta
                        theta = extra_angle_data(nang+1,nn)
                        x = cos(theta)
                        jx = nint((x+1.0d0)/delta_jx_10)
                        if(jx < 0)jx = 0
                        if(jx > max_jx_10)jx = max_jx_10
                        Inelastic_Ang_Dist(jx,jstate,in) = Inelastic_Ang_Dist(jx,jstate,in) +      &
                                                           tally_weight/real(num_theta,kind=8)
                     end do
                     do L = 0, Ang_L_max
                        Inelastic_Ang_L_avg(L,jstate,in) = Inelastic_Ang_L_avg(L,jstate,in) +      &
                                                           part_Ang_data(L,nn)*tally_weight
                     end do
                  else                                          !  Inelastic in some other fashion
                     Inelastic_cs(0,in) = Inelastic_cs(0,in) + tally_weight
                     Inelastic_count(0,in) = Inelastic_count(0,in) + 1
                     Inelastic_total(in) = Inelastic_total(in) + tally_weight
                     event_g_mult = 0.0d0
                     Exit_channel(ichann)%Channel_cs(in,ictype) =                                  &
                         Exit_channel(ichann)%Channel_cs(in,ictype) + tally_weight
                     Exit_Channel(ichann)%num_event = Exit_Channel(ichann)%num_event + 1
                     do nn = 1, num_part                           !  loop over particles and collect spectra
                        jstate = nint(part_data(5,nn))
                        idb = nint(part_data(6,nn))
                        jproj = nint(part_data(2,nn))
                        k = jproj
!------    If it is a discrete state, add a bit of smear to avoid accidental collision with boundary of spectrum
!------    bin
                        e_shift = 0.0d0
                        if(idb == 1)e_shift = (2.0d0*random_64(iseed) - 1.0d0)*de_spec*0.5d0
                        icc = int((part_data(12,nn) + e_shift)/de_spec) + 1
!                        theta = part_data(16,nn)
                        if(k >= 0 .and. icc >= 0 .and. icc <= num_e)then
                           Exit_Channel(ichann)%part_mult(k,ictype,in) =                           &
                                Exit_Channel(ichann)%part_mult(k,ictype,in) + tally_weight
                           Exit_Channel(ichann)%Spect(k,ictype,in)%E_spec(icc) =                   &
                               Exit_Channel(ichann)%Spect(k,ictype,in)%E_spec(icc) +               &
                               tally_weight/de_spec
                           Exit_Channel(ichann)%Spect(k,ictype,in)%E_count(icc) =                  &
                               Exit_Channel(ichann)%Spect(k,ictype,in)%E_count(icc) + 1
                           do nang = 1, num_theta
                              theta = extra_angle_data(num_theta+1,nn)
                              x = cos(theta)
                              jx = nint((x+1.0d0)/delta_jx_10)
                              if(jx < 0)jx = 0
                              if(jx > max_jx_10)jx = max_jx_10
                              Exit_Channel(ichann)%Spect(k,ictype,in)%E_Ang_Dist(jx,icc) =            &
                                  Exit_Channel(ichann)%Spect(k,ictype,in)%E_Ang_Dist(jx,icc) +        &
                                  tally_weight/(delta_jx_10*de_spec)/real(num_theta,kind=8)
                              Exit_Channel(ichann)%Spect(k,ictype,in)%Ang_Dist(jx) =                  &
                                  Exit_Channel(ichann)%Spect(k,ictype,in)%Ang_Dist(jx) +              &
                                  tally_weight/delta_jx_10/real(num_theta,kind=8)
                           end do
                        end if
                     end do

                     Inelastic_g_mult(in) = Inelastic_g_mult(in) + event_g_mult
                     Inelastic_g_mult_var(in) = Inelastic_g_mult_var(in) + event_g_mult**2
                  end if
               else

                  Exit_channel(ichann)%Channel_cs(in,ictype) =                                     &
                      Exit_channel(ichann)%Channel_cs(in,ictype) + tally_weight
                  Exit_Channel(ichann)%num_event = Exit_Channel(ichann)%num_event + 1
                  do nn = 1, num_part
                     k = nint(part_data(2,nn))
!------    If it is a discrete state, add a bit of smear to avoid accidental collision with boundary of spectrum
!------    bin
                     icc = int(part_data(12,nn)/de_spec) + 1
                     if(k >= 0 .and. icc >= 0 .and. icc <= num_e)then
                        Exit_Channel(ichann)%part_mult(k,ictype,in) =                              &
                             Exit_Channel(ichann)%part_mult(k,ictype,in) + tally_weight
                        Exit_Channel(ichann)%Spect(k,ictype,in)%E_spec(icc) =                      &
                            Exit_Channel(ichann)%Spect(k,ictype,in)%E_spec(icc) +                  &
                            tally_weight/de_spec
                        Exit_Channel(ichann)%Spect(k,ictype,in)%E_count(icc) =                     &
                            Exit_Channel(ichann)%Spect(k,ictype,in)%E_count(icc) + 1
                        do nang = 1, num_theta
                           theta = extra_angle_data(num_theta + 1,nn)
                           x = cos(theta)
                           jx = nint((x+1.0d0)/delta_jx_10)
                           if(jx < 0)jx = 0
                           if(jx > max_jx_10)jx = max_jx_10
                           Exit_Channel(ichann)%Spect(k,ictype,in)%E_Ang_Dist(jx,icc) =               &
                               Exit_Channel(ichann)%Spect(k,ictype,in)%E_Ang_Dist(jx,icc) +           &
                               tally_weight/(delta_jx_10*de_spec)/real(num_theta,kind=8)
                           Exit_Channel(ichann)%Spect(k,ictype,in)%Ang_Dist(jx) =                     &
                               Exit_Channel(ichann)%Spect(k,ictype,in)%Ang_Dist(jx) +                 &
                               tally_weight/delta_jx_10/real(num_theta,kind=8)
                        end do
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
                        write(20,*)'num_part = ', num_part
                        write(20,*)'Ix_I > j_max before decay loop, Ix_i =', Ix_i
                        write(20,*)'Check: pop_calc = ',pop_calc
                        write(20,*)'Check: compound = ',compound
                        write(20,*)'Check: preeq_decay = ', preeq_decay
                        write(20,*)'Check: direct = ',direct
                        write(20,*)'Check: cc_decay = ',cc_decay
                        write(20,*)'Check: dwba_decay = ',dwba_decay
                        write(20,*)'idb = ',idb
                        write(6,*)'Issue with part_Lab_spectrum inde > num_spect_e'
                     end if      
                  end if
               end do
            end if


            if(preeq_decay)then
               k = nint(part_data(2,1))
               if(k >= 0 .and. k <= 7)then
                  preeq_css(k,in) = preeq_css(k,in) + tally_weight
                  icc = int(part_data(15,1)/de_spec2)
                  if(icc >= 0 .and. icc <= num_e)then
                     preeq_spect_full(k,icc) = preeq_spect_full(k,icc) + tally_weight
                     if(.not. fission_decay)                                                          &
                        preeq_spect(k,icc) = preeq_spect(k,icc) + tally_weight
                        preeq_test_cs = preeq_test_cs + tally_weight
                        test_count = test_count + tally_weight
                  end if
               end if
            end if
            if(direct .and. cc_decay)then
               icc = int(part_data(15,1)/de_spec2)
               if(icc >= 0 .and. icc <= num_e)direct_spectrum(icc) = direct_spectrum(icc) + tally_weight
            end if
            if(direct .and. dwba_decay)then
               icc = int(part_data(15,1)/de_spec2)
               if(icc >= 0 .and. icc <= num_e)dwba_spectrum(icc) = dwba_spectrum(icc) + tally_weight
            end if
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Finished Monte Carlo Sampling Loop     -----------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
         end do                      !   end do nsamp loop 

!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Finished Monte Carlo Sampling Loop     -----------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
         do k = 0, 6
            smear(k) = 0.2d0
         end do

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
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------  Collect particle spectra ssociated with various reaction types, pre-equilibirum   ----
!-------  direct, dwba, etc.                                                                ----
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Collect spectra for each particle type (iproj,x*jproj)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(.not. pop_calc)then
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
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Pre-equilibrium decays 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
            if(PREEQ_Model > 0)then
               preeq_test_cs = preeq_test_cs*reaction_cs(in)/tally_norm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Fold in a resolution smear
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               do k = 0,6
                  resp(k,0:num_e) = 0.0d0 
                  resp2(k,0:num_e) = 0.0d0 
                  do i = 0, num_e
                     e1 = real(i,kind=8)*de_spec2
                     do j = 0, num_e
                        e2 = real(j,kind=8)*de_spec2
                        resp(k,j) = resp(k,j) +                                                    &
                           Preeq_spect(k,i)*exp(-((e1-e2)/smear(1))**2/2.0d0)/sqrt(2.0d0*pi*smear(1)**2)
                        resp2(k,j) = resp2(k,j) +                                                  &
                           Preeq_spect_full(k,i)*exp(-((e1-e2)/smear(1))**2/2.0d0)/sqrt(2.0d0*pi*smear(1)**2)
                     end do
                  end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Normalize for distribution
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  xnorm1 = 0.0d0
                  xnorm2 = 0.0d0
                  do i = 0, num_e
                     e1 = real(i,kind=8)*de_spec2
                     xnorm1 = xnorm1 + resp(k,i)*de_spec2
                     xnorm2 = xnorm2 + Preeq_Spect(k,i)*de_spec2
                  end do
                  if(xnorm1 > 0.0d0)then
                     do i = 0, num_e
                        e1 = real(i,kind=8)*de_spec2
                        resp(k,i) = resp(k,i)*xnorm2/xnorm1
                        if(resp(k,i) < 1.0d-20)resp(k,i) = 1.0d-20
                     end do
                  end if
                  xnorm1 = 0.0d0
                  xnorm2 = 0.0d0
                  do i = 0, num_e
                     e1 = real(i,kind=8)*de_spec2
                     xnorm1 = xnorm1 + resp2(k,i)*de_spec2
                     xnorm2 = xnorm2 + Preeq_Spect_full(k,i)*de_spec2
                  end do
                  if(xnorm1 > 0.0d0)then
                     do i = 0, num_e
                        e1 = real(i,kind=8)*de_spec2
                        resp2(k,i) = resp2(k,i)*xnorm2/xnorm1
                        if(resp2(k,i) < 1.0d-20)resp2(k,i) = 1.0d-20
                     end do
                  end if
               end do


               directory(1:ilib_dir) = lib_dir(1:ilib_dir)
               idir = ilib_dir + 1
               directory(idir:idir) = '/'
               idir = idir + 1
               directory(idir:idir+6) = 'Spectra'
               idir = idir + 7
               directory(idir:idir) = '/'

               open(unit=24,file=                                                                        &
                   directory(1:idir)//file_name(1:ifile)//'-Preeq-Spectrum.dat',status='unknown')
               write(24,'(''#Pre-equilibirum spectrum for each particle type '')')
               write(24,'(''#Calculated as barns/MeV '')')
               write(24,'(''#Incident energy = '',f10.4)')e_in
               if(fission)then
                  write(24,'(''#First set includes fission events, second set excludes fission events'')')
                  write(24,'(''#'',10x,14(1x,13x,i2))')(k, k = 0, nucleus(1)%num_decay-1),         &
                                                       (k, k = 0, nucleus(1)%num_decay-1)
               else
                  write(24,'(''#'',10x,14(1x,13x,i2))')(k, k = 0, nucleus(1)%num_decay-1)
               end if
               do i = 0, num_e
                  energy = real(i,kind=8)*de_spec2
                  if(fission)then
                     write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,                                  &
                                (resp2(k,i),k = 0, nucleus(1)%num_decay-1),                        &
                                (resp(k,i),k = 0, nucleus(1)%num_decay-1)
                     energy = energy + de_spec2
                     write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,                                  &
                                (resp2(k,i),k = 0, nucleus(1)%num_decay-1),                        &
                                (resp(k,i),k = 0, nucleus(1)%num_decay-1)
                  else
                     write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,(resp2(k,i),k = 0, nucleus(1)%num_decay-1)
                     energy = energy + de_spec2
                     write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,(resp2(k,i),k = 0, nucleus(1)%num_decay-1)
                  end if
               end do
               close(unit=24)
            end if

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Direct and DWBA states decays 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
            if(OpticalCS%numcc > 1)then
               res1(0:num_e) = 0.0d0 
               res2(0:num_e) = 0.0d0 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Fold in a resolution smear
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               do i = 0, num_e
                  e1 = real(i,kind=8)*de_spec2
                  do j = 0, num_e
                     e2 = real(j,kind=8)*de_spec2
                     res1(j) = res1(j) +                                                           &
                        direct_Spectrum(i)*exp(-((e1-e2)/smear(1))**2/2.0d0)/sqrt(2.0d0*pi*smear(1)**2)
                     res2(j) = res2(j) +                                                           &
                        dwba_Spectrum(i)*exp(-((e1-e2)/smear(1))**2/2.0d0)/sqrt(2.0d0*pi*smear(1)**2)
                  end do
               end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Normalize for distribution
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               xnorm1 = 0.0d0
               xnorm2 = 0.0d0
               xnorm3 = 0.0d0
               xnorm4 = 0.0d0
               do i = 0, num_e
                  e1 = real(i,kind=8)*de_spec2
                  xnorm1 = xnorm1 + res1(i)*de_spec2
                  xnorm2 = xnorm2 + res2(i)*de_spec2
                  xnorm3 = xnorm3 + direct_Spectrum(i)*de_spec2
                  xnorm4 = xnorm4 + dwba_Spectrum(i)*de_spec2
               end do
               if(xnorm1 > 0.0d0)then
                  do i = 0, num_e
                     e1 = real(i,kind=8)*de_spec2
                     res1(i) = res1(i)*xnorm3/xnorm1
                     if(res1(i) < 1.0d-20)res1(i) = 1.0d-20
                  end do
               end if
               if(xnorm2 > 0.0d0)then
                  do i = 0, num_e
                     e1 = real(i,kind=8)*de_spec2
                     res2(i) = res2(i)*xnorm4/xnorm2
                     if(res2(i) < 1.0d-20)res2(i) = 1.0d-20
                  end do
               end if

               directory(1:ilib_dir) = lib_dir(1:ilib_dir)
               idir = ilib_dir + 1
               directory(idir:idir) = '/'
               idir = idir + 1
               directory(idir:idir+6) = 'Spectra'
               idir = idir + 7
               directory(idir:idir) = '/'

               open(unit=24,file=                                                                        &
                   directory(1:idir)//file_name(1:ifile)//'-Direct-Spectrum.dat',status='unknown')
               write(24,'(''#Outgoing spectra for direct and DWBA reactions '')')
               write(24,'(''#Calculated as barns/MeV '')')
               write(24,'(''#Incident energy = '',f10.4)')e_in
               write(24,'(''#'',10x,10x,''Direct'',12x,''DWBA'')')
               do i = 0, max_num
                  energy = real(i,kind=8)*de_spec2
                  write(24,'(1x,f10.4,14(1x,1pe15.7))')energy, res1(i), res2(i)
                  energy = energy + de_spec2
                  write(24,'(1x,f10.4,14(1x,1pe15.7))')energy, res1(i), res2(i)
               end do
               close(unit=24)
            end if
         end if       ! ---- finished with spectra associated with reactions only

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

!     write(6,*)'num ',xnum_eff, xnum_elastic, num_elastic

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Here, collect the spectra for particles emitted and do not distinguish by channel
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
         do k = 0, 6
            if(k == iproj .and. iproj == 1 .and. .not. pop_calc)then
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
                  ran = random_64(iseed)
                  do j = 1, ixx_max
                     x = real(j,kind=8)*delta_ix - 1.0d0 -                                         &
                         random_64(iseed)*delta_ix*0.9999999d0
                    if(SE_prob(j) > ran)exit
                  end do
                  theta_0 = acos(x)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------   Now phi
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  ran = random_64(iseed)
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
                  if(icc <= num_e)x_particle_Spectrum(icc,k) = x_particle_Spectrum(icc,k) + 1.0d0
                  x_particle_cs(in,k) = x_particle_cs(in,k) + 1.0d0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------     Write out events to file     -----------------------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                  if(dump_events .and. .not. binary_event_file)then                           !   Formatted output
                     write(88,'(2(1x,i10),3(1x,f10.5))')nsamp, num_part, E_in,                     &
                         nucleus(1)%Kinetic_Energy, tally_weight
                     do j = 1, num_part
                        write(88,'(1x,i5,21(2x,f12.5))')j,(part_data(m,j), m = 1, n_dat)
                     end do
                     flush(88)
                  end if
                  if(dump_events .and. binary_event_file)then                           !   Unformatted output
                     write(88)nsamp, num_part, E_in, nucleus(1)%Kinetic_Energy,tally_weight,       &
                             ((part_data(m,j), m = 1, n_dat), j =1, num_part)
                     flush(88)
                  end if

               end do
            end if
            x_particle_cs(in,k) = x_particle_cs(in,k)*temp_cs/xnum_eff
            xnorm = 0.0d0
            do i = 0, num_e
               x_particle_Spectrum(i,k) = x_particle_Spectrum(i,k)*temp_cs/xnum_eff/de_spec2
               xnorm = xnorm + x_particle_Spectrum(i,k)*de_spec2
            end do

            if(xnorm > 0.0d0)then
               do i = 0, num_e
                  x_particle_Spectrum_dist(i,k) = x_particle_Spectrum(i,k)/xnorm
               end do
            end if
         end do

         directory(1:ilib_dir) = lib_dir(1:ilib_dir)
         idir = ilib_dir + 1
         directory(idir:idir) = '/'
         idir = idir + 1
         directory(idir:idir+6) = 'Spectra'
         idir = idir + 7
         directory(idir:idir) = '/'

         open(unit=24,file=                                                                        &
              directory(1:idir)//file_name(1:ifile)//'-xParticle-Spectrum.dat',status='unknown')
         write(24,'(''#Outgoing spectra for each particle type '')')
         write(24,'(''#Incident energy = '',f10.4)')e_in
         write(24,'(''#Calculated as barns/MeV '')')
         write(24,'(''#k='',8x,7(10x,i2,4x))')(k, k = 0, 6)
         write(24,'(''#'',10x,7(3x,''----------------''))')
         x_particle_Spectrum_res(0:num_e,0:6) = 0.0d0 
         x_particle_Spectrum_res_dist(0:num_e,0:6) = 0.0d0 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Fold in a resolution smear
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do k = 0, 6
            do i = 0, num_e
               e1 = real(i,kind=8)*de_spec2
               do j = 0, num_e
                  e2 = real(j,kind=8)*de_spec2
                  x_particle_spectrum_res(j,k) = x_particle_spectrum_res(j,k) +                    &
                       x_particle_Spectrum(i,k)*exp(-((e1-e2)/smear(k))**2/2.0d0)/                 &
                        sqrt(2.0d0*pi*smear(k)**2)
               end do
            end do
         end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Normalize the distribution
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do k = 0, 6
            xnorm1 = 0.0d0
            xnorm2 = 0.0d0
            do i = 0, num_e
               e1 = real(i,kind=8)*de_spec2
               xnorm1 = xnorm1 + x_particle_Spectrum_res(i,k)*de_spec2
               xnorm2 = xnorm2 + x_particle_Spectrum(i,k)*de_spec2
            end do
            if(xnorm1 > 0.0d0)then
               do i = 0, num_e
                  e1 = real(i,kind=8)*de_spec2
                  x_particle_Spectrum_res_dist(i,k) = x_particle_Spectrum_res(i,k)/xnorm1
                  if(x_particle_Spectrum_res_dist(i,k) < 1.0d-20)                                  &
                     x_particle_Spectrum_res_dist(i,k) = 1.0d-20
               end do
            end if
            if(xnorm2 > 0.0d0)then
               do i = 0, num_e
                  e1 = real(i,kind=8)*de_spec2
                  x_particle_Spectrum_res(i,k) = x_particle_Spectrum_res(i,k)*xnorm2/xnorm1
                  if(x_particle_Spectrum_res(i,k) < 1.0d-20)x_particle_Spectrum_res(i,k)=1.0d-20
               end do
            end if
         end do

         do i = 0, max_num
            energy = real(i,kind=8)*de_spec2
            write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,                                           &
               (x_particle_Spectrum_res(i,k)*cs_scale,k=0, 6),                                     &
               (x_particle_Spectrum_res_dist(i,k),k=0, 6)    
            energy = energy + de_spec2
            write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,                                           &
               (x_particle_Spectrum_res(i,k)*cs_scale,k=0, 6),                                     &
               (x_particle_Spectrum_res_dist(i,k),k=0, 6)    
         end do
         close(unit=24)


         if(dump_events)close(unit=88)


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------------   Distribute cs hung in discrete states to final states             +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(biased_sampling)then
            write(6,*)'Statistics on continuous bin with hung decays'
            write(6,'(''                                         #hung              #total      fraction'')')
            write(6,'(''                              ----------------    ----------------    ----------'')')
            write(13,*)'Statistics on continuous bin with hung decays'
            write(13,'(''                                         #hung              #total      fraction'')')
            write(13,'(''                              ----------------    ----------------    ----------'')')
         else
            write(6,*)'Statistics on continuous bin with hung decays'
            write(6,'(''                                         #hung              #total      fraction       Num events'')')
            write(6,'(''                              ----------------    ----------------    ----------     ------------'')')
            write(13,*)'Statistics on continuous bin with hung decays'
            write(13,'(''                                         #hung              #total      fraction       Num events'')')
            write(13,'(''                              ----------------    ----------------    ----------     ------------'')')
         end if
         do ichann = 1, num_channels
            hang = Exit_Channel(ichann)%Channel_cs(in,-1)
            sum = 0.0d0
            do j = 1, Exit_Channel(ichann)%num_cs
               sum = sum + Exit_Channel(ichann)%Channel_cs(in,j)
            end do
            ratio = 0.0d0
            ilast = index(Exit_Channel(ichann)%Channel_Label,' ') - 1
            if(sum > 0.0d0)ratio = Exit_Channel(ichann)%Channel_cs(in,-1)/sum
            if(biased_sampling)then
               write(6,'(''Channel: '',i4,1x,a12,2(4x,f16.3),4x,f10.4''%'')')                       &
               ichann, Exit_Channel(ichann)%Channel_Label(1:ilast), hang, sum, ratio*100.
               write(13,'(''Channel: '',i4,1x,a12,2(4x,f16.3),4x,f10.4''%'')')                      &
               ichann, Exit_Channel(ichann)%Channel_Label(1:ilast), hang, sum, ratio*100.
            else
               write(6,'(''Channel: '',i4,1x,a12,2(4x,f16.3),4x,f10.4''%'',4x,i12)')                &
               ichann, Exit_Channel(ichann)%Channel_Label(1:ilast),                                &
               hang, sum, ratio*100., Exit_Channel(ichann)%num_event
               write(13,'(''Channel: '',i4,1x,a12,2(4x,f16.3),4x,f10.4''%'',4x,i12)')               &
               ichann, Exit_Channel(ichann)%Channel_Label(1:ilast),                                &
               hang, sum, ratio*100., Exit_Channel(ichann)%num_event
            end if
         end do
         write(6,*)'Hung bins are forced to decay to discrete states randomly'
         write(13,*)'Hung bins are forced to decay to discrete states randomly'

         icharr = 32
         out_buff(1:icharr) = 'For more details consult file : '
         il1 = index(out_file,' ')-1
         out_buff(icharr+1:icharr+il1)=out_file(1:il1)
         icharr = icharr + il1
         out_buff(icharr+1:icharr+15) = '.sample_trouble'
         icharr = icharr + 15
         if(re_baseline > 0)then
           write(6,*)
           write(6,'(''*********************************************************************************'')')
           write(6,'(i10,'' Decays out of '',i10,'' had to be re-baselined due to energy'')')re_baseline,num_mc_samp
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
           do k = 1, 6
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
           end do
        end if
        if(.not. pop_calc .and. .not. event_generator)then
           if(Inelastic_cs(0,in) >= 1.0d-8)then
              Inelastic_g_mult(in) = Inelastic_g_mult(in)/Inelastic_cs(0,in)
              Inelastic_g_mult_var(in) = Inelastic_g_mult_var(in)/Inelastic_cs(0,in)
              Inelastic_g_mult_var(in) = sqrt(dabs(Inelastic_g_mult_var(in) - Inelastic_g_mult(in)**2))
           end if
           if(.not.allocated(xvalue))allocate(xvalue(0:max_jx_10))
           do jx = 0, max_jx_10
              xvalue(jx) = real(jx,kind=8)*delta_jx_10 - 1.0d0
           end do
           LL_max = max_jx_10 - 2
           do j = 1, nucleus(itarget)%num_discrete
!-----   Normalize parameters for angular distributions for inelastic scattering to discrete states 
              if(inelastic_cs(j,in) > 1.0d-6)then
                 do L = 0, Ang_L_max
                    Inelastic_Ang_L_avg(L,j,in) = Inelastic_Ang_L_avg(L,j,in)/Inelastic_cs(j,in)
                 end do
              end if
!-----   Inelastic cross sections
              Inelastic_cs(j,in) = Inelastic_cs(j,in)*reaction_cs(in)/tally_norm
              Inelastic_Ang_Dist(0,j,in) = Inelastic_Ang_Dist(0,j,in)*2.0d0
              Inelastic_Ang_Dist(max_jx_10,j,in) = Inelastic_Ang_Dist(max_jx_10,j,in)*2.0d0

              sum = 0.0d0
              do jx = 0, max_jx_10 - 1
                 sum = sum + (Inelastic_Ang_Dist(jx,j,in) + Inelastic_Ang_Dist(jx+1,j,in))*         &
                              delta_jx_10*0.5d0
              end do
             
              Inelastic_L_max(j,in) = 0
              Inelastic_Ang_L(0:Ang_L_max,j,in) = 0.0d0
              Inelastic_Ang_L(0,j,in) = 0.5d0
             
              if(sum/tally_norm >= 1.0d-6)then

                 do jx = 0, max_jx_10
                    Inelastic_Ang_Dist(jx,j,in) = Inelastic_Ang_Dist(jx,j,in)/sum
                 end do

                 Inelastic_L_max(j,in) = 8
                 if(Inelastic_count(j,in) < 4000)Inelastic_L_max(j,in) = 6
                 if(Inelastic_count(j,in) < 2000)Inelastic_L_max(j,in) = 4
                 if(Inelastic_count(j,in) < 1000)Inelastic_L_max(j,in) = 2
                 if(Inelastic_count(j,in) < 500)Inelastic_L_max(j,in) = 0
                 call Legendre_expand(max_jx_10+1,xvalue(0),Inelastic_Ang_Dist(0,j,in),            &
                                      Inelastic_L_max(j,in),Inelastic_Ang_L(0,j,in))
              end if
           end do
           if(allocated(xvalue))deallocate(xvalue)
        end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------    Now the other Channels    -------------------------------------------------
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        check_sum = 0.0d0
!
        if( .not. event_generator)then
           if(.not.allocated(xvalue))allocate(xvalue(0:max_jx_10))
           do jx = 0, max_jx_10
              xvalue(jx) = real(jx,kind=8)*delta_jx_10 - 1.0d0
           end do
!
           do i = 1, num_channels
              inuc = Exit_Channel(i)%Final_nucleus
              n_min = 1
              if(fission)n_min = 0
              do n = n_min, Exit_Channel(i)%num_cs
                 do k = 0, 6                                  !   Loop over particle types
                    if(Exit_Channel(i)%Channel_cs(in,n) >= 1.0d-8)then
                       Exit_Channel(i)%part_mult(k,n,in) =                                          &
                          Exit_Channel(i)%part_mult(k,n,in)/Exit_Channel(i)%Channel_cs(in,n) 
                    end if
!-----   New approach to Spectrum and Angular Distributions
                    do icc = 0, num_e
                       if(Exit_Channel(i)%Spect(k,n,in)%E_count(icc) == 0)cycle
                       Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(0,icc) =                            &
                          Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(0,icc)*2.0d0
                       Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(max_jx_10,icc) =                    &
                          Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(max_jx_10,icc)*2.0d0
                       sum = 0.0d0
!------    Normalize for each energy
                       do jx = 1, max_jx_10
                          sum = sum + (Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(jx,icc) +           &
                                       Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(jx-1,icc))*         &
                                      delta_jx_10*0.5d0
                       end do
                       if(sum > 1.0d-8)then
                          do jx = 0, max_jx_10
                             Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(jx,icc) =                     &
                                  Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(jx,icc)/sum
                          end do
                       end if

                       sum = Exit_Channel(i)%Spect(k,n,in)%E_count(icc)
!-----   Set max L based on number of counts in energy bin
                       Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_max(icc) = 8
                       if(sum < 8000)Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_max(icc) = 6
                       if(sum < 6000)Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_max(icc) = 4
                       if(sum < 4000)Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_max(icc) = 2
                       if(sum < 2000)Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_max(icc) = 0
!-----   Calculate average, avgerage difference from average, and expected difference
!-----   Check if data have enough statistics to warrant more Legendre coefficients
                       if(Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_max(icc) > 0)then
                          avg = 0.0d0
                          do jx = 0, max_jx_10
                             avg = avg + Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(jx,icc)
                          end do
                          avg = avg/real(max_jx_10 + 1,kind=8)
                          avg_diff = 0.0d0
                          do jx = 0, max_jx_10
                             avg_diff = avg_diff + abs(Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(jx,icc) - avg)
                          end do
                          avg_diff = avg_diff/real(max_jx_10+1,kind=8)
                          avg_diff = avg_diff/avg
                          expected_diff = real(Exit_Channel(i)%Spect(k,n,in)%E_count(icc),kind=8)/    &
                                          real(max_jx_10 + 1,kind=8)
                          expected_diff = sqrt(expected_diff)/expected_diff/sqrt(real(num_theta,kind=8))
                          if(avg_diff < 2.5d0*expected_diff)Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_max(icc) = 4
                          if(avg_diff < 1.75d0*expected_diff)Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_max(icc) = 2
                          if(avg_diff < expected_diff)Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_max(icc) = 0
                       end if

                       call Legendre_expand(max_jx_10+1,xvalue(0),                                    &
                                            Exit_Channel(i)%Spect(k,n,in)%E_Ang_Dist(0,icc),          &
                                            Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_max(icc),           &
                                            Exit_Channel(i)%Spect(k,n,in)%E_Ang_L(0,icc))

!----------------------------------------------------------------------------------------------------
                    end do                       !  End of Loop over icc (outgoing energies
                 end do                          !  End of Loop over k, particles

              end do
              if(track_gammas)then
                 do n = 1, nucleus(inuc)%num_discrete
                    do j = 1, nucleus(inuc)%state(n)%nbranch
                       Exit_channel(i)%state(n)%cs(j,in) =                                            &
                          Exit_channel(i)%state(n)%cs(j,in)/tally_norm
                   end do
                 end do
             end if                               
           end do
!------    Deallocate temporary array with x-values
           if(allocated(xvalue))deallocate(xvalue)
        end if
!

!-----   Normalize the Channel cross sections and collect check_sum to make sure all
!-----   events are accounted for
        do i = 1, num_channels
           n_min = 1
           if(fission)n_min = 0
           do n = n_min, Exit_Channel(i)%num_cs
              Exit_Channel(i)%Channel_cs(in,n) =                                                   &
                   Exit_Channel(i)%Channel_cs(in,n)/tally_norm
              if(n > 0)check_sum = check_sum +                                                     &
                           Exit_Channel(i)%Channel_cs(in,n)
           end do
        end do


        check_sum = check_sum + sum_inelastic
        if(fission)check_sum = check_sum + fission_cs(in)
        write(6,*)
        write(6,*)'Check that all events are accounted for'
        write(6,'('' Check_sum = '',f12.8)') check_sum
        if(dabs(check_sum - 1.0d0) <= 1.0d-8)write(6,*)'All events are accounted for'

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------  if we are tracking gammas, print them out now, rather than with other channels 
!-------  after we finish all the incident energies. Print to a separate file for each 
!-------  incident energy
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ifile = core_file
        file_name(ifile+1:ifile+1) = '-'
        ifile =ifile +1
        if(e_in < 10.0)then
           file_name(ifile+1:ifile+1)='-'
           ifile=ifile+1
           write(file_name(ifile+1:ifile+6),'(f6.4)')e_in
           ifile=ifile+6
        else
           write(file_name(ifile+1:ifile+7),'(f7.4)')e_in
           ifile=ifile+7
        end if


         if(track_gammas .and. .not. event_generator)then
            ilab = ntar
            file_lab(1:ntar) = target_label(1:ntar)
            ilab = ilab + 1
            file_lab(ilab:ilab) = '('
            ilab = ilab + 1
            file_lab(ilab:ilab) = particle(iproj)%label
            ilab = ilab + 1
            file_lab(ilab:ilab) = ','
            do i = 1, num_channels
               inuc = Exit_Channel(i)%Final_nucleus
               call nucleus_label(inuc,nres,resid_label)
               directory(1:ilib_dir) = lib_dir(1:ilib_dir)
               idir = ilib_dir + 1
               directory(idir:idir) = '/'
               ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
               directory(idir+1:idir+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
               idir = idir + ilast
               idir = idir + 1
               directory(idir:idir) = '/'

               Q = Exit_Channel(i)%Q_value
               Q = -Q

               ifile = 8
               outfile(1:ifile) = 'Channel_'
               outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
               ifile = ifile + ilast
               outfile(ifile+1:ifile+5) = '_Ein_'
               ifile = ifile + 5
               if(e_in < 10.0)then
                  outfile(ifile+1:ifile+1) ='0'
                  ifile = ifile + 1
                  write(outfile(ifile+1:ifile+6),'(f6.4)')e_in
                  ifile = ifile + 6
               else
                  write(outfile(ifile+1:ifile+7),'(f7.4)')e_in
                  ifile = ifile + 7
               end if
               open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'_cs_gammas.dat',status = 'unknown')
               file_lab2(1:20) = ' '
               ilab2 = ilab
               file_lab2(1:ilab2) = file_lab(1:ilab)
               ilab2 = ilab2 + 1
               ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
               file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
               ilab2 = ilab2 + ilast
               file_lab2(ilab2:ilab2) = ')'
               write(100,'(''# '',a20)')file_lab2
!-----------------------------------
               write(100,'(''#'')')
               write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
               write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
               if(iproj > 0 .and. iproj < 7)then
                  write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
               else
                  write(100,'(''# Population decay, no projectile'')')
               end if
               write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
               write(100,'(''#'')')
               if(Exit_channel(i)%num_particles > 0)then
                  if(explicit_channels)then
                     write(100,'(''# Mass of emitted particles '')')
                     do m = 1, Exit_channel(i)%num_particles
                        k = Exit_channel(i)%decay_particles(m)
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m, particle(k)%mass/mass_u
                     end do
                  else
                     write(100,'(''# Mass of emitted particles '')')
                     m = 0
                     do k = 1, 6
                        do nn = 1, Exit_channel(i)%num_part(k)
                           m = m + 1
                           write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                           write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m,particle(k)%mass/mass_u
                        end do
                     end do
                  end if
               end if

               write(100,'(''#'')')
               write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
               write(100,'(''#'')')
               write(100,'(''# Listing discrete gammas observed in this particle decay channel'')')
               write(100,'(''# Note these are not exclusive to a specific discrete state'')')
               write(100,'(''# Incident Energy = '',f10.5)')E_in
               sum = 0.0d0

               do nn = 1, Exit_Channel(i)%num_cs
                  sum = sum + Exit_channel(i)%Channel_cs(in,nn)
               end do

               if(.not. pop_calc)then
                  sum = sum*reaction_cs(in)*cs_scale
                  write(100,'(''#'')')
                  write(100,'(''# Total cross section to this exit decay channel = '',1pe15.7,1x,a2)')sum,cs_units
                  write(100,'(''#'')')
               else
                  sum = sum*reaction_cs(in)
                  write(100,'(''#'')')
                  write(100,'(''# Fraction of decay strength to this exit decay channel = '',1pe15.7,1x,a2)')sum
                  write(100,'(''#'')')
               end if


               write(100,'(''# Cross section for discrete gamma transitions in '',a5)')resid_label(1:nres)
               write(100,'(''#'')')
               write(100,1900)cs_units
1900           format('#   i    Jp     f    Jp     e-gamma (MeV)        sigma (',a2,')     sigma(g)/channel')
               write(100,'(''#  --- ----    --- ----     -------------     -------------     -------------'')')

               check_sum = 0.0d0
               do k = 1, nucleus(inuc)%num_discrete
                  do m = 1, nucleus(inuc)%state(k)%nbranch
                     nf  = nucleus(inuc)%state(k)%ibranch(m)
                     ipi = nint((nucleus(inuc)%state(k)%parity+1.0d0)/2.0d0)
                     ipf = nint((nucleus(inuc)%state(nf)%parity+1.0d0)/2.0d0)
                     cs = Exit_channel(i)%state(k)%cs(m,in)*reaction_cs(in)*cs_scale
                     ratio = 0.0d0
                     if(sum > 1.0d-10)ratio = Exit_channel(i)%state(k)%cs(m,in)/sum
                     check_sum = check_sum + ratio 
                     write(100,'(1x,i4,1x,f4.1,a1,''  '',i4,1x,f4.1,a1,8x,f10.5,3x,e15.7,3x,e15.7,3x,e15.7)')  &
                                                           k,nucleus(inuc)%state(k)%spin,ch_par(ipi),          &
                                                           nf,nucleus(inuc)%state(nf)%spin,ch_par(ipf),        &
                                                           nucleus(inuc)%state(k)%egamma(m),                   &
                                                           cs,ratio
                  end do
               end do
               close(unit = 100)
            end do                               
         end if

 1901    continue

         if(.not.pop_calc)then
            if(associated(channel_prob))nullify(channel_prob)
            if(associated(ichannel))nullify(ichannel)
            spin_target=nucleus(itarget)%state(istate)%spin
            spin_proj=particle(iproj)%spin
            sp1=abs(spin_target-spin_proj)
            sp2=spin_target+spin_proj
            isp=nint(sp2-sp1)
            isp_max = nint(2.0d0*spin_proj)
            Ix = 0
            do l = 0, particle(iproj)%lmax
               xj_min = abs(dfloat(l) - spin_proj)
               xj_max = abs(dfloat(l) + spin_proj)
               isp = nint(xj_max - xj_min)
!               do is = 0, isp
               xj = real(l,kind=8) - spin_proj
               do is = 0, isp_max
                  xj = xj + dfloat(is)
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
!-------   Notify user about events that failed energy conservation limits

     if(num_bad_e > 0)then
        ifile=index(out_file,' ') - 1
        write(6,*)'There were ',num_bad_e,' events that violated energy conservation by 0.1 keV or more'
        write(6,*)'Events written to ',out_file(1:ifile)//'.bad_energy'
     else
        write(6,*)'There were no events that violated energy conservation at 0.1 keV level'
     end if


!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------  Output data to Library Files     -----------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------  Arrays to sum emission spectrum and angular distributions  -+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--rem      if(.not.allocated(E_spec))allocate(E_spec(0:num_e,0:6))
!--rem      if(.not.allocated(E_norm))allocate(E_norm(0:6))
!--rem      if(.not.allocated(Ang_norm))allocate(Ang_norm(0:6))
!---------------------------------------------------------------------+

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Core Directory     -------------------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      directory(1:ilib_dir) = lib_dir(1:ilib_dir)
      idir = ilib_dir + 1
      directory(idir:idir) = '/'

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Decay pattern for each for each nucleus   --------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      do i = 1,num_comp
         iA = nucleus(i)%A
         if(iA < 10)then
            ifile = 1
            write(outfile(1:ifile),'(i1)')iA
         elseif(iA < 99)then
            ifile = 2
            write(outfile(1:ifile),'(i2)')iA
         elseif(iA < 1000)then
            ifile = 3
            write(outfile(1:ifile),'(i3)')iA
         end if
         if(nucleus(i)%atomic_symbol(1:1) == ' ')then
            outfile(ifile+1:ifile+1) = nucleus(i)%atomic_symbol(2:2)
            ifile = ifile + 1
         else
            outfile(ifile+1:ifile+2) = nucleus(i)%atomic_symbol(1:2)
            ifile = ifile + 2
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'-decay-properties.dat',status = 'unknown')            
         write(100,'(''#Decay properties of discrete states for'',a5)')outfile(1:ifile)
         write(100,'(''#Mass = '',1pe23.16,'' MeV'')')nucleus(i)%mass
         write(100,'(''#Mass = '',1pe23.16,'' amu'')')nucleus(i)%mass/mass_u
         write(100,'(''#AMU = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''#'')')
         do k = 1, nucleus(i)%num_discrete, 1
            write(100,'(''    #       Ex     J   par #decay     Isomer'')')
            write(100,'(''  ---     ------  --- ----  -----  --------------'')')
            iso_label = 'Isomer = False'
            if(nucleus(i)%state(k)%isomer)iso_label = 'Isomer = True '
            write(100,'(1x,i4,1x,f10.4,2(1x,f4.1),1x,i5,3x,a14)')k-1,nucleus(i)%state(k)%energy,    &
                        nucleus(i)%state(k)%spin,nucleus(i)%state(k)%parity,                        &
                        nucleus(i)%state(k)%nbranch, iso_label
            write(100,'(''#      i --->   f            branch      prob_gamma         prob_ic'')')
            write(100,'(''#    ---      ---     -------------   -------------   -------------'')')
            do m = 1, nucleus(i)%state(k)%nbranch
               write(100,'(4x,i4,'' --->'',i4,4(3x,1pe15.7))')k-1,nucleus(i)%state(k)%ibranch(m)-1, &
                                                              nucleus(i)%state(k)%branch(m),        &
                                                              nucleus(i)%state(k)%p_gamma(m),       &
                                                              nucleus(i)%state(k)%p_ic(m)
            end do
            write(100,'(''#'')')
         end do
      end do


      if(event_generator)then
         write(6,*)'****************************************************'
         write(6,*)'----   Finished simulating and writing events.   ---'
         write(6,*)'----   All finished.                             ---'
         write(6,*)'****************************************************'
         call exit
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


      if(.not. pop_calc)then
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Total cross section          ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(iproj == 1)then
            ifile = 8
            outfile(1:ifile) = 'Total_cs'
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = 'X'
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                        &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''# Total Cross Section data '')')
            write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
            write(100,'(''#'',2(''   ----------------''))')
            do in = 1, num_energies
               e_in=projectile%energy(in)
               cs = (reaction_cs(in) + SE_cs(in))*cs_scale
               write(100,'(1x,4(3x,1pe16.7))')e_in, cs
            end do
            close(unit=100)
         end if
         
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Reaction cross section          ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ifile = 11
         outfile(1:ifile) = 'Reaction_cs'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = 'X'
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''# Reaction Cross Section data = Total - Elastic'')')
         write(100,'(''# Reaction Cross Section data = Absorption + direct'')')
         write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
         write(100,'(''#'',2(''   ----------------''))')
         do in = 1, num_energies
            e_in=projectile%energy(in)
            cs = reaction_cs(in)*cs_scale
            write(100,'(1x,4(3x,1pe16.7))')e_in, cs
         end do
         close(unit=100)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Absorption cross section          ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ifile = 13
         outfile(1:ifile) = 'Absorption_cs'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = 'X'
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''# Absorption Cross Section data'')')
         write(100,'(''# Absorption Cross Section data = Reaction - Direct - Compound Elastic'')')
         write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
         write(100,'(''#'',2(''   ----------------''))')
         do in = 1, num_energies
            e_in=projectile%energy(in)
            cs = absorption_cs(in)*cs_scale
            write(100,'(1x,4(3x,1pe16.7))')e_in, cs
         end do
         close(unit=100)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Preequilibrium cross section          ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
         if(PREEQ_Model > 0)then
            ifile = 17
            outfile(1:ifile) = 'Preequilibrium_cs'
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = 'X'
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''# Preequilibrium Cross section. Units = '',a2)')
            write(100,'(''#'',3x,6x,''E_in'',6x,3x,4x,''Tot('',a2,'')'',7(3x,6x,a1,''('',a2,'')'',5x))')  &
               cs_units,(particle(k)%label,cs_units,k = 0,6)
            write(100,'(''#'',9(''   ----------------''))')
            do in = 1, num_energies
               e_in=projectile%energy(in)
               cs = 0.0d0
               do k = 0, 6
                  cs = cs + preeq_css(k,in)
               end do
               cs = cs*reaction_cs(in)*cs_scale
               write(100,'(1x,9(3x,1pe16.7))')e_in, cs, (preeq_css(k,in)*reaction_cs(in)*cs_scale, k = 0, 6)
            end do
            close(unit=100)
         end if
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Direct cross sections           ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
         ifile = 9
         outfile(1:ifile) = 'Direct_cs'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = 'X'
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''# Direct Cross Section data'')')
         write(100,'(''# Reaction Cross Section data'')')
         write(100,1100)cs_units,cs_units,cs_units
 1100    format('#         E_in      ',10x,'cc(',a2,')',3x,8x,'dwba(',a2,')',3x,                    &
                      9x,'tot(',a2,')')
         write(100,'(''#'',4(''   ----------------''))')
         do in = 1, num_energies
            e_in=projectile%energy(in)
            cs = absorption_cs(in)*cs_scale
            write(100,'(1x,4(3x,1pe16.7))')e_in, direct_cc(in)*cs_scale, direct_dwba(in)*cs_scale,  &
                    direct_tot(in)*cs_scale
         end do
         close(unit=100)

      end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Fission cross section           ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(fission)then
            directory(1:ilib_dir) = lib_dir(1:ilib_dir)
            idir = ilib_dir + 1
            directory(idir:idir) = '/'
            idir = idir + 1
            directory(idir:idir) = 'f'
            idir = idir + 1
            directory(idir:idir) = '/'
            ifile = 10
            outfile(1:ifile) = 'fission_cs'
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = 'f'
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                        &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''# Fission cross section data '')')
            write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
            write(100,'(''#'',2(''   ----------------''))')
            do in = 1, num_energies
               e_in=projectile%energy(in)
               write(100,'(1x,4(3x,1pe16.7))')e_in, fission_cs(in)*reaction_cs(in)*cs_scale
            end do
            close(unit=100)
!--------------   Output data on fission <J> <J**2>
            open(unit=24,file='Fission-Stats.dat',status='unknown')
            do in = 1, num_energies
               E_in = projectile%energy(in)
               do icomp = 1, num_comp
                  if(Fiss_tally(icomp,in) > 0.0d0)then
                     Fiss_J_avg(icomp,in) = Fiss_J_avg(icomp,in)/Fiss_tally(icomp,in)
                     Fiss_J_var(icomp,in) = Fiss_J_var(icomp,in)/Fiss_tally(icomp,in)
                     Fiss_J_var(icomp,in) = sqrt(abs(Fiss_J_var(icomp,in) - Fiss_J_avg(icomp,in)**2))
                  end if
               end do
               write(24,'(25(1x,f10.5))')E_in,(Fiss_J_avg(icomp,in),Fiss_J_var(icomp,in),icomp = 1, num_comp)
            end do
            close(unit=24)
         end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------  Elastic and inelastic channels             ---------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

      if(.not. pop_calc)then          ! Information for elastic and inelastic channels - 
         directory(1:ilib_dir) = lib_dir(1:ilib_dir)
         idir = ilib_dir + 1
         directory(idir:idir) = '/'
         idir = idir + 1
         directory(idir:idir) = particle(iproj)%label
         idir = idir + 1
         directory(idir:idir) = '/'
                                      ! don't have with population calculation
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Elastic scattering      --------------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         if(iproj == 1)then     !   Doesn't make sense for charged particles
            ifile = 16
            outfile(1:ifile) = 'Shape_Elastic_cs'
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
!            file_lab2(ilab2:ilab2) = 'n'
            file_lab2(ilab2:ilab2) = particle(iproj)%label
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            j = istate
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                        &
                  nucleus(itarget)%state(istate)%energy
            if(iproj == 1)write(100,'(''# Shape Elastic cross section data '')')
            if(iproj > 1)write(100,'(''# <Shape Elastic/Rutherford> '')')
!            write(100,'(''# To Level # '',i3)')j
            write(100,'(''#'')')
            write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
            write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            if(iproj > 0 .and. iproj < 7)then
               write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
            else
               write(100,'(''# Population decay, no projectile'')')
            end if
            write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
            write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
            write(100,'(''#'')')
            if(iproj == 1)then
               write(100,'(''#          E_in              xs'',''('',a2,'')'')')cs_units
            elseif(iproj > 1)then
               write(100,'(''#          E_in         <OM/Rutherford>'')')
            end if
            write(100,'(''#'',2(''   ----------------''))')
            do in = 1, num_energies
               cs = SE_cs(in)
               e_in = projectile%energy(in)
               if(iproj == 1)then
                  write(100,'(1x,4(3x,1pe16.7))')e_in, cs*cs_scale
               elseif(iproj > 1)then
                  write(100,'(1x,4(3x,1pe16.7))')e_in, cs*0.5d0
               end if
            end do
            close(unit=100)
!         end if

         if(.not.allocated(Ang_Dist))allocate(Ang_Dist(0:max_jx_100))
         ifile = 17
         outfile(1:ifile) = 'Shape_Elastic_Ang'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = particle(iproj)%label
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         j = istate
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
               nucleus(itarget)%state(istate)%energy
         if(iproj == 1)write(100,'(''# Shape Elastic Angular Distribution data '')')
         if(iproj == 1)write(100,'(''# Angular Distribution normalized to unity '')')
         if(iproj > 1)write(100,'(''# Elastic Scattering Angular Distribution data '')')
         if(iproj > 1)write(100,'(''# Calculated as ratio to Rutherford: OM/Rutherford '')')
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
         write(100,'(''#'')')
         write(100,'(''# Frame = COM'')')
         alf = 0.0d0
         bet = 0.0d0
         do in = 1, num_energies
            e_in = projectile%energy(in)
            if(iproj == 1)then
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                          e_in,SE_cs(in),cs_units
            elseif(iproj > 1)then
               write(100,'(''# E_in = '',1pe16.7,3x,''OM/Rutherford = '',1pe16.7)')            &
                          e_in,SE_cs(in)
            end if
            if(iproj == 1)then
               write(100,'(''#         E_in            cos(theta)            Prob'')')
            elseif(iproj > 1)then
               write(100,'(''#         E_in            cos(theta)        OM/Rutherford'')')
            end if
            write(100,'(''#'',3(''   ----------------''))')
            do jx = 0, max_jx_100
               x = real(jx,kind=8)*delta_jx_100 - 1.0d0
               if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.999845
               alf = 0.0d0
               bet = 0.0d0
               sum = 0.0d0
               do L = 0, Ang_L_max
                  sum = sum + SE_Ang(L,in)*poly(L,1,alf,bet,x)
               end do
               Ang_Dist(jx) = sum
            end do
            xnorm = 0.0d0
            do jx = 1, max_jx_100
               xnorm = xnorm + (Ang_Dist(jx-1) + Ang_Dist(jx))*delta_jx_100*0.5d0
            end do
            do jx = 0, max_jx_100
               x = real(jx,kind=8)*delta_jx_100 - 1.0d0
               Temp = Ang_Dist(jx)
               if(iproj == 1 .and. xnorm > 1.0d-20)Temp = Temp/xnorm
               write(100,'(1x,3(3x,1pe16.7))')e_in, x, Temp
            end do
         end do
         close(unit=100)
         if(allocated(Ang_Dist))deallocate(Ang_Dist)
!
!---   Also print out Legendre Coefficients
!
         ifile = 21
         outfile(1:ifile) = 'Shape_Elastic_Ang_Leg'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = particle(iproj)%label
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         j = istate
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
               nucleus(itarget)%state(istate)%energy
         if(iproj == 1)write(100,'(''# Shape Elastic Angular Distribution data '')')
         if(iproj == 1)write(100,'(''# Angular Distribution normalized to unity '')')
         if(iproj == 1)write(100,'(''# Angular Distribution normalized to unity '')')
         if(iproj > 1)write(100,'(''# Elastic Scattering Angular Distribution data '')')
         if(iproj > 1)write(100,'(''# Calculated as ratio to Rutherford: OM/Rutherford '')')
         write(100,'(''# Coefficients of Legendre expansion'')')
         write(100,'(''# Ang_dist(x) = Sum_L a(L)*P(L,x)'')')
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
         write(100,'(''#'')')
         write(100,'(''# Frame = COM'')')

         do in = 1, num_energies
            e_in = projectile%energy(in)
            if(iproj == 1)then
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                          e_in,SE_cs(in),cs_units
            elseif(iproj > 1)then
               write(100,'(''# E_in = '',1pe16.7,3x,''OM/Rutherford = '',1pe16.7)')                 &
                          e_in,SE_cs(in)
            end if
            write(100,'(''#         E_in            L            a(L)'')')
            write(100,'(''#   ----------------     ---     ----------------'')')
            xnorm = 2.0d0*SE_Ang(0,in)
            do L = 0, Ang_L_max
               Temp = SE_Ang(L,in)
               if(iproj == 1)Temp = SE_Ang(L,in)/xnorm
               write(100,'(1x,3x,1pe16.7,3x,i5,3x,1pe16.7)')e_in, L, Temp
            end do 
         end do
         close(unit=100)


         ifile = 19
         outfile(1:ifile) = 'Compound_Elastic_cs'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = particle(iproj)%label
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         j = istate
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''# Compound Elastic cross section data '')')
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
         write(100,'(''#'')')
         write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
         write(100,'(''#'',2(''   ----------------''))')
         do in = 1, num_energies
            cs = Inelastic_cs(istate,in)
            e_in=projectile%energy(in)
            write(100,'(1x,4(3x,1pe16.7))')e_in, cs*cs_scale
         end do
         close(unit=100)


         if(.not.allocated(Ang_Dist))allocate(Ang_Dist(0:max_jx_100))
         ifile = 20
         outfile(1:ifile) = 'Compound_Elastic_Ang'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = particle(iproj)%label
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         j = istate
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''# Compound Elastic Angular Distribution data '')')
         write(100,'(''# Distribution normalized to unity'')')
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')   &
                  j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
         write(100,'(''#'')')
         write(100,'(''# Frame = COM'')')
         alf = 0.0d0
         bet = 0.0d0
         do in = 1, num_energies
            e_in = projectile%energy(in)
            write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                       e_in,Inelastic_cs(istate,in),cs_units
            if(Inelastic_cs(istate,in) < cs_threshold)cycle
               write(100,'(''#         E_in            cos(theta)            Prob'')')
            write(100,'(''#'',3(''   ----------------''))')
            xnorm = 0.0d0
            do jx = 0, max_jx_100
               x = real(jx,kind=8)*delta_jx_100 - 1.0d0
               if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.999845
               alf = 0.0d0
               bet = 0.0d0
               sum = 0.0d0
               do L = 0, Inelastic_L_max(istate,in)
                  sum = sum + Inelastic_Ang_L(L,istate,in)*poly(L,1,alf,bet,x)
               end do
               Ang_Dist(jx) = sum
               if(jx >= 1)xnorm = (Ang_Dist(jx-1) + Ang_Dist(jx))*delta_jx_100*0.5d0
            end do
            do jx = 0, max_jx_100
               x = real(jx,kind=8)*delta_jx_100 - 1.0d0
               Temp = Ang_Dist(jx)
               if(xnorm > 1.0d-7)Temp = Temp/xnorm
               write(100,'(1x,3(3x,1pe16.7))')e_in, x, Temp
            end do
         end do
         close(unit=100)
         if(allocated(Ang_Dist))deallocate(Ang_Dist)
!
!---   Also print out Legendre Coefficients
!
         ifile = 24
         outfile(1:ifile) = 'Compound_Elastic_Ang_Leg'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = particle(iproj)%label
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         j = istate
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
               nucleus(itarget)%state(istate)%energy
         if(iproj == 1)write(100,'(''# Compound Elastic Angular Distribution data '')')
         if(iproj == 1)write(100,'(''# Angular Distribution normalized to unity '')')
         write(100,'(''# Coefficients of Legendre expansion'')')
         write(100,'(''# Ang_dist(x) = Sum_L a(L)*P(L,x)'')')
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
         write(100,'(''#'')')
         write(100,'(''# Frame = COM'')')

         do in = 1, num_energies
            e_in = projectile%energy(in)
            write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                       e_in,Inelastic_cs(istate,in),cs_units
            if(Inelastic_cs(istate,in) < cs_threshold)cycle
            if(iproj == 1)then
               write(100,'(''#         E_in            L            a(L)'')')
            elseif(iproj > 1)then
               write(100,'(''#         E_in            L         a(L)('',a2,''/sr)'')')cs_units
            end if
            write(100,'(''#   ----------------     ---     ----------------'')')
            xnorm = 2.0d0*Inelastic_Ang_L(0,istate,in)
            do L = 0, Inelastic_L_max(istate,in)
               if(xnorm > 1.0d-8)Temp = Inelastic_Ang_L(L,istate,in)/xnorm
               if(iproj > 1)Temp = Temp*Inelastic_cs(istate,in)
               write(100,'(1x,3x,1pe16.7,3x,i5,3x,1pe16.7)')e_in, L, Temp
            end do 
         end do
         close(unit=100)

!
!----   Angular Distribution for Elastic
!
         if(.not.allocated(Ang_Dist))allocate(Ang_Dist(0:max_jx_100))
         ifile = 11
         outfile(1:ifile) = 'Elastic_Ang'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = particle(iproj)%label
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         j = istate
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                        &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''# Elastic Angular Distribution data '')')
         if(iproj ==1)write(100,'(''# Distribution normalized to unity'')')
         if(iproj > 1)write(100,'(''# Distribution given as Elastic/Rutherford'')')
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
         write(100,'(''#'')')
         write(100,'(''# Frame = COM'')')

         xA = (nucleus(itarget)%mass+nucleus(itarget)%state(istate)%energy)/particle(iproj)%mass
         Coul = (particle(iproj)%Z*nucleus(itarget)%Z*                                           &
                (1.0d0+xA)/xA*fine_structure*hbar_c)**2*0.25d0
         if(.not.allocated(xvalue))allocate(xvalue(0:max_jx_100))
         do jx = 0, max_jx_100
            x = real(jx,kind=8)*delta_jx_100 - 1.0d0
            if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.99999d0
            xvalue(jx) = x
         end do

         Elastic_cs(1:num_energies) = 0.0d0
         sig_C = 0.0d0  
         do in = 1, num_energies
            e_in = projectile%energy(in)
!-------   Shape + Compound Elastic
            xnorm = 0.0d0
            do jx = 0, max_jx_100
               x = real(jx,kind=8)*delta_jx_100 - 1.0d0
               if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.99999d0
               if(iproj > 1)sig_C = Coul/e_in**2/(1.0d0 - x)**2*fmsq_eq_barn                ! Rutherford in b/sr
               alf = 0.0d0
               bet = 0.0d0
               shape = 0.0d0
               comp = 0.0d0
               do L = 0, Ang_L_max
                  P_L = poly(L,1,alf,bet,x)
                  shape = shape + SE_Ang(L,in)*P_L
                  comp = comp + Inelastic_Ang_L(L,istate,in)*P_L
               end do
               shape = shape
               comp = comp*Inelastic_cs(istate,in)
               if(iproj > 1)comp = comp/Sig_C
               Ang_Dist(jx) = shape + comp
               if(jx > 0)xnorm = xnorm + (Ang_Dist(jx-1) + Ang_Dist(jx))*delta_jx_100*0.5d0
            end do
            if(iproj == 1)then                   !  neutrons - cross section
               Elastic_Ang(0:Ang_L_max,in) = 0.0d0
               do L = 0, Ang_L_max
                  Temp = SE_Ang(L,in) + Inelastic_Ang_L(L,istate,in)*Inelastic_cs(istate,in)
                  Elastic_Ang(L,in) = Temp/(SE_cs(in) + Inelastic_cs(istate,in))
               end do
               Elastic_cs(in) = (SE_cs(in) + Inelastic_cs(istate,in))
            elseif(iproj> 1)then                 !  Charged particles - ratio to Rutherford
               call Legendre_expand(max_jx_100+1,xvalue(0),Ang_Dist(0),                          &
                                    Ang_L_max,Elastic_Ang(0,in))
               do L = 0, Ang_L_max
                  Elastic_Ang(L,in) = Elastic_Ang(L,in)*xnorm
               end do
               Elastic_cs(in) = xnorm*0.5d0
            end if                  
            if(iproj == 1)then
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                          e_in,(SE_cs(in)+Inelastic_cs(istate,in)),cs_units
            elseif(iproj > 1)then
               write(100,'(''# E_in = '',1pe16.7,3x,'' <EL/Rutherford> = '',1pe16.7)')           &
                          e_in,xnorm
            end if
            if(iproj == 1)then
               write(100,'(''#         E_in            cos(theta)            Prob'')')
            elseif(iproj > 1)then
               write(100,'(''#         E_in            cos(theta)        EL/Rutherford'')')
            end if
            write(100,'(''#'',3(''   ----------------''))')
            do jx = 0, max_jx_100
               x = real(jx,kind=8)*delta_jx_100 - 1.0d0
               if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.99999d0
               Temp = Ang_Dist(jx)
               if(iproj == 1 .and. xnorm > 1.0d-8)Temp = Temp/xnorm
               if(iproj > 1 )Temp = Temp
               write(100,'(1x,3(3x,1pe16.7))')e_in, x, Temp
            end do
         end do
         close(unit=100)
         if(allocated(Ang_Dist))deallocate(Ang_Dist)
         if(allocated(xvalue))deallocate(xvalue)
!----
!----   Moved this below angular distribution to collect ratio for charged particles
!----
         ifile = 10
         outfile(1:ifile) = 'Elastic_cs'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = particle(iproj)%label
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         j = istate
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
               nucleus(itarget)%state(istate)%energy
         write(100,'(''# Total Elastic cross section data '')')
         if(iproj == 1)write(100,'(''# Total Elastic cross section data '')')
         if(iproj > 1)write(100,'(''# <Total Elastic/Rutherford> '')')
         if(iproj > 1)write(100,'(''# Integrated over d(cos(theta)) '')')
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
               j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
         write(100,'(''#'')')
         if(iproj == 1)then
            write(100,'(''#          E_in              xs'',''('',a2,'')'')')cs_units
         elseif(iproj > 1)then
            write(100,'(''#          E_in         <EL/Rutherford>'')')
         end if
         write(100,'(''#'',2(''   ----------------''))')
         do in = 1, num_energies
            cs = Elastic_cs(in)
            if(iproj == 1) cs = cs*cs_scale
            e_in=projectile%energy(in)
            write(100,'(1x,4(3x,1pe16.7))')e_in, cs
         end do
         close(unit=100)




!
!---   Also print out Legendre Coefficients
!
!         if(iproj == 1)then 
         ifile = 15
         outfile(1:ifile) = 'Elastic_Ang_Leg'
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = particle(iproj)%label
!         file_lab2(ilab2:ilab2) = 'n'
         ilab2 = ilab2 + 1
         file_lab2(ilab2:ilab2) = ')'
         j = istate
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
               j-1,nucleus(itarget)%state(j)%spin, ch_par(ipi),                                    &
               nucleus(itarget)%state(j)%energy
         if(iproj == 1)write(100,'(''# Elastic Angular Distribution data '')')
         if(iproj == 1)write(100,'(''# Angular Distribution normalized to unity '')')
         if(iproj > 1)write(100,'(''# Elastic Scattering Angular Distribution data '')')
         if(iproj > 1)write(100,'(''# Calculated as ratio to Rutherford: OM/Rutherford '')')
         write(100,'(''# Coefficients of Legendre expansion'')')
         write(100,'(''# Ang_dist(x) = Sum_L a(L)*P(L,x)'')')
         write(100,'(''#'')')
         write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
         write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         if(iproj > 0 .and. iproj < 7)then
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
         else
            write(100,'(''# Population decay, no projectile'')')
         end if
         write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
         ipf = nint((nucleus(itarget)%state(j)%parity + 1)/2.0d0)
         write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
         write(100,'(''#'')')
         write(100,'(''# Frame = COM'')')

         if(iproj == 1)then
            write(100,'(''#         E_in            L            a(L)'')')
         elseif(iproj > 1)then
            write(100,'(''#         E_in            L            a(L)'')')
         end if
         write(100,'(''#   ----------------     ---     ----------------'')')
         do in = 1, num_energies
            e_in = projectile%energy(in)
            write(100,'(''#'')')
            do L = 0, Ang_L_max
               write(100,'(1x,3x,1pe16.7,3x,i5,3x,1pe16.7)')e_in, L, Elastic_Ang(L,in)
            end do 
         end do
         close(unit=100)
!         end if
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Inelastic to discrete states    ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

         if(.not.allocated(Ang_Dist))allocate(Ang_Dist(0:max_jx_50))
         do j = 1, nucleus(itarget)%num_discrete
            if(j == istate)cycle
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    cross section v. E_in           ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            ifile = 17
            outfile(1:ifile) = 'Inelastic_cs_L000'
            l = j - 1
            if(l < 10)then
               write(outfile(ifile:ifile),'(i1)')l
            elseif(l >= 10 .and. l < 100)then
               write(outfile(ifile-1:ifile),'(i2)')l
            elseif(l >= 100 .and. l < 1000)then
               write(outfile(ifile-2:ifile),'(i3)')l
            else
               write(6,*)'too many discrete states inelastic_cs ',l
               stop
            end if
!**********************************************************************************
!-----   Check if this state has appreciable cross section. If not, skip     -----*
!-----   Threshold is based on the number of samples, and the smallest       -----*
!-----   cross section observable with one sample, i.e., abs_cs/num_mc_samp  -----*
!-----   We use 10 times this as a cut off. Cross ection below this will be  -----*
!-----   unreliable due to limitations in Monte Carlo sampling               -----*
!**********************************************************************************
            print_cs = .false.
            do in = 1, num_energies
               thresh = 20.0d0*absorption_cs(in)/real(num_mc_samp,kind=8)
               if(Inelastic_cs(j,in) > thresh)print_cs = .true.
            end do
            if(.not. print_cs)cycle
!**********************************************************************************

            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = particle(iproj)%label
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = quote
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''# Inelastic cross section'')')
            write(100,'(''#'')')
            write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
            write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
            inuc = itarget
            write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
            ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
            write(100,'(''# Final state  = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
            write(100,'(''#'')')
            write(100,'(''#         E_in              xs'',''('',a2,'')'')')cs_units
            write(100,'(''#'',2(''   ----------------''))')
            do in = 1, num_energies
               e_in = projectile%energy(in)
               write(100,'(1x,4(3x,1pe16.7))')e_in, Inelastic_cs(j,in)*cs_scale
            end do
            close(unit=100)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Angular distributions           ------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ifile = 18
            outfile(1:ifile) = 'Inelastic_Ang_L000'
            l = j - 1
            if(l < 10)then
               write(outfile(ifile:ifile),'(i1)')l
            elseif(l >= 10 .and. l < 100)then
               write(outfile(ifile-1:ifile),'(i2)')l
            elseif(l >= 100 .and. l < 1000)then
               write(outfile(ifile-2:ifile),'(i3)')l
            else
               write(6,*)'too many discrete states inelastic_ANG ',l
               stop
            end if
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = particle(iproj)%label
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = quote
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''# Inelastic Angular Distribution data '')')
            write(100,'(''#'')')
            write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
            write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
            inuc = itarget
            write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
            ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
            write(100,'(''# Final state  = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
            write(100,'(''# Frame = COM'')')
            write(100,'(''#'')')
            do in = 1, num_energies
               e_in = projectile%energy(in)
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                       e_in,Inelastic_cs(j,in),cs_units
               if(Inelastic_cs(j,in) < cs_threshold)cycle
               write(100,'(''#         E_in            cos(theta)            Prob'')')
               write(100,'(''#'',3(''   ----------------''))')
               xnorm = 0.0d0
               do jx = 0, max_jx_50
                  x = real(jx,kind=8)*delta_jx_50 - 1.0d0
                  if(iproj > 1 .and. abs(x - 1.0d0) <= 1.0d-5) x = 0.99999d0
                  alf = 0.0d0
                  bet = 0.0d0
                  sum = 0.0d0
                  do L = 0, Inelastic_L_max(j,in)
                     sum = sum + Inelastic_Ang_L(L,j,in)*poly(L,1,alf,bet,x)
                  end do
                  Ang_Dist(jx) = sum
                  if(jx >= 1)xnorm = xnorm + (Ang_Dist(jx-1) + Ang_Dist(jx))*delta_jx_50*0.5d0
               end do
               do jx = 0, max_jx_50
                  x = real(jx,kind=8)*delta_jx_50 - 1.0d0
                  if(xnorm > 1.0d-8)Ang_Dist(jx) = Ang_Dist(jx)/xnorm
                  write(100,'(1x,3(3x,1pe16.7))')e_in, x, Ang_Dist(jx)
               end do
            end do
            close(unit=100)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Angular distributions Legendre Coefficients        -----+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ifile = 22
            outfile(1:ifile) = 'Inelastic_Ang_Leg_L000'
            l = j - 1
            if(l < 10)then
               write(outfile(ifile:ifile),'(i1)')l
            elseif(l >= 10 .and. l < 100)then
               write(outfile(ifile-1:ifile),'(i2)')l
            elseif(l >= 100 .and. l < 1000)then
               write(outfile(ifile-2:ifile),'(i3)')l
            else
               write(6,*)'too many discrete states Inelastic_ANG_Leg ',l
               stop
            end if
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = particle(iproj)%label
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = quote
            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''# Inelastic Angular Distribution data '')')
            write(100,'(''# Coefficients of Legendre expansion'')')
            write(100,'(''# Ang_dist(x) = Sum_L a(L)*P(L,x)'')')
            write(100,'(''#'')')
            write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
            write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
            inuc = itarget
            write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
            ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
            write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
            write(100,'(''# Frame = COM'')')
            write(100,'(''#'')')

            do in = 1, num_energies
               e_in = projectile%energy(in)
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                       e_in,Inelastic_cs(j,in),cs_units
               if(Inelastic_cs(j,in) < cs_threshold)cycle
               write(100,'(''#         E_in            L            a(L)'')')
               write(100,'(''#   ----------------     ---     ----------------'')')
!               write(100,'(''#'')')
               xnorm = 2.0d0*Inelastic_Ang_L(0,j,in)
               do L = 0, Inelastic_L_max(j,in)
                  if(xnorm > 1.0d-8)Temp = Inelastic_Ang_L(L,j,in)/xnorm
                  write(100,'(1x,3x,1pe16.7,3x,i5,3x,1pe16.7)')e_in, L, Temp
               end do 
            end do
            close(unit=100)
!------------------------------------------------------------------------------------
         end do
         if(allocated(Ang_Dist))deallocate(Ang_Dist)
      end if                    !!!!   Finished with printing out elastic & inelastic info

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Data for exit channels       ---------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      if(.not.allocated(Temp_Ang_Dist))allocate(Temp_Ang_Dist(0:max_jx_20,0:num_e,0:6))
      do i = 1, num_channels
         inuc = Exit_Channel(i)%Final_nucleus
         directory(1:ilib_dir) = lib_dir(1:ilib_dir)
         idir = ilib_dir + 1
         directory(idir:idir) = '/'
         ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
         directory(idir+1:idir+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         idir = idir + ilast
         idir = idir + 1
         directory(idir:idir) = '/'

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------    Channel reaction data   -------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         do n = 1, Exit_Channel(i)%num_cs
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------    Channel cross section data   --------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            j = Exit_Channel(i)%StateLabel(n)
            Q = Exit_Channel(i)%Q_value - nucleus(inuc)%state(j)%energy +                           &
                nucleus(itarget)%state(istate)%energy

            ifile = 8
            outfile(1:ifile) = 'Channel_'
            outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ifile = ifile + ilast
            outfile(ifile+1:ifile+8) = '_cs_L000'
            ifile = ifile + 8
            l = j - 1
            if(l < 10)then
               write(outfile(ifile:ifile),'(i1)')l
            elseif(l >= 10 .and. l < 100)then
               write(outfile(ifile-1:ifile),'(i2)')l
            elseif(l >= 100 .and. l < 1000)then
               write(outfile(ifile-2:ifile),'(i3)')l
            else
               write(6,*)'too many discrete states Channel_cs ',l
               stop
            end if
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
!            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''# To level #'',i3)')l
!-----------------------------------
            write(100,'(''#'')')
            write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
            write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            if(iproj > 0 .and. iproj < 7)then
               write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
            else
               write(100,'(''# Population decay, no projectile'')')
            end if
            inuc = Exit_Channel(i)%Final_nucleus
            write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
            ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
            write(100,'(''# Xs collected in Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
            write(100,'(''#'')')
            if(Exit_channel(i)%num_particles > 0)then
               if(explicit_channels)then
                  write(100,'(''# Mass of emitted particles '')')
                  do m = 1, Exit_channel(i)%num_particles
                     k = Exit_channel(i)%decay_particles(m)
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m, particle(k)%mass/mass_u
                  end do
               else
                  write(100,'(''# Mass of emitted particles '')')
                  m = 0
                  do k = 1, 6
                     do nn = 1, Exit_channel(i)%num_part(k)
                        m = m + 1
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m,particle(k)%mass/mass_u
                     end do
                  end do
               end if
            end if
            write(100,'(''#'')')
            write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
            write(100,'(''#'')')
            write(100,'(''# Cross section data '')')
            write(100,1902)cs_units,(particle(k)%label, k = 0, 6)
 1902       format('#         E_in              xs(',a2,')',5x,7(8x,'Mult(',a1,')',4x))
            write(100,'(''#'',9(''   ----------------''))')
            do in = 1, num_energies
               e_in=projectile%energy(in)
               cs = Exit_Channel(i)%Channel_cs(in,n)*reaction_cs(in)*cs_scale
               write(100,'(1x,9(3x,1pe16.7))')e_in, cs, (Exit_Channel(i)%part_mult(k,n,in),k=0,6)
            end do
            close(unit=100)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Angular distribution and spectrum  ---------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------   Output P(Eout,x|Ein)   ------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
            ifile = 8
            outfile(1:ifile) = 'Channel_'
            outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ifile = ifile + ilast
            outfile(ifile+1:ifile+20) = '_Eout_Ang_G_Ein_L000'
            ifile = ifile + 20
            l = j - 1
            if(l < 10)then
               write(outfile(ifile:ifile),'(i1)')l
            elseif(l >= 10 .and. l < 100)then
               write(outfile(ifile-1:ifile),'(i2)')l
            elseif(l >= 100 .and. l < 1000)then
               write(outfile(ifile-2:ifile),'(i3)')l
            else
               write(6,*)'too many discrete states Channel_Eout_ANG_G_Ein ',l
               stop
            end if
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
!            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                        &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''#'')')
            write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
            write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            if(iproj > 0 .and. iproj < 7)then
               write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
            else
               write(100,'(''# Population decay, no projectile'')')
            end if
            inuc = Exit_Channel(i)%Final_nucleus
            write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
            ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
            write(100,'(''# Xs collected in Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
            write(100,'(''#'')')

            if(Exit_channel(i)%num_particles > 0)then
               if(explicit_channels)then
                  write(100,'(''# Mass of emitted particles '')')
                  do m = 1, Exit_channel(i)%num_particles
                     k = Exit_channel(i)%decay_particles(m)
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m, particle(k)%mass/mass_u
                  end do
               else
                  write(100,'(''# Mass of emitted particles '')')
                  m = 0
                  do k = 1, 6
                     do nn = 1, Exit_channel(i)%num_part(k)
                        m = m + 1
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m,particle(k)%mass/mass_u
                     end do
                   end do
               end if
            end if
            write(100,'(''#'')')
            write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
            write(100,'(''#'')')
            write(100,'(''# Angular Distribution data  P(x,Eout|Ein) '')')
            write(100,'(''# Normalized to unity over integral Int(-1,1) int(0,E_max)  P(x,Eout|Ein)dx,dEout '')')
            write(100,'(''# dx = '',1pe16.7,'' dEout = '',1pe16.7,'' MeV'')')delta_jx_20, de_spec
            write(100,'(''# Frame = COM'')')
            do in = 1, num_energies
               e_in = projectile%energy(in)
               if(e_in < -Q)cycle
               cs = Exit_Channel(i)%Channel_cs(in,n)*reaction_cs(in)
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                       e_in,cs,cs_units
               if(cs < cs_threshold)cycle
               write(100,'(''#'')')
               write(100,1910)(particle(k)%label, k = 0,6)
 1910       format('#         E_in               E_out           cos(theta)   ',                    &
                   7(3x,(5x,'Prob(',a1,')',4x)))
               write(100,'(''#'',10(''   ----------------''))')
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------    Find max value of icc to print out
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Normalize the particle spectra and find max icc to print out
               Temp_Ang_Dist(0:max_jx_20,0:num_e,0:6) = 0.0d0
               do k = 0, 6
                  xnorm = 0.0d0
                  do icc = 0, num_e
                     xnorm = xnorm + Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)*de_spec
                  end do
                  if(xnorm > 1.0d-8)then
                     do icc = 0, num_e
                        Exit_Channel(i)%Spect(k,n,in)%E_spec(icc) =                                 &
                          Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)/xnorm
                     end do
                  else
                     do icc = 0, num_e
                        Exit_Channel(i)%Spect(k,n,in)%E_spec(icc) = 0.0d0
                     end do
                  end if
               end do
               icc_max = num_e
               do icc = num_e - 1, 0, -1
                  sum = 0.0d0
                  do k = 0, 6
                     sum = sum + Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)
                  end do
                  if(sum >= 1.0d-8)then
                     icc_max = icc + 1
                     exit
                  end if
               end do
               if(icc_max < 1)icc_max = 1
!------   Create angular distribution array - normalize over angle and multiply by prob 
               alf = 0.0d0
               bet = 0.0d0
               do k = 0, 6
                  do icc = 0, icc_max
                     do jx = 0, max_jx_20
                        x = real(jx,kind=8)*delta_jx_20 - 1.0d0
                        sum = 0.0d0
                        do L = 0, Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_Max(icc)
                           sum = sum + Exit_Channel(i)%Spect(k,n,in)%E_Ang_L(L,icc)*                &
                                       poly(L,1,alf,bet,x)
                        end do
                        Temp_Ang_Dist(jx,icc,k) = sum
                     end do
                     xnorm = 0.0d0
                     do jx = 0, max_jx_20 - 1
                        xnorm = xnorm + (Temp_Ang_Dist(jx,icc,k) + Temp_Ang_Dist(jx+1,icc,k))*      &
                                delta_jx_20*0.5d0
                     end do
                     if(xnorm > 0.0d0)then
                        do jx = 0, max_jx_20
                           Temp_Ang_Dist(jx,icc,k) = Temp_Ang_Dist(jx,icc,k)/xnorm*                 &
                                                     Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)
                        end do
                     end if
                  end do
                  xnorm = 0.0d0
                  temp = 0.0d0
                  do icc = 0, icc_max
                     sum = 0.0d0
                     do jx = 0, max_jx_20 - 1
                        sum = sum + (Temp_Ang_Dist(jx,icc,k) + Temp_Ang_Dist(jx+1,icc,k))*          &
                                     0.5d0*delta_jx_20
                     end do
                     if(Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)> 1.0d-8)                          &
                     xnorm = xnorm + sum*de_spec
                     temp = temp + Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)*de_spec
                  end do
               end do

!----------   Print out data 
               do icc = 0, icc_max
                  e_out = real(icc,kind=8)*de_spec
                  if(icc >0)e_out = e_out - de_spec/2.0d0
                  do jx = 0, max_jx_20
                     x = real(jx,kind=8)*delta_jx_20 - 1.0d0
                     write(100,'(1x,10(3x,1pe16.7))')e_in, e_out, x,                                &
                        (Temp_Ang_Dist(jx,icc,k),k = 0,6)
                  end do
               end do
             end do
            close(unit=100)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Angular distributions Legendre Coefficients        -----+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ifile = 8
            outfile(1:ifile) = 'Channel_'
            outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ifile = ifile + ilast
            outfile(ifile+1:ifile+24) = '_Eout_Ang_Leg_G_Ein_L000'
            ifile = ifile + 24
            l = j - 1
            if(l < 10)then
               write(outfile(ifile:ifile),'(i1)')l
            elseif(l >= 10 .and. l < 100)then
               write(outfile(ifile-1:ifile),'(i2)')l
            elseif(l >= 100 .and. l < 1000)then
               write(outfile(ifile-2:ifile),'(i3)')l
            else
               write(6,*)'too many discrete states Channel_Eout_ANG_Leg_G_Ein ',l
               stop
            end if
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
!            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                        &
                  nucleus(itarget)%state(istate)%energy
!            write(100,'(''# To level # '',i3)')l
            write(100,'(''#'')')
            write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
            write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            if(iproj > 0 .and. iproj < 7)then
               write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
            else
               write(100,'(''# Population decay, no projectile'')')
            end if
            inuc = Exit_Channel(i)%Final_nucleus
            write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
            ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
            write(100,'(''# Xs collected in Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
            write(100,'(''#'')')
            if(Exit_channel(i)%num_particles > 0)then
               if(explicit_channels)then
                  write(100,'(''# Mass of emitted particles '')')
                  do m = 1, Exit_channel(i)%num_particles
                     k = Exit_channel(i)%decay_particles(m)
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m, particle(k)%mass/mass_u
                  end do
               else
                  write(100,'(''# Mass of emitted particles '')')
                  m = 0
                  do k = 1, 6
                     do nn = 1, Exit_channel(i)%num_part(k)
                        m = m + 1
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m,particle(k)%mass/mass_u
                     end do
                   end do
               end if
            end if
            write(100,'(''#'')')
            write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
            write(100,'(''#'')')
            write(100,'(''# Angular Distribution data  P(x,Eout|Ein) '')')
            write(100,'(''# Normalized to unity over integral Int(-1,1) int(0,E_max)  P(x,Eout|Ein)dx,dEout '')')
            write(100,'(''# dx = '',1pe16.7,'' dEout = '',1pe16.7,'' MeV'')')delta_jx_20, de_spec
            write(100,'(''# Frame = COM'')')

            do in = 1, num_energies
               e_in = projectile%energy(in)
               icc_max = num_e
               cs = Exit_Channel(i)%Channel_cs(in,n)*reaction_cs(in)
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                       e_in,cs,cs_units
               if(cs < cs_threshold)cycle
               do icc = num_e - 1, 0, -1
                  sum = 0.0d0
                  do k = 0, 6
                     sum = sum + Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)
                  end do
                  if(sum >= 1.0d-8)then
                     icc_max = icc + 1
                     exit
                  end if
               end do
               if(icc_max < 1)icc_max = 1
               do icc = 0, icc_max
                  E_out = real(icc,kind=8)*de_spec
                  if(icc > 0)E_out = E_out - 0.5d0*de_spec
                  write(100,'(''#'')')
                  write(100,'(''#'',36x,''Prob(Eout)'',7(3x,1pe16.7))')                             &
                     (Exit_Channel(i)%Spect(k,n,in)%E_spec(icc),k=0,6)
                  write(100,'(''#'')')
            write(100,2110)(particle(k)%label, k = 0,6)
 2110       format('#         E_in               E_out           L',7(3x,(5x,'a(L)(',a1,')',4x)))
            write(100,'(''#   ----------------   ----------------     ---'',7(3x,''----------------''))')
                  L_max = 0
                  do k = 0, 6
                     if(Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_Max(icc) > L_max)                     &
                        L_max = Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_Max(icc)
                  end do
                  do L = 0, L_max
                    write(100,'(1x,3x,1pe16.7,3x,1pe16.7,3x,i5,7(3x,1pe16.7))')e_in, e_out, L,     &
                       (Exit_Channel(i)%Spect(k,n,in)%E_Ang_L(L,icc),k=0,6)
                  end do
               end do 
            end do

            close(unit=100)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------   Output P(x|Ein)   ------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(.not.allocated(Ang_dis))allocate(Ang_dis(0:max_jx_20,0:6))
            ifile = 8
            outfile(1:ifile) = 'Channel_'
            outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ifile = ifile + ilast
            outfile(ifile+1:ifile+15) = '_Ang_G_Ein_L000'
            ifile = ifile + 15
            l = j - 1
            if(l < 10)then
               write(outfile(ifile:ifile),'(i1)')l
            elseif(l >= 10 .and. l < 100)then
               write(outfile(ifile-1:ifile),'(i2)')l
            elseif(l >= 100 .and. l < 1000)then
               write(outfile(ifile-2:ifile),'(i3)')l
            else
               write(6,*)'too many discrete states Channel_ANG_G_Ein ',l
               stop
            end if
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
!            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
!            write(100,'(''# To level # '',i3)')l
            write(100,'(''#'')')
            write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
            write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            if(iproj > 0 .and. iproj < 7)then
               write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
            else
               write(100,'(''# Population decay, no projectile'')')
            end if
            inuc = Exit_Channel(i)%Final_nucleus
            write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
            ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
            write(100,'(''# Xs collected in Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
            write(100,'(''#'')')
            if(Exit_channel(i)%num_particles > 0)then
               if(explicit_channels)then
                  write(100,'(''# Mass of emitted particles '')')
                  do m = 1, Exit_channel(i)%num_particles
                     k = Exit_channel(i)%decay_particles(m)
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m, particle(k)%mass/mass_u
                  end do
               else
                  write(100,'(''# Mass of emitted particles '')')
                  m = 0
                  do k = 1, 6
                     do nn = 1, Exit_channel(i)%num_part(k)
                        m = m + 1
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m,particle(k)%mass/mass_u
                     end do
                   end do
               end if
            end if
            write(100,'(''#'')')
            write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
            write(100,'(''#'')')
            write(100,'(''# Angular Distribution data  P(x|Ein) '')')
            write(100,'(''# Normalized to unity over integral Int(-1,1) P(x|Ein)dx '')')
            write(100,'(''# dx = '',1pe16.7)')delta_jx_20
            write(100,'(''# Frame = COM'')')
            Ang_dis(0:max_jx_10,0:6) = 0.0d0
            do in = 1, num_energies
               e_in=projectile%energy(in)
               if(e_in < -Q)cycle
               cs = Exit_Channel(i)%Channel_cs(in,n)*reaction_cs(in)
               write(100,'(''#'')')
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                       e_in,cs,cs_units
               if(cs < cs_threshold)cycle
                write(100,1903)(particle(k)%label, k = 0,6)
 1903       format('#         E_in             cos(theta)   ',3x,7(3x,(5x,'Prob(',a1,')',4x)))
            write(100,'(''#'',9(''   ----------------''))')
!----    Calculate Angular distribution again and then integrate over energy
               alf = 0.0d0
               bet = 0.0d0
               do k = 0, 6
                  do jx = 0, max_jx_20
                     Ang_dis(jx,k) = 0.0d0
                     sum = 0.0d0
                     do icc = 0, num_e
                        x = real(jx,kind=8)*delta_jx_20 - 1.0d0
                        do L = 0, Exit_Channel(i)%Spect(k,n,in)%E_Ang_L_Max(icc)
                           sum = sum + Exit_Channel(i)%Spect(k,n,in)%E_Ang_L(L,icc)*                &
                                       poly(L,1,alf,bet,x)*                                         &
                                       Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)*de_spec
                        end do
                     end do
                     Ang_dis(jx,k) = sum
                  end do
                  xnorm = 0.0d0
                  do jx = 1, max_jx_20
                     xnorm = xnorm + (Ang_dis(jx-1,k)+Ang_dis(jx,k))*delta_jx_20*0.5d0
                  end do
                  if(xnorm > 0.0d0)then
                     do jx = 0, max_jx_20
                        Ang_dis(jx,k) = Ang_dis(jx,k)/xnorm
                     end do
                  end if
               end do
               do jx = 0, max_jx_20
                  x = real(jx,kind=8)*delta_jx_20 - 1.0d0
                  write(100,'(1x,9(3x,1pe16.7))')e_in, x, (Ang_dis(jx,k),k = 0,6)
               end do
            end do
            close(unit=100)
            if(allocated(Ang_dis))deallocate(Ang_dis)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------   Output P(Eout|Ein)   ------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ifile = 8
            outfile(1:ifile) = 'Channel_'
            outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ifile = ifile + ilast
            outfile(ifile+1:ifile+16) = '_Eout_G_Ein_L000'
            ifile = ifile + 16
            l = j - 1
            if(l < 10)then
               write(outfile(ifile:ifile),'(i1)')l
            elseif(l >= 10 .and. l < 100)then
               write(outfile(ifile-1:ifile),'(i2)')l
            elseif(l >= 100 .and. l < 1000)then
               write(outfile(ifile-2:ifile),'(i3)')l
            else
               write(6,*)'too many discrete states Channel_G_Ein ',l
               stop
            end if
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
!            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                        &
                  nucleus(itarget)%state(istate)%energy
!            write(100,'(''# To level # '',i3)')l
            write(100,'(''#'')')
            write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
            write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            if(iproj > 0 .and. iproj < 7)then
               write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
            else
               write(100,'(''# Population decay, no projectile'')')
            end if
            inuc = Exit_Channel(i)%Final_nucleus
            write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
            ipf = nint((nucleus(inuc)%state(j)%parity + 1)/2.0d0)
            write(100,'(''# Xs collected in Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  j-1,nucleus(inuc)%state(j)%spin, ch_par(ipf), nucleus(inuc)%state(j)%energy
            write(100,'(''#'')')
            if(Exit_channel(i)%num_particles > 0)then
               if(explicit_channels)then
                  write(100,'(''# Mass of emitted particles '')')
                  do m = 1, Exit_channel(i)%num_particles
                     k = Exit_channel(i)%decay_particles(m)
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m, particle(k)%mass/mass_u
                  end do
               else
                  write(100,'(''# Mass of emitted particles '')')
                  m = 0
                  do k = 1, 6
                     do nn = 1, Exit_channel(i)%num_part(k)
                        m = m + 1
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m,particle(k)%mass/mass_u
                     end do
                   end do
               end if
            end if
            write(100,'(''#'')')
            write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
            write(100,'(''#'')')
            write(100,'(''# Emission Spectrum data  P(Eout|Ein) '')')
            write(100,'(''# Normalized to unity over integral Int(0,Emax) P(Eout|Ein)dEout '')')
            write(100,'(''# dEout = '',1pe16.7)')de_spec
            write(100,'(''# Frame = COM'')')
!--rem            E_spec(0:num_e,0:6) = 0.0d0
!--rem           E_norm(0:6) = 0.0d0
            do in = 1, num_energies
               e_in=projectile%energy(in)
               if(e_in < -Q)cycle
               write(100,'(''#'')')
               cs = Exit_Channel(i)%Channel_cs(in,n)*reaction_cs(in)
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')  &
                       e_in,cs,cs_units
               if(cs < cs_threshold)cycle
            write(100,1904)(particle(k)%label, k = 0,6)
  1904      format('#         E_in               E_out      ',7(3x,(5x,'Prob(',a1,')',4x)))
            write(100,'(''#'',9(''   ----------------''))')
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------    Find max value of icc to print out
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----
               icc_max = num_e
               do icc = num_e - 1, 0, -1
                  sum = 0.0d0
                  do k = 0, 6
                     sum = sum + Exit_Channel(i)%Spect(k,n,in)%E_spec(icc)
                  end do
                  if(sum >= 1.0d-8)then
                     icc_max = icc + 1
                     exit
                  end if
               end do
               if(icc_max < 1)icc_max = 1
               do icc = 0, icc_max
                  e_out = real(icc,kind=8)*de_spec
                  if(icc > 0)e_out = e_out - de_spec/2
                  write(100,'(1x,9(3x,1pe16.7))')e_in, e_out,                                      &
                     (Exit_Channel(i)%Spect(k,n,in)%E_spec(icc),k = 0,6)
               end do
            end do
            close(unit=100)
         end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------    Channel fission cross section data   ------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         if(fission)then
            ifile = 8
            outfile(1:ifile) = 'Channel_'
            outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ifile = ifile + ilast
            outfile(ifile+1:ifile+11) = '_fission_cs'
            ifile = ifile + 11
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
!            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''# To level = N/A'')')
            write(100,'(''# Cross section data '')')
            write(100,'(''#         E_in              xs('',a2,'')'')')cs_units
            write(100,'(''#'',2(''   ----------------''))')
            do in = 1, num_energies
               e_in=projectile%energy(in)
               cs = Exit_Channel(i)%Channel_cs(in,0)*reaction_cs(in)*cs_scale
               write(100,'(1x,9(3x,1pe16.7))')e_in,cs
            end do
            close(unit=100)
         end if
      end do
      if(allocated(Temp_Ang_Dist))deallocate(Temp_Ang_Dist)

!
!-------------------------------------------------------------------------------------------
!-----   Print out gammas for each channel if track_gammas == .true. for each energy
!-------------------------------------------------------------------------------------------
!
      if(track_gammas)then
         do i = 1, num_channels
            inuc = Exit_Channel(i)%Final_nucleus
            directory(1:ilib_dir) = lib_dir(1:ilib_dir)
            idir = ilib_dir + 1
            directory(idir:idir) = '/'
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            directory(idir+1:idir+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            idir = idir + ilast
            idir = idir + 1
            directory(idir:idir) = '/'
            Q = Exit_Channel(i)%Q_value
            Q = -Q
            ifile = 8
            outfile(1:ifile) = 'Channel_'
            outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ifile = ifile + ilast
            outfile(ifile+1:ifile+13) = '_Gammas_G_Ein'
            ifile = ifile + 13
            open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
!            ilab2 = ilab2 + 1
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')  &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                         &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''# Mass amu = '',1pe15.7,'' MeV'')')mass_u
            write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
            if(iproj < 7)then
               write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')particle(iproj)%mass/mass_u
            else
               write(100,'(''# Population decay, no projectile'')')
            end if
            write(100,'(''# Mass of residual = '',1pe23.16,'' amu'')')nucleus(inuc)%mass/mass_u
            write(100,'(''#'')')
            if(Exit_channel(i)%num_particles > 0)then
               if(explicit_channels)then
                  write(100,'(''# Mass of emitted particles '')')
                  do m = 1, Exit_channel(i)%num_particles
                     k = Exit_channel(i)%decay_particles(m)
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m, particle(k)%mass/mass_u
                  end do
               else
                  write(100,'(''# Mass of emitted particles '')')
                  m = 0
                  do k = 1, 6
                     do nn = 1, Exit_channel(i)%num_part(k)
                        m = m + 1
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')m,particle(k)%mass/mass_u
                     end do
                   end do
               end if
            end if
            write(100,'(''#'')')
            write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
            write(100,'(''#'')')
            if(.not. pop_calc)then
               write(100,'(''# Gamma-ray probabilities for each incident energy '')')
            else
               write(100,'(''# Gamma-ray probabilities for each initial excitation energy '')')
            end if


            if(.not. pop_calc)then
               write(100,'(''# Cross section for discrete gamma transitions in '',a5)')resid_label(1:nres)
               write(100,'(''#'')')
            else
               write(100,'(''# Fraction of decay for discrete gamma transitions in '',a5)')resid_label(1:nres)
               write(100,'(''#'')')
            end if

            do k = 1, nucleus(inuc)%num_discrete
               do m = 1, nucleus(inuc)%state(k)%nbranch
                  nf  = nucleus(inuc)%state(k)%ibranch(m)
                  ipi = nint((nucleus(inuc)%state(k)%parity+1.0d0)/2.0d0)
                  ipf = nint((nucleus(inuc)%state(nf)%parity+1.0d0)/2.0d0)

                  write(100,'(''# Transitions involving gamma ray #'',i4,'' from state #'',i4)')m,k

                  write(100,'(''#   i    Jp     f    Jp     e-gamma (MeV)'')')
                  write(100,'(''#  --- ----    --- ----     ------------- '')')
                  write(100,1911)k-1,nucleus(inuc)%state(k)%spin,ch_par(ipi),                      &
                                    nf-1,nucleus(inuc)%state(nf)%spin,ch_par(ipf),                 &
                                    nucleus(inuc)%state(k)%egamma(m)
1911        format('#',i4,1x,f4.1,a1,2x,i4,1x,f4.1,a1,8x,f10.5) 

                  if(.not. pop_calc)then
                     write(100,1906)cs_units
                  else
                     write(100,1920)
                  end if
1906        format('#        E_in    i    Jp     f    Jp     e-gamma (MeV)        sigma (',a2,    &
                   ')     sigma(g)/channel')
1920        format('#        E_in    i    Jp     f    Jp     e-gamma (MeV)        Fraction ',     &
                   '      Fraction/channel')
            write(100,1907)
1907        format('#  ----------   --- ----    --- ----     -------------     -------------     -------------')

                  do in = 1, num_energies
                     e_in = projectile%energy(in)
                     cs = Exit_channel(i)%state(k)%cs(m,in)*reaction_cs(in)*cs_scale
                     ratio = 0.0d0
                     sum = 0.0d0
                     do nn = 1, Exit_Channel(i)%num_cs
                        sum = sum + Exit_channel(i)%Channel_cs(in,nn)
                     end do

                     if(sum > 1.0d-10)ratio = cs/sum
1908   format(1x,f12.6,1x,i4,1x,f4.1,a1,'  ',i4,1x,f4.1,a1,8x,f10.5,3x,e15.7,3x,e15.7,3x,e15.7) 
                     write(100,1908)e_in,k-1,nucleus(inuc)%state(k)%spin,ch_par(ipi),              &
                                    nf-1,nucleus(inuc)%state(nf)%spin,ch_par(ipf),                 &
                                    nucleus(inuc)%state(k)%egamma(m),                              &
                                    cs,ratio
                  end do
               end do
            end do
            close(unit = 100)
         end do                               
      end if


      close(unit = 13)

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------    Finished Library output    -----------------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!-------  Arrays to sum emission spectrum and angular distributions  -+
!--rem     if(allocated(E_spec))deallocate(E_spec)
!--rem     if(allocated(E_norm))deallocate(E_norm)
!--rem     if(allocated(Ang_norm))deallocate(Ang_norm)
!---------------------------------------------------------------------+

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
!*****************************************************************************80
!
   use nuclei
   implicit none
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
!*****************************************************************************80
!
   implicit none
   real(kind=8) xj
   jhat=2.0d0*xj+1.0d0
   return
end function jhat
!
!*****************************************************************************80
!
real(kind=8) function interpolate(itype, x_in, num, grid, vec)
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
!*****************************************************************************80
!
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
!------------------------------------------------------
   real(kind=8) :: det3_3
!------------------------------------------------------

  i0 = 0
  if(x_in < grid(1))stop 'E_in smaller than input grid'
  if(x_in > grid(num))then
     write(6,*)num
     write(6,*)x_in,grid(num)
     stop 'E_in larger than input grid'
  end if
  do i = 1, num - 1                                  !  start from the bottom of the grid
     if(x_in == grid(i))then                  !  it is exactly on a grid point
        interpolate = vec(i)
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
  interpolate = 0.0d0
  if(itype == 0)then
     interpolate = y
  elseif(itype == 1)then
     interpolate = exp(y)
  end if
end function interpolate


real(kind=8) function det3_3(mat)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This functions computes the detemrinant of a 3x3 matrix
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
!*****************************************************************************80
!
   implicit none
   real(kind=8) :: mat(3,3)
   det3_3 = mat(1,1)*(mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)) -     &
            mat(1,2)*(mat(2,1)*mat(3,3) - mat(2,3)*mat(3,1)) +     &
            mat(1,3)*(mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1))
   return
end function det3_3


real(kind=8) function Gauss_var(iseed)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This function returns Gaussian distributed variables
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
!*****************************************************************************80
!
    use variable_kinds
    use constants
    implicit none
    integer(kind=int_64), intent(inout) :: iseed
    real(kind=8) :: u, v
!--------   External functions
    real(kind=8) :: random_64
!------------------------------------------------------------------------------
    u = random_64(iseed)
    v = random_64(iseed)
    if(random_64(iseed) < 0.5) then
       Gauss_var = sqrt(-2.0*log(u))*cos(two_pi*v)
    else
       Gauss_var = sqrt(-2.0*log(u))*sin(two_pi*v)
    end if
    return
end function Gauss_var


integer(kind=4) function find_channel(n_dat, dim_part, num_part, part_data, num_part_type)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This function returns the channel number given a set of decay particles
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
!*****************************************************************************80
!

   use options
   use Channel_info
   implicit none
   integer(kind=4), intent(in) :: n_dat, dim_part, num_part
   real(kind=8), intent(in) :: part_data(n_dat, dim_part)
   integer(kind=4), intent(in) :: num_part_type(6)
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

