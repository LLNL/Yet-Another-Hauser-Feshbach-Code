module variable_kinds
   integer(kind=4), parameter :: int_32 = selected_int_kind(9)
   integer(kind=4), parameter :: int_64 = selected_int_kind(18)
   integer(kind=4), parameter :: real_32 = selected_real_kind(6,37)
   integer(kind=4), parameter :: real_64 = selected_real_kind(15,307)
end module variable_kinds

!
!*******************************************************************************
!
module options
!
!*******************************************************************************
!
!  Discussion:
!
!------  Modules for controlling options
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
   character(len=132) version
   parameter (version = 'MC-3.12')
   integer(kind=int_64) :: iseed
   logical ran_setup
   integer(kind=4) PREEQ_Model
   logical analytic_preeq
   logical fission
   logical channels
   logical dump_events
   logical binary_event_file
   logical track_gammas, track_primary_gammas
   integer(kind=4) :: output_mode
   integer(kind=4) :: preeq_pair_model
   real(kind=8) :: preeq_delta
   integer(kind=4) :: preeq_fwell
   integer(kind=4) :: WF_model
   real(kind=8) :: trans_p_cut, trans_e_cut
   integer(kind=4) :: write_me
   real(kind=8) :: prob_cut
   real(kind=8) :: pop_check
   real(kind=8) sig_sum, sig_sumg, sig_sumn, sig_sep_e, sig_in
   real(kind=8) sigb(3)
   integer(kind=4) :: numn1, numn2, numn3
   integer(kind=4) :: part_lmax
   integer(kind=4) :: E1_model
   integer(kind=4)  e_l_max,m_l_max
   character(len=50) ex_pop_file
   logical pop_calc
   logical j_pop_calc
   logical fit_Gamma_gamma
   logical All_gammas
   logical Out_gammas_vs_E
   integer(kind=4) :: num_pop_e,num_pop
   real(kind=8) :: rho_cut
   type Pop_type
      real(kind=8) :: Ex_pop
      real(kind=8) :: dEx_pop
      integer(kind=4) :: num_pop
      real(kind=8), allocatable, dimension(:) :: j_pop
      real(kind=8), allocatable, dimension(:) :: par_pop
      real(kind=8), allocatable, dimension(:) :: bin_pop
   end type Pop_type

   type(Pop_type), allocatable, dimension(:) :: Pop_data

   integer(kind=4) :: pair_model
   real(kind=8) :: Preeq_V, Preeq_V1, Preeq_K
!   real(kind=8) :: Preeq_V, Preeq_V_part(2)
   real(kind=8) :: Preeq_g_div
   logical Preeq_g_a
!------   Scale factor for cross sections - b vs mb
   real(kind=8) :: cs_scale
   character(len = 2) :: cs_units
   
!------   Renormalize Transmission coefficients
   real(kind=8) :: T_norm

   logical quasi_elastic
   real(kind=8) :: quasi_e

   integer(kind=4) :: num_dec(0:6)
   logical :: ex_set
   integer(kind=4) :: Max_J_allowed
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Optical Model
   character(len=6)  optical
!   character(len=20) :: OM_pot
!   logical :: OM_pot_set
!   integer(kind=4) :: OM_option
   integer(kind=4) :: ifresco_shape
   character(len=132) :: local_cc_file
   logical exist_cc_file
   real(kind=8) :: cc_scale
   logical :: scale_elastic
   logical :: do_dwba

   real(kind=8) :: elastic_scale, elastic_shift, elastic_damp
   real(kind=8) :: Fiss_Max_J

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Monte Carlo sampling 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer(kind=4) :: nsamp, num_mc_samp
   logical :: biased_sampling

   logical :: trans_avg_l

   integer(kind=4) :: tries, total_tries, passes, num_force_decay
   integer(kind=4) :: re_baseline

   real(kind=8) :: Init_Kinetic_Energy, dInit_Kinetic_Energy


end module options
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
module useful_data
!
!*******************************************************************************
!
!  Discussion:
!
!------  Modules containing useful data across the code system
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
   character(len=132) out_file
   logical compound_setup
   real(kind=8) :: de
   real(kind=8) :: de_spec
   character(len=2) mass_file
   integer(kind=4) :: lev_option
   integer(kind=4) :: max_particle(0:6)         !  maximum number of each particle to be tracked
   logical primary_setup
   real(kind=8), allocatable, dimension (:,:) :: clb_l
   integer(kind=4) :: i_bind

   real(kind=8) :: test_spectrum(0:10000)

!.. names of stable isotopes, H-Es
   character(len=2) symb(99)                 !   moved to module useful_data
!.. mixed case
   data symb/    &
          'H ','He','Li','Be',' B',' C',' N',' O',' F','Ne','Na', &
          'Mg','Al','Si','P ',' S','Cl','Ar',' K','Ca','Sc','Ti', &
          ' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As', &
          'Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru', &
          'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs', &
          'Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy', &
          'Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir', &
          'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra', &
          'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es'/
end module useful_data
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
module Channel_info
!
!*******************************************************************************
!
!  Discussion:
!
!------  Modules containing information of decay channels
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

   type first_decay
      integer(kind=4) :: num_decay
      real(kind=8), allocatable, dimension (:) :: decay_prob
      integer(kind=4), allocatable, dimension (:) :: decay_list
   end type first_decay

   type channel_data
      real(kind=8) :: Channel_HFden
      integer(kind=4) :: num_decay
      integer(kind=4), allocatable, dimension (:) :: decay_to
      integer(kind=4), allocatable, dimension (:) :: decay_particle
      real(kind=8), allocatable, dimension (:) :: Channel_prob
      real(kind=8), allocatable, dimension (:) :: Channel_trans
      type(first_decay), allocatable, dimension (:) :: Channel_decay
   end type channel_data
      
   type(channel_data), allocatable, dimension(:,:,:) :: Channel


   integer(kind=4) :: num_channels
   logical :: Apply_Coulomb_Barrier
   real(kind=8) :: Coulomb_Barrier(6)

   type discrete_state
      integer(kind=4) :: nbranch
      real(kind=8), allocatable, dimension (:,:) :: cs
   end type discrete_state

   type E_A_Spectrum
       real(kind=8), allocatable, dimension(:) :: E_spec
       integer(kind=4), allocatable, dimension(:) :: E_count
       real(kind=8), allocatable, dimension(:) :: Ang_Dist
       real(kind=8), allocatable, dimension(:) :: Ang_L
       integer(kind=4) :: Ang_L_Max
       real(kind=8), allocatable, dimension(:,:) :: E_Ang_Dist
       real(kind=8), allocatable, dimension(:,:) :: E_Ang_L
       integer(kind=4), allocatable, dimension(:) :: E_Ang_L_Max
   end type E_A_Spectrum

   type ExitChannel
      real(kind=8) :: Q_value
      integer(kind=4) :: num_particles
      integer(kind=4), allocatable, dimension(:) :: decay_particles
      integer(kind=4) :: Channel_code
      character(len=20) :: Channel_Label
      integer(kind=4) :: Final_nucleus
      integer(kind=4) :: num_cs
      integer(kind=4) :: num_event
      real(kind=8), allocatable, dimension(:,:) :: Channel_cs
      integer(kind=4), allocatable, dimension(:) :: StateLabel
      type(discrete_state), allocatable, dimension(:) :: state
      real(kind=8), allocatable, dimension(:,:,:) :: part_mult
!      real(kind=8), allocatable, dimension(:,:,:,:,:) :: Spectrum
      real(kind=8), allocatable, dimension(:,:,:,:,:) :: Ang_L
      real(kind=8), allocatable, dimension(:,:,:,:) :: Max_L
      type(E_A_Spectrum), allocatable, dimension(:,:,:) :: Spect
      
   end type ExitChannel

   type(ExitChannel), allocatable, dimension(:) :: Exit_Channel

   type cc_state
      integer(kind=4) :: istate
      real(kind=8) :: spin
      real(kind=8) :: K
      real(kind=8) :: energy
      real(kind=8) :: parity
      integer(kind=4) :: state_type
      integer(kind=4) :: Ix_min
      integer(kind=4) :: Ix_max
      real(kind=8), allocatable :: spin_prob(:)
      real(kind=8) :: E_min
      real(kind=8) :: Delta_E
   end type cc_state

   type OpticalData
      integer(kind=4) :: nume
      integer(kind=4) :: numcc
      integer(kind=4) :: Max_L
      integer(kind=4) :: ielastic
      real(kind=8), allocatable :: energy(:)
      real(kind=8), allocatable :: optical_cs(:,:)
      real(kind=8), allocatable :: optical_leg(:,:,:)
      type(cc_state), allocatable :: state(:)
   end type OpticalData

   type(OpticalData) :: OpticalCS

end module Channel_info
!
!*******************************************************************************
!
module nuclei
!
!*******************************************************************************
!
!  Discussion:
!
!------  Modules for Derived type for data associated with bins, including
!------  decay information
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
   type bin_decay
      integer(kind=4) :: num_decay
      real(kind=8), allocatable, dimension (:) :: decay_prob
      integer(kind=4), allocatable, dimension (:) :: decay_list
   end type bin_decay

   type bin_data
      real(kind=8) :: rho
      real(kind=8) :: HF_den
      integer(kind=4) :: num_decay
      integer(kind=4), allocatable, dimension (:) :: decay_to
      integer(kind=4), allocatable, dimension (:) :: decay_particle
      real(kind=8), allocatable, dimension (:) :: HF_prob
      real(kind=8), allocatable, dimension (:) :: HF_prob2
      real(kind=8), allocatable, dimension (:) :: HF_trans
      type(bin_decay), allocatable, dimension (:) :: nuke_decay
   end type bin_data

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!---------    derived type for discrete states
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   type discrete_state
      logical evaluated
      integer(kind=4) :: ensdf_index
      real(kind=8) :: energy                                          !  Excitation energy for discrete state
      real(kind=8) :: spin                                            !  angular momentum
      real(kind=8) :: parity                                          !  parity
      real(kind=8) :: t12                                             !  half life
      logical isomer                                          !  logical identifying it as an isomer or not
      integer(kind=4) :: exit_lab
      real(kind=8) :: pop
      integer(kind=4) :: nbranch                                      !  # of states it gamma decays to
      integer(kind=4), allocatable, dimension (:) :: ibranch            !  array identifying states it gamma decays to
      real(kind=8), allocatable, dimension (:) :: branch                !  array with branching ratios for discrete gamma decays
      real(kind=8), allocatable, dimension (:) :: p_gamma                !  array with branching ratios for discrete gamma decays
      real(kind=8), allocatable, dimension (:) :: egamma                !  array with gamma-ray energies for decay
      real(kind=8), allocatable, dimension (:) :: p_ic                   !  array with internal conversion coefficients
      real(kind=8), allocatable, dimension (:) :: cs                    !  array with gamma-ray energies for decay
   end type discrete_state
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Derived type for Fission-barrier properties
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   type Fission_Barrier
      real(kind=8) :: barrier
      real(kind=8) :: hbw
      real(kind=8) :: barrier_damp(3)
      real(kind=8) :: Max_J
      integer(kind=4) :: symmetry
      real(kind=8) :: beta_2
      real(kind=8) :: level_param(20)
!      integer(kind=4) :: vib_mode
      real(kind=8) :: vib_enh(3)
      real(kind=8) :: rot_enh(5)
      real(kind=8) :: ecut
      integer(kind=4) :: num_discrete
      integer(kind=4) :: ncut
      real(kind=8), allocatable, dimension(:) :: state_e
      real(kind=8), allocatable, dimension(:) :: state_j
      real(kind=8), allocatable, dimension(:) :: state_pi
   end type Fission_Barrier

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   derived type describing each nucleus exisiting in the calculation
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   type define_nuclei                     !  data for each compound nucleus in the chain
      integer(kind=4) :: Z                                            !  Number of protons
      integer(kind=4) :: A                                            !  Number of nucleons
      character(len=2) atomic_symbol                               !  Character for atomic symbol
      real(kind=8) :: BE                                              !  mass excess
      real(kind=8) :: ME                                              !  mass excess
      real(kind=8) :: Mass
      real(kind=8) :: sep_e(0:6)                                        !  separation energies for n,p,d,t,he3,alpha
      real(kind=8) :: Delta_exp
      real(kind=8) :: Delta_exp_p
      real(kind=8) :: Delta_exp_n
      integer(kind=4) :: lev_option
      integer(kind=4) :: pair_model
      real(kind=8) :: Ex_max                                          !  maximum excitation energy
      real(kind=8) :: Kinetic_energy                                  !  Lab kinetic energy of nucleus
      real(kind=8) :: dKinetic_energy                                  !  Lab kinetic energy of nucleus
      logical PREEQ
!-----  Info definging bin structure
      integer(kind=4) :: nbin                                         !  number of excitation energy bins
      real(kind=8) :: jshift                                          !  = 0.5 of odd A = 0 for even A
      integer(kind=4) :: j_max                                        !  number of angular momentum bins
      real(kind=8), allocatable, dimension (:) :: e_grid               !  Energies for the level density grid
      real(kind=8), allocatable, dimension (:,:,:) :: pop              !  Populations (E,J,pi)
      real(kind=8), allocatable, dimension (:,:) :: PREEQ_pop          !  Populations removed from compound due to pree-equilibirum for J,pi states
      real(kind=8), allocatable, dimension (:,:,:) :: rho              !  level densities (E,J,pi)
      real(kind=8), allocatable, dimension (:,:,:) :: HF_den           !  Hauser-Feshbach denominators for each bin
      real(kind=8), allocatable, dimension (:,:,:,:) :: HF_prob           !  Hauser-Feshbach decay prob for particle type
      type(bin_data), allocatable, dimension(:,:,:) :: bins

!------  Data structures holding cross section data
      integer(kind=4) :: num_cs
      integer(kind=4), allocatable, dimension (:) :: state_lab
      real(kind=8), allocatable, dimension (:,:) :: channel_cs
      real(kind=8), allocatable, dimension (:) :: PREEQ_cs
      real(kind=8), allocatable, dimension (:,:) :: PREEQ_part_cs
      real(kind=8), allocatable, dimension (:) :: fission_cs
      real(kind=8), allocatable, dimension (:) :: hang_cs

!-----  Discrete state data
      integer(kind=4) :: num_discrete                                 !  # of discrete states
      integer(kind=4) :: ncut                                 !  # of discrete states
      type(discrete_state), allocatable, dimension (:) :: state    !  variable describing properties of discrete states
      real(kind=8), allocatable, dimension (:) :: HF_den_state         !  HF denominator for discrete states

!-----  How this nucleus decays
      integer(kind=4) :: num_decay
      integer(kind=4), allocatable, dimension (:) :: decay_to
      integer(kind=4), allocatable, dimension (:) :: decay_to_CN
      integer(kind=4), allocatable, dimension (:) :: decay_particle
      integer(kind=4) :: nump(6)                                      !  number of emitted particles from first compound
      integer(kind=4) :: num_part_decay
      integer(kind=4), allocatable, dimension (:) :: channel
      integer(kind=4) :: nbin_part
!-----  Spectra of particles emitted from this compound nucleus
      real(kind=8), allocatable, dimension (:,:) :: part_spectrum
      real(kind=8), allocatable, dimension (:,:) :: Lab_part_spectrum
      real(kind=8), allocatable, dimension (:,:) :: PREEQ_part_spectrum


!-----  Electromagnetic strength functions
      integer(kind=4) :: num_res
      real(kind=8) :: er_E1(4)                       !   Giant-dipole data
      real(kind=8) :: gr_E1(4)
      real(kind=8) :: sr_E1(4)
      logical :: E1_default
      real(kind=8), allocatable, dimension (:) :: er_M  ! magnetic properties
      real(kind=8), allocatable, dimension (:) :: gr_M
      real(kind=8), allocatable, dimension (:) :: sr_M
      real(kind=8), allocatable, dimension (:) :: er_E  ! electric propoerties for L > 1
      real(kind=8), allocatable, dimension (:) :: gr_E
      real(kind=8), allocatable, dimension (:) :: sr_E
      integer(kind=4) :: nbin_em
      integer(kind=4) :: lmax_E
      integer(kind=4) :: lmax_M
      integer(kind=4) :: e1_model
      real(kind=8), allocatable, dimension (:,:) :: f_E
      real(kind=8), allocatable, dimension (:,:) :: f_M
      real(kind=8), allocatable, dimension (:,:) :: str_E
      real(kind=8), allocatable, dimension (:,:) :: str_M
!------   Fission information
      logical :: fission
      integer(kind=4) :: F_n_barr
      type(Fission_Barrier), allocatable, dimension (:) :: F_Barrier    !  variable describing properties Fission Barriers
      real(kind=8) :: Fiss_cs
!------  Level Density information
      real(kind=8) :: beta(2:6)
      real(kind=8) :: sig2
      real(kind=8) :: sig2_perp
      real(kind=8) :: sig2_ax
      real(kind=8) :: vib_enh(3)
      real(kind=8) :: rot_enh(5)
      real(kind=8) :: target_spin
      integer(kind=4) :: target_ipar
      real(kind=8) :: s1                                              !  spin #1 for D0 J_gs(z,a-1)-0.5
      real(kind=8) :: s2                                              !  spin #1 for D0 J_gs(z,a-1)+0.5
      real(kind=8) :: p1                                              !  spin #1 for D0 J_gs(z,a-1)-0.5
      real(kind=8) :: p2                                              !  spin #1 for D0 J_gs(z,a-1)+0.5
      real(kind=8) :: p3                                              !  spin #1 for D0 J_gs(z,a-1)-0.5
      integer(kind=4) ipar
      real(kind=8) :: D0
      real(kind=8) :: D0exp
      real(kind=8) :: dD0exp
      real(kind=8) :: D1
      real(kind=8) :: D1exp
      real(kind=8) :: dD1exp
      real(kind=8) :: Gamma_g
      real(kind=8) :: Gamma_g_exp
      real(kind=8) :: dGamma_g_exp
      real(kind=8) :: Gamma_g_1
      real(kind=8) :: Gamma_g_1_exp
      real(kind=8) :: dGamma_g_1_exp
      character(len=10) :: level_model
      real(kind=8) :: level_ecut
      real(kind=8) :: level_param(20)                                 !  array for level density parameters: (1-8)= aparam,spin_cut,del,shell,gamma,ematch,ecut,sg2cut
      real(kind=8) :: a_Sn
      real(kind=8) :: sig2_Sn
      logical fit_D0
      logical fit_aparam
      logical fit_ematch
      logical fission_read
   end type define_nuclei
!
!----------  Array for each nucleus in the calculation
!
   type(define_nuclei), allocatable, dimension (:) :: nucleus
!
!----------  Type defining target nucleus
!
   type target_nucleus
      logical specified
      integer(kind=4) :: Z
      integer(kind=4) :: A
      integer(kind=4) :: icomp
      integer(kind=4) :: istate
      real(kind=8), allocatable, dimension(:,:) :: pop_xjpi
      real(kind=8), allocatable, dimension(:,:,:,:) :: pop_channel
   end type target_nucleus
!---------  Target variable
   type(target_nucleus) target
!
!---------  Type defining projectile
!
   type projectile_nucleus
      logical specified
      integer(kind=4) :: Z
      integer(kind=4) :: A
      integer(kind=4) :: particle_type
      real(kind=8), allocatable, dimension (:) :: energy
      integer(kind=4) :: num_e
      real(kind=8) :: e_max
      real(kind=8) :: e_min
      real(kind=8) :: e_step
   end type projectile_nucleus
!
!---------  projectile variable
!
   type(projectile_nucleus) projectile
!
!---------   Arrays needed for pre-equilibrium models
!
   integer(kind=4) :: pp_max,pn_max
   real(kind=8), allocatable, dimension (:,:,:,:) :: dWk
   real(kind=8), allocatable, dimension (:,:,:) :: Wk
   real(kind=8), allocatable, dimension (:,:) :: W
   real(kind=8), allocatable, dimension (:,:) :: tau
   real(kind=8), allocatable, dimension (:,:) :: taup
   real(kind=8), allocatable, dimension (:,:) :: Pre_Prob
   real(kind=8), allocatable, dimension (:,:) :: Gam_p
   real(kind=8), allocatable, dimension (:,:) :: Gam_n
   real(kind=8), allocatable, dimension (:,:) :: Gam_pp
   real(kind=8), allocatable, dimension (:,:) :: Gam_pn
   real(kind=8), allocatable, dimension (:,:) :: Gam_0_pn
   real(kind=8), allocatable, dimension (:,:) :: Gam_0_np
   real(kind=8), allocatable, dimension (:,:) :: Lamb_p_p
   real(kind=8), allocatable, dimension (:,:) :: Lamb_p_n
   real(kind=8), allocatable, dimension (:,:) :: Lamb_0_pn
   real(kind=8), allocatable, dimension (:,:) :: Lamb_0_np
   integer(kind=4) :: p0(2)
end module nuclei
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
module particles_def
!
!*******************************************************************************
!
!  Discussion:
!
!------   Module holding data on particle properties
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
   type particles
      integer(kind=4) :: Z
      integer(kind=4) :: A
      real(kind=8) :: spin
      real(kind=8) :: ME
      real(kind=8) :: mass
      real(kind=8) :: par
      logical in_decay
      real(kind=8) :: min_e
      real(kind=8) :: max_e
      character(len=1) label
      character(len=8) name
      integer(kind=4) :: lmax
      integer(kind=4) :: nume
      integer(kind=4) :: opt_pot
      integer(kind=4) :: max_opt_pot
      integer(kind=4) :: om_option
      logical :: opt_pot_set
      real(kind=8), allocatable:: e_grid(:)
      real(kind=8), allocatable, dimension (:,:,:) :: trans_read
      integer(kind=4) :: nbin
      real(kind=8), allocatable, dimension (:,:,:) :: trans
      integer(kind=4) :: nbin_spectrum
      real(kind=8), allocatable, dimension (:) :: Lab_spectrum
      real(kind=8), allocatable, dimension (:) :: COM_spectrum
      real(kind=8), allocatable, dimension (:) :: PREEQ_spectrum
   end type particles
!----------------------------------------------------------------------
   type(particles), allocatable, dimension (:) :: particle
end module particles_def
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
module directory_structure
!
!*******************************************************************************
!
!  Discussion:
!
!------   Module holding data on the directory structure
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
   character(len=200) home_path
   integer(kind=4) len_home
   character(len=200) data_path      ! path where data files are kept
   integer(kind=4) len_path             ! length of data_path
end module directory_structure        
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
module constants
!
!*******************************************************************************
!
!  Discussion:
!
!------   Module holding useful constants
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
   real(kind=8) :: pi,two_pi
   real(kind=8) :: hbar
   real(kind=8) :: hbar_c
   real(kind=8) :: fmsq_eq_barn
   real(kind=8) :: barn_eq_fmsq
   real(kind=8) :: fmsq_eq_mbarn
   real(kind=8) :: e_sq
   real(kind=8) :: mass_u
   real(kind=8) :: fine_structure
   real(kind=8) :: d_zero,d_half,d_one,d_two,d_three,d_four
   real(kind=8) :: factorial(0:50)
   real(kind=8) :: lfactorial(0:101)
   real(kind=8) :: g_metric(0:3,0:3)
end module constants
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
module Gauss_integration
!
!*******************************************************************************
!
!  Discussion:
!
!------   Module holding data used for Gauss-quadrature integeration
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
   real(kind=8) :: alf, bet
   integer(kind=4) :: n_glag
   real(kind=8), allocatable :: x_glag(:), w_glag(:) 
   integer(kind=4) :: n_gleg
   real(kind=8), allocatable :: x_gleg(:), w_gleg(:) 
end module Gauss_integration
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
module lev_fit
!
!*******************************************************************************
!
!  Discussion:
!
!------   Module holding data used for fitting level densities
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
   real(kind=8) :: sep_e_l,s1_lev,s2_lev
   integer(kind=4) :: ipar_lev
   integer(kind=4) :: AA_lev
end module lev_fit
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
module nodeinfo
!
!*******************************************************************************
!
!  Discussion:
!
!------   Module holding data for node information when MPI is implemented
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
   integer(kind=4) :: iproc,nproc,icomm
   integer(kind=4) :: group
end module nodeinfo
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
module print_control
!
!*******************************************************************************
!
!  Discussion:
!
!  Module for MPI node information
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
   logical print_output
   logical print_spectra
end module print_control
!
!*******************************************************************************
!
module pre_equilibrium_no_1
!
!*******************************************************************************
!
!  Discussion:
!
!  Module for Pre-equilibrium model #1
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
   real(kind=8) :: M2_C1
   real(kind=8) :: M2_C2
   real(kind=8) :: M2_C3
   real(kind=8) :: M2_Rpp
   real(kind=8) :: M2_Rnn
   real(kind=8) :: M2_Rpn
   real(kind=8) :: M2_Rnp
end module pre_equilibrium_no_1
