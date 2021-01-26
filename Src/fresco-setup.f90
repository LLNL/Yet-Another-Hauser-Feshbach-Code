!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!subroutine make_fresco_tco(data_path, len_path, tco_file, len_tco,       &
!                           de, pindex, iproj, itarget, istate, Ang_L_max)
subroutine make_fresco_tco(data_path, len_path, tco_file, len_tco,       &
                           pindex, iproj, itarget, istate, Ang_L_max)
!
!*******************************************************************************
!
!  Discussion:
!
!    This routine sets up and runs the optical model program
!    FRESCO. Transmission coefficients, elastic scattering, and
!    coupled-channels data are collected and written to files for
!    processing during the Hauser-Feshbach calculation.
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
  use nodeinfo
  use particles_def
  use nuclei
  use Channel_info
  use options
  implicit none
  character(len=200), intent(in) :: data_path
  integer(kind=4), intent(in) :: len_path
  character(len=80), intent(in) :: tco_file
  integer(kind=4), intent(in) ::  len_tco
!  real(kind=8), intent(in) :: de
  integer(kind=4), intent(in) :: pindex, iproj, itarget, istate, Ang_L_max
!-----------------------------------------------------------------------
  integer(kind=4) iztarget,iatarget

  real(kind=8) :: target_spin, target_parity
  real(kind=8) :: beta_def(6)
  character(len=2) symb
  character(len=1) label
  character(len=8) partname
  real(kind=8) :: e_min, e_max
  integer(kind=4) :: lmax
  integer(kind=4) :: zpart, apart
  real(kind=8) :: mass
  real(kind=8) :: spin
  integer(kind=4) :: ifirst, ilast
!---------------------------------------------------------------------
  integer(kind=4)iaaa
  character(len=5) nuke
  character(len=6) fname
  character(len=100) fresco_name
  character(len=100) out_name
  character(len=100) err_name
  character(len=100) tco_name
  integer(kind=4) :: i, k, l, n
  integer(kind=4) :: istart, istop
  integer(kind=4) :: iend, iendf, iend_out, iend_tco, iend_err
  integer(kind=4) :: iA, iZ
  real(kind=8) :: A, Z
  real(kind=8) :: Ap, xA
!---------------------------------------------------------------------
  integer(kind=4) :: num_temp
  integer(kind=4) :: iztemp, iatemp
  integer(kind=4), allocatable :: cc_index(:)
  real(kind=8), allocatable :: cc_state_e(:), cc_state_j(:)
  integer(kind=4), allocatable :: cc_state_par(:), cc_state_type(:)
  integer(kind=4), allocatable :: cc_state_k(:), cc_state_kpp(:)
  real(kind=8), allocatable :: cc_state_str(:)
!-rem  integer(kind=4) :: ibin
  integer(kind=4) :: if_state
  real(kind=8) :: e_cut

  integer(kind=4) :: num_states

  real(kind=8) :: K_band, J_gs, K_state, J_state
  integer(kind=4) :: par_band
!-rem  integer(kind=4) :: Ix_min, Ix_max, 
  integer(kind=4) :: Ix
!-rem  real(kind=8) :: J_min, J_max
  real(kind=8) :: zzero
  real(kind=8) :: xk_factor

  real(kind=8) :: xKK, xJmin, xJmax, xJ, STR
  integer(kind=4) :: jndex, jjndex, KK, num_K, ncount, par_state

  character(len=132) cc_file
  integer(kind=4) :: ilen_cc_file

  logical remove

  real(kind=8) :: beta(6), def_l(6)
  real(kind=8) :: radius

  real(kind=8) :: ener, ebin_min
  real(kind=8) :: d3

  real(kind=8) :: V_pot(2,3)
  real(kind=8) :: R_pot(2,3)
  real(kind=8) :: A_pot(2,3)
  real(kind=8) :: RC
  integer(kind=4) :: ipot

  logical deformed
  logical cc_found


  real(kind=8) :: ac

  real(kind=8) :: mass_proj, mass_target

  character(len=15) char_energy
  character(len=8) namet
  character(len=132) command, line

  integer(kind=4)ii
  integer(kind=4) :: ipot_end
  integer(kind=4) :: isp_max

!--------------   coupled channels data   -------------------
  integer(kind=4) :: itar, ichan, num_th, num_ang
  real(kind=8) :: th_inc, th_min, th_max
!--------------   coupled channels cross sections, etc.
  integer(kind=4)ie
  real(kind=8) :: tco_data
  real(kind=8), allocatable :: energy(:)
  real(kind=8) :: e_lab
  real(kind=8), allocatable :: optical_cs(:,:)
  real(kind=8), allocatable :: optical_leg(:,:,:)

  real(kind=8), allocatable :: theta(:)
  real(kind=8), allocatable :: ang_dist(:,:,:)
!
!------   Warning, using diffent Legendre weights here than in main
!------   these are separate here, and names are changed to avoid
!------   potential conflict in module Gauss_integration
!
!  integer(kind=4) :: Ang_L_max
  integer(kind=4) :: n_gleg2
  real(kind=8), allocatable :: x_gleg2(:), w_gleg2(:)

  real(kind=8) :: ang_theta, ang_rad
  real(kind=8) :: sum, xnorm, value

  real(kind=8) :: cross

  logical read_file
  integer(kind=4)isp, it
  real(kind=8) diff
  logical flap
  
  character(len=2) rela

  integer(kind=4) :: kp, ncc, nccba, nex, jtmax
  integer(kind=4) :: num_cc, num_dwba, num_read, num_tot, num_cc_read
  integer(kind=4) :: iter
  real(kind=8) :: absend
  real(kind=8) :: emin, emax
  real(kind=8) :: x, y
  integer(kind=4) :: nume
  real(kind=8) :: xnume
  real(kind=8) :: factor

  integer(kind=4) :: ll
  integer(kind=4) :: minl, maxl, num_l_lines, num_e_lines

  real(kind=8) :: Coul, xmu
  real(kind=8) :: sig_C


  character(len=100) title

  real(kind=8) :: alf, bet

  character(len=2) :: opt_label

!--------------------------------------------------------------------
  integer(kind=4) :: numw
  integer(kind=4) :: startw(66), stopw(66)
!--------------------------------------------------------------------
!--------    External functions   -----------------------------------

!  real(kind=8) :: Legendre
  real(kind=8) :: poly
  real(kind=8) :: interp
  real(kind=8) :: clebr

!-----   set up templates for protons and neutrons so that we can run 
!-----   RunTemplate
!-----   neutrons first

!   open(unit=88, file = 'test.dat',status = 'unknown')

!   write(6,*)'Target state = ',istate

  zzero = 0.0d0
  rela = 'bg'
  kp = 1

  deformed = .false.

  par_band = -100

  if(iproc == 0)then
     write(6,*)'Transmission coefficient file not found'
     write(6,*)'Will create using default optical model parameters'
     write(6,*)'Creating with Fresco'
  end if
  iaaa=2
  nuke='     '
!------------------------------------------------------------
  iztarget = nucleus(itarget)%Z
  iatarget = nucleus(itarget)%A
  target_spin = nucleus(itarget)%state(istate)%spin
  target_parity = nucleus(itarget)%state(istate)%parity
  beta_def(1:6) = 0.0d0
  do l = 2, 6
     beta_def(l) = nucleus(itarget)%beta(l)
  end do
  symb = nucleus(itarget)%atomic_symbol
  label = particle(pindex)%label
  partname = particle(pindex)%name
!---------    e_min and e_max need to be in the COM frame
  e_min = particle(pindex)%min_e
  e_max = particle(pindex)%max_e
  zpart = particle(pindex)%Z
  apart = particle(pindex)%A
  mass = particle(pindex)%mass
  spin = particle(pindex)%spin
  lmax = particle(pindex)%lmax

  isp_max = nint(2.0d0*spin)

!----------------------------------------------------------
  iend = 0
  if(iatarget < 10)then
     write(nuke(iend+1:iend+1),'(i1)')iatarget
     iend = iend + 1
  elseif(iatarget >= 10 .and. iatarget < 100)then
     write(nuke(iend+1:iend+2),'(i2)')iatarget
     iend = iend + 2
  elseif(iatarget >= 100 .and. iatarget < 1000)then
     write(nuke(iend+1:iend+3),'(i3)')iatarget
     iend = iend + 3
  endif
  if(symb(1:1).ne.' ')then
     nuke(iend+1:iend+2) = symb
     iend = iend + 2
  else
     nuke(iend+1:iend+1) = symb(2:2)
     iend = iend + 1
  end if
!--------------------------------------------------------------
  y = 0.0
  fname = label//nuke(1:iend)
  iendf = iend+1

  beta(1:6) = 0.0d0
  def_l(1:6) = 0.0d0
  flap = .false.
  radius = 1.1957*real(iatarget,kind=8)**(1./3.)
  ac = real(iatarget,kind=8)**(1./3.)

  iA = iatarget
  iZ = iztarget
  A = real(iatarget,kind=8)
  Z = real(iztarget,kind=8)
  mass_target = nucleus(itarget)%mass + nucleus(itarget)%state(istate)%energy
  mass_target = mass_target/mass_u
  mass_proj = particle(pindex)%mass/mass_u

  xA = mass_target/mass_proj
  Coul = (particle(pindex)%Z*nucleus(itarget)%Z*                                                    & 
         (1.0d0 + xA)/xA*fine_structure*hbar_c)**2*0.25d0



  Ap = Apart
  if(pindex <= 2) Ap = 0.0d0

  ebin_min = nucleus(itarget)%e_grid(1) - 0.5d0*nucleus(itarget)%delta_e(1)

  nex = 0
  nccba = 0
  ncc = 0

  write(6,*)'Check ',particle(pindex)%do_dwba


  cc_found = .false.
!----   check for coupled channels ONLY in protons and neutrons
  if(pindex <= 2)then
     if(exist_cc_file)then
        ilen_cc_file = index(local_cc_file,' ')-1
        cc_file(1:ilen_cc_file) = local_cc_file(1:ilen_cc_file)
        write(6,*)'Using local coupled-channels file: ', cc_file(1:ilen_cc_file)
     else
        cc_file = data_path(1:len_path)//'Coupled-Channels.txt'
        ilen_cc_file = index(cc_file,' ')-1
        write(6,*)'Using system coupled-channels file: ', cc_file(1:ilen_cc_file)
     end if

     open(unit = 21,file = cc_file(1:ilen_cc_file),status = 'old')

     num_states = nucleus(itarget)%ncut
     if(All_gammas)num_states = nucleus(itarget)%num_discrete

  do while(.not. cc_found)      !   - loop is on iz and ia nuclei not on lines in file
     read(21,'(a)')line
     if(line(1:3) == 'END' .or. line(1:3) == 'End' .or. line(1:3) =='end')exit
     read(21,*)iztemp,iatemp
     if(iztemp == iztarget .and. iatemp == iatarget)then
        deformed = .true.
        cc_found = .true.
        read(21,*)beta(2), beta(4), beta(6)
        if(nint(beta(2)) == -1)then
           beta(2) = beta_def(2)
           beta(4) = beta_def(4)
           beta(6) = 1.0d-5
        else
           beta_def(2) = beta(2)
           beta_def(4) = beta(4)
           beta_def(6) = beta(6)
        end if
        if(abs(beta(6)) < 1.0d-5)beta(6) = 1.0d-5
        def_l(2) = radius*beta(2)
        def_l(4) = radius*beta(4)
        def_l(6) = radius*beta(6)
        read(21,*)num_cc_read, num_dwba
!----   Read in and compute total DWBA states if num_read > num_cc
!------    Don't do DWBA for particles different from projectile
        if(pindex /= iproj)then
           num_read = num_cc_read
           num_tot = num_cc_read
        end if
!------    Don't do DWBA if do_dwba = .false.
        if(.not. particle(pindex)%do_dwba)then
           num_read = num_cc_read
           num_tot = num_cc_read
        end if
        K_band = 200.0d0
        J_gs = real(nucleus(itarget)%state(istate)%spin,kind=8)
!----   First, get value for K for the ground-state, coupled-channels band
!----   It is the lowest spin in the band, not necessarily the bandhead
        num_tot = 0
        num_cc = 0
        do i = 1, num_cc_read
           read(21,'(a)')line
           call parse_string(line, numw, startw, stopw)
           read(line(startw(1):stopw(1)),*)jndex
!----   Make sure that the discrete state in the coupled channels calculation
!----   is in the calculation, the endsf index is stored in state(i)%ensdf_index
           do n = 1, num_states
              jjndex = nucleus(itarget)%state(n)%ensdf_index
              if(jndex == jjndex)then
                 if(nucleus(itarget)%state(n)%energy > nucleus(itarget)%level_param(7) &
                    .and. .not. All_gammas)exit
                 xJ = nucleus(itarget)%state(n)%spin
                 if(xJ < K_band)then
                    K_band = xJ
                    par_band = nint(nucleus(itarget)%state(n)%parity)
                 end if
                 num_cc = num_cc + 1
                 num_tot = num_tot + 1
                 exit
             end if
           end do
        end do
!----  Now we have K, read in DWBA states and work out couplings to get total
!----  number of states that will be in full FRESCO calculation
        do i = 1, num_dwba
           read(21,'(a)')line
           call parse_string(line, numw, startw, stopw)
           read(line(startw(1):stopw(1)),*)jndex
           read(line(startw(2):stopw(2)),*)ener
           read(line(startw(3):stopw(3)-1),*)xJ
           if(jndex > 0)then
              do n = 1, num_states
                 jjndex = nucleus(itarget)%state(n)%ensdf_index
                 if(jndex == jjndex)then
                    if(nucleus(itarget)%state(n)%energy > nucleus(itarget)%level_param(7) &
                       .and. .not. All_gammas)exit
                    num_tot = num_tot + 1
                    exit
                 end if
              end do
           else
!-------    Make sure that energy of DWBA state is above lowest continuous
!-------    energy bin
              if(ener >= ebin_min)then
                 if(xJ >= 0.0d0)then
                    num_tot = num_tot + 1
                 else
                    read(line(startw(4):stopw(4)),*)KK
                    xKK = real(KK,kind=8)
                    xJmin = max(abs(xKK - K_band),K_band)
                    xJmax = xKK + K_band
                    num_K = nint(xJmax - xJmin) + 1
                    num_tot = num_tot + num_K
                 end if
              end if
           end if
        end do
!------    Allocate important arrays
        if(.not.allocated(cc_index))allocate(cc_index(num_tot))
        if(.not.allocated(cc_state_j))allocate(cc_state_j(num_tot))
        if(.not.allocated(cc_state_e))allocate(cc_state_e(num_tot))
        if(.not.allocated(cc_state_par))allocate(cc_state_par(num_tot))
        if(.not.allocated(cc_state_type))allocate(cc_state_type(num_tot))
        if(.not.allocated(cc_state_k))allocate(cc_state_k(num_tot))
        if(.not.allocated(cc_state_kpp))allocate(cc_state_kpp(num_tot))
        if(.not.allocated(cc_state_str))allocate(cc_state_str(num_tot))
        cc_state_j(1) = nucleus(itarget)%state(1)%spin
        cc_state_e(1) = nucleus(itarget)%state(1)%energy
        cc_state_par(1) = nint(nucleus(itarget)%state(1)%parity)
        cc_state_type(1) = 1
!-------   Read in coupled-channels information
!-------   First rewind file and get back to the point for reading in data

        do i = 1, num_cc_read + num_dwba
           backspace(21)
        end do
        ncount = 0
        do i = 1, num_cc_read
           read(21,'(a)')line
           call parse_string(line, numw, startw, stopw)
           read(line(startw(1):stopw(1)),*)jndex
           if(jndex <= 0)then
              write(6,*)'Error in coupled channels, this state has an improper cc_index'
              write(6,'(a)')line(1:stopw(numw))
              stop
           end if
!-----  See if state in cc list is in calculation
!-----  Check against all discrete states for target - jndex is index in ensdf evaluated file
!-----  This is stored in %state(n)%ensdf_index
!-----  The loop n is over all descrete states
           do n = 1, num_states
              jjndex = nucleus(itarget)%state(n)%ensdf_index
              if(jndex == jjndex)then
                 if(nucleus(itarget)%state(n)%energy > nucleus(itarget)%level_param(7) &
                    .and. .not. All_gammas)exit
                 cc_index(i) = n
                 ii = cc_index(i)
                 if(ii == istate)then
                    if_state = i
                 end if
                 cc_state_j(i) = nucleus(itarget)%state(ii)%spin
                 cc_state_e(i) = nucleus(itarget)%state(ii)%energy
                 cc_state_par(i) = nint(nucleus(itarget)%state(ii)%parity)
                 cc_state_type(i) = 1
                 cc_state_kpp(i) = 1
                 cc_state_k(i) = 0
                 cc_state_str(i) = 0.0d0
                 ncount = ncount + 1
                 exit
              end if
           end do
       end do
!------    read in information for DWBA states
        do i = 1, num_dwba
           read(21,'(a)')line
           call parse_string(line, numw, startw, stopw)
           read(line(startw(1):stopw(1)),*)jndex
           if(numw < 5)then
              write(6,*)'Error reading in coupled channels'
              write(6,*)'Not enough entries for DWBA state'
              write(6,*)'This is entry #',i,' in list'
              write(6,'(a)')line(1:stopw(numw))
           end if
           if(jndex > 0)then
              do n = 1, num_states
                 if(jndex == nucleus(itarget)%state(n)%ensdf_index)then
                    if(nucleus(itarget)%state(n)%energy > nucleus(itarget)%level_param(7) &
                       .and. .not. All_gammas)exit
                    ncount = ncount + 1
                    cc_index(ncount) = n
                    ii = cc_index(ncount)
                    cc_state_j(ncount) = nucleus(itarget)%state(ii)%spin
                    cc_state_e(ncount) = nucleus(itarget)%state(ii)%energy
                    cc_state_par(ncount) = nint(nucleus(itarget)%state(ii)%parity)
                    cc_state_type(ncount) = 1
                    cc_state_kpp(ncount) = 2
                    read(line(startw(4):stopw(4)),*)cc_state_k(ncount)
                    read(line(startw(5):stopw(5)),*)cc_state_str(ncount)
                    exit
                 end if
              end do
           elseif(jndex == 0)then
              read(line(startw(2):stopw(2)),*)ener
              read(line(startw(3):stopw(3)-1),*)xJ
              read(line(startw(4):stopw(4)),*)KK
              read(line(startw(5):stopw(5)),*)STR
!-------    Make sure that energy of DWBA state is above lowest continuous
!-------    energy bin
              if(ener >= ebin_min)then
                 if(xJ >= 0.0d0)then
                    ncount = ncount + 1
                    cc_index(ncount) = jndex
                    cc_state_e(ncount) = ener
                    cc_state_j(ncount) = xJ
                    cc_state_par(ncount) = 1
                    if(line(stopw(3):stopw(3)) == '-')cc_state_par(ncount) = -1
                    cc_state_type(ncount) = 0
                    cc_state_k(ncount) = KK
                    cc_state_str(ncount) = STR
                    cc_state_kpp(ncount) = 2
                 else
                    xKK = real(KK,kind=8)
                    xJmin = max(abs(xKK - K_band),K_band)
                    xJmax = xKK + K_band
                    num_K = nint(xJmax - xJmin) + 1
                    xJmin = xJmin - 1.0d0
                    if(par_band == -100)stop 'par_band not set correctly'
                    par_state = par_band*(-1)**KK
                    do k = 1, num_K
                       ncount = ncount + 1
                       cc_index(ncount) = jndex
                       cc_state_e(ncount) = ener
                       cc_state_j(ncount) = xJmin + real(k,kind=8)
                       cc_state_par(ncount) = 1
                       cc_state_par(ncount) = par_state
                       cc_state_type(ncount) = 0
                       cc_state_k(ncount) = KK
                       cc_state_str(ncount) = STR
                       cc_state_kpp(ncount) = 2
                    end do
                  end if
              end if
           elseif(jndex == -1)then
              read(line(startw(2):stopw(2)),*)ener
              read(line(startw(3):stopw(3)-1),*)xJ
              read(line(startw(4):stopw(4)),*)KK
              read(line(startw(5):stopw(5)),*)STR
!-------    Make sure that energy of DWBA state is above lowest continuous
!-------    energy bin
              if(ener >= ebin_min)then
                 ncount = ncount + 1
                 cc_index(ncount) = jndex
                 cc_state_e(ncount) = ener
                 cc_state_j(ncount) = XJ
                 cc_state_par(ncount) = 1
                 if(line(stopw(3):stopw(3)) == '-')cc_state_par(ncount) = -1
                 cc_state_type(ncount) = -1
                 cc_state_k(ncount) = kk
                 cc_state_str(ncount) = STR
                 cc_state_kpp(ncount) = 2
              end if
           end if
        end do
        ncc = num_cc
        nex = num_tot
     else
!-------    Nucleus not found yet, so keep reading
        read(21,*)
        read(21,*)num_cc,num_temp
        do i = 1, num_cc + num_temp
           read(21,*)
        end do
     end if
     end do
     close(unit=21)
  end if
!------     Nucleus not found in coupled channels file, so perform spherical optical model
!------     set do_dwba = .false.
  if(.not. cc_found)then
     particle(pindex)%do_dwba = .false.
     num_cc = 1
     num_tot = num_cc
     nex = num_tot
     ncc = num_cc
     if(.not.allocated(cc_index))allocate(cc_index(num_tot))
     if(.not.allocated(cc_state_j))allocate(cc_state_j(num_tot))
     if(.not.allocated(cc_state_e))allocate(cc_state_e(num_tot))
     if(.not.allocated(cc_state_par))allocate(cc_state_par(num_tot))
     if(.not.allocated(cc_state_type))allocate(cc_state_type(num_tot))
     if(.not.allocated(cc_state_k))allocate(cc_state_k(num_tot))
     if(.not.allocated(cc_state_kpp))allocate(cc_state_kpp(num_tot))
     if(.not.allocated(cc_state_str))allocate(cc_state_str(num_tot))
     cc_index(1) = istate
     cc_state_j(1) = nucleus(itarget)%state(istate)%spin
     cc_state_e(1) = nucleus(itarget)%state(istate)%energy
     cc_state_par(1) = nint(nucleus(itarget)%state(1)%parity)
     cc_state_type(1) = 1
     cc_state_kpp(1) = 1
     cc_state_k(1) = 0
     cc_state_str(1) = 0.0d0
     K_band = nucleus(itarget)%state(istate)%spin
     J_gs = 0.0d0
  end if


!---------    Sort through coupled channels states, remove if above e_cut and not 
!---------    tracking all discretes. Also, remove if CCBA states are degenerate

  e_cut = nucleus(itarget)%level_param(7)
  n = 2
  do while(n <= nex)
     remove = .false.
     if(cc_state_type(n) == 1)then
        if(cc_state_e(n) > e_cut .and. .not. All_gammas)then    ! remove from the list, push others down
           do i = n + 1,nex
              cc_index(i-1) = cc_index(i)
              cc_state_e(i-1) = cc_state_e(i)
              cc_state_j(i-1) = cc_state_j(i)
              cc_state_par(i-1) = cc_state_par(i)
              cc_state_type(i-1) = cc_state_type(i)
              cc_state_kpp(i-1) = cc_state_kpp(i)
              cc_state_k(i-1) = cc_state_k(i)
              cc_state_str(i-1) = cc_state_str(i)
           end do
           cc_index(nex) = 0
           cc_state_e(nex) = 0.0d0
           cc_state_j(nex) = 0.0d0
           cc_state_par(nex) = 0
           cc_state_type(nex) = -1
           cc_state_kpp(nex) = 0
           cc_state_k(nex) = 0
           cc_state_str(nex) = 0.0d0
           ncc = ncc - 1
           nex = nex - 1
           remove = .true.
        end if
     end if
     if(.not. remove) n = n + 1
  end do

!---   We have all the states read in and set up to do DWBA, but check if DWBA is
!---   wanted. If not, set nex = ncc
  if(.not. particle(pindex)%do_dwba)nex = ncc


!---    Find K of the coupled-channels band, smallest J 

!---------    Coupled channels control
!---------    ncc = total number of coupled channles
!---------    nex = total number of excited states in calculaiton including CCBA states
!---------    if ncc /= nex run one iteration after coupled channels for DWBA states.
  iter = 0
  if(ncc /= nex)iter = 1

!
!-------   Set up energy grid. multiplication, to make log-log grid
!

  nume = 200
  xnume = real(nume,kind=8)
  emin = e_min
  emax = 1.5d0*e_max
  factor = (log(emax)-log(emin))/xnume
  factor = exp(factor)
!  if(pindex == 2)then
!     emin = 0.05d0
!  elseif(pindex == 3)then
!     emin = 0.020d0
!  elseif(pindex == 4)then
!     emin = 0.020d0
!  elseif(pindex == 5)then
!     emin = 0.100d0
!  elseif(pindex == 6)then
!     emin = 0.100d0
!  end if             
  factor = (log(emax)-log(emin))/xnume
  factor = exp(factor)

  ener = emin/factor
  nume = 0

!  write(6,*)'e_max = ',e_max

  emax = max(e_max + 5.0d0, 25.0d0)

!  write(6,*)emin,emax,factor

  do while(ener < emax)
     nume = nume + 1
     ener = ener*factor
  end do


  particle(pindex)%nume = nume
  if(.not.allocated(particle(pindex)%e_grid))allocate(particle(pindex)%e_grid(nume))
  particle(pindex)%e_grid(1:nume) = 0.0d0
  if(.not.allocated(particle(pindex)%trans_read))allocate(particle(pindex)%trans_read(nume,0:isp_max,0:lmax))
  particle(pindex)%trans_read(1:nume,0:isp_max,0:lmax) = 0.0d0
!------------   temporary coupled channels cross section data 
  allocate(energy(nume))
  energy(1:nume) = 0.0d0
  ener = emin/factor
  do ie = 1, nume
     ener = ener*factor
!     energy(ie) = ener
     particle(pindex)%e_grid(ie) = ener*mass_target/(mass_target + mass_proj)                  !   COM frame
!     e_rel = ener*mass_target/(mass_target + mass_proj)
!     e_lab = ener*(mass_target + mass_proj)/mass_target
     e_lab = ener
     energy(ie) = e_lab                                  !   Lab frame
     write(6,*)ie, energy(ie), particle(pindex)%e_grid(ie)
  end do

  jtmax=20
  absend=0.001d0

  num_ang = 181
  allocate(theta(num_ang))
  theta(1:num_ang) = 0.0d0
  allocate(ang_dist(num_ang, nex, nume))
  ang_dist(1:num_ang,1:nex,1:nume) = 0.0d0
!-------    Gauss-Legendre integration
!  Ang_L_max = 30
  n_gleg2 = Ang_L_max + 1
  allocate(x_gleg2(n_gleg2), w_gleg2(n_gleg2))
  alf = 0.0d0
  bet = 0.0d0
  call gauss_quad(n_gleg2,1,alf,bet,x_gleg2,w_gleg2)
!------------------------------------------------------------------------------------------
!--------   Direct, elastic, etc is used ONLY for the projectile    -----------------------
!------------------------------------------------------------------------------------------
  if(.not.allocated(optical_cs))allocate(optical_cs(nume, nex))                   !  Ecis cross sections
   optical_cs(1:nume, 1:nex) = 0.0d0
  if(pindex == iproj .and. .not.allocated(optical_leg))then
     allocate(optical_leg(nume, 0:Ang_L_max, nex)) ! Legendre coeeficients
     optical_leg(1:nume, 0:Ang_L_max, 1:nex) = 0.0d0
  end if

  th_min = 0.0d0
  th_max = 180.0d0
  th_inc = 1.0d0
!  do i = 1, num_ang
!     theta(i) = th_min + real(i-1,kind=8)*th_inc
!  end do

!  ipot_end = index(OM_pot,' ')-1

  namet(1:8) = ' '
  namet(9-iend:8) = nuke(1:iend)

  kp = 1                            !  Potential number, here only one potential is used
  ie = 0

  out_name(1:iendf) = tco_file(1:len_tco)
  iend_out = len_tco + 11
  out_name(iendf+1:iend_out) = '_fresco.out'

  err_name(1:iendf) = tco_file(1:len_tco)
  iend_err = len_tco + 15
  out_name(iendf+1:iend_err) = '_fresco_err.out'

  write(6,*)'particle = ',pindex
  write(6,*)'Optical potential # = ',particle(pindex)%opt_pot


  do ie = 1, nume
     ener = energy(ie)
     char_energy(1:15) = ' '
     write(char_energy,'(e15.7)')ener
     if(pindex == 1)then
        RC = 1.2d0
        if(particle(pindex)%opt_pot == 1)                         &
           call KD_potential(pindex, iA, iZ, ener, V_pot, R_pot, a_pot, RC, D3)
        if(particle(pindex)%opt_pot == 2)                         &
           call soukhovitskii_potential(pindex, iA, iZ, particle(pindex)%om_option,  &
                                        ener, V_pot, R_pot, a_pot, RC)
        if(particle(pindex)%opt_pot == 3)                         &
           call maslov_03_potential(ener, V_pot, R_pot, a_pot, RC)
     elseif(pindex == 2)then
        RC = 1.2d0
        if(particle(pindex)%opt_pot == 1)                         &
           call KD_potential(pindex, iA, iZ, ener, V_pot, R_pot, a_pot, RC, D3)
        if(particle(pindex)%opt_pot == 2)                         &
           call soukhovitskii_potential(pindex, iA, iZ, particle(pindex)%om_option,  &
                                        ener, V_pot, R_pot, a_pot, RC)
        if(particle(pindex)%opt_pot == 3)                         &
           call maslov_03_potential(ener, V_pot, R_pot, a_pot, RC)
     elseif(pindex == 3)then
        if(particle(pindex)%opt_pot == 1)                         &
        call perey_d_potential(pindex, iA, iZ, ener, V_pot, R_pot, a_pot, RC)
     elseif(pindex == 4)then
        if(particle(pindex)%opt_pot == 1)                         &
        call becchetti_t_potential(pindex, iA, iZ, ener, V_pot, R_pot, a_pot, RC)
     elseif(pindex == 5)then
        if(particle(pindex)%opt_pot == 1)                         &
        call becchetti_h_potential(pindex, iA, iZ, ener, V_pot, R_pot, a_pot, RC)
     elseif(pindex == 6)then
        if(particle(pindex)%opt_pot == 1)                         &
        call avrigeanu_a_potential(pindex, iA, iZ, ener, V_pot, R_pot, a_pot, RC)
     end if

     fresco_name(1:iendf) = tco_file(1:len_tco)
     istart = len_tco + 1
     fresco_name(istart:istart+3) = '-ie_'
     istart = istart + 4
     iend = istart + 3
     fresco_name(istart:iend) = '0000'
     if(ie < 10)then
        write(fresco_name(iend:iend),'(i1)')ie
     elseif(ie < 100)then
        write(fresco_name(iend-1:iend),'(i2)')ie
     elseif(ie < 1000)then
        write(fresco_name(iend-2:iend),'(i3)')ie
     elseif(ie < 10000)then
        write(fresco_name(iend-3:iend),'(i4)')ie
     else
        stop 'ie > 9999 cannot run fresco'
     end if

!----    Open and create fresco input file for this energy
     open(unit=20,file = fresco_name(1:iend)//'.in', status='unknown')
!----   write fresco input
     if(particle(pindex)%opt_pot < 10)then
        write(opt_label(1:1),'(i1)')particle(pindex)%opt_pot
        ipot_end = 1
     else
        write(opt_label(1:2),'(i2)')particle(pindex)%opt_pot
        ipot_end = 2
     end if

     write(20,'(a)') fname(1:iendf)//' with potential #'//opt_label(1:ipot_end)//', at E_lab ='//char_energy
     write(20,'(a)') 'NAMELIST'
     write(20,'(a)') ' &Fresco  hcm= 0.1 rmatch=  20.000'
     absend = 1.0d-4
     if(ener <= 0.5d0)absend = 1.0d-6
     write(20,'(a,i4,a,f8.6)') '    jtmin=   0.0 jtmax=',jtmax,' absend= ',absend
     write(20,14) th_min, th_inc, th_max, ncc
14      format('    thmin=',f3.1,' thinc=',f3.1,' thmax=',f5.1,' iblock=',i3)
     write(20,'(''    chans= 1 smats= 2 xstabl= 1 tcfile= 3 iter='',i2)')iter
     write(20,15) ener
15      format('    elab=',e15.7,'  pel=1 exl=1 lab=1 lin=1 lex=1 /')
     write(20,*) 
     if(ncc == nex)then
        write(20,16) label, mass_proj, zpart, nex
16         format('&Partition  namep=''',a1,'       '' massp= ',f12.8,' zp= ',i3,' nex=',i3)
     else
        write(20,166) label, mass_proj, zpart, nex
166        format('&Partition  namep=''',a1,'       '' massp= ',f12.8,' zp= ',i3,' nex=',i3,4x,'mixpot=2')
     end if
     write(20,17) namet,mass_target,nint(Z)
17      format('            namet=''',a8,''' masst= ',f12.8,' zt= ',i3,' qval=  0.000/')
     write(20,18)spin,1,cc_state_kpp(1),cc_state_j(1),cc_state_par(1),cc_state_e(1), K_band
18      format('&States jp= ',f3.1,' ptyp=',i2,' ep=  0.000000  cpot='i3,' jt=',f4.1,' ptyt=',i2,' et=',f8.4,' kkt = ',f4.1,'/')

     do i = 2, ncc
        write(20,21)cc_state_kpp(i),cc_state_j(i),cc_state_par(i),cc_state_e(i),K_band
21      format('&States copyp= 1                       cpot=',i3,' jt=',f4.1,' ptyt=',i2,' et=',f8.4,' kkt = ',f4.1,'/')
     end do
     write(20,*)
     do i = ncc+1, nex
        write(20,21)cc_state_kpp(i),cc_state_j(i),cc_state_par(i),cc_state_e(i),K_band
     end do

     write(20,'(''&Partition /'')')
     write(20,*)

30   format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:3)=',3f10.5,'/')
31   format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:6)=',6f10.5,'/')
32   format('&POT /'/)
!---------    Write data to generate potential
     kp = 1                           !   potential for coupled-channels states
!---------    First print Coulomb
     write(20,30) kp,0,0,A,Ap,RC
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Now loop over indvidual components of nuclear potentials
!---------   ipot = 1, 6
!---------   ipot = 1 : Volume
!---------   ipot = 2 : derivative - surface
!---------   ipot = 3 : Spin-orbit
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     do ipot = 1, 3
        if(abs(V_pot(1,ipot)) > 1.0d-6)then
           write(20,31) kp,ipot,0,V_pot(1,ipot),R_pot(1,ipot),a_pot(1,ipot),0.,0.,0.
           if(ipot /= 3 .and. deformed)write(20,31)kp,11,ifresco_shape,(beta(k)*ac*R_pot(1,ipot),k=1,6)
        end if     
        if(abs(V_pot(2,ipot)) > 1.0d-6)then
           write(20,31) kp,ipot,0,0.,0.,0.,V_pot(2,ipot),R_pot(2,ipot),a_pot(2,ipot)
           if(ipot /= 3 .and. deformed)write(20,31)kp,11,ifresco_shape,(beta(k)*ac*R_pot(2,ipot),k=1,6)
        end if     
     end do
     if(ncc /= nex)then
        write(20,*)
        kp = 2
        write(20,30) kp,0,0,A,Ap,RC
        ipot = 1
        write(20,31) kp,ipot,0,V_pot(1,ipot),R_pot(1,ipot),a_pot(1,ipot),V_pot(2,ipot),R_pot(2,ipot),a_pot(2,ipot)
        ipot = 2
        write(20,31) kp,-ipot,0,V_pot(1,ipot),R_pot(1,ipot),a_pot(1,ipot),V_pot(2,ipot),R_pot(2,ipot),a_pot(2,ipot)
        write(20,31)kp,13,10,0.0,0.5,0.5,0.5,0.0,0.0

 33  format('  &step ib=',i3,1x,'ia= ',i2,1x,'k=',i2,' str=',f10.6,1x'/')
 34  format('  &step /')
        do i = ncc + 1, nex
           xk_factor = sqrt(2.0d0*J_gs + 1.0d0)
           if(cc_state_type(i) == 0)then
              K_state = real(cc_state_k(i),kind=8)
              J_state = real(cc_state_j(i),kind=8)
              xk_factor = sqrt(2.0d0*J_gs + 1.0d0)*clebr(J_gs,K_band,K_state,zzero,J_state,K_band)
           end if
           write(20,33)i,if_state,cc_state_k(i),cc_state_str(i)*ac*R_pot(1,1)*xk_factor*cc_scale
        end do
        write(20,34)
        ipot = 3
        write(20,31) kp,ipot,0,V_pot(1,ipot),R_pot(1,ipot),a_pot(1,ipot),V_pot(2,ipot),R_pot(2,ipot),a_pot(2,ipot)
     end if
     write(20,*)     
     write(20,32)
     write(20,'(''&Overlap /'')')
     write(20,'(''&Coupling /'')')

     close(unit=20)

!-------     Now run fresco
     write(6,*)'**********************************************'
     write(6,*)'Calling unix system command to execute fresco'

     line(1:132) = ' '
     ifirst = 1
     ilast = 12
     line(ifirst:ilast) = 'Calculating '
     ifirst = 13
     ilast = 13
     line(ifirst:ilast) = particle(pindex)%label
     ifirst = 14
     ilast = 16
     line(ifirst:ilast) = ' + '
     ifirst = 17
     if(iA < 10)then
        ilast = 18
        write(line(ifirst:ilast),'(i1)')iA
     elseif(iA < 100)then
        ilast = 19
        write(line(ifirst:ilast),'(i2)')iA
     elseif(iA < 1000)then
        ilast = 20
        write(line(ifirst:ilast),'(i3)')iA
     end if
     ifirst = ilast + 1
     if(symb(1:1).ne.' ')then
        ilast = ifirst + 1
        line(ifirst:ilast) = symb
     else
        ilast = ifirst
        line(ifirst:ilast) = symb(2:2)
     end if
  

     write(6,'(a)')line(1:ilast)//', at E_lab ='//char_energy

     command(1:132) = ' '
     istart = 1
     istop = 10
     command(istart:istop) = 'frescox < '
     istart = istop + 1
     istop = istop + iend + 3
     command(istart:istop) = fresco_name(1:iend)//'.in'
     istart = istop + 1
     istop = istart + 2
     command(istart:istop) = ' > '
     istart = istop + 1
     istop = istart + iend + 4
     command(istart:istop) = fresco_name(1:iend)//'.out'
     istart = istop + 1
     istop = istart + 4
     command(istart:istop) = ' 2>> '
     istart = istop + 1
     istop = istart + iend_err - 1
     command(istart:istop) = out_name(1:iend_err)

     write(6,'(a)')command(1:istop)
     call system(command(1:istop))

!--------------------------------------------------------------
!----------   Output is in fort.5420



     open(unit=20,file = 'fort.5420', status = 'old')

!----------     Read in transmission coefficients   -----------
     do ii = 0, isp_max
        read(20,'(a)')line
        call parse_string(line,numw,startw,stopw)
        if(line(1:1) == '#')then
           read(line(startw(7):stopw(7)),*)diff
           isp = nint(diff + spin)
        end if
        read_file = .true.
        do while(read_file)
           read(20,'(a)')line
           call parse_string(line,numw,startw,stopw)
           if(numw == 1)then
              if(line(startw(1):stopw(1)) == '&')read_file = .false.
           else
              read(line,*)l,tco_data
              if(tco_data < 1.0d-17)tco_data = 1.0d-17
              if(l <= particle(pindex)%lmax)particle(pindex)%trans_read(ie,ii,l) = tco_data
           end if
        end do
     end do


!---------    Read line, has number of coupled channels in it
     read(20,'(a)')line
!---------    Read in angular distributions and cross sections to 
!---------    coupled states

!------------------------------------------------------------------------------------------
!------------   Angular distributions for elastic and direct parts    -------------------
!------------------------------------------------------------------------------------------
     if(pindex == iproj)then
        do n = 1, nex
           read(20,*)itar, ichan, num_th, th_inc, th_min, cross
              if(cross < 1.0d-9)cross = 0.0d0
              if(isnan(cross))cross = 0.0d0
              if(ichan .ne. itar .and. cross > 1.0d4)cross = 0.0d0
              optical_cs(ie,n) = cross/1000.0d0
              do it = 1, num_th
                 read(20,*)x, y
                 theta(it) = x
                 if(cross < 1.0d-9)y = 0.0d0
                 ang_dist(it,n,ie) = y/1000.0d0
                 if(pindex > 1)then
                    sig_C = 0.0d0
                    ang_rad = theta(it)*pi/180.0d0
                    xmu = cos(ang_rad)
                    sig_C = Coul/energy(ie)**2/(1.0d0 - xmu)**2*fmsq_eq_barn   ! Rutherford in b/sr
                    ang_dist(it,n,ie) = ang_dist(it,n,ie)/sig_C
                 end if
              end do
              
           xnorm = 0.0d0
           do ix = 1, n_gleg2
              ang_theta = acos(x_gleg2(ix))*180.0d0/pi
              value = interp(ang_theta, num_ang, theta, ang_dist(1,n,ie))
              xnorm = xnorm + w_gleg2(ix)*value
           end do
           if(iproj > 1)optical_cs(ie,n) = xnorm

           do L = 0, Ang_L_max
              optical_leg(ie,l,n) = 0.0d0
              alf = 0.0d0
              bet = 0.0d0
              sum = 0.0d0
              do ix = 1, n_gleg2
                 ang_theta = acos(x_gleg2(ix))*180.0d0/pi
                 value = interp(ang_theta, num_ang, theta, ang_dist(1,n,ie))
                 sum = sum + w_gleg2(ix)*value*poly(L,1,alf,bet,x_gleg2(ix))
              end do
              if(xnorm > 1.0d-20)optical_leg(ie,L,n) = sum*0.5d0*(2.0d0*real(L,kind=8)+1.0d0)/xnorm
              optical_leg(ie,L,n) = sum*0.5d0*(2.0d0*real(L,kind=8)+1.0d0)
           end do
           do it = 1, num_th
              ang_rad = theta(it)*pi/180.0d0
              x = cos(ang_rad)
              sum = 0.0d0
              do L = 0, Ang_L_max
                 sum = sum + optical_leg(ie,L,n)*poly(L,1,alf,bet,x)
              end do
              sum = sum
           end do
        end do
        close(unit=20)
     end if
!-----    Delete fort.* files, which will ensure that it is not possible
!-----    to continue on if fresco crashes
     call system('rm fort.*')
  end do

  tco_name(1:iendf) = fname(1:iendf)
  iend_tco = iendf + 4

  write(6,*)'Writing transmission coefficients'

!-------    Write transmission coefficients to file for general use later
  open(unit = 51, file = tco_file(1:len_tco)//'.tcoef', status = 'unknown')
  title(1:1) = particle(pindex)%label
  title(2:5) = ' on '
  iend = 5 + iendf - 1
  title(6:iend) = tco_file(2:iendf)
  write(51,'(''# '',a50)')title(1:iend)
  write(51,'(''# Transmission coefficients computed using FRESCO'')')
  write(51,'(''# nume = '',i6)')nume
  write(51,'(''# lmax = '',i6)')lmax
  write(51,'(''# Energy list '')')
  num_e_lines = nume/10
  if(10*num_e_lines /= nume)num_e_lines = num_e_lines + 1
  do ll = 1, num_e_lines
     minl = (ll-1)*10 + 1
     maxl = min(nume,minl+9)
     write(51,'(10(1x,1pe15.9))')(particle(pindex)%e_grid(ie),ie = minl, maxl)
  end do
  num_l_lines = (lmax+1)/10
  if(10*num_l_lines /= lmax)num_l_lines = num_l_lines + 1
  do ii = 0, isp_max
     write(51,'(''# J-L = '',f4.1)')real(ii,kind=8) - spin
     do ll = 1, num_l_lines
        minl = (ll-1)*10
        maxl = min(lmax,minl+9)
        write(51,'(''#         energy'',10(1x,8x,''l = '',i3))')(l,l=minl,maxl)
        write(51,'(''#---------------'',10(1x,''---------------''))')
        do ie = 1, nume
           write(51,'(11(1x,1pe15.9))')particle(pindex)%e_grid(ie),(particle(pindex)%trans_read(ie,ii,l),l=minl,maxl)
        end do
     end do
  end do
  close(unit=51)


  if(pindex == iproj .and. .not. pop_calc)then
     open(unit = 12, file = tco_file(1:len_tco)//'-CC.data', status = 'unknown')

     write(12,'(3(1x,i5),1x,f10.4)')nume, nex, Ang_L_max

     do n = 1, nex
        write(12,'(2(1x,i5),1x,f4.1,1x,f4.1,1x,f4.1,1x,f8.4,1x,i6)')                &
              n, cc_index(n), cc_state_j(n),                                        &
              real(cc_state_par(n),kind=8), real(cc_state_K(n),kind=8),             &
              cc_state_e(n), cc_state_type(n)
     end do

     num_l_lines = (nex+1)/10
     if(10*num_l_lines /= nex)num_l_lines = num_l_lines + 1
     do ll = 1, num_l_lines
        minl = (ll-1)*10+1
        maxl = min(nex,minl+9)
        write(12,'(''#         energy'',10(1x,8x,''n = '',i3))')(l,l=minl,maxl)
        do ie = 1, nume
           write(12,'(1x,e15.7,10(1x,(e15.7)))')energy(ie),(optical_cs(ie,n),n = minl,maxl)
        end do
     end do

     do n = 1, nex
        write(12,'(1x,i5,1x,i5)')n,Ang_L_max
        do ie = 1, nume
           write(12,'(1x,e15.7,50(1x,e15.7))')energy(ie),(optical_leg(ie,l,n),l = 0, Ang_L_max)
        end do
     end do
     close(unit=12)
  end if

!------------------------------------------------------------------------------------------
!---------    Now finish up and clean up uneeded variables      ---------------------------
!------------------------------------------------------------------------------------------

  if(allocated(energy))deallocate(energy)
  if(allocated(optical_cs))deallocate(optical_cs)
  if(allocated(optical_leg))deallocate(optical_leg)
  if(allocated(theta))deallocate(theta)
  if(allocated(ang_dist))deallocate(ang_dist)
  if(allocated(x_gleg2))deallocate(x_gleg2)
  if(allocated(w_gleg2))deallocate(w_gleg2)
  if(allocated(cc_index))deallocate(cc_index)
  if(allocated(cc_state_j))deallocate(cc_state_j)
  if(allocated(cc_state_e))deallocate(cc_state_e)
  if(allocated(cc_state_par))deallocate(cc_state_par)
  if(allocated(cc_state_type))deallocate(cc_state_type)
  if(allocated(cc_state_k))deallocate(cc_state_k)
  if(allocated(cc_state_kpp))deallocate(cc_state_kpp)
  if(allocated(cc_state_str))deallocate(cc_state_str)

  return
end subroutine make_fresco_tco
