!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine get_spectrum(data_path, len_path, overide, symb, iz, ia, inuc)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine reads in spectroscopic data and stores it in arrays
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
   use nodeinfo
   use nuclei
   implicit none
   character(len=200), intent(in) :: data_path
   integer(kind=4), intent(in) :: len_path
   character(len=1), intent(in) :: overide
   character(len=2), intent(in) :: symb
   integer(kind=4), intent(in) :: iz,ia
   integer(kind=4), intent(in) :: inuc
!-------------   Data for reading levels 
   real(kind=8) :: ex_energy,spin
   integer(kind=4) :: ilab
   integer(kind=4) :: par
!--------------
   character(len=1) :: toveride
   character(len=80) :: fname
   character(len=300) :: line
   integer(kind=4) :: len_fname
   character(len=5) :: char
   integer(kind=4) :: ipar
   integer(kind=4) :: i, j, m, n
   integer(kind=4) :: num_s
   character(len=2) :: symbb
   integer(kind=4) :: izp, iap
   integer(kind=4) :: nol, nog, nmax, nc
   integer(kind=4) :: n_cut
   real(kind=8) :: sn,sp
   integer(kind=4) :: ng, ngr, nd
   integer(kind=4) :: num, numd
   real(kind=8) :: t12
   integer(kind=4) :: nf(99)
   real(kind=8) :: eg(99), pg(99), pe(99), icc(99)   
   character(len=1) :: nf_mod(99), eg_mod(99), pg_mod(99), pe_mod(99), icc_mod(99)   
   character(len=1) :: unc
   real(kind=8) :: shift
   integer(kind=4) :: iband
   real(kind=8) :: s1, s2
   real(kind=8) :: sum
   logical fexist
   logical evaluated
   real(kind=8) :: target_spin
   integer(kind=4) :: target_ipar
   real(kind=8) :: alpha
   real(kind=8) :: decay_prob(10)
   character(len=7) :: decay_type(10)
   integer(kind=4) :: isymb, ilast
   logical :: check

   character(len=1) :: e_mod, j_mod, par_mod, t12_mod, ngr_mod
  
   integer(kind=4), allocatable :: complete_decay(:)
   logical :: complete
   integer(kind=4), allocatable :: state_map(:)
   real(kind=8) :: ecut
   logical :: found
   integer(kind=4) :: io_error
!-----------------------------------------------------------------------  
   evaluated=.false.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------ find file name with data
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   fname='z000.dat'
   if(iz < 10)write(fname(4:4),'(i1)')iz
   if(iz >= 10.and.iz < 100)write(fname(3:4),'(i2)')iz
   if(iz >= 100)write(fname(2:4),'(i3)')iz
   len_fname=index(fname,' ')-1

   isymb = 1
   if(symb(1:1) == ' ')isymb = 2
   ilast = 5
   if(nucleus(inuc)%Label(5:5) == ' ')ilast = 4
!
!----------------------------------------------------------------------
   target_spin = -1.0d0
   ipar = -1
   open(unit=51,                                                      &
        file=data_path(1:len_path)//'levels/'//fname(1:len_fname),    &
        status='old')
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------   First look for Z,A-1 element to get its ground state spin
!----------   Used to compute D0, i.e., the level spacing at the 
!----------   neutron-separation energy 
!------  Search for element in data file
!------  in A and element name if same get data
!------   Use a default value to signal if it is known or not
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   s1 = -1.5d0
   s2 = -0.5d0
   found = .false.
   io_error = 0
   do while(.not. found .and. io_error == 0)
      read(51,'(i3,a2)',iostat=io_error)iap,symbb !  find element, end -. end of subroutine abort
      if(symbb(2:2) == ' ')then
         symbb(2:2) = symbb(1:1)
         symbb(1:1) = ' '
      end if
      if(ia-1 == iap .and. symb == symbb)then
         read(51,'(i3,1x,f10.6,1x,f5.1,i3,1x,(e10.2),i3,1x,a1,32x,f9.4)')   &
         ilab, ex_energy, spin, par, t12, ng, unc
         found = .true.
         if(par == -1)then           !  par =-1,1 (0 if indeterminant)
            ipar = 1                 !  but we'll use 0,1 for positive and negative parity
         elseif(par == 1)then                        !  i.e., parity =(-1)**ipar
            ipar = 0                 !  this lets us use 0,1 for an index
         elseif(par == 0)then        !  Parity not defined, make even
            ipar = 0
         end if
         if(spin < 0.0d0)then        !   ground-state spin is not known, and undefined, make 0 or 0.5
            spin = 0.0d0
            if(iand(ia-1,1) == 1)spin = 0.5d0
         end if
         if(spin <= 1.0d-4)then
            s1 = spin + 0.5d0
            s2 = spin + 0.5d0
         else
            s1 = spin - 0.5d0
           s2 = spin + 0.5d0
          end if
         target_spin = spin
         target_ipar = ipar
      end if
   end do
   if(.not. found)then
      spin = 0.0d0
      ipar = 1
      if(iand(ia-1,1) == 1)spin = 0.5d0
      s1 = spin + 0.5d0
      s2 = spin + 0.5d0
      target_spin = spin
      target_ipar = ipar
   end if

   close(unit=51)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------  Open data file
!------   First check to see if it exists in the evaluated area
!------   set internal overide switch to that passed in
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   evaluated = .false.
   toveride = overide 
   fexist = .false.
   found = .false.
   inquire(file=data_path(1:len_path)//'levels-eval/'//fname(1:len_fname),exist=fexist)
   if(fexist)then
      if(iproc == 0)write(6,'(''Looking in evaluated file for '',a5)')nucleus(inuc)%Label(1:ilast)
      open(unit=51,                                                                   &
         file=data_path(1:len_path)//'levels-eval/'//fname(1:len_fname),              &
         status='old')
      io_error = 0
      do while(.not. found .and. io_error == 0)
         read(51,'(i3,a2)',iostat=io_error)iap,symbb !  find element, end -. end of subroutine abort
         if(symbb(2:2) == ' ')then
            symbb(2:2) = symbb(1:1)
            symbb(1:1) = ' '
         end if
!  write(6,*)ia,iap,symb, symbb
         if(ia == iap .and. symb == symbb)then
            backspace(51)
            read(51,'(a5,6i5,2f12.6)',iostat=io_error)char, iap, izp, nol, nog, nmax, nc, sn, sp
            found = .true.
            evaluated = .true.               ! skip optional additional input
         end if
      end do
      if(.not. found)then
          if(iproc == 0)write(6,'(a5,'' Not found in file '',a70)')                           &
                                     nucleus(inuc)%Label(1:ilast),                            &
                                     data_path(1:len_path)//'levels-eval/'//fname(1:len_fname)
          close(unit=51)
      end if
   end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----    It isn't in the evaluated file, so try standard RIPL file
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(.not. found)then
      if(iproc == 0)write(6,'(''Looking in standard RIPL-3 file '',a70)')                  &
            data_path(1:len_path)//'levels/'//fname(1:len_fname)
      open(unit=51,file=data_path(1:len_path)//'levels/'//fname(1:len_fname),              &
           status='old')
      io_error = 0
      do while(.not. found .and. io_error == 0)
         read(51,'(i3,a2)',iostat=io_error)iap,symbb !  find element, end -. end of subroutine abort
         if(symbb(2:2) == ' ')then
            symbb(2:2) = symbb(1:1)
            symbb(1:1) = ' '
         end if
         if(ia == iap .and. symb == symbb)then
            backspace(51)
            read(51,'(a5,6i5,2f12.6)',iostat=io_error)char, iap, izp, nol, nog, nmax, nc, sn, sp
            found = .true.
         end if
      end do
   end if
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----   Found nucleus, read data and set up transitions
!-----   First, we read everything and then trim
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(found)then
      if(.not. allocated(complete_decay))allocate(complete_decay(nol))
      if(.not. allocated(state_map))allocate(state_map(nol))
      nucleus(inuc)%num_discrete = nol
      if(.not. allocated(nucleus(inuc)%state))allocate(nucleus(inuc)%state(nol))
      complete_decay(1:nol) = 1
!      complete_decay(1) = 1
      state_map(1:nol) = 1
      do i = 1, nol
         line(1:300) = ' '
         read(51,'(a)')line
!         read(line,'(i3,1x,f10.6,1x,f5.1,i3,1x,e10.3,i3,1x,a1,24x,i3,220x,f10.6,1x,i2)')      &
!              ilab, ex_energy, spin, par, t12, ngr, unc, nd, shift, iband
         read(line,'(i3,1x,f10.6,a1,f5.1,a1,i2,a1,e10.3,a1,i2,a1,a1,24x,i3,220x,f10.6,1x,i2)')      &
              ilab, ex_energy, e_mod, spin, j_mod, par, par_mod, t12, t12_mod, ngr, ngr_mod,        &
              unc, nd, shift, iband
         nucleus(inuc)%state(i)%energy = ex_energy
         nucleus(inuc)%state(i)%spin = spin
         nucleus(inuc)%state(i)%parity = par
         nucleus(inuc)%state(i)%shift = shift
         nucleus(inuc)%state(i)%level_float = .false.
         if(shift > 0.0d0)nucleus(inuc)%state(i)%level_float = .true.
         nucleus(inuc)%state(i)%iband = iband
         nucleus(inuc)%state(i)%isomer = .false.
         nucleus(inuc)%state(i)%t12 = t12
         nucleus(inuc)%state(i)%state_modified = .false.
         if(e_mod /= ' '.or. j_mod /= ' ' .or. par_mod /= ' ' .or. t12_mod /= ' ' .or. ngr_mod /= ' ')&
            nucleus(inuc)%state(i)%state_modified = .true.
         if(spin < 0.0d0 .or. par == 0 .or. shift > 0.0d0)state_map(i) = -1
         if(ngr > 99)then
            if(iproc == 0)write(6,*)'Too many gammas, increase ngr in spectrum to at least ',ngr
            call MPI_Abort(icomm,101,ierr)
         end if
         ng = 0
         do m = 1, ngr
            line(1:300) = ' '
            read(51,'(a)')line
!            read(line,'(39x,i4,1x,f10.4,3(1x,e10.3))')nf(m),eg(m),pg(m),pe(m),icc(m)
            read(line,'(39x,i4,a1,f10.4,a1,e10.3,a1,e10.3,a1,e10.3,a1)')                               &
                                                      nf(m), nf_mod(m), eg(m), eg_mod(m),           &
                                                      pg(m), pg_mod(m), pe(m), pe_mod(m),           &
                                                      icc(m), icc_mod(m)
            if(pe(m) < pg(m))pe = pg(m)*(1.0d0 + icc(m))           !   In case there is an error in the RIPL file
            if(pe(m) > 1.0d-30)ng = ng + 1                            !   Check and remove transtions that are zero
         end do
         nucleus(inuc)%state(i)%nbranch = ng
         if(ng > 0)then
            allocate(nucleus(inuc)%state(i)%ibranch(ng))       !  states
            allocate(nucleus(inuc)%state(i)%branch(ng))        !  branching ratio
            allocate(nucleus(inuc)%state(i)%p_gamma(ng))       !  branching ratio
            allocate(nucleus(inuc)%state(i)%egamma(ng))        !  gamma-ray energy
            allocate(nucleus(inuc)%state(i)%p_ic(ng))          !  internal conversion coefficient
            allocate(nucleus(inuc)%state(i)%cs(ng))            !  Decay cross section
            allocate(nucleus(inuc)%state(i)%branch_modified(ng))            !  Decay cross section
            nucleus(inuc)%state(i)%branch_modified(1:ng) = .false.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Get the gamma decay branches
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            sum = 0.0d0                                         !   Norm of em decay if other branches are open
            complete = .true.
            n= 0
            do m = 1, ngr
               if(pe(m) > 1.0d-30)then                                    !   Include only non-zero transitions
                  n = n + 1
                  sum = sum + pe(m)
                  nucleus(inuc)%state(i)%p_gamma(n) = pg(m)
                  nucleus(inuc)%state(i)%ibranch(n) = nf(m)
                  nucleus(inuc)%state(i)%egamma(n) = eg(m)
                  nucleus(inuc)%state(i)%branch(n) = pe(m)
                  nucleus(inuc)%state(i)%p_ic(n) = icc(m)
                  nucleus(inuc)%state(i)%cs(n) = 0.0d0
                  if(complete_decay(nf(m)) == 0 .or. state_map(nf(m)) == -1)complete = .false.         !  Check if final state has a complrete decay chain
                  nucleus(inuc)%state(i)%branch_modified(n) = .false.
                  if(nf_mod(m) /= ' ' .or. eg_mod(m) /= ' ' .or. pg_mod(m) /= ' ' .or.           &
                     pe_mod(m) /= ' ' .or. icc_mod(m) /= ' ')then
                     nucleus(inuc)%state(i)%branch_modified(n) = .true.
                     nucleus(inuc)%state(i)%state_modified = .true.
                  end if
               end if
            end do
            if(.not. complete .and. .not. nucleus(inuc)%state(i)%isomer)complete_decay(i) = 0   ! isomers are compelte
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   They need to be properly normalized since they are < 1 if other decays are present
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(sum > 0.0d0)then
               do m = 1, ng
                  nucleus(inuc)%state(i)%branch(m) = nucleus(inuc)%state(i)%branch(m)/sum
                  alpha = (1.0d0 + nucleus(inuc)%state(i)%p_ic(m))
                  nucleus(inuc)%state(i)%p_gamma(m) = nucleus(inuc)%state(i)%branch(m)/alpha
                  if(nucleus(inuc)%state(i)%branch(m) > 0.0d0)then
                     nucleus(inuc)%state(i)%p_gamma(m) = nucleus(inuc)%state(i)%p_gamma(m)/   &
                                                          nucleus(inuc)%state(i)%branch(m)  
                     nucleus(inuc)%state(i)%p_ic(m) = 1.0d0 - nucleus(inuc)%state(i)%p_gamma(m)
                  end if
               end do
            end if
         elseif(i > 1)then
            complete_decay(i) = 0
         end if
         if(i/= 1 .and. t12 > t12_isomer)then
            nucleus(inuc)%state(i)%isomer = .true.
            complete_decay(i) = 1
         end if
         if(nd == 0)then
            nd = 1
            nucleus(inuc)%state(i)%num_decay = nd
            allocate(nucleus(inuc)%state(i)%decay_prob(nd))
            nucleus(inuc)%state(i)%decay_prob(nd) = 1.0d0
            allocate(nucleus(inuc)%state(i)%decay_type(nd))
            nucleus(inuc)%state(i)%decay_type(nd) = 0
            allocate(nucleus(inuc)%state(i)%decay_to(nd))
            nucleus(inuc)%state(i)%decay_to(nd) = inuc
         else
            nucleus(inuc)%state(i)%num_decay = nd
            allocate(nucleus(inuc)%state(i)%decay_prob(nd))
            allocate(nucleus(inuc)%state(i)%decay_type(nd))
            allocate(nucleus(inuc)%state(i)%decay_to(nd))
            nucleus(inuc)%state(i)%decay_to(1:nd) = 0
            line(1:300) = ' '
            read(line,'(66x,10(4x,e10.4,1x,a7))')(decay_prob(m),decay_type(m),m=1,nd)
            sum = 0.0d0
            do m = 1, nd
               sum = sum + decay_prob(m)
               if(decay_type(m)(1:3) == '%IT')nucleus(inuc)%state(i)%decay_type(m) = 0
               if(decay_type(m)(1:2) == '%G')nucleus(inuc)%state(i)%decay_type(m) = 0
               if(decay_type(m)(1:2) == '%N')nucleus(inuc)%state(i)%decay_type(m) = 1
               if(decay_type(m)(1:2) == '%P')nucleus(inuc)%state(i)%decay_type(m) = 2
               if(decay_type(m)(1:2) == '%A')nucleus(inuc)%state(i)%decay_type(m) = 6
               if(decay_type(m)(1:5) == '%3HE')nucleus(inuc)%state(i)%decay_type(m) = 5
               if(decay_type(m)(1:2) == '%B')nucleus(inuc)%state(i)%decay_type(m) = -1
               if(decay_type(m)(1:2) == '%E')nucleus(inuc)%state(i)%decay_type(m) = -1
            end do
            do m = 1, nd
               nucleus(inuc)%state(i)%decay_prob(m) = decay_prob(m)/sum
            end do
         end if
         if(nucleus(inuc)%state(i)%isomer)complete_decay(i) = 1
      end do
      close(unit=51)
      par = nint(nucleus(inuc)%state(1)%parity)
      nucleus(inuc)%ipar = (par+1)/2
      nucleus(inuc)%s1 = s1
      nucleus(inuc)%s2 = s2
      nucleus(inuc)%target_spin = target_spin
      nucleus(inuc)%target_ipar = target_ipar
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------   Rare case where ground-state spin is not known
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(nucleus(inuc)%state(1)%spin < 0.0d0)then
         if(iand(ia,1) == 1)then             !   Assume if A odd J=0.5 or J=0 if A even
            spin = 0.5d0
         else
            spin = 0.0d0
         end if
         nucleus(inuc)%ipar = 0
         nucleus(inuc)%state(1)%spin = spin
         nucleus(inuc)%state(1)%parity = 1.0d0
         nucleus(inuc)%state(1)%t12 = 100000.d0 
         nucleus(inuc)%state(1)%nbranch = 0              !  # of branches       
         nucleus(inuc)%state(1)%pop = 0.0d0
         nucleus(inuc)%state(1)%exit_lab = 1
      end if
!--------     Nucleus wasn't found, so set up with dummy values
   else
      close(unit=51)
      if(iproc == 0)then
         write(6,'(i3,a2,'' Not found in file '',a70)')ia, symb(1:isymb),                   &
                                     data_path(1:len_path)//'levels-eval/'//fname(1:len_fname)
         write(6,*)'Will set up with no levels with "dummy" spins and parities to continue'
      end if
      nucleus(inuc)%num_discrete = 1
      allocate(nucleus(inuc)%state(nucleus(inuc)%num_discrete))
      nucleus(inuc)%state(1)%energy = 0.0d0
      if(iand(ia,1) == 1)then             !   Assume if A odd J=0.5 or J=0 if A even
         spin = 0.5d0
      else
         spin = 0.0d0
      end if
      nucleus(inuc)%ipar = 0
      nucleus(inuc)%state(1)%spin = spin
      nucleus(inuc)%s1 = s1
      nucleus(inuc)%s2 = s2
      nucleus(inuc)%state(1)%parity = 0
      nucleus(inuc)%state(1)%t12 = 1000.
      nucleus(inuc)%state(1)%nbranch = 0              !  # of branches       
      nucleus(inuc)%state(1)%pop = 0.0
      nucleus(inuc)%state(1)%exit_lab = 1
      nucleus(inuc)%target_spin = target_spin
      nucleus(inuc)%target_ipar = target_ipar
      nucleus(inuc)%ncut = 1
      nucleus(inuc)%level_ecut = nucleus(inuc)%state(nucleus(inuc)%ncut)%energy + 0.001d0
      return
   end if

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------    The arrays are filled, now it is time to cull states
!------    Remove unassigned J and parity, states that are floating,
!------    states without complete decay paths, unless an isomer
!------    First find Ecut and ncut=number of levels below Ecut. Define Ecut as first level
!------    without definitive spin, parity, a complete gamma decay chain or no gammas. 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   n_cut = 1
   do i = 2, nucleus(inuc)%num_discrete
!     write(6,*)i, nucleus(inuc)%state(i)%level_float
!     write(6,*)i, nucleus(inuc)%state(i)%spin
!     write(6,*)i, nucleus(inuc)%state(i)%parity
!     write(6,*)i, complete_decay(i)
!     write(6,*)i, nucleus(inuc)%state(i)%nbranch
      if(nucleus(inuc)%state(i)%level_float)cycle
      if(nucleus(inuc)%state(i)%spin < 0.0d0)exit
      if(nint(nucleus(inuc)%state(i)%parity) == 0)exit
      if(complete_decay(i) == 0)exit
      if(nucleus(inuc)%state(i)%nbranch == 0 .and.              &
         .not. nucleus(inuc)%state(i)%isomer)exit
      n_cut = i
   end do
   nucleus(inuc)%ncut = n_cut
   ecut = nucleus(inuc)%state(n_cut)%energy + 0.010d0
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Now set state_map for the states to be removed
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   write(6,*)'state_map'
!----   Always have the ground state
   numd = 1
   state_map(1) = numd
   do i = 2, nucleus(inuc)%num_discrete
      state_map(i) = -1
      check = .true.
      if(nucleus(inuc)%state(i)%spin < 0.0d0)check = .false.
      if(nint(nucleus(inuc)%state(i)%parity) == 0)check = .false.
      if(complete_decay(i) == 0)check = .false.
      if(nucleus(inuc)%state(i)%nbranch == 0 .and. .not. nucleus(inuc)%state(i)%isomer )check = .false. 
      if(nucleus(inuc)%state(i)%level_float)check = .false.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Failsafe - make sure that all transitions go to a safe, i.e., 
!----   not removed state. If it does, also mark this state for removal.
!----   Not likely, but if it does, it would create a seg-fault later
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do m = 1, nucleus(inuc)%state(i)%nbranch
         if(state_map(nucleus(inuc)%state(i)%ibranch(m)) == -1)check = .false.
      end do
      if(check)then              !  it passed, increment kept states and update map
         numd = numd + 1
         state_map(i) = numd
      end if
!      write(6,*)i, state_map(i), nucleus(inuc)%state(i)%spin, nucleus(inuc)%state(i)%parity,      &
!                nucleus(inuc)%state(i)%energy, nucleus(inuc)%state(i)%nbranch, complete_decay(i), &
!                nucleus(inuc)%state(i)%level_float
   end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Now remove the unwanted states
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   num = nucleus(inuc)%num_discrete
   call remove_states(inuc, num, state_map)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Set flag that identifies state to be collected as
!------   decay channel
!------   Also, find ncut again, after levels are removed
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   num_s = 1
   nucleus(inuc)%state(1)%exit_lab = num_s
   do j = 2, nucleus(inuc)%num_discrete
      nucleus(inuc)%state(j)%exit_lab = 0
      if(nucleus(inuc)%state(j)%isomer)then
         num_s = num_s + 1
         nucleus(inuc)%state(j)%exit_lab = num_s
      end if
      if(nucleus(inuc)%state(j)%energy < ecut)n_cut = j
   end do
   nucleus(inuc)%ncut = n_cut
   nucleus(inuc)%level_ecut = nucleus(inuc)%state(nucleus(inuc)%ncut)%energy + 0.001d0
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Should now be done
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   write(6,*)'inuc = ',inuc
!   write(6,*)'num_discrete = ',nucleus(inuc)%num_discrete
!   write(6,*)'ncut = ',nucleus(inuc)%ncut
!   do i = 1, nucleus(inuc)%num_discrete
!      write(6,*)nucleus(inuc)%state(i)%spin, nucleus(inuc)%state(i)%parity,                       &
!                nucleus(inuc)%state(i)%energy, nucleus(inuc)%state(i)%isomer
!   end do
!   do i = 1, nucleus(inuc)%num_discrete
!      write(6,*)i, nucleus(inuc)%state(i)%energy, nucleus(inuc)%state(i)%nbranch
!      do n = 1, nucleus(inuc)%state(i)%nbranch
!         write(6,*)nucleus(inuc)%state(i)%ibranch(n),nucleus(inuc)%state(i)%egamma(n),            &
!                   nucleus(inuc)%state(i)%branch(n),nucleus(inuc)%state(i)%p_gamma(n),            &
!                   nucleus(inuc)%state(i)%p_ic(n)
!      end do
!   end do

   if(allocated(complete_decay))deallocate(complete_decay)
   if(allocated(state_map))deallocate(state_map)

   return
end subroutine get_spectrum
