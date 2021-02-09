!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine get_spectrum(data_path,len_path,overide,       &
                        t12_isomer,symb,iz,ia,icomp)
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
  real(kind=8), intent(in) :: t12_isomer
  character(len=2), intent(in) :: symb
  integer(kind=4), intent(in) :: iz,ia
  integer(kind=4), intent(in) :: icomp
!-------------  Data for reading levels 
  real(kind=8) :: spectrum,spin
  integer(kind=4) :: ilab
  integer(kind=4) :: par
!--------------
  character(len=1) :: toveride
  character(len=80) :: fname
  character(len=120) :: line
  integer(kind=4) :: len_fname
  character(len=5) :: char
  integer(kind=4) :: ipar
  integer(kind=4) :: i,j
  integer(kind=4) :: num_s
  character(len=2) :: symbb
  integer(kind=4) :: izp,iap
  integer(kind=4) :: nol,nog,nmax,nc
  integer(kind=4) :: ncc,ncount
  real(kind=8) :: sn,sp
  integer(kind=4) :: ng
  character(len=1) :: jf
  character(len=2) :: char1
  real(kind=8) :: t12
  integer(kind=4) :: nf, nff
  real(kind=8) :: eg,pg,pe,icc   
  real(kind=8) :: btot, bicc
  real(kind=8) :: s1,s2
  real(kind=8) :: spec_old
  logical ifail
  logical isomer
  logical fexist
  logical evaluated
  real(kind=8) :: target_spin
  integer(kind=4) :: target_ipar
  real(kind=8) :: alpha, xnorm
  
  integer(kind=4), allocatable :: complete_decay(:)
  integer(kind=4) :: complete
  integer(kind=4) :: num_states, n_cut, ii
  integer(kind=4), allocatable :: state_map(:)
  logical E_cut_set
!-----------------------------------------------------------------------  
  evaluated=.false.
!------ find file name with data
  fname='z000.dat'
  if(iz < 10)write(fname(4:4),'(i1)')iz
  if(iz >= 10.and.iz < 100)write(fname(3:4),'(i2)')iz
  if(iz >= 100)write(fname(2:4),'(i3)')iz
  len_fname=index(fname,' ')-1
!
!----------------------------------------------------------------------
  target_spin = -1.0d0
  ipar = -1
  num_states = -1
  ncount = -1
  open(unit=51,                                                      &
       file=data_path(1:len_path)//'levels/'//fname(1:len_fname),    &
       status='old')
!--------------------------------------------------------------------------
!----------   First look for Z,A-1 element to get its ground state spin
!----------   Used to compute D0, i.e., the level spacing at the 
!----------   neutron-separation energy 
!------  Search for element in data file
!------  in A and element name if same get data
!------  Use a default value to signal if it is known or not
  s1 = -1.5d0
  s2 = -0.5d0
 21   read(51,'(i3,a2)',end=11)iap,symbb !  find element, end -. end of subroutine abort
  if(symbb(2:2) == ' ')then
     symbb(2:2) = symbb(1:1)
     symbb(1:1) = ' '
  end if
  if(.not.(ia-1 == iap .and. symb == symbb))goto 21    ! not same read another line
     read(51,'(i3,1x,f10.6,1x,f5.1,i3,1x,(e10.2),i3,1x,a1)',err=120)   &
           ilab,spectrum,spin,par,t12,ng,jf
  if(par == -1)then           !  par =-1,1 (0 if indeterimant)
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
  target_spin = spin
  target_ipar = ipar
  if(spin <= 1.0d-4)then
     s1 = spin + 0.5d0
     s2 = spin + 0.5d0
  else
     s1 = spin - 0.5d0
     s2 = spin + 0.5d0
  end if
  goto 12
 11   continue                                         !   nucleus not there
 12   continue
  close(unit=51)

!------  Open data file
!------   First check to see if it exists in the evaluated area
!------   set internal overide switch to that passed in
  toveride=overide 
  fexist=.false.
  inquire(file=data_path(1:len_path)//'levels-eval/'//fname(1:len_fname),exist=fexist)
  if(fexist)then
     if(iproc == 0)write(6,'(''Looking in evaluated file for '',i3,a2)')ia,symb
     open(unit=51,                                                                   &
        file=data_path(1:len_path)//'levels-eval/'//fname(1:len_fname),              &
        status='old')
 100     read(51,'(i3,a2)',end=101)iap,symbb
     if(symbb(2:2) == ' ')then
        symbb(2:2)=symbb(1:1)
        symbb(1:1)=' '
     end if
     if(.not.(ia == iap .and. symb == symbb))goto 100    ! not same read another line
     toveride = ' '                   ! only states marked with u
     evaluated = .true.               ! skip optional additional input
     goto 102
 101     continue
         if(iproc == 0)write(6,'(i3,a2,'' Not found in file '',a70)')                &
         ia,symb,data_path(1:len_path)//'levels-eval/'//fname(1:len_fname)
         close(unit=51)
     else
         if(iproc == 0)write(6,'(''Using standard RIPL-3 file '',a70)')              &
            data_path(1:len_path)//'levels/'//fname(1:len_fname)
  end if
!-----   Found data, open file
  open(unit=51,                                                                      &
       file=data_path(1:len_path)//'levels/'//fname(1:len_fname),                    &
       status='old')
 1    read(51,'(i3,a2)',end=10)iap,symbb        !  find element, end -. end of subroutine abort
  if(symbb(2:2) == ' ')then
     symbb(2:2)=symbb(1:1)
     symbb(1:1)=' '
  end if
  if(.not.(ia == iap .and. symb == symbb))goto 1    ! not same read another line
!-------  Found the nucleus we want
 102  ifail=.false.                                   ! when ifail = true stop reading
  backspace(51)
  read(51,'(a5,6i5,2f12.6)')char,iap,izp,nol,nog,nmax,nc,sn,sp
!------   Count how many states we can keep
!------   Slight modification to limit states based on observation for Ce132 where
!------   higher super-deformed bands are intruduced, but the exact band-head energy isn't 
!------   known, hence it resets to zero. Therefore, stop counting if excitation energy
!------   is less than previous state. Possible trap if all prior states have known
!------   gamma transitions.
  spec_old = -1.0d0
  ncount = 0                                        !  total number of reads
  ncc = 0                                           !  number of states we'll keep

  if(.not. allocated(complete_decay))allocate(complete_decay(nol))
  if(.not. allocated(state_map))allocate(state_map(nol))

  complete_decay(1:nol) = 0
  complete_decay(1) = 1
  num_states = 0

  do i = 1, nol
     ifail = .false.
     char1='  '
     line(1:120) = ' '
     read(51,'(a)',end=120)line

     ncount = ncount + 1                            !  increment read counter
     if(line(4:5) /= '  ')goto 120
     read(line,'(i3,1x,f10.6,1x,f5.1,i3,1x,(e10.2),i3,1x,a1)')                       &
          ilab,spectrum,spin,par,t12,ng,jf

     ifail = .false.
     if(spin < 0.0d0)ifail = .true.                 !  spin=-1.0 indicates unknown spin assignment in ENSDF file
     if(par == 0)ifail = .true.                     !  par=0 also indicates unknwon parity assignment.

     if(spectrum < spec_old)ifail = .true.                 !  energy of next decreases - probably an upper level super-deformed band
     if(ng == 0)ifail = .true.
     if(ng > 0)then
        complete = 1
        do j = 1, ng
           read(51,'(a)')line
           read(line,'(39x,i4,1x,f10.3,3(1x,e10.3))')nf,eg,pg,pe,icc
           if(complete_decay(nf) /= 1)complete = 0
           ncount = ncount + 1                           !  increment read counter
        end do
        if(complete == 1 .and. .not. ifail)complete_decay(i) = 1
     end if
     if(t12 > t12_isomer .and. .not. ifail)complete_decay(i) = 1
     spec_old = spectrum
  end do
 120  continue
!------------   Counted how many good states, so go back to beginning of this nucleus


  do i=1,ncount
     backspace(51)
  end do
  
!-----    Went through all discrete states. No check for first noncomplete decay, 
!-----    which specifies E_cut
!-----    Also, if All_gammas = .true. set up discrete states above E_cut if there
!-----    are any


  E_cut_set = .false.
  n_cut = 1
  ncc = 1
  do i = 2, nol
     if(complete_decay(i) == 1 )then
        ncc = ncc + 1
     end if
     if(complete_decay(i) == 0 .and. .not. E_cut_set)then
        n_cut = ncc
        E_cut_set = .true.
        if(.not. All_gammas) exit
     end if
  end do

  nucleus(icomp)%num_discrete = ncc
  nucleus(icomp)%ncut = n_cut
  allocate(nucleus(icomp)%state(ncc))
  
  ii = 0
  do i = 1, nol
     if(complete_decay(i) == 1)then
        ii = ii + 1
        if(ii <= ncc)then
           state_map(i) = ii
        end if
     end if
  end do
  

  if(ncc == 0) goto 10                               !   No known levels, treat as if it doesn't exist


  ii = 0
  do i = 1, nol
     char1='  '
     line(1:120) = ' '
     read(51,'(a)',err=220)line
     if(line(4:5) /= '  ')goto 220

     read(line,'(i3,1x,f10.6,1x,f5.1,i3,1x,(e10.2),i3,1x,a1)')     &
             ilab,spectrum,spin,par,t12,ng,jf
     if(t12 < 1.0d-15)t12=1.0d-15               !   min value for lifetime
!--------     Put state information into data array

     if(complete_decay(i) == 1 )then

        ii = ii + 1  
        if(.not. All_gammas .and. ii > n_cut)exit
        if(ii == 1)then
           nucleus(icomp)%ipar = (par+1)/2
           nucleus(icomp)%s1 = s1
           nucleus(icomp)%s2 = s2
           nucleus(icomp)%target_spin = target_spin
           nucleus(icomp)%target_ipar = ipar
        end if
        nucleus(icomp)%state(ii)%ensdf_index = i
        nucleus(icomp)%state(ii)%energy = spectrum
        nucleus(icomp)%state(ii)%spin = spin
        nucleus(icomp)%state(ii)%parity = par
        nucleus(icomp)%state(ii)%t12 = t12
!--------   Check if state is an isomer
        isomer = .false.
        if(t12 > t12_isomer)isomer = .true.
        nucleus(icomp)%state(ii)%isomer = .false.
        nucleus(icomp)%state(ii)%pop = 0.0d0
        if(isomer.and.ii > 1)nucleus(icomp)%state(ii)%isomer = .true.
!--------  Set up data for gamma branching
        nucleus(icomp)%state(ii)%nbranch = ng              !  # of branches
        if(ii > 1 .and. ng > 0)then
           allocate(nucleus(icomp)%state(ii)%ibranch(ng))   !  states
           allocate(nucleus(icomp)%state(ii)%branch(ng))    !  branching ratio
           allocate(nucleus(icomp)%state(ii)%p_gamma(ng))    !  branching ratio
           allocate(nucleus(icomp)%state(ii)%egamma(ng))    !  gamma-ray energy
           allocate(nucleus(icomp)%state(ii)%p_ic(ng))       !  internal conversion coefficient
           allocate(nucleus(icomp)%state(ii)%cs(ng))        !  Decay cross section
           btot = 0.0d0                                     !  norm for branching ratio
           bicc = 0.0d0
           do j = 1, ng                                     !  loop over branches
              read(51,'(a)')line
              read(line,'(39x,i4,1x,f10.3,3(1x,e10.3))')nf,eg,pg,pe,icc
              nff = state_map(nf)
              nucleus(icomp)%state(ii)%ibranch(j) = nff        !  state lable for final state
              nucleus(icomp)%state(ii)%branch(j) = pe          !  branching probablity
              nucleus(icomp)%state(ii)%p_ic(j) = icc            !  internal conversion coefficient
              nucleus(icomp)%state(ii)%egamma(j) = nucleus(icomp)%state(ii)%energy -  &
                                                   nucleus(icomp)%state(nff)%energy
              btot= btot + pe

           end do
           if(btot < 1.0d-4)then            !  compute based internal conversion coefficients
              xnorm = 0.0d0
              do j = 1, ng                                       !  Guess branches based on phase space arguments, rate ~ E_gamma**3
                 xnorm = xnorm + nucleus(icomp)%state(ii)%egamma(j)**3 
              end do
              do j = 1, ng                                       !  normalize gamma branches to unity
                 nff = nucleus(icomp)%state(ii)%ibranch(j) 
                 nucleus(icomp)%state(ii)%branch(j) = nucleus(icomp)%state(ii)%egamma(j)**3/xnorm
                 alpha = nucleus(icomp)%state(ii)%p_ic(j) 
                 nucleus(icomp)%state(ii)%p_gamma(j) = 1.0d0/(1.0d0 + alpha)
                 nucleus(icomp)%state(ii)%p_ic(j) = alpha/(1.0d0 + alpha)
              end do
           else
              xnorm = 0.0d0
              do j = 1, ng                                       !  normalize gamma branches to unity
                 xnorm = xnorm + nucleus(icomp)%state(ii)%branch(j)
              end do
              do j = 1, ng                                       !  normalize gamma branches to unity
                 nff = nucleus(icomp)%state(ii)%ibranch(j) 
                 nucleus(icomp)%state(ii)%branch(j) = nucleus(icomp)%state(ii)%branch(j)/xnorm
                 alpha = nucleus(icomp)%state(ii)%p_ic(j) 
                 nucleus(icomp)%state(ii)%p_gamma(j) = 1.0d0/(1.0d0 + alpha)
                 nucleus(icomp)%state(ii)%p_ic(j) = alpha/(1.0d0 + alpha)
              end do
           end if
        end if
        if(ii >= ncc)exit
     else
        if(ng > 0)then
           do j = 1, ng
              read(51,*)
           end do
        end if
     end if
  end do
 220  continue
 

  nucleus(icomp)%level_ecut = nucleus(icomp)%state(n_cut)%energy + 0.001d0
  close(51)
  nucleus(icomp)%state(1)%exit_lab = 1
  num_s = 1
!
!-----    Rare occasion where ground state spin is not known spin = -1
!-----    Make a guess and assume J = 0 if A=even or 0.5 if A= odd
!-----    Assume parity is even
!
  if(nucleus(icomp)%state(1)%spin < 0.0d0)then
     if(iand(ia,1) == 1)then             !   Assume if A odd J=0.5 or J=0 if A even
        spin = 0.5d0
         else
        spin = 0.0d0
     end if
     if(spin <= 1.0d-4)then
        s1 = spin + 0.5d0
        s2 = spin + 0.5d0
     else
        s1 = spin - 0.5
        s2 = spin + 0.5
     end if
     nucleus(icomp)%ipar = 0
     nucleus(icomp)%state(1)%spin = spin
     nucleus(icomp)%s1 = s1
     nucleus(icomp)%s2 = s2
     nucleus(icomp)%state(1)%parity = 1
     nucleus(icomp)%state(1)%t12 = 1000.
     nucleus(icomp)%state(1)%nbranch = 0              !  # of branches       
     nucleus(icomp)%state(1)%pop = 0.0
     nucleus(icomp)%state(1)%exit_lab = 1
     nucleus(icomp)%target_spin = spin
     nucleus(icomp)%target_ipar = 0
  end if



  do j = 2, nucleus(icomp)%num_discrete
     nucleus(icomp)%state(j)%exit_lab = 0
     if(nucleus(icomp)%state(j)%isomer)then
        num_s = num_s + 1
        nucleus(icomp)%state(j)%exit_lab = num_s
     end if
  end do

  if(allocated(complete_decay))deallocate(complete_decay)
  if(allocated(state_map))deallocate(state_map)

  return

 10   continue                                        !  nucleus not found set up with one level for ground state (unkown)
  if(iproc == 0)then
     write(6,'(''Nucleus '',i3,a2,'' not found'')')ia,symb
     write(6,*)'Setting up with no levels'
  end if
  nucleus(icomp)%num_discrete = 1
  allocate(nucleus(icomp)%state(nucleus(icomp)%num_discrete))
  nucleus(icomp)%state(1)%energy = 0.0d0
  if(iand(ia,1) == 1)then             !   Assume if A odd J=0.5 or J=0 if A even
     spin = 0.5d0
      else
     spin = 0.0d0
  end if
  if(spin <= 1.0d-4)then
     s1 = spin + 0.5d0
     s2 = spin + 0.5d0
  else
     s1 = spin - 0.5
     s2 = spin + 0.5
  end if
  nucleus(icomp)%ipar = 0
  nucleus(icomp)%state(1)%spin = spin
  nucleus(icomp)%s1 = s1
  nucleus(icomp)%s2 = s2
  nucleus(icomp)%state(1)%parity = 1
  nucleus(icomp)%state(1)%t12 = 1000.
  nucleus(icomp)%state(1)%nbranch = 0              !  # of branches       
  nucleus(icomp)%state(1)%pop = 0.0
  nucleus(icomp)%state(1)%exit_lab = 1
  nucleus(icomp)%target_spin = spin
  nucleus(icomp)%target_ipar = 0

  
  if(allocated(complete_decay))deallocate(complete_decay)
  if(allocated(state_map))deallocate(state_map)

  return
end subroutine get_spectrum
