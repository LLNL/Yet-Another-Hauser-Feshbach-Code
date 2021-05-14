!
!*******************************************************************************
!
subroutine parse_command(icommand, command, finish)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine parses commands controlling the calculation
!
!  Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        constants
!        nodeinfo
!        print_control
!        useful_data
!        nuclei
!        particles_def
!        directory_structure
!        pre_equilibrium_no_1
!
!     Subroutines:
!
!        parse_string
!        lower_case_word
!        char_logical
!        print_command_error
!        find_ZA
!        extract_ZA_data
!        incident_energy_setup
!        init_barrier_data
!        get_lev_den
!
!     External functions:
!
!        integer(kind=4) :: particle_index
!
!     MPI routines
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
!*******************************************************************************
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
   use particles_def
   use directory_structure
   use pre_equilibrium_no_1
!---------------------------------------------------------------------
   implicit none
!---------------------------------------------------------------------
   integer(kind=4), intent(inout) :: icommand
   character(len=132), intent(inout) :: command
!---------------------------------------------------------------------
   integer(kind=4) :: numw
   integer(kind=4) :: startw(66), stopw(66)
   logical :: finish
   character(len=50) :: read_file
!---------------------------------------------------------------------
   integer(kind=4) :: i, j, k, n, num, itemp_read
   integer(kind=4) :: nread
   integer(kind=4) :: istart, istop, ilast
   integer(kind=4) :: ibegin, iend
   real(kind=8) :: e_min,e_max
   character(len=132) :: temp_char

   integer(kind=4) :: ix, lx, gsf_type
   integer(kind=4) :: num_data, n_gsf

   integer(kind=4) :: io_error

   real(kind=8) :: e_temp

   logical :: interact
   logical :: logic_char
   logical :: read_error
   integer(kind=4) :: ndat
   integer(kind=4) :: iZ,IA
   real(kind=8) :: x(66)
   real(kind=8) :: emin, emax, estep
   integer(kind=4) :: nchar
   integer(kind=4) :: nw
   real(kind=8) :: Max_J
   real(kind=8) :: beta_2
   real(kind=8) :: sig2_perp, sig2_ax
   real(kind=8) :: xnorm
   integer(kind=int_64) :: one_int_64
   integer(kind=4) :: iseed
!-----------------   External functions ------------------------------
   integer(kind=4) :: particle_index
   real(kind=8) :: Kxl
!---------------------------------------------------------------------
   interact=.false.
   if(command(1:1) == '#' .or. command(1:1) == '!')return
   istart = 1
   istop = index(command,' ')-1

   call parse_string(command,numw,startw,stopw)

   if(numw <= 0)then
      finish = .true.
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   istart = startw(1)
   istop = stopw(1)
   if(command(startw(1):stopw(1)) == 'end')then
      finish=.true.
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'file')then
      
      out_file(1:132) = ' '

      nchar = stopw(2) - startw(2) + 1

      out_file(1:nchar) = command(startw(2):stopw(2))

      temp_char(1:132) = ' '
      temp_char(1:nchar) = out_file(1:nchar)
      call lower_case_word(nchar,temp_char(1:nchar))

      if(out_file(1:nchar) == 'none')print_output = .false.

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'fission')then
      icommand = icommand + 1
      if(numw < 2)then
!         write(6,*)'Error in input for option "fission"'
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      fission = logic_char

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

   if(command(startw(1):stopw(1)) == 'event_generator')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      event_generator = logic_char

      if(numw == 3)then
         call char_logical(command(startw(2):stopw(2)), logic_char, read_error)
         if(read_error)then
            call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
            return
         end if
         file_energy_index = logic_char
      end if

      if(event_generator)then
         dump_events = .true.
         binary_event_file = .true.
         print_libraries = .false.
         if(iproc == 0)then
            write(6,*)'***********************************************************'
            write(6,*)'----   YAHFC is being run in event-generator mode       ---'
            write(6,*)'----   Events will be printed to binary file            ---'
            write(6,*)'----   For an ascii file, use command                   ---'
            write(6,*)'----   "dump_events y a".                               ---'
            write(6,*)'----   Note that an answer "n/f/0" will be overridden   ---'
            write(6,*)'----   to "y/t/1" as events must be written to a file.  ---'
            write(6,*)'----   with this option.                                ---'
            write(6,*)'***********************************************************'
         end if
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'use_unequal_bins')then
      icommand = icommand + 1
      if(numw /= 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      use_unequal_bins = logic_char

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'xs_only')then
      icommand = icommand + 1
      if(numw /= 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      xs_only = logic_char

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'target')then
      icommand = icommand + 1

      nchar = stopw(2) - startw(2) + 1
      call find_ZA(nchar,command(startw(2):stopw(2)),iZ,iA)
      if(iZ == -1 .or. iA == -1)then
         if(numw /= 3)then
            if(iproc ==0)write(6,*)'Error specifying target - not enough data'
            call MPI_Abort(icomm,101,ierr)
         end if
         read(command(startw(2):stopw(2)),*)target%Z
         read(command(startw(3):stopw(3)),*)target%A
      else
         target%Z = iZ
         target%A = iA
      end if
      target%specified = .true.
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'target_state')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)target%istate
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'projectile')then
      icommand = icommand + 1
      if(numw == 2)then
         i = particle_index(command(startw(2):stopw(2)))
         if(i < 0)then
            if(iproc ==0)write(6,*)'Particle misidentified in command "projectile"'
            call MPI_Abort(icomm,101,ierr)
         end if
         projectile%particle_type = i
         projectile%Z = particle(i)%Z
         projectile%A = particle(i)%A
         projectile%specified = .true.
         particle(i)%do_dwba = .false.
!-----    Set values for pree-equilibrium model for incident neutrons and protons
         if(i == 1)then
            Preeq_V = 38.0d0
            Preeq_V1 = 12.0d0
            Preeq_K = 245.0d0
         elseif(i == 2)then
            Preeq_V = 38.0d0
            Preeq_V1 = 22.0d0
            Preeq_K = 450.0d0
         end if
!------   If user didn't specify either proj_e_file or proj_eminmax
!------   set up default incident energy grid

         if(.not. energy_input)call incident_energy_setup
!
!------   Make sure max_particle(iproj) >= 1
!
!              if(max_particle(i) >= 0)max_particle(i) = max(max_particle(i),1)
         if(max_particle(i) >= 0)max_particle(i) = -1
         return
      end if
      read(command(startw(2):stopw(2)),*)projectile%Z
      read(command(startw(3):stopw(3)),*)projectile%A
      projectile%particle_type = -1
      do i = 0, 6
         if(projectile%Z == particle(i)%Z.and.  &
            projectile%A == particle(i)%A)then
              projectile%particle_type = i
              projectile%specified = .true.
              particle(i)%do_dwba = .false.
!-----    Set values for pree-equilibrium model for incident neutrons and protons
              if(i == 1)then
                 Preeq_V = 38.0d0
                 Preeq_V1 = 12.0d0
                 Preeq_K = 245.0d0
              elseif(i == 2)then
                 Preeq_V = 38.0d0
                 Preeq_V1 = 22.0d0
                 Preeq_K = 450.0d0
              end if
!------   If user didn't specify either proj_e_file or proj_eminmax
!------   set up default incident energy grid
              if(.not. energy_input)call incident_energy_setup
!
!------   Make sure max_particle(iproj) >= 1
!
              if(max_particle(i) >= 0)max_particle(i) = -1
              return
         end if
      end do
!-------    Reserved to define a calculation populations input
!-------    Populations for a given energy are read in for a set of J^\pi values
      if(projectile%Z == -1 .and. projectile%A == 0)then
         if(numw < 4)then
            if(iproc ==0)write(6,*)'Not enough input in command projectile for population calculation'
            call MPI_Abort(icomm,101,ierr)
         end if
         projectile%particle_type = -2
         projectile%Z = 0
         projectile%A = 0
         j_pop_calc = .false.
         ex_pop_file(1:50)=' '
         ilast = stopw(4) - startw(4) + 1
         ex_pop_file(1:ilast) = command(startw(4):stopw(4))
         open(unit=8, file = ex_pop_file(1:ilast), status='old')
!----   First do some counting to make sure there is proper input
!----   Remove energies with no populations - just as a fail safe
         Pop_max_J = 0.0d0
         k = 0
         read(8,*)nread
         do j = 1, nread
            read(8,*)x(1),x(2)
            read(8,*)n
            if(n > 0) k = k + 1
            do i = 1, n
               read(8,*)
            end do
         end do
         num_pop_e = k
         rewind(8)
         if(num_pop_e > 0)then
            if(.not.allocated(Pop_data))allocate(Pop_data(1:num_pop_e))
            read(8,*)nread
            k = 0
            do j = 1, nread
!-----   Read line with energy, energy spread, and possibly the overall normalization          
               read(8,'(a)')temp_char
               call parse_string(temp_char,numw,startw,stopw)
               if(k == 1 .and. numw == 2)then
                  pop_calc_prob = .true.
               elseif(k == 1 .and. numw == 3)then
                  pop_calc_prob = .false.
               end if
!----  Read in number of J^pi populations to be read in
               read(8,*)n
               if(n > 0)then
!----   First check if the proper amount of data is on line for the two types of calculation
                  if(numw < 2)then
                     if(iproc == 0)then
                        write(6,*)'Error in setting energy for population calculation, not enough data for entry #',n
                        call MPI_Abort(icomm, 101, ierr)
                     end if
                  end if
!----   Error for failure to enter a normalization for any subsequent entry if
!----   first one is set, and therefore pop_calc_prob = .false.
                  if(numw == 2 .and. .not. pop_calc_prob)then
                     if(iproc == 0)then
                        write(6,*)'Error: No input for population normalization for entry #',n
                        call MPI_Abort(icomm, 101, ierr)
                     end if
                  end if
!----   Increment population counter and fill arrays 
                  k = k + 1
                  read(temp_char(startw(1):stopw(1)),*)Pop_data(k)%Ex_pop
                  read(temp_char(startw(2):stopw(2)),*)Pop_data(k)%dEx_pop
                  if(numw == 3)then
                     read(temp_char(startw(3):stopw(3)),*)Pop_data(k)%Norm_pop
                  else
                     Pop_data(k)%Norm_pop = 1.0d0
                  end if
!                  Pop_data(k)%Ex_pop = x(1)
!                  Pop_data(k)%dEx_pop = x(2)
                  Pop_data(k)%num_pop = n
                  if(.not.allocated(Pop_data(k)%j_pop))allocate(Pop_data(k)%j_pop(n))
                  if(.not.allocated(Pop_data(k)%par_pop))allocate(Pop_data(k)%par_pop(n))
                  if(.not.allocated(Pop_data(k)%bin_pop))allocate(Pop_data(k)%bin_pop(n))
                  xnorm = 0.0d0
                  do i = 1, Pop_data(k)%num_pop
                     read(8,*)Pop_data(k)%j_pop(i), Pop_data(k)%par_pop(i), Pop_data(k)%bin_pop(i)
                     xnorm = xnorm + Pop_data(k)%bin_pop(i)
                     if(Pop_data(k)%j_pop(i) > Pop_max_J)Pop_max_J = Pop_data(k)%j_pop(i)
                  end do
                  if(xnorm < 1.0d-6)then
                     if(iproc == 0)write(6,*)'Error!! The total population is too small < 1.0d-6'
                     call MPI_Abort(icomm,101,ierr)
                  end if
                  do i = 1, Pop_data(k)%num_pop
                     Pop_data(k)%bin_pop(i) = Pop_data(k)%bin_pop(i)/xnorm
                  end do     
               end if
            end do
            close(unit=8)
         else
            if(iproc ==0)write(6,*)'Something wrong with Population option, no populations are specified'
            call MPI_Abort(icomm,101,ierr)
         end if
         ex_set = .true.
         pop_calc = .true.

         projectile%specified=.true.
         projectile%num_e = num_pop_e
         if(.not.allocated(projectile%energy))allocate(projectile%energy(num_pop_e))
         do i = 1, num_pop_e
            projectile%energy(i) = Pop_data(i)%Ex_pop
         end do
         return
      end if
!
!----   Population with a single energy value and populations set up based on level density
!
      if(projectile%Z == -1 .and. projectile%A == -1)then
         if(numw < 4)then
            if(iproc ==0)write(6,*)'Not enough input in command projectile for LD population calculation'
            call MPI_Abort(icomm,101,ierr)
         end if
         projectile%particle_type = -2
         projectile%Z = 0
         projectile%A = 0
         pop_calc = .true.
         j_pop_calc = .true.
         num_pop_e = 1
         if(.not.allocated(Pop_data))allocate(Pop_data(1:num_pop_e))
         read(command(startw(4):stopw(4)),*)Pop_data(1)%Ex_pop
         Pop_data(1)%dEx_pop = 0.0d0
         if(numw == 5)read(command(startw(5):stopw(5)),*)Pop_data(1)%dEx_pop
         Pop_data(1)%Norm_pop = 1
         pop_calc_prob = .true. 
         if(numw == 6)then
            read(command(startw(6):stopw(6)),*)Pop_data(1)%Norm_pop
            pop_calc_prob = .false.
         end if

         Pop_data(1)%num_pop = 0

         ex_set = .true.

         projectile%specified=.true.
         projectile%num_e = num_pop_e

         if(.not.allocated(projectile%energy))allocate(projectile%energy(num_pop_e))
         projectile%energy(1) = Pop_data(1)%Ex_pop
         return
      end if
      if(iproc == 0)write(6,*)'Error in specifying projectile'
      call MPI_Abort(icomm,101,ierr)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'max_particle')then
      icommand = icommand + 1
      if(numw < 3)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
!      read(command(startw(2):stopw(2)),*)k
      k = particle_index(command(startw(2):stopw(2)))
      if(k < 1)then
          if(iproc == 0)write(6,*)'Particle missidentified in command "max_particle"'
          call MPI_Abort(icomm,101,ierr)
      end if
      read(command(startw(3):stopw(3)),*)i
!---   set max_particle(k). Note if k == iproj, max_particle >=1
!---   namely, it can't be zero, so prevent user from making an error
!---   that prevents anything from running.
      if(k == 1 .and. i >= 0 .and. iproc == 0)then
         write(6,*)'*************************************************'
         write(6,*)'*----  WARNING!!!! WARNING!!!! WARNINING!!! ----*'
         write(6,*)'*----  It is strongly recommended that      ----*'
         write(6,*)'*----  neutrons not be limited in the       ----*'
         write(6,*)'*----  number that can be emitted. It is    ----*'
         write(6,*)'*----  recommend to use "max_particle n -1" ----*'
         write(6,*)'*----  or the default, which limits the     ----*'
         write(6,*)'*----  number of neutrons by energy only.   ----*'  
         write(6,*)'*----  WARNING!!!! WARNING!!!! WARNINING!!! ----*'
         write(6,*)'*************************************************'
      end if
      if(k == projectile%particle_type)then
         if(i < 0)then
            max_particle(k) = -1
         elseif(i == 0)then
            max_particle(k) = 1
         else
            max_particle(k) = i
         end if
      else
         max_particle(k) = i
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'proj_e_file' .and. .not. ex_set)then
      icommand=icommand+1
      if(projectile%particle_type == -1)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read_file(1:50)=' '
      iend = 101
      ibegin = 101
      do n = 100, istop + 1, -1
         if(command(n:n) /= ' ')then
            iend = n
            exit
         end if
      end do
      do n = iend, istop + 1, -1
         if(command(n:n) == ' ')then
            ibegin = n + 1
            exit
         end if
      end do
      read_file = command(ibegin:iend)
      ilast = iend - ibegin + 1
      open(unit=8, file = read_file(1:ilast), status='old')
      num = 0
 7    read(8,*,end=8)
      num = num + 1
      goto 7
 8    rewind(8)
      allocate(projectile%energy(num))
      projectile%num_e = num
      do i = 1, num
         read(8,*)projectile%energy(i)
      end do
      close(unit=8)
      e_min = 1000000.
      e_max = -1000000.
      do i = 1, num
         if(projectile%energy(i) < e_min)    &
            e_min = projectile%energy(i)
            if(projectile%energy(i) > e_max) &
            e_max = projectile%energy(i)
      end do
      projectile%e_min = e_min
      projectile%e_max = e_max
      ex_set = .true.    
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'proj_eminmax' .and. .not. ex_set)then
      icommand = icommand + 1
      if(projectile%particle_type == -1)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      if(numw < 4)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)emin
      read(command(startw(3):stopw(3)),*)emax
      read(command(startw(4):stopw(4)),*)estep

      if(emin <= de)then
         emin = de/2.0d0
      else
         emin = emin - de/2.0d0
         num = int(emin/de)
         emin = real(num,kind=8)*de + de/2.0d0
         emax = emax - de/2.0d0
         num = int(emax/de)
         emax = real(num,kind=8)*de + de/2.0d0
      end if
      projectile%e_min = emin
      projectile%e_max = emax
      projectile%e_step = estep
      num = nint((projectile%e_max - projectile%e_min)/projectile%e_step) + 1
      allocate(projectile%energy(1:num))
      projectile%num_e = num
      do i = 1, num
         projectile%energy(i)=projectile%e_min + dfloat(i-1)*projectile%e_step
      end do
      projectile%e_max = projectile%energy(num)
      ex_set = .true.
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'cs_scale')then
      icommand = icommand + 1
      if(iproc == 0)then
         write(6,*)'*******************************************************'
         write(6,*)'*  WARNING!! WARNING!! WARNING!! WARNING!! WARNING!!  *' 
         write(6,*)'*       Option "cs_units" is no longer valid          *'
         write(6,*)'*       please use cs_units instead                   *'
         write(6,*)'*       "cs_units b" for calculations in barns        *'
         write(6,*)'*       "cs_units mb" for calculations in milibarns   *'
         write(6,*)'*******************************************************'
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'cs_units')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      if(command(startw(2):stopw(2)) == 'b')then
         cs_scale = 1.0d0
         cs_units = ' b'
      elseif(command(startw(2):stopw(2)) == 'mb')then
         cs_scale = 1.0d3
         cs_units = 'mb'
      else
         if(iproc == 0)then
            write(6,*)'Error in cross section scale factor'
            write(6,*)'Only barns (b) or millibarns (mb) are allowed'
            write(6,*)'Defaulting to barns'
         end if
         cs_scale = 1.0d0
         cs_units = ' b'
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_aparam')then
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%level_param(1)=x(1)
            nucleus(i)%fit_aparam = .false.
            if(nucleus(i)%fission .and..not. nucleus(i)%fission_read)then
               do j = 1, nucleus(i)%F_n_barr
                  nucleus(i)%F_Barrier(j)%level_param(1) = nucleus(i)%level_param(1)
               end do
            end if
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_d0')then
      icommand = icommand + 1

      ndat = 2
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%D0exp = x(1)
            nucleus(i)%dD0exp = x(2)
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_spin_cut')then
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(x(1) > 0.0)nucleus(i)%level_param(2) = x(1)
            if(nucleus(i)%fission .and..not. nucleus(i)%fission_read)then
               do j = 1, nucleus(i)%F_n_barr
                  nucleus(i)%F_Barrier(j)%level_param(2) = nucleus(i)%level_param(2)
               end do
            end if
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_sig_model')then
      icommand = icommand + 1


      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(x(1) > 0.0)nucleus(i)%level_param(13)=x(1)
            if(nucleus(i)%fission .and..not. nucleus(i)%fission_read)then
               do j = 1, nucleus(i)%F_n_barr
                  nucleus(i)%F_Barrier(j)%level_param(13) = nucleus(i)%level_param(13)
               end do
            end if
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_delta')then
      icommand = icommand + 1


      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%level_param(3) = x(1)
            nucleus(i)%pair_model = 3
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_ecut')then
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(x(1) < nucleus(i)%level_ecut)then
               nucleus(i)%level_ecut = x(1)
               nucleus(i)%level_param(7)=x(1)
               do j = 1, nucleus(i)%num_discrete
                  if(nucleus(i)%state(j)%energy > x(1))then
                     nucleus(i)%ncut = j
                     exit
                  end if
               end do
               if(.not. all_discrete_states)then
                  nucleus(i)%num_discrete = j - 1
               end if
               return
            end if
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_parity_fac')then
      icommand = icommand + 1

      ndat = 3
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(nint(x(1)) == 1) nucleus(i)%level_param(16) = 1.0d0
            if(nint(x(1)) == -1) nucleus(i)%level_param(16) = 2.0d0
            nucleus(i)%level_param(17) = x(2)
            nucleus(i)%level_param(18) = x(3)
            if(nucleus(i)%fission .and..not. nucleus(i)%fission_read)then
               do j = 1, nucleus(i)%F_n_barr
                  nucleus(i)%F_Barrier(j)%level_param(17) = nucleus(i)%level_param(17)
                  nucleus(i)%F_Barrier(j)%level_param(18) = nucleus(i)%level_param(18)
               end do
            end if
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_shell')then
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%level_param(4) = x(1)
            if(nucleus(i)%fission .and. .not. nucleus(i)%fission_read)then
               do j = 1, nucleus(i)%F_n_barr
                  nucleus(i)%F_Barrier(j)%level_param(4) = nucleus(i)%level_param(4)
               end do
            end if
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_gamma')then
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%level_param(5) = x(1)
            if(nucleus(i)%fission .and. .not. nucleus(i)%fission_read)then
               do j = 1, nucleus(i)%F_n_barr
                  nucleus(i)%F_Barrier(j)%level_param(5) = nucleus(i)%level_param(5)
               end do
            end if
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_ematch')then
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%level_param(6) = x(1)
            nucleus(i)%fit_ematch = .false.
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_fit_d0')then
      icommand = icommand + 1
      if(numw /= 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         nucleus(i)%fit_D0 = logic_char
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_fit_aparam')then
      icommand = icommand + 1
      if(numw /= 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      call char_logical(command(startw(2):stopw(2)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         nucleus(i)%fit_aparam = logic_char
      end do

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 't12_isomer')then
      icommand = icommand + 1
      if(numw /= 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)t12_isomer

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'fit_aparam')then
      icommand = icommand + 1
      if(iproc ==0)then
         write(6,*)
         write(6,*)'---------------------------------------------'
         write(6,*)'---  The command "fit_aparam" is obsolete ---'
         write(6,*)'---- Use "lev_fit_d0" or "lev_fit_aparam" ---'
         write(6,*)'---- See manual for input guidance        ---'
         write(6,*)'---------------------------------------------'
         write(6,*)'*********************************************'
         write(6,*)'**** Continuing as if:                    ***'
         write(6,*)'**** "lev_fit_d0 y"                       ***'
         write(6,*)'**** "lev_fit_aparam n"                   ***'
         write(6,*)'**** Fit to D0 for all nuclei by          ***'
         write(6,*)'**** adjusting the shell correction       ***'
         write(6,*)'**** parameter                            ***'
         write(6,*)'*********************************************'
         write(6,*)
      end if
      do i = 1, num_comp
         nucleus(i)%fit_D0 = .true.
         nucleus(i)%fit_aparam = .false.
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_fit_ematch')then
      icommand = icommand + 1
      if(numw /= 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         nucleus(i)%fit_ematch = logic_char
      end do

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_rot_enhance')then
      icommand = icommand + 1

      ndat = 3
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(nucleus(i)%lev_option < 2) return
            do j = 1, 3
               if(x(j) >= 0.0d0)nucleus(i)%rot_enh(j) = x(j)
            end do
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_vib_enhance')then
      icommand = icommand + 1

      ndat = 3
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(nucleus(i)%lev_option < 2) return
            nucleus(i)%level_param(11) = real(k,kind=8)
            do j = 1, 3
               if(x(j) >= 0.0d0)nucleus(i)%vib_enh(j) = x(j)
            end do
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_vib_enhance_mode')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)k
      do i = 1, num_comp
         if(nucleus(i)%lev_option < 2) return
         nucleus(i)%level_param(11) = real(k,kind=8)
         if(nucleus(i)%fission)call init_barrier_data(i)
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'initial_ke')then
      icommand = icommand+1
      if(.not. pop_calc)then
         if(iproc == 0)write(6,*)'This option is only compatible with a pupulation calculation'
         if(iproc == 0)write(6,*)'and is being ignored'
         return
      end if
      if(numw < 3)then
         if(iproc == 0)write(6,*)'Error in input for option "initial_ke"'
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      read(command(startw(3):stopw(3)),*)x(2)
      Init_Kinetic_Energy = x(1)
      dInit_Kinetic_Energy = x(2)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'quasi_elastic')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)quasi_e
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'beta_2')then
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%beta(2) = x(1)
            nucleus(i)%rot_enh(4) = (1.0d0 + x(1)/3.0d0)
            nucleus(i)%rot_enh(5) = sqrt(pi/2.0d0)*(1.0d0 - 2.0d0*x(1)/3.0d0)
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_option')then          !   global setting of this parameter
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)j

      if(j > 2)then
         if(iproc == 0)write(6,*)'Bad input for "lev_option"'
         call MPI_Abort(icomm,101,ierr)
      end if

      do i = 1, num_comp
         iZ = nucleus(i)%Z
         iA = nucleus(i)%A
         if(iA > 20)then
            nucleus(i)%lev_option = j
            nucleus(i)%level_param(9) = real(j,kind=8)
            if(j == 2) then
               nucleus(i)%level_param(10) = 1.0d0
               nucleus(i)%level_param(11) = 1.0d0
            end if
            call get_lev_den(data_path,len_path,                         &
                             symb(iZ),iZ,iA,i)
            if(nucleus(i)%fission .and. .not. nucleus(i)%fission_read)   &
               call init_barrier_data(i)

         elseif(j > 1 .and. iproc == 0)then
            if(iproc == 0)write(6,*)'Warning, A is too small and lev_option > 1 is dangerous'
            if(iproc == 0)write(6,*)'Keeping default option 0'
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!

   if(command(startw(1):stopw(1)) == 'lev_nuc_option')then          !   setting of this parameter for nucleus iZ,iA
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      j = nint(x(1),kind=4)

      if(j > 2)then
         if(iproc == 0)write(6,*)'Bad input for "lev_nuc_option"'
         call MPI_Abort(icomm,101,ierr)
      end if
      if(iA <= 20 .and. j > 1)then
         if(iproc == 0)write(6,*)'Warning, A is too small and lev_option > 1 is dangerous'
         if(iproc == 0)write(6,*)'Keeping default option 0'
         return
      end if
      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%lev_option = j
            nucleus(i)%level_param(9) = X(1)
            if(j == 2) then
               nucleus(i)%level_param(10) = 1.0d0
               nucleus(i)%level_param(11) = 1.0d0
            end if
            call get_lev_den(data_path,len_path,                         &
                             symb(iZ),iZ,iA,i)
            if(nucleus(i)%fission .and. .not. nucleus(i)%fission_read)   &
               call init_barrier_data(i)

            return
         end if         
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_nuc_fit_D0')then
      icommand = icommand + 1

      ndat = 0
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      nw = nw + 1
      if(numw < nw)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(nw):stopw(nw)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%fit_D0 = logic_char
            return
         end if
      end do

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_nuc_fit_aparam')then
      icommand = icommand + 1


      ndat = 0
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      nw = nw + 1
      if(numw < nw)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(nw):stopw(nw)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%fit_aparam = logic_char
            return
         end if
      end do

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_nuc_fit_ematch')then
      icommand = icommand + 1

      ndat = 0
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      nw = nw + 1
      call char_logical(command(startw(nw):stopw(nw)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%fit_ematch = logic_char
            return
         end if
      end do

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'e1_param')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 4
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      j = nint(X(1),kind=4)
      X(1) = X(2)
      X(2) = X(3)
      X(3) = X(4)

      if(j > max_num_gsf)then
         if(iproc ==0)write(6,*)'Bad input for e1_param, requesting more than max components =',max_num_gsf
         call MPI_Abort(icomm,101,ierr)
      end if
      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(nucleus(i)%E1_default)nucleus(i)%E1_default = .false.
            nucleus(i)%EL_mode(1)%gsf(j)%er = X(1)
            nucleus(i)%EL_mode(1)%gsf(j)%gr = X(2)
            nucleus(i)%EL_mode(1)%gsf(j)%sr = X(3)
            nucleus(i)%fit_gamma_gamma = .false.
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'el_param')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 5
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      lx = nint(X(1),kind=4)
      j = nint(X(2),kind=4)
      X(1) = X(3)
      X(2) = X(4)
      X(3) = X(5)

      if(j > max_num_gsf)then
         if(iproc ==0)write(6,*)'Bad input for el_param, requesting more than max components =',max_num_gsf
         call MPI_Abort(icomm,101,ierr)
      end if
      if(lx > e_l_max)then
         if(iproc ==0)write(6,*)'Bad input for el_param, requesting more than max components =',max_num_gsf
         call MPI_Abort(icomm,101,ierr)
      end if
      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(nucleus(i)%E1_default)nucleus(i)%E1_default = .false.
!            nucleus(i)%er_E1(j) = x(1)
!            nucleus(i)%gr_E1(j) = x(2)
!            nucleus(i)%sr_E1(j) = x(3)
            nucleus(i)%EL_mode(lx)%gsf(j)%er = X(1)
            nucleus(i)%EL_mode(lx)%gsf(j)%gr = X(2)
            nucleus(i)%EL_mode(lx)%gsf(j)%sr = X(3)
            nucleus(i)%fit_gamma_gamma = .false.
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'ml_param')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 5
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      lx = nint(X(1),kind=4)
      j = nint(X(2),kind=4)
      X(1) = X(3)
      X(2) = X(4)
      X(3) = X(5)

      if(j > max_num_gsf)then
         if(iproc ==0)write(6,*)'Bad input for ml_param, requesting more than max components =',max_num_gsf
         call MPI_Abort(icomm,101,ierr)
      end if
      if(lx > m_l_max)then
         if(iproc ==0)write(6,*)'Bad input for el_param, requesting more than max components =',max_num_gsf
         call MPI_Abort(icomm,101,ierr)
      end if
      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(nucleus(i)%E1_default)nucleus(i)%E1_default = .false.
            nucleus(i)%ML_mode(lx)%gsf(j)%er = X(1)
            nucleus(i)%ML_mode(lx)%gsf(j)%gr = X(2)
            nucleus(i)%ML_mode(lx)%gsf(j)%sr = X(3)
            nucleus(i)%fit_gamma_gamma = .false.
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_num_barrier')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      j = nint(X(1),kind=4)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(iproc ==0)then
               write(6,*)'Default Fission Barriers are overridden for ', nucleus(i)%Label
               write(6,*)'New defaults established.'
               write(6,*)'Barrier heights = 6.0 MeV, hbw = 0.6 MeV, symmetric level densities'
               write(6,*)'Use input commands to specify all Fission parameters!!!'
            end if
            nucleus(i)%fission = .true.
            nucleus(i)%F_n_barr = j
            if(allocated(nucleus(i)%F_Barrier))deallocate(nucleus(i)%F_Barrier)
            allocate(nucleus(i)%F_Barrier(nucleus(i)%F_n_barr))
!-----    Reset barrier data. Barrier symmetries are unknown and will be reset in 
!-----    sybroutine init_barrier_data. However, the user should manually reset 
            nucleus(i)%F_Barrier(1:nucleus(i)%F_n_barr)%symmetry = 0
            call init_barrier_data(i)
!-----
            return
         end if
      end do
      if(iproc ==0)write(6,*)'Nucleus not found: f_Num_Barrier',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_barrier')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 3
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      j = nint(X(1),kind=4)
      X(1) = X(2)
      X(2) = X(3)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'too many barriers for F_Barrier',' Z = ',iZ,' A = ',iA
               call MPI_Abort(icomm,101,ierr)
            end if
            nucleus(i)%F_Barrier(j)%barrier = x(1)
            nucleus(i)%F_Barrier(j)%hbw = x(2)
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: f_Barrier',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    damping factor for fission barriers
!    F_barrier = F_barrier*x(1)*exp(-x(3)**2*(Ex-x())**2)
!   
   if(command(startw(1):stopw(1)) == 'f_barrier_damp')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 3
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      j = nint(X(1),kind=4)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc ==0)write(6,*)'too many barriers for F_Barrier',' Z = ',iZ,' A = ',iA
               call MPI_Abort(icomm,101,ierr)
            end if
            nucleus(i)%F_Barrier(j)%barrier_damp(2) = x(2)
            nucleus(i)%F_Barrier(j)%barrier_damp(3) = x(3)
!-------   Make damping factor = 1.0 at Ex = 0.0
            x(1) = exp((x(3)*x(2))**2)
            nucleus(i)%F_Barrier(j)%barrier_damp(1) = x(1)
            return
         end if
      end do
      if(iproc ==0)write(6,*)'Nucleus not found: F_Barrier_damp',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    Set symmetry for barrier
!   
   if(command(startw(1):stopw(1)) == 'f_barrier_symmetry')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      j = nint(X(1),kind=4)
      nw = nw + 1
      if(numw < nw)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'Error: requesting too many barriers for F_Barrier for nucleus Z = ',&
                                       iZ,' A = ',iA
               call MPI_Abort(icomm,101,ierr)
            end if
            if(command(startw(nw):stopw(nw)) == 's' .or. command(startw(nw):stopw(nw)) == '1')then
               nucleus(i)%F_Barrier(j)%symmetry = 1
               nucleus(i)%F_Barrier(j)%level_param(10) = real(nucleus(i)%F_Barrier(j)%symmetry,kind=8)
               return
            end if
            if(command(startw(nw):stopw(nw)) == 'lr-a' .or. command(startw(nw):stopw(nw)) == '2')then
               nucleus(i)%F_Barrier(j)%symmetry = 2
               nucleus(i)%F_Barrier(j)%level_param(10) = real(nucleus(i)%F_Barrier(j)%symmetry,kind=8)
               if(j == 1 .and. iproc == 0)write(6,*)'WARNING!!!! ----  Setting first barrier to left-right asymmetric'
               return
            end if
            if(command(startw(nw):stopw(nw)) == 'ta-lr' .or. command(startw(nw):stopw(nw)) == '3')then
               nucleus(i)%F_Barrier(j)%symmetry = 3
               nucleus(i)%F_Barrier(j)%level_param(10) = real(nucleus(i)%F_Barrier(j)%symmetry,kind=8)
               return
            end if
            if(command(startw(nw):stopw(nw)) == 'ta-nlr' .or. command(startw(nw):stopw(nw)) == '4')then
               nucleus(i)%F_Barrier(j)%symmetry = 4
               nucleus(i)%F_Barrier(j)%level_param(10) = real(nucleus(i)%F_Barrier(j)%symmetry,kind=8)
               if(j == 1 .and. iproc == 0)write(6,*)'WARNING!!!! ----  Setting first barrier to triaxial no left-right asymmetry'
               return
            end if
         end if
      end do
      if(iproc ==0)then
         write(6,*)'Nucleus not found: f_barrier_symmetry',' Z = ',iZ,' A = ',iA
         write(6,*)'Keeping the default value for this nucleus and barrier'
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_ecut')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 2
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      j = nint(X(1),kind=4)
      X(1) = X(2)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'too many barriers for F_ecut',' Z = ',iZ,' A = ',iA
               call MPI_Abort(icomm,101,ierr)
            end if
            nucleus(i)%F_Barrier(j)%ecut = x(1)
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: F_ecut',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_aparam')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 2
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      j = nint(X(1),kind=4)
      X(1) = X(2)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'too many barriers for f_lev_aparam',' Z = ',iZ,' A = ',iA
               call MPI_Abort(icomm,101,ierr)
            end if
            nucleus(i)%F_Barrier(j)%level_param(1) = x(1)
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: f_lev_aparam',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_spin')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 2
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      j = nint(X(1),kind=4)
      X(1) = X(2)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'too many barriers for F_lev_spin',' Z = ',iZ,' A = ',iA
               call MPI_Abort(icomm,101,ierr)
            end if
            nucleus(i)%F_Barrier(j)%level_param(2) = x(1)
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_spin',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(command(startw(1):stopw(1)) == 'f_lev_delta')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 2
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      j = nint(X(1),kind=4)
      X(1) = X(2)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'too many barriers for F_lev_delta',' Z = ',iZ,' A = ',iA
               call MPI_Abort(icomm,101,ierr)
            end if
            nucleus(i)%F_Barrier(j)%level_param(3) = x(1)
            nucleus(i)%F_Barrier(j)%level_param(6) = 2.5 + 150./real(iA,kind=8) + x(1)
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_delta',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_shell')then
      icommand = icommand+1

      ndat = 2
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      j = nint(X(1),kind=4)
      X(1) = X(2)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'too many barriers for F_lev_shell',' Z = ',iZ,' A = ',iA
               call MPI_Abort(icomm,101,ierr)
            end if
            nucleus(i)%F_Barrier(j)%level_param(4) = x(1)
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_shell',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_gamma')then
      icommand = icommand + 1

      ndat = 2
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      j = nint(X(1),kind=4)
      X(1) = X(2)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'too many barriers for F_lev_gamma'
               call MPI_Abort(icomm,101,ierr)
            end if
            nucleus(i)%F_Barrier(j)%level_param(5) = x(1)
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_gamma',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_rot_enhance')then
      icommand = icommand + 1

      ndat = 4
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      j = nint(X(1),kind=4)
      X(1) = X(2)
      X(2) = X(3)
      X(3) = X(4)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'Error in F_lev_rot_enhance index > # of barriers'
               call MPI_Abort(icomm,101,ierr)
            end if
            do k = 1,3
               if(x(k) >= 0.0d0)nucleus(i)%F_Barrier(j)%rot_enh(k) = x(k)
            end do
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_rot_enhance',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_vib_enhance')then
      icommand = icommand + 1

      ndat = 4
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      j = nint(X(1),kind=4)
      X(1) = X(2)
      X(2) = X(3)
      X(3) = X(4)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'Error in F_vib_vib_enhance index > # of barriers'
               call MPI_Abort(icomm,101,ierr)
            end if
            do k = 1, 3
               if(x(k) >= 0.0d0)nucleus(i)%F_Barrier(j)%vib_enh(k) = x(k)
            end do
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_vib_enhance',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_ematch')then
      icommand = icommand + 1

      ndat = 2
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      j = nint(X(1),kind=4)
      X(1) = X(2)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*)'too many barriers for F_lev_ematch'
               call MPI_Abort(icomm,101,ierr)
            end if
            if(x(1) <= nucleus(i)%F_Barrier(j)%level_param(3) + 0.25)then
               if(iproc == 0)then
                  write(6,*)'****************************************************'
                  write(6,*)'*           WARNING!!!                             *'
                  write(6,*)'*  Attempting to set ematch less than delta + 0.25 *'
                  write(6,*)'*  in nucleus = ',i,' barrier # ',j
                  write(6,*)'*  A value this small can lead to numerical issues *'
                  write(6,*)'*  Resetting to delta + 0.25 MeV                   *'
                  write(6,*)'****************************************************'
               end if
               nucleus(i)%F_Barrier(j)%level_param(6) =        &
                  nucleus(i)%F_Barrier(j)%level_param(6)  + 0.25d0
            else
               nucleus(i)%F_Barrier(j)%level_param(6) = x(1)
            end if
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_ematch',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_beta2')then          !   global setting of this parameter
      icommand = icommand + 1

      ndat = 2
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      j = nint(X(1),kind=4)
      X(1) = X(2)

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(j > nucleus(i)%F_n_barr)then
               if(iproc == 0)write(6,*) 'WARNING -- too many barriers for F_beta_2 in nucleus ',i
               call MPI_Abort(icomm,101,ierr)
            end if
            beta_2 = x(1)
            nucleus(i)%F_Barrier(j)%beta_2 = beta_2
            sig2_perp = (1.0d0 + beta_2/3.0d0)
            sig2_ax = sqrt(pi/2.0d0)*(1.0d0 - 2.0d0*beta_2/3.0d0)
            nucleus(i)%F_Barrier(j)%rot_enh(4) = sig2_perp
            nucleus(i)%F_Barrier(j)%rot_enh(5) = sig2_ax
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: f_beta2',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_barr_levels')then
      icommand = icommand + 1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      j = nint(X(1),kind=4)
      nw = nw + 1

      read_file(1:50) = ' '
      ilast = stopw(nw) - startw(nw) + 1

      read_file(1:ilast) = command(startw(nw):stopw(nw))

      do i = 1, num_comp
         if(nucleus(i)%Z == iZ .and. nucleus(i)%A == iA)then
            open(unit=8, file = read_file(1:ilast), status='old')
            num = 0
 9          read(8,*,end=10)
            num = num + 1
            goto 9
 10         rewind(8)
            nucleus(i)%F_Barrier(j)%num_discrete = num
            if(.not.allocated(nucleus(i)%F_barrier(j)%state_e))then
               allocate(nucleus(i)%F_barrier(j)%state_e(num))
            else
               deallocate(nucleus(i)%F_barrier(j)%state_e)
               allocate(nucleus(i)%F_barrier(j)%state_e(num))
            end if
            if(.not.allocated(nucleus(i)%F_barrier(j)%state_j))then
               allocate(nucleus(i)%F_barrier(j)%state_j(num))
            else
               deallocate(nucleus(i)%F_barrier(j)%state_j)
               allocate(nucleus(i)%F_barrier(j)%state_j(num))
            end if
            if(.not.allocated(nucleus(i)%F_barrier(j)%state_pi))then
               allocate(nucleus(i)%F_barrier(j)%state_pi(num))
            else
               deallocate(nucleus(i)%F_barrier(j)%state_pi)
               allocate(nucleus(i)%F_barrier(j)%state_pi(num))
            end if
            do k = 1, num
               read(8,*)nucleus(i)%F_Barrier(j)%state_e(k),             &
                         nucleus(i)%F_Barrier(j)%state_j(k),            &
                         nucleus(i)%F_Barrier(j)%state_pi(k)
            end do
            close(unit=8)
         end if
      end do

      if(iproc == 0)write(6,*)'Nucleus not found: F_Barr_levels',' Z = ',iZ,' A = ',iA
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'ran_seed')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)iseed
      iseed_32 = iseed
      iseed_64 = int(iseed,kind=int_64)
      one_int_64 = 1_int_64
      if(iand(iseed_64,one_int_64) /= one_int_64)iseed_64 = iseed_64 + one_int_64
      if(iseed_64 > 0)iseed_64 = -iseed_64
      if(iand(iseed_32,1) /= 1)iseed_32 = iseed_32 + 1
      if(iseed_32 > 0)iseed_32 = -iseed_32      
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'num_mc_samp')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)j
      if(j < 1)then
         if(iproc == 0)write(6,*)'Invalid value using default - num_mc_samp = ',num_mc_samp
      end if
      num_mc_samp = j
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_model')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)j
      if(j > 1)then
         if(iproc == 0)write(6,*)'Invalid choice using default'
         return
      end if
      PREEQ_Model = j
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_m2_c1')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      M2_C1 = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_m2_c2')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      M2_C2 = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_m2_c3')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      M2_C3 = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_m2_rpp')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      M2_Rpp = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_m2_rnn')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      M2_Rnn = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_m2_rpn')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      M2_Rpn = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_m2_rnp')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      M2_Rnp = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_well_v')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      Preeq_V = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_gam_fact')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      Preeq_gam_fact = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_well_v1')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      Preeq_V1 = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!----   changed 12/15/20 to accomdate change in well function
!----   don't need values for both protons and neutrons since 
!----   it depends on the projectile, and not emitted particle
!----   Before didn't have correct values for both protons or neutrons
   if(command(startw(1):stopw(1)) == 'preeq_well_k')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      Preeq_K = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_g_div')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      Preeq_g_div = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_g_a')then
      icommand = icommand+1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)i
      if(i == 0)then
         Preeq_g_a = .false.
      elseif(i == 1)then
         Preeq_g_a = .true.
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_pair_model')then
      icommand = icommand+1
      if(numw < 3)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)preeq_pair_model
      if(preeq_pair_model > 2)preeq_pair_model = 1
      if(preeq_pair_model == 2)then
         read(command(startw(3):stopw(3)),*)preeq_delta
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_fwell')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)preeq_fwell
      if(preeq_fwell < 0 .or. preeq_fwell > 2)preeq_fwell=1
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_analytic')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      analytic_preeq = logic_char
      return 

   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'wf_model')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)WF_model
      if(WF_model > 1)WF_model = 1
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'delta_e')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)de
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'e1_model')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)E1_model
      do i = 1, num_comp
         nucleus(i)%e1_model = E1_model
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'pair_model')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)j
      pair_model = j
      if(j <= 2)then
         do i = 1, num_comp
            iZ = nucleus(i)%Z
            iA = nucleus(i)%A
            nucleus(i)%pair_model = j
            call get_lev_den(data_path,len_path,               &
                             symb(iZ),iZ,iA,i)
         end do
      else
         if(iproc == 0)write(6,*)'Invalid input for option: pair_model'
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'output_mode')then
      icommand = icommand + 1
      if(iproc == 0)then
          write(6,*)'Warning - Ignoring "output_mode" as it is longer a valid input command'
      end if
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'max_j_allowed')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)itemp_read
      if(pop_calc)then
         if(itemp_read > max_J_allowed)max_J_allowed = itemp_read
      else
         max_J_allowed = itemp_read
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'out_gammas_vs_e')then
      icommand = icommand + 1
      if(iproc == 0)then
          write(6,*)'Warning - "out_gammas_vs_e" is longer a valid input command'
          write(6,*)'Use "track_gammas y" instead. Continuing as if "track_gammas" was called'
      end if
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      track_gammas = logic_char         
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'track_gammas')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      track_gammas = logic_char         

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'prob_cut')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      prob_cut = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'trans_p_cut')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      trans_p_cut = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'trans_e_cut')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      trans_e_cut = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'rho_cut')then
      icommand=icommand+1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      rho_cut = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'optical_potential')then
      icommand = icommand + 1
      if(numw < 3)then
         if(iproc == 0)write(6,*)'Error in input for option "optical_potential"'
         if(iproc == 0)write(6,*)'Not enough input to set optical potential'
         return
      end if
!----   particle type is in 2nd word
      istart = startw(2)
      istop = stopw(2)
      k = particle_index(command(startw(2):stopw(2)))
      if(k < 1)then
         if(iproc == 0)write(6,*)'Particle misidentified in command "optical_potential"'
         call MPI_Abort(icomm,101,ierr)
      end if
!----   optical potential type is in 3rd word
      istart = startw(3)
      istop = stopw(3)
      read(command(istart:istop),*)j
      if(j <= particle(k)%max_opt_pot)then
         particle(k)%opt_pot = j
         particle(k)%opt_pot_set = .true.
         return
      else
         if(iproc == 0)write(6,*)'Error in input for option "optical_potential"'
         if(iproc == 0)write(6,*)'Requesting optical potential that is not available'
         if(iproc == 0)write(6,'(1x,''for particle '',a8)')particle(k)%name
         return
      end if
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'om_option')then
      icommand = icommand + 1
      if(numw < 3)then
         if(iproc == 0)write(6,*)'Error in input for option "om_option"'
         if(iproc == 0)write(6,*)'Not enough input to set optical potential'
         return
      end if
!----   particle type is in 2nd word
      istart = startw(2)
      istop = stopw(2)
      k = particle_index(command(startw(2):stopw(2)))
      if(k <= 0 .or. k > 6)then
         if(iproc == 0)write(6,*)'Error in input for option "om_option"'
         if(iproc == 0)write(6,*)'Illegal particle type to set optical potential'
         return
      end if
!----   optical potential type is in 3rd word
      read(command(startw(3):stopw(3)),*)particle(k)%om_option
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'explicit_channels')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      explicit_channels = .false.

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      explicit_channels = logic_char

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'dump_events')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      if(.not. event_generator)then
         dump_events = logic_char
      else
         if(iproc == 0)then
            write(6,*)'*****************************************************************'
            write(6,*)'----   Event generator mode has been selected, therefore,    ----'
            write(6,*)'----   continuing with "dump_events = .true.".               ----'
            write(6,*)'----   This option cannot be changed.                        ----'
            write(6,*)'*****************************************************************'
         end if
      end if


      if(dump_events)then
         nchar = stopw(3) - startw(3) + 1
         call lower_case_word(nchar,command(startw(3):stopw(3)))
         if(command(startw(3):stopw(3)) == 'b')then
            binary_event_file = .true.
         elseif(command(startw(3):stopw(3)) == 'a')then
            binary_event_file = .false.
         else
            binary_event_file = .true.
            if(iproc == 0)write(6,*)'error specifying event file type. Will write to binary file'
         end if
         return
      end if

      if(iproc == 0)write(6,*)'dump_events = ',dump_events
      return
   end if

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'e_l_max')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)j
      e_l_max = j
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'm_l_max')then
      icommand=icommand+1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)j
      m_l_max = j
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'fit_gamma_gamma')then
      icommand=icommand+1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      fit_gamma_gamma = logic_char
      do i = 1, num_comp
         nucleus(i)%fit_gamma_gamma = fit_gamma_gamma
      end do   

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'nuc_fit_gamma_gamma')then
      icommand=icommand+1

      ndat = 0
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      nw = nw + 1
      if(numw < nw)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(nw):stopw(nw)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%fit_gamma_gamma = logic_char
            return
         end if
      end do


      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'do_dwba')then
      icommand=icommand+1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      k = projectile%particle_type

      particle(k)%do_dwba = logic_char

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'set_gamma_gamma')then
      icommand=icommand+1

      ndat = 1
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            nucleus(i)%Gamma_g_exp = x(1)
            nucleus(i)%dGamma_g_exp = -1.0d0
            return
         end if
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_den_read')then
      icommand=icommand+1

      ndat = 0
      call extract_ZA_data(command, numw, startw, stopw, ndat,                 &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      if(numw - nw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(nucleus(i)%lev_den_both)then
               if(iproc == 0)then
                  write(6,*)'Error with command "lev_den_read"'
                  write(6,*)'Reading both parities has already specified'
               end if
               call MPI_Abort(icomm, 101, ierr)
            end if              
            nucleus(i)%lev_den_read = .true.
            nucleus(i)%fit_D0 = .false.
            nucleus(i)%fit_aparam = .false.
            ilast = stopw(nw+2) - startw(nw+2) + 1
            if(command(startw(nw+1):stopw(nw+1)) == 'b' .or.                    &
               command(startw(nw+1):stopw(nw+1)) == 'B')then
               nucleus(i)%lev_den_both = .true.
                   nucleus(i)%lev_den_file(0)(1:ilast) = command(startw(nw+2):stopw(nw+2))
            elseif(command(startw(nw+1):stopw(nw+1)) == '0' .or.                &
                   command(startw(nw+1):stopw(nw+1)) == 'n' .or.                &
                   command(startw(nw+1):stopw(nw+1)) == 'N' .or.                &
                   command(startw(nw+1):stopw(nw+1)) == '-')then
                   nucleus(i)%lev_den_file(0)(1:ilast) = command(startw(nw+2):stopw(nw+2))
            elseif(command(startw(nw+1):stopw(nw+1)) == '1' .or.                &
                   command(startw(nw+1):stopw(nw+1)) == 'p' .or.                &
                   command(startw(nw+1):stopw(nw+1)) == 'P' .or.                &
                   command(startw(nw+1):stopw(nw+1)) == '+')then
               nucleus(i)%lev_den_file(1)(1:ilast) = command(startw(nw+2):stopw(nw+2))
            end if
            return
         end if
      end do
      if(iproc == 0)write(6,*)'Nucleus not found: lev_den_read',' Z = ',iZ,' A = ',iA
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'max_j_allowed')then
      icommand=icommand+1
      read(command(startw(2):stopw(2)),*)Max_J_allowed
      if(Max_J_allowed < 20)then
         Max_J_allowed = 20
         if(iproc == 0)write(6,*)'Attempt to set max_J_allowed < 20, set to 20'
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'all_gammas')then
      icommand=icommand+1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      if(iproc == 0)then
         write(6,*)
         write(6,*)'---------------------------------------------'
         write(6,*)'---  The command "all_gammas" is obsolete ---'
         write(6,*)'---- Use "all_discrete_states"            ---'
         write(6,*)'---- See manual for input guidance        ---'
         write(6,*)'---------------------------------------------'
         write(6,*)'*********************************************'
         write(6,*)'**** Continuing as if:                    ***'
         write(6,*)'**** "all_discrete_states y"              ***'
         write(6,*)'*********************************************'
         write(6,*)
      end if

      all_discrete_states = .false.

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      all_discrete_states = logic_char

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'all_discrete_states')then
      icommand=icommand+1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      all_discrete_states = .false.

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      all_discrete_states = logic_char

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'primary_decay')then
      icommand=icommand+1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      primary_decay = .false.

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      primary_decay = logic_char

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'optical_code')then
      icommand = icommand + 1
      if(iproc == 0)then
         write(6,*)'*******************************************************'
         write(6,*)'*  WARNING!! WARNING!! WARNING!! WARNING!! WARNING!!  *' 
         write(6,*)'*  The option: "optical_code" is no longer available  *'
         write(6,*)'*  The default optical model code is FRESCO           *'
         write(6,*)'*  This calculation will proceed using FRESCO         *'
         write(6,*)'*******************************************************'
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'cc_file')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      istart = startw(2)
      istop = stopw(2)
      local_cc_file(1:istop-istart+1) = command(istart:istop)
      inquire(file = local_cc_file(1:istop-istart+1), exist = exist_cc_file)
      if(.not. exist_cc_file)then
         if(iproc == 0)write(6,*)'Input file with coupled-channels input data does not exist'
         if(iproc == 0)write(6,*)'Input different name than',local_cc_file(1:istop-istart+1)
         stop
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'cc_scale')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)cc_scale
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'scale_elastic')then
      icommand = icommand + 1
      if(projectile%particle_type /= 1)then
         scale_elastic = .false.
         if(iproc == 0)then
            write(6,*)'*********************************************************'
            write(6,*)'* ---- WARNING!!! WARNING!!! WARNING!!! WARNING!!! ---- *'
            write(6,*)'* Attempting to scale elastic cross section for an      *'
            write(6,*)'* incident particle with electric charge. This is ill   *'
            write(6,*)'* defined and will automatically be overriden.          *'
            write(6,*)'* Elastic scattering will not be rescaled.              *'
            write(6,*)'* You should remove command "scale_elastic" from run    *'
            write(6,*)'* file for charged particles to remove this warning     *'
            write(6,*)'*********************************************************'
         end if
         return
      end if
      if(numw < 4)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      elastic_scale = 1.0d0
      elastic_shift = 1.0d0
      elastic_damp = 0.0d0
      if(numw < 4)then
         if(iproc == 0)write(6,*)'Wrong input for command scale_elastic'
         if(iproc == 0)write(6,*)'command: scale_elastic  elastic_scale, elastic_shift, elastic_damp'
         stop
      end if 
      scale_elastic = .true.   
      istart = startw(2)
      istop = stopw(2)
      read(command(istart:istop),*)elastic_scale
      istart = startw(3)
      istop = stopw(3)
      read(command(istart:istop),*)elastic_shift
      istart = startw(4)
      istop = stopw(4)
      read(command(istart:istop),*)elastic_damp
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'fiss_max_j')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)Fiss_Max_J
      do i = 1, num_comp
          Max_J = Fiss_Max_J
          iA = nucleus(i)%A
          if(iand(iA,1) == 1)Max_J = Max_J + 0.5
          do j = 1, nucleus(i)%F_n_barr
             nucleus(i)%F_Barrier(i)%Max_J = Max_J
          end do
      end do
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'trans_avg_l')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      if(numw /= 2)then
         if(iproc == 0)write(6,*)'Wrong input for command trans_avg_l'
         call MPI_Abort(icomm,101,ierr)
      end if
      trans_avg_l = .false.

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      trans_avg_l = logic_char

      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'num_theta_angles' )then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)num_theta_angles
      if(num_theta_angles < 0)num_theta_angles = 1
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'trans_norm' )then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      T_norm = x(1)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'max_num_gsf' )then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)x(1)
      max_num_gsf = nint(x(1),kind=4)
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'fresco_shape')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      read(command(startw(2):stopw(2)),*)ifresco_shape
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'biased_sampling')then
      icommand = icommand + 1
      if(numw < 2)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      biased_sampling = .false.

      call char_logical(command(startw(2):stopw(2)),logic_char,read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      biased_sampling = logic_char

      return 
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'start')then
      icommand = icommand + 1
      if(iproc == 0)write(6,*)'Command copied from from Fission Barrier file, ignoring'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'verbose_output')then
      icommand = icommand + 1
      if(numw < 2)then
!         write(6,*)'Error in input for option "fission"'
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      verbose = logic_char

      return      !   if it gets here a proper input wasn't given so keep default
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'write_outfile')then
      icommand = icommand + 1
      if(numw < 2)then
!         write(6,*)'Error in input for option "fission"'
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      call char_logical(command(startw(2):stopw(2)), logic_char, read_error)

      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if

      print_output = logic_char

      return      !   if it gets here a proper input wasn't given so keep default
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! read_el_gsf  Z  A  L  f/s  file  norm
   if(command(startw(1):stopw(1)) == 'read_el_gsf')then
      icommand = icommand + 1

      ndat = 0
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      if(numw < nw + 3)then
         if(iproc == 0)write(6,*)'Error, not enough input in call to "read_e_l_gsf"'
         call MPI_Abort(icomm,101,ierr)
      end if
      ix = 0
!-----   Find the multipolarity 
      read(command(startw(nw+1):stopw(nw+1)),*)lx
      if(lx < 1 .or. lx > e_l_max)then               !   error in multipolarity
         if(iproc == 0)then
            write(6,*)'Error in multipolarity in "read_e_l_gsf" L = ',lx
            write(6,*)'Maximum L value allowed for Electric transition is ',e_l_max
         end if
         call MPI_Abort(icomm,101,ierr)
      end if
      gsf_type = -1
!-----   Find out if we are reading a strength function: f/0 
!-----   or is it a photoabsorption cross section: s/1
      if(command(startw(nw+2):stopw(nw+2)) == 'F' .or.                      &
         command(startw(nw+2):stopw(nw+2)) == 'f' .or.                      &
         command(startw(nw+2):stopw(nw+2)) == '0')then
           gsf_type = 0
      elseif(command(startw(nw+2):stopw(nw+2)) == 'S' .or.                  &
         command(startw(nw+2):stopw(nw+2)) == 's' .or.                      &
         command(startw(nw+2):stopw(nw+2)) == '1')then
           gsf_type = 1
      end if
      if(gsf_type == -1)then
         if(iproc == 0)write(6,*)'Error specifying if reading strength function or cross section'
         call MPI_Abort(icomm,101,ierr)
      end if

!-----   The file name
      ilast = stopw(nw+3) - startw(nw+3) + 1
      read_file(1:50) = ' '
      read_file(1:ilast) = command(startw(nw+3):stopw(nw+3))

!-----   check if normalization is to be read in
      xnorm = 1.0d0
      if(numw == nw + 4)read(command(startw(nw+4):stopw(nw+4)),*)xnorm

!-----   Find the nucleus for data to be read in
      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(.not. nucleus(i)%EL_mode(lx)%gsf_read)                       &
               nucleus(i)%EL_mode(lx)%num_gsf = 0
            nucleus(i)%EL_mode(lx)%gsf_read = .true.
            nucleus(i)%EL_mode(lx)%num_gsf =                                           &
               nucleus(i)%EL_mode(lx)%num_gsf + 1
            nucleus(i)%fit_gamma_gamma = .false.    !  turn off fit to gamma_gamma
            n_gsf = nucleus(i)%EL_mode(lx)%num_gsf
            nucleus(i)%EL_mode(lx)%gsf(n_gsf)%gsf_type = gsf_type
!------   Check that there aren't too many strength functions - maximum is 3
            if(n_gsf > max_num_gsf)then
               if(iproc == 0)write(6,*)'Error, attempting to read in too many strength functions, max = ', &
               max_num_gsf
               call MPI_Abort(icomm,101,ierr)
            end if
            open(unit=10,file = read_file(1:ilast), status = 'old')
!------   Find out how much data we will be reading in
            num_data = 0
            io_error = 0
            do while(io_error == 0)
               temp_char(1:132) = ' '
               read(10,'(a)',iostat = io_error)temp_char
               if(temp_char(1:1) == '#' .or. temp_char(1:1) == '!')cycle
               num_data = num_data + 1
            end do
            num_data = num_data - 1
!-----   Set number of data for this gsf
            nucleus(i)%EL_mode(lx)%gsf(n_gsf)%num_data = num_data
!-----   Allocate energy and gsf arrays
            if(.not. allocated(nucleus(i)%EL_mode(lx)%gsf(n_gsf)%e_gsf_r))             &
               allocate(nucleus(i)%EL_mode(lx)%gsf(n_gsf)%e_gsf_r(num_data))
            nucleus(i)%EL_mode(lx)%gsf(n_gsf)%e_gsf_r(1:num_data) = 0.0d0
            if(.not. allocated(nucleus(i)%EL_mode(lx)%gsf(n_gsf)%gsf_r))               &
               allocate(nucleus(i)%EL_mode(lx)%gsf(n_gsf)%gsf_r(num_data))
            nucleus(i)%EL_mode(lx)%gsf(n_gsf)%gsf_r(1:num_data) = 0.0d0
!-----   read in energy and gsf data from file
            rewind(10)
            j = 1
            io_error = 0
            do while(j <= num_data .and. io_error == 0)
               temp_char(1:132) = ' '
               read(10,'(a)',iostat = io_error)temp_char
               if(temp_char(1:1) == '#' .or. temp_char(1:1) == '!')cycle
               read(temp_char,*)                                                    &
                  nucleus(i)%EL_mode(lx)%gsf(n_gsf)%e_gsf_r(j),                     &
                  nucleus(i)%EL_mode(lx)%gsf(n_gsf)%gsf_r(j)
!----   If it is photo-absorption, convert to strength function
               if(gsf_type == 1)then
                  e_temp = nucleus(i)%EL_mode(lx)%gsf(n_gsf)%e_gsf_r(j)
                  if(e_temp > 0.0d0)then
                     nucleus(i)%EL_mode(lx)%gsf(n_gsf)%gsf_r(j) =                   &
                         nucleus(i)%EL_mode(lx)%gsf(n_gsf)%gsf_r(j)*                &
                         Kxl(lx)/e_temp**(2*lx+1)
                  else
                     nucleus(i)%EL_mode(lx)%gsf(n_gsf)%gsf_r(j) = 0.0d0
                  end if
               end if
               j = j + 1
            end do
!-----   If io_error == 1, then something went wrong with the read
            if(io_error == 1)then
               if(iproc == 0)write(6,*)'Error reading in energy and strength function data'
               call MPI_Abort(icomm,101,ierr)
            end if
!-----   Now check if the grid energy is >= maximum excitation energy
            if(nucleus(i)%EL_mode(lx)%gsf(n_gsf)%e_gsf_r(num_data) <                  &
               nucleus(i)%Ex_max)then
               write(6,*)'Error in gamma strength function read in file ',read_file(1:ilast) 
               write(6,*)'Excitation energy is less than maximum excitation energy for this nucleus'
               write(6,*)'It must be greater than or equal to ', nucleus(i)%Ex_max
               call MPI_Abort(icomm,101,ierr)
            end if
            nucleus(i)%EL_mode(lx)%gsf(n_gsf)%gsf_norm = xnorm
            return
         end if
      end do
!----    If we get here, the nucleus wasn't found
      if(iproc == 0)then
          write(6,*)'************************************************************************************'
          write(6,*)'Warning nucleus not found for command "read_e_l_gsf". Looking for Z,A = ',iZ,iA
          write(6,*)'Possibly the incident energy was changed in this run, and the nucleus is no longer'
          write(6,*)'in the decay chain. Check'
          write(6,*)'************************************************************************************'
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! read_ml_gsf  Z  A  L  f/s  file  norm
   if(command(startw(1):stopw(1)) == 'read_ml_gsf')then
      icommand = icommand + 1

      ndat = 0
      call extract_ZA_data(command, numw, startw, stopw, ndat,         &
                           iZ, iA, X, nw, read_error)
      if(read_error)then
         call print_command_error(stopw(1)-startw(1)+1,command(startw(1):stopw(1)))
         return
      end if
      if(numw < nw + 3)then
         if(iproc == 0)write(6,*)'Error, not enough input in call to "read_m_l_gsf"'
         call MPI_Abort(icomm,101,ierr)
      end if
      ix = 0
!-----   Find the multipolarity 
      read(command(startw(nw+1):stopw(nw+1)),*)lx
      if(lx < 1 .or. lx > e_l_max)then               !   error in multipolarity
         if(iproc == 0)then
            write(6,*)'Error in multipolarity in "read_m_l_gsf" L = ',lx
            write(6,*)'Maximum L value allowed for Electric transition is ',m_l_max
         end if
         call MPI_Abort(icomm,101,ierr)
      end if
      gsf_type = -1
!-----   Find out if we are reading a strength function: f/0 
!-----   or is it a photoabsorption cross section: s/1
      if(command(startw(nw+2):stopw(nw+2)) == 'F' .or.                  &
         command(startw(nw+2):stopw(nw+2)) == 'f' .or.                  &
         command(startw(nw+2):stopw(nw+2)) == '0')then
           gsf_type = 0
      elseif(command(startw(nw+2):stopw(nw+2)) == 'S' .or.              &
         command(startw(nw+2):stopw(nw+2)) == 's' .or.                  &
         command(startw(nw+2):stopw(nw+2)) == '1')then
           gsf_type = 1
      end if
      if(gsf_type == -1)then
         if(iproc == 0)write(6,*)'Error specifying if reading strength function or cross section'
         call MPI_Abort(icomm,101,ierr)
      end if

!-----   The file name
      ilast = stopw(nw+3) - startw(nw+3) + 1
      read_file(1:50) = ' '
      read_file(1:ilast) = command(startw(nw+3):stopw(nw+3))

!-----   check if normalization is to be read in
      xnorm = 1.0d0
      if(numw == nw + 4)read(command(startw(nw+4):stopw(nw+4)),*)xnorm

!-----   Find the nucleus for data to be read in
      do i = 1, num_comp
         if(iZ == nucleus(i)%Z .and. iA == nucleus(i)%A)then
            if(.not. nucleus(i)%ML_mode(lx)%gsf_read)                       &
               nucleus(i)%ML_mode(lx)%num_gsf = 0
            nucleus(i)%ML_mode(lx)%gsf_read = .true.
            nucleus(i)%fit_gamma_gamma = .false.    !  turn off fit to gamma_gamma
            nucleus(i)%ML_mode(lx)%num_gsf =                                           &
               nucleus(i)%ML_mode(lx)%num_gsf + 1
            n_gsf = nucleus(i)%ML_mode(lx)%num_gsf
            nucleus(i)%ML_mode(lx)%gsf(n_gsf)%gsf_type = gsf_type
!------   Check that there aren't too many strength functions - maximum is 3
            if(n_gsf > max_num_gsf)then
               if(iproc == 0)write(6,*)'Error, attempting to read in too many strength functions, max = ', &
               max_num_gsf
               call MPI_Abort(icomm,101,ierr)
            end if
            open(unit=10,file = read_file(1:ilast), status = 'old')
!------   Find out how much data we will be reading in
            num_data = 0
            io_error = 0
            do while(io_error == 0)
               temp_char(1:132) = ' '
               read(10,'(a)',iostat = io_error)temp_char
               if(temp_char(1:1) == '#' .or. temp_char(1:1) == '!')cycle
               num_data = num_data + 1
            end do
            num_data = num_data - 1
!-----   Set number of data for this gsf
            nucleus(i)%ML_mode(lx)%gsf(n_gsf)%num_data = num_data
!-----   Allocate energy and gsf arrays
            if(.not. allocated(nucleus(i)%ML_mode(lx)%gsf(n_gsf)%e_gsf_r))             &
               allocate(nucleus(i)%ML_mode(lx)%gsf(n_gsf)%e_gsf_r(num_data))
            nucleus(i)%ML_mode(lx)%gsf(n_gsf)%e_gsf_r(1:num_data) = 0.0d0
            if(.not. allocated(nucleus(i)%ML_mode(lx)%gsf(n_gsf)%gsf_r))               &
               allocate(nucleus(i)%ML_mode(lx)%gsf(n_gsf)%gsf_r(num_data))
            nucleus(i)%ML_mode(lx)%gsf(n_gsf)%gsf_r(1:num_data) = 0.0d0
!-----   read in energy and gsf data from file
            rewind(10)
            j = 1
            io_error = 0
            do while(j <= num_data .and. io_error == 0)
               temp_char(1:132) = ' '
               read(10,'(a)',iostat = io_error)temp_char
               if(temp_char(1:1) == '#' .or. temp_char(1:1) == '!')cycle
               read(temp_char,*)                                                    &
                  nucleus(i)%ML_mode(lx)%gsf(n_gsf)%e_gsf_r(j),                     &
                  nucleus(i)%ML_mode(lx)%gsf(n_gsf)%gsf_r(j)
               if(gsf_type == 1)then
                  e_temp = nucleus(i)%ML_mode(lx)%gsf(n_gsf)%e_gsf_r(j)
                  if(e_temp > 0.0d0)then
                     nucleus(i)%ML_mode(lx)%gsf(n_gsf)%gsf_r(j) =                   &
                         nucleus(i)%ML_mode(lx)%gsf(n_gsf)%gsf_r(j)*                &
                         Kxl(lx)/e_temp**(2*lx+1)
                  else
                     nucleus(i)%ML_mode(lx)%gsf(n_gsf)%gsf_r(j) = 0.0d0
                  end if
               end if
               j = j + 1
            end do
!-----   If io_error == 1, then something went wrong with the read
            if(io_error == 1)then
               if(iproc == 0)write(6,*)'Error reading in energy and strength function data'
               call MPI_Abort(icomm,101,ierr)
            end if
!-----   Now check if the grid energy is >= maximum excitation energy
            if(nucleus(i)%ML_mode(lx)%gsf(n_gsf)%e_gsf_r(num_data) <                  &
               nucleus(i)%Ex_max)then
               write(6,*)'Error in gamma strength function read in file ',read_file(1:ilast) 
               write(6,*)'Excitation energy is less than maximum excitation energy for this nucleus'
               write(6,*)'It must be greater than or equal to ', nucleus(i)%Ex_max
               call MPI_Abort(icomm,101,ierr)
            end if
            nucleus(i)%ML_mode(lx)%gsf(n_gsf)%gsf_norm = xnorm
            return
         end if
      end do
!----    If we get here, the nucleus wasn't found
      if(iproc == 0)then
          write(6,*)'************************************************************************************'
          write(6,*)'Warning nucleus not found for command "read_e_l_gsf". Looking for Z,A = ',iZ,iA
          write(6,*)'Possibly the incident energy was changed in this run, and the nucleus is no longer'
          write(6,*)'in the decay chain. Check'
          write(6,*)'************************************************************************************'
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!-----   If we get here, then the command wasn't recognized
   if(iproc == 0)then
      write(6,*)'*******************************************************'
      write(6,*)'*  WARNING!! WARNING!! WARNING!! WARNING!! WARNING!!  *' 
      write(6,*)'*          COMMAND NOT RECOGNIZED  --   CHECK         *'
      write(6,*)command(1:stopw(numw))
      write(6,*)'*******************************************************'
   end if

   return
end subroutine parse_command
!
!*******************************************************************************
!
subroutine parse_string(input,numw,startw,stopw)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine breaks up a character(len=132) string into separate
!    words by specifying start and stop positions in the string for each word
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
   character(len=132), intent(in) :: input
   integer(kind=4), intent(out) :: numw
   integer(kind=4), intent(out) :: startw(66), stopw(66)
   integer(kind=4) k
   logical word
   startw(1:66) = 0
   stopw(1:66) = 0
   numw = 0
   word = .false.
   k = 1
 1 if(k > 1 .and. (input(k:k) == '!' .or. input(k:k) == '#'))then
      goto 4
   elseif(input(k:k) /= ' ')then
      if(.not.word)then
         numw = numw + 1
         word = .true.
         startw(numw) = k
         if(k == 132)then
            stopw(numw) = k
            goto 4
         end if
      end if
   elseif(input(k:k) == ' ')then
      if(word)then
         word = .false.
         stopw(numw) = k - 1
      end if
      if(k == 132) goto 4
   end if
   k = k + 1
   if(k > 132)goto 4
   goto 1
 4 continue
   return
end subroutine parse_string
!
!*******************************************************************************
!
subroutine lower_case_word(nchar,word)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine changes all lower case characters in a string 
!    to upper case
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
   integer(kind=4), intent(in) :: nchar
   character(len=nchar), intent(inout) :: word
!-------------------------------------------------------------
   integer(kind=4) :: i, j
   do i = 1, nchar
      j = iachar(word(i:i))
      if(j >= 65 .and. j <= 90)then       !  Upper case: 65-90
         j = j + 32                       !  Lower case: 97-122
         word(i:i) = achar(j)
      end if
   end do
   return
end subroutine lower_case_word
!
!*******************************************************************************
!
subroutine upper_case_word(nchar,word)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine changes all lower case characters in a string 
!    to upper case
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
!*******************************************************************************
!
   implicit none
   integer(kind=4), intent(in) :: nchar
   character(len=nchar), intent(inout) :: word
!-------------------------------------------------------------
   integer(kind=4) :: i, j
   do i = 1, nchar
      j = iachar(word(i:i))
      if(j >= 97 .and. j <= 122)then      !  Lower case: 97-122
         j = j - 32                       !  Upper case: 65-90
         word(i:i) = achar(j)
      end if
   end do
   return
end subroutine upper_case_word
!
!*******************************************************************************
!
logical function is_char_number(char)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function checks if a character variable is a number
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
!*******************************************************************************
!
   implicit none
   character(len=1), intent(in) :: char
!------------------------------------------------
   integer(kind=4) :: i
   i = iachar(char)
   is_char_number = .false.
   if(i >= 48 .and. i <= 57)is_char_number = .true.
   return
end function is_char_number
!
!*******************************************************************************
!
logical function is_char_letter(char)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function checks if a character variable is a letter - upper or 
!    lower case
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
!*******************************************************************************
!
   implicit none
   character(len=1), intent(in) :: char
!------------------------------------------------
   integer(kind=4) :: i
   i = iachar(char)
   is_char_letter = .false.
   if(i >= 60 .and. i <= 90)is_char_letter = .true.
   if(i >= 97 .and. i <= 122)is_char_letter = .true.
   return
end function is_char_letter
!
!*******************************************************************************
!
integer(kind=4) function rank_commands(command)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function assigns a rank to each command to define its
!    relative importance. Zero being the "highest"
!    rank. This adds flexibility to easily change
!    the relative ranking, i.e., add more granularity
!    if needed later.
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
   character(len=132), intent(in) :: command
!----------------------------------------------------------------
   integer(kind=4) :: numw
   integer(kind=4) :: startw(66), stopw(66)
!-----   Start with base assumption of rank of 20 for all commands
!-----   and set to different value based individual importance
   rank_commands = 20
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   call parse_string(command,numw,startw,stopw)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'all_gammas')then
      rank_commands = 0 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'all_discrete_states')then
      rank_commands = 0 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 't12_isomer')then
      rank_commands = 0 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'explicit_channels')then
      rank_commands = 0 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'event_generator')then
      rank_commands = 0 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'e_l_max')then
      rank_commands = 0 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'm_l_max')then
      rank_commands = 0 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'max_num_gsf')then
      rank_commands = 0 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'projectile')then
      rank_commands = 1
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'max_particle')then
      rank_commands = 2
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'target')then
      rank_commands = 3
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'target_state')then
      rank_commands = 5
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'delta_e')then
      rank_commands = 5
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'max_j_allowed')then
      rank_commands = 5
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'file')then
      rank_commands = 6
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'proj_e_file')then
      rank_commands = 6
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'proj_eminmax')then
      rank_commands = 6
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_option')then
      rank_commands = 7 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_aparam')then
      rank_commands = 7
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'fit_gamma_gamma')then
      rank_commands = 7
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_fit_d0')then
      rank_commands = 8
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_fit_aparam')then
      rank_commands = 8
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'pair_model')then
      rank_commands = 8 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'pair_vib_enhance_mode')then
      rank_commands = 8 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'nuc_fit_gamma_gamma')then
      rank_commands = 9
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_nuc_option')then
      rank_commands = 9 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_nuc_fit_d0')then
      rank_commands = 10
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_nuc_fit_aparam')then
      rank_commands = 10
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_delta')then
      rank_commands = 11
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_shell')then
      rank_commands = 11
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_gamma')then
      rank_commands = 11
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_ematch')then
      rank_commands = 11
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_ecut')then
      rank_commands = 11
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_d0')then
      rank_commands = 11
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_spin_cut')then
      rank_commands = 11
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_sig_model')then
      rank_commands = 11
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_parity_fac')then
      rank_commands = 11
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'lev_vib_enhance_mode')then
      rank_commands = 12 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):startw(1)+3) == 'num_theta_angles')then
      rank_commands = 13
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):startw(1)+3) == 'optical_potential')then
      rank_commands = 13
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'preeq_model')then
      rank_commands = 13 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'wf_model')then
      rank_commands = 13 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'do_dwba')then
      rank_commands = 13 
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):startw(1)+3) == 'om_option')then
      rank_commands = 11
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_num_barrier')then
      rank_commands = 22
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_barrier')then
      rank_commands = 23
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_delta')then
      rank_commands = 24
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_shell')then
      rank_commands = 25
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_gamma')then
      rank_commands = 26
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_spin')then
      rank_commands = 27
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_ematch')then
      rank_commands = 28
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_barrier_symmetry')then
      rank_commands = 29
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_rot_enhance')then
      rank_commands = 30
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_vib_enhance')then
      rank_commands = 30
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):startw(1)+1) == 'f_')then
      rank_commands = 31
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):startw(1)+1) == 'lev_den_read')then
      rank_commands = 35
      return
   end if
!
!------   Didn't get caught in command filter so it keeps default rank
!
   return
end function rank_commands
!
!*******************************************************************************
!
integer function particle_index(char)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns particle index 'k' given the particle index
!    or particle label
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
   implicit none
   character(len=1) char
   particle_index = -1
   if(char == '0' .or. char == 'g')particle_index = 0 
   if(char == '1' .or. char == 'n')particle_index = 1 
   if(char == '2' .or. char == 'p')particle_index = 2 
   if(char == '3' .or. char == 'd')particle_index = 3 
   if(char == '4' .or. char == 't')particle_index = 4 
   if(char == '5' .or. char == 'h')particle_index = 5 
   if(char == '6' .or. char == 'a')particle_index = 6 
   return
end function particle_index
!
!*******************************************************************************
!
subroutine find_ZA(nchar, word, iZ, iA)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns the Z and A of a nucleus input with mass and
!    element symbol, e.g., for word = 238U, iZ = 92 and iA = 238
!    If the element is not found, either iZ = -1 or iA = -1 is returned
!    signifying an error
!
!   Dependencies:
!
!     Modules:
!
!        useful_data
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
   use useful_data
   implicit none
   integer(kind=4), intent(in) :: nchar
   character, intent(in) :: word(nchar)
   integer(kind=4), intent(out) :: iZ, iA
!-----------------------------------------------------
   character(len=3) :: number
   character(len=2) :: Element
   integer(kind=4) :: i, n, nn, nc
   logical is_char_number
   logical is_char_letter
!-------------------------------------------------------------------------
!-----   Element symbols are found in module useful_data in the 
!-----   file modules.f90
!-------------------------------------------------------------------------

   iZ = -1
   iA = -1

!-------------------------------------------------------------------------
!----   Check if mass number is at the front or the end
!-------------------------------------------------------------------------
   n = 0
   do i = 1, nchar
      if(is_char_number(word(i)))then
         n = n + 1
         number(n:n) = word(i)
      end if
   end do

   if(n == 0)return

   read(number(1:n),*)iA

   nc = 0
   do i = 1, nchar
      if(is_char_letter(word(i)))nc = nc + 1
   end do

   Element(1:2) = '  '

   nn = 2 - nc
   nc = 0
   do i = 1, nchar
      if(is_char_letter(word(i)))then
         nn = nn + 1
         nc = nc + 1
         Element(nn:nn) = word(i)
         if(nc == 1)then
           call upper_case_word(1,Element(nn:nn))
         else
           call lower_case_word(1,Element(nn:nn))
         end if
       end if
   end do

   do i = 1, num_elements
      if(Element == symb(i))then
         iZ = i
         exit
      end if
   end do

   return

end subroutine find_ZA
!
!*******************************************************************************
!
subroutine extract_ZA_data(command, numw, startw, stopw, ndat,     &
                           iZ, iA, X, nw, read_error)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine reads the commands to extract iZ, iA, and ndat elements
!    of data stored in X(ndat). if insufficient data is provided to be read
!    read_error is returned as .false.
!    The integer nw is the number of words in the line that were read
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
!*******************************************************************************
!
   implicit none
   character(len=132), intent(in) :: command
   integer(kind=4), intent(in) :: numw
   integer(kind=4), dimension(numw), intent(in) :: startw, stopw
   integer(kind=4), intent(in) :: ndat
   integer(kind=4), intent(out) :: iZ, iA
   real(kind=8), dimension(ndat), intent(out) :: X
   integer(kind=4), intent(out) :: nw
   logical, intent(out) :: read_error
!--------------------------------------------------------------------
   integer(kind=4) :: nchar
   integer(kind=4) :: n
!--------------------------------------------------------------------
   read_error = .false.
   nchar = stopw(2) - startw(2) + 1

   call find_ZA(nchar,command(startw(2):stopw(2)),iZ,iA)
   if(iZ == -1 .or. iA == -1)then
      nw = 3
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
   else
      nw = 2
   end if
   if(ndat == 0)then
      return
   end if      
   if(numw < ndat + nw)then
      read_error = .true.
      return
   end if
   do n = 1, ndat
      nw = nw + 1
      read(command(startw(nw):stopw(nw)),*)X(n)
   end do

   return

end subroutine extract_ZA_data
!
!*******************************************************************************
!
subroutine char_logical(char, logic_char, read_error)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function checks a character and returns .true. if char = 'y', 't', or '1';
!    .false. if char = 'n', 'f', or '0'. Otherwise, it returns read_error = .true.
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
!*******************************************************************************
!
   implicit none
   character(len=1), intent(inout) :: char
   logical, intent(out) :: logic_char
   logical, intent(out) :: read_error
!-------------------------------------------------------------------
   logic_char = .false.
   read_error = .false.
!-------------------------------------------------------------------
   call lower_case_word(1, char) 
   if(char == 'y' .or. char == 't' .or. char == '1')then
      logic_char = .true.
   elseif(char == 'n' .or. char == 'f' .or. char == '0')then
      logic_char = .false.
   else
      read_error = .true.
   end if
   return
end subroutine char_logical
!
!*******************************************************************************
!
subroutine print_command_error(nchar,word)
!
!*******************************************************************************
!
!  Discussion:
!
!    This Subroutine prints that an error in the input for a given command has 
!    has been detected
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
   use nodeinfo
   integer(kind=4), intent(in) :: nchar
   character(len=132), intent(in) :: word
!-------------------------------------------------------------------------------
   character(len=132) :: error_line
!-------------------------------------------------------------------------------
   error_line(1:132) = ' '
   error_line(1:27) = 'Error in input for option "'
   ierr_last = 27
   if(iproc == 0)write(6,*)error_line(1:ierr_last)//word(1:nchar)//'"'
   return
end subroutine print_command_error
!
!*****************************************************************************80
!
subroutine incident_energy_setup
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to setup the incident energy grid in the event that the
!    user doesnot specify it with either proj_e_file or proj_eminmax 
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
   use nuclei
!-----------------------------------------------------
   integer(kind=4) :: i, num
   num = 221
   if(.not. allocated(projectile%energy))allocate(projectile%energy(num))
   projectile%energy(1) = 1.0d-4
   projectile%energy(2) = 1.0d-3
   do i = 3, 22
      projectile%energy(i) = projectile%energy(i-1) + 5.0d-3
   end do
   do i = 23, 221
      projectile%energy(i) = projectile%energy(i-1) + 1.0d-1
   end do
   projectile%num_e = num
   projectile%e_min = projectile%energy(1)
   projectile%e_max = projectile%energy(221)
   ex_set = .true.
   return
end subroutine incident_energy_setup



