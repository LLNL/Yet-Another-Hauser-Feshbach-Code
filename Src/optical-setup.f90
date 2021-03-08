!
!*******************************************************************************
!
subroutine optical_setup(data_path, len_path, iproj, itarget,                  &
                         istate, de, num_comp, Ang_L_max)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine sets up and runs the optical model program
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
   use nuclei
   use options
   use Channel_info
   use particles_def
   implicit none
!-------------------------------  External data
   character(len=200), intent(in) :: data_path
   integer(kind=4), intent(in) :: len_path
   integer(kind=4), intent(in) :: iproj, itarget, istate
   real(kind=8), intent(in) :: de
   integer(kind=4), intent(in) :: num_comp
   integer(kind=4), intent(in) :: Ang_L_max

!-------------------------------   Internal data
   character(len=80) :: bigblank        ! 80 blank character
   character(len=80) :: tco_file
   integer(kind=4) :: len_tco
   character(len=80) :: prn_file        ! file.tco Transmission coefficient file
   integer(kind=4) :: len_prn             ! length of tc_file 
   logical tco_exist, Ang_exist
   integer(kind=4), allocatable :: cc_index(:)
   logical optical_run 
!-------------------------
   integer(kind=4) :: i, j, k, l, n, idummy, in, isp
   integer(kind=4) :: ie
   integer(kind=4) :: ii
   integer(kind=4) :: ll
   integer(kind=4) :: maxl, minl
   integer(kind=4) :: num_l_lines
   integer(kind=4) :: num_e_lines
   integer(kind=4) :: lmax, nume

   integer(kind=4) :: iZ, iA

   integer(kind=4) :: jtarget, jstate
   real(kind=8) :: energy
   logical match
   integer(kind=4) :: max_nume
   integer(kind=4) :: nbin
   real(kind=8) :: xdummy
   real(kind=8) :: xj, spin

   real(kind=8) :: ee

   real(kind=8) :: J_gs, J_state, K_state
   real(kind=8) :: J_min, J_max
   real(kind=8) :: K_band
   integer(kind=4) :: Ix_min, Ix_max, Ix
   real(kind=8) :: zzero, sum, xk_factor
   logical :: check_dwba

   integer(kind=4) :: n_elastic

   real(kind=8) :: E_in, E_rel, mass_target, mass_proj

!   real(kind=8) :: de_read

   real(kind=8) :: total_cs, new_total_cs
   real(kind=8) :: sum_direct
   real(kind=8) :: elastic_cs, new_elastic_cs
   real(kind=8) :: absorption_cs, new_absorption_cs
   real(kind=8) :: T_factor

   integer(kind=4) :: ibin
   real(kind=8) :: ex_min

   real(kind=8) :: dx, emin, emax


!--------------------External functions--------------------------------
   real(kind=8) :: tco_interpolate
   real(kind=8) :: clebr
   real(kind=8) :: comp_cs
!
!-------------------------------------------------------------------------------+
!-------                   Execute subroutine                                   +
!-------------------------------------------------------------------------------+
!
   zzero = 0.0d0
   len_prn = 0
   bigblank(1:80) = ' ' 


   do k = 1, 6                  !   loop over particles
      isp = nint(2.0d0*particle(k)%spin)
      if(.not.particle(k)%in_decay .and. k /= iproj)cycle        !   if not in decay chain cycle over
      if(print_me)write(6,*)'Setting up transmission coefficients for ',particle(k)%name
!
!----------   Setup grids energy grids for transmission coefficients 
!----------   for each particle type
!
!----   Find target for this decay particle
      iZ = nucleus(1)%Z - particle(k)%Z
      iA = nucleus(1)%A - particle(k)%A
      do jtarget = 1, num_comp
         if(iZ == nucleus(jtarget)%Z .and. iA == nucleus(jtarget)%A)then
            jstate = 1
            if(jtarget == itarget)jstate = istate
            exit
         end if
      end do
 
      tco_file = bigblank
      tco_file(1:1) = particle(k)%label
      len_tco = 1
      if(nucleus(jtarget)%atomic_symbol(1:1) == ' ')then
         tco_file(len_tco+1:len_tco+1) = nucleus(jtarget)%atomic_symbol(2:2)
         len_tco=len_tco+1
      else
         tco_file(len_tco+1:len_tco+2) = nucleus(jtarget)%atomic_symbol(1:2)
         len_tco=len_tco+2
      end if
      if(nucleus(jtarget)%A < 10)then
         write(tco_file(len_tco+1:len_tco+1),'(i1)')nucleus(jtarget)%A
         len_tco=len_tco+1
      elseif(nucleus(jtarget)%A < 100)then
         write(tco_file(len_tco+1:len_tco+2),'(i2)')nucleus(jtarget)%A
         len_tco=len_tco+2
      elseif(nucleus(jtarget)%A < 1000)then
         write(tco_file(len_tco+1:len_tco+3),'(i3)')nucleus(jtarget)%A
         len_tco=len_tco+3
      end if

!   write(6,*)jtarget
!   write(6,*)k,nucleus(jtarget)%Z, nucleus(jtarget)%A
!   write(6,*)tco_file(1:len_tco)
      tco_exist = .false.
      inquire(file=tco_file(1:len_tco)//'.tcoef',exist=tco_exist)
!----    Only need Ang_exist for k=iproj
      Ang_exist = .true.
      if(k == iproj)then 
         prn_file = tco_file(1:len_tco)//'-CC'
         len_prn = len_tco + 3
         inquire(file=prn_file(1:len_prn)//'.data',exist=Ang_exist)
      end if
      optical_run = .false.
 19   if(.not.tco_exist .or. .not. Ang_exist)then   !  we need to run the optical model
!-------------------------------------------------------------------------+
!------  Set up optical model calculation                                 +
!-------------------------------------------------------------------------+
!         call make_fresco_tco(data_path,len_path,tco_file,len_tco,      &
!                              de,k,iproj,jtarget,jstate,Ang_L_max)
         call make_fresco_tco(data_path,len_path,tco_file,len_tco,      &
                              k,iproj,jtarget,jstate,Ang_L_max)
      end if
!***********************************************************************
!-----   Transmission coefficient file exists, read in needed info
!***********************************************************************
      open(unit=50,file=tco_file(1:len_tco)//'.tcoef',status='old')
      max_nume = -10000
      read(50,*)
      read(50,*)
      read(50,'(9x,i6)')nume
      read(50,'(9x,i6)')lmax
      particle(k)%lmax = lmax
      particle(k)%nume = nume
      if(allocated(particle(k)%e_grid))                                &
          deallocate(particle(k)%e_grid)
      if(allocated(particle(k)%trans_read))                            &
          deallocate(particle(k)%trans_read)
      if(.not.allocated(particle(k)%e_grid))                           &
          allocate(particle(k)%e_grid(nume))
      if(.not.allocated(particle(k)%trans_read))                       &
          allocate(particle(k)%trans_read(1:nume,0:isp,0:lmax))
      read(50,*)
      num_e_lines = nume/10
      if(10*num_e_lines < nume)num_e_lines = num_e_lines + 1
      do ll = 1, num_e_lines
         minl = (ll-1)*10 + 1
         maxl = min(nume,minl+9)
!         read(50,'(10(1x,1pe16.9))')(particle(k)%e_grid(i),i = minl, maxl)
         read(50,*)(particle(k)%e_grid(i),i = minl, maxl)
      end do


      if(particle(k)%e_grid(nume) < particle(k)%max_e)then        ! Needs in calculation are incompatible 
         if(iproc == 0)then
            write(6,*)'Warning requested maximum energy is ',     &  ! with data file, need to remake
                      'greater than the maximum energy in'
            write(6,*)'the transmission coefficient file'
            write(6,*)'Will try to call tco subroutine to build'
         end if
         tco_exist=.false.
         Ang_exist = .false.
         deallocate(particle(k)%e_grid)
         deallocate(particle(k)%trans_read)
         close(unit=50)
         goto 19
      end if

      num_l_lines = (lmax+1)/10
      if(10*num_l_lines < lmax+1)num_l_lines = num_l_lines + 1
      do ii = 0, isp
         read(50,*)
         do ll = 1, num_l_lines
            minl = (ll-1)*10
            maxl = min(lmax,minl+9)
            read(50,*)
            read(50,*)
            do ie = 1, nume
               read(50,*)xdummy,(particle(k)%trans_read(ie,ii,l),l = minl, maxl)
               do l = minl, maxl
                  if(particle(k)%trans_read(ie,ii,l) < 1.0d-9)particle(k)%trans_read(ie,ii,l) = 1.0d-9
               end do
!   write(80,'(2(1x,i4),10(1x,1pe16.7))')isp,ie,(particle(k)%trans_read(ie,ii,l),l= minl, maxl)
            end do
         end do
      end do
      close(unit=50)

      do ie = 1, nume 
         do l = 0, lmax
            do j = 0, isp
               particle(k)%trans_read(ie,j,l) = particle(k)%trans_read(ie,j,l)*T_norm
            end do
         end do
      end do

      spin = particle(k)%spin

      if(abs(spin - 0.5d0) < 1.0d-5 .and. trans_avg_l)then
         do ie = 1, nume 
            do l = 0, lmax
               xdummy = 0.0d0
               xj = real(l,kind=8) - spin
               do j = 0, isp
                  xj = xj + real(j,kind=8)
                  if(xj < 0.0d0)cycle
                  xdummy = xdummy + real(l+j,kind=8)*particle(k)%trans_read(ie,j,l)
               end do
               xdummy = xdummy/real(2*l+1,kind=8)
               xj = real(l,kind=8) - spin
               do j = 0, isp
                  xj = xj + real(j,kind=8)
                  if(xj < 0.0d0)cycle
                  particle(k)%trans_read(ie,j,l) = xdummy
               end do
            end do
         end do
      end if


!--------   Read in data associated with Optical-Model calculation
!--------   Elastic and coupled-channels cross sections
!--------   And Legendre coefficients for the angular distributions
!--------   Data is read in on Ecis energy grid, later we will interpolate
!--------   on this grid to fill in output data

!   write(6,*)'k = ', ' iproj = ',iproj,' scale_elastic', scale_elastic

      if(k == iproj)then
         open(unit=50, file=prn_file(1:len_prn)//'.data',status='old')
         read(50,*)OpticalCS%nume, OpticalCS%numcc, OpticalCS%Max_L           ! data dimensions
!----   Now OpticalCS%max_L is tied to Ang_L_max and should always be the same
!----   But, if one had increased Ang_L_max, ran FRESCO, and then went back 
!----   and decreased Ang_L_max, there would be a mismatch and a possible 
!----   array out of bounds if OpticalCS%max_L > Ang_L_max
!----   A dynamic fix is difficult since many arrays are allocated with
!----   Ang_L_max before optical_setup is run.
         if(OpticalCS%max_L > Ang_L_max)then
            if(iproc == 0)then
               write(6,*)'---- WARNING ---- WARNING ---- WARNING ---- WARNING ---- WARNING ----'
               write(6,*)'Issue in optical-setup.f90'
               write(6,*)'OpticalCS%max_L > Ang_L_max for Legendre expansion of elastic scattering'
               write(6,*)'This is indicative of a mismatch in versions as they should be the same'
               write(6,*)'For the time being, setting OpticalCS%max_L = Ang_L_max'
               write(6,*)'But, you should rerun FRESCO and start over'
               write(6,*)'---- WARNING ---- WARNING ---- WARNING ---- WARNING ---- WARNING ----'
            end if
            OpticalCS%max_L = Ang_L_max
         end if
!
!----   Changed 5-21-2020 not to tie the DWBA calculation to the explicit energy grid.
!----
         allocate(OpticalCS%state(OpticalCS%numcc))                                 ! coupled-channels states
         allocate(OpticalCS%energy(OpticalCS%nume))                                 !  Ecis cross sections
         allocate(OpticalCS%optical_cs(OpticalCS%nume, OpticalCS%numcc))                  !  Ecis cross sections
         allocate(OpticalCS%optical_leg(OpticalCS%nume, 0:OpticalCS%Max_L, OpticalCS%numcc)) ! Legendre coeeficients
         allocate(cc_index(OpticalCS%numcc))
         K_band = 1000.0
!----   Put in check to see if DWBA calculation was done. This is denoted by 
!----   OpticalCS%state(i)%state_type < 1 (generally 0). If no state has OpticalCS%state(i)%state_type < 1
!----   then DWBA calculation was not performed, and there is no need to check if excitation energy 
!----   grid was the same.

         ex_min = nucleus(jtarget)%e_grid(1) - nucleus(jtarget)%delta_e(1)/2.0d0

         check_dwba = .false.
         do i = 1, OpticalCS%numcc
            read(50,*)idummy, cc_index(i), OpticalCS%state(i)%spin,OpticalCS%state(i)%parity,    &
                      OpticalCS%state(i)%K, OpticalCS%state(i)%energy, OpticalCS%state(i)%state_type
!            if(OpticalCS%state(i)%state_type < 1)particle(k)%do_dwba = .true.
            if(OpticalCS%state(i)%state_type < 1)check_dwba = .true.
         end do

         if(particle(iproj)%do_dwba .and. .not. check_dwba)then
            if(iproc == 0)then
               write(6,*)'ERROR!!!  -  do_dwba = .true. but Optical Model calculation was performed without DWBA states'
               write(6,*)'Edit command file and set do_dwba = .false. with the command "do_dwba n"'
            end if
            call MPI_Abort(icomm,101,ierr)
         else if(.not. particle(iproj)%do_dwba .and. check_dwba)then
            if(iproc == 0)then
               write(6,*)'ERROR!!!  -  do_dwba = .false. but Optical Model calculation was performed with DWBA states'
               write(6,*)'Edit command file and set do_dwba = .true. with the command "do_dwba y"'
            end if
            call MPI_Abort(icomm,101,ierr)
         end if


         rewind(50)
         read(50,*)
         do i = 1, OpticalCS%numcc
            read(50,*)idummy, cc_index(i), OpticalCS%state(i)%spin,OpticalCS%state(i)%parity,    &
                      OpticalCS%state(i)%K, OpticalCS%state(i)%energy, OpticalCS%state(i)%state_type
!
!----   Changed 5-21-2020 not to tie the DWBA calculation to the xplicit energy grid.
!----
!            if(do_dwba .and. OpticalCS%state(i)%state_type > 0 .and. abs(de_read - de) > 0.000999999)then
!                write(6,*)'+++++++++++++++++   ERROR!!   +++++++++++++++++++++++++++++++++++++++'
!                write(6,*)'Energy bin width in DWBA calculation is not the same as requested   +'
!                write(6,*)'in this calculation. This is an error. Stop, rerun fresco with      +'
!                write(6,*)'the same delta_e option                                             +'
!                write(6,*)'+++++++++++++++++   ERROR!!   +++++++++++++++++++++++++++++++++++++++'
!                stop
!            end if
!----   Before cc_index was setup in fresco-setup and was defined as an energy bin. Now, the state is computed
!----   as a state with a defined energy. We need the energy bin for this state in this calculation.
            if(cc_index(i) <= 0)then
               do ibin = 1, nucleus(jtarget)%nbin
                  ee = OpticalCS%state(i)%energy
                  if(ee >= nucleus(jtarget)%e_grid(ibin) - 0.5d0*nucleus(jtarget)%delta_e(ibin) .and.   &
                     ee < nucleus(jtarget)%e_grid(ibin) + 0.5d0*nucleus(jtarget)%delta_e(ibin))exit
               end do              
!               ibin = int((OpticalCS%state(i)%energy - ex_min)/de) + 1
               cc_index(i) = ibin
            end if
            if(OpticalCS%state(i)%state_type > 0 .and. OpticalCS%state(i)%spin < K_band)     &
               K_band = OpticalCS%state(i)%spin
         end do
         J_gs = OpticalCS%state(1)%spin
         do i = 1, OpticalCS%numcc
            if(OpticalCS%state(i)%state_type == -1)then
               K_state = real(OpticalCS%state(i)%K,kind=8)
               J_min = max(K_band,abs(K_band - K_state))
               J_max = K_band + K_state
               Ix_min = nint(J_min - nucleus(jtarget)%jshift)
               Ix_max = nint(J_max - nucleus(jtarget)%jshift)
               OpticalCS%state(i)%Ix_min = Ix_min
               OpticalCS%state(i)%Ix_max = Ix_max
               allocate(OpticalCS%state(i)%spin_prob(Ix_min:Ix_max))
               sum = 0.0d0
               do Ix = Ix_min, Ix_max
                  J_state = real(Ix_min,kind=8) + nucleus(jtarget)%jshift
                  xk_factor = clebr(J_gs,K_band,K_state,zzero,J_state,K_band)
                  OpticalCS%state(i)%spin_prob(Ix) = xk_factor
                  sum = sum + xk_factor
               end do
               OpticalCS%state(i)%spin_prob(Ix_min) = OpticalCS%state(i)%spin_prob(Ix_min)/sum
               do Ix = Ix_min + 1, Ix_max
                  OpticalCS%state(i)%spin_prob(Ix) = OpticalCS%state(i)%spin_prob(Ix)/sum + &
                                                  OpticalCS%state(i)%spin_prob(Ix-1)
               end do
            end if
         end do

!-----   Before proceeding, check that these coupled channels match with existing 
!-----   states in the target nucleus
!-----   Assume that they do, and prove that they don't
         match = .true.
         do in = 1, OpticalCS%numcc
            j = cc_index(in)
            if(OpticalCS%state(in)%state_type == 1)then
               if(abs(nucleus(jtarget)%state(j)%spin - OpticalCS%state(in)%spin) > 1.0d-3)match = .false.
               if(abs(nucleus(jtarget)%state(j)%energy - OpticalCS%state(in)%energy) > 1.0d-4)match = .false.
               if(.not. match)then
                  if(iproc ==0)write(6,*)'Error in coupled channels, state = ',in,' does not match with a state in target nucleus'
                  call MPI_Abort(icomm,101,ierr)
               end if
             end if
         end do
!
!-----   Also, regarding DWBA states, set up excitation energy window for sampling them. This addresses an issue 
!-----   arising when the spacing between the stated DWBA states is greater than thespacing between the bins
!-----   This leads to a funny series of peaks. Overall, these DWBA states to bins, are meant to be
!-----   distributed "evenly" in the energy window between adjacent DWBA states
         if(OpticalCS%numcc > 1)then
            do in = 1, OpticalCS%numcc
               OpticalCS%state(in)%E_min = 0.0d0
               OpticalCS%state(in)%Delta_E = 0.0d0
               if(OpticalCS%state(in)%state_type /= 1)then
                  if(in == OpticalCS%numcc)then              !   last one, no states above
                     if(OpticalCS%state(in-1)%state_type /= 1)then
                        OpticalCS%state(in)%Delta_E = OpticalCS%state(in)%energy -             &
                                                      OpticalCS%state(in-1)%energy
                        OpticalCS%state(in)%E_min = OpticalCS%state(in)%energy -               &
                                                    OpticalCS%state(in)%Delta_E/2.0d0
                     else
                        OpticalCS%state(in)%Delta_E = OpticalCS%state(in+1)%energy -             &
                                                      OpticalCS%state(in)%energy
!                        OpticalCS%state(in)%Delta_E = de
                        OpticalCS%state(in)%E_min = OpticalCS%state(in)%energy -               &
                                                    OpticalCS%state(in)%Delta_E/2.0d0
                     end if
                     cycle
                  end if
                  if(OpticalCS%state(in-1)%state_type /= 1)then
                     OpticalCS%state(in)%Delta_E = (OpticalCS%state(in)%energy -             &
                                                    OpticalCS%state(in-1)%energy)/2.0d0 +    &
                                                   (OpticalCS%state(in+1)%energy -           &
                                                    OpticalCS%state(in)%energy)/2.0d0
                     OpticalCS%state(in)%E_min = OpticalCS%state(in)%energy -                &
                                                 (OpticalCS%state(in)%energy -               &
                                                  OpticalCS%state(in-1)%energy)/2.0d0
                  else                                !  First DWBA state, no states below
                     dx = (OpticalCS%state(in+1)%energy - OpticalCS%state(in)%energy)
                     emax = OpticalCS%state(in+1)%energy - dx/2.0
                     emin = max(OpticalCS%state(in)%energy - dx/2.0,Ex_min + 0.0005)
                     OpticalCS%state(in)%Delta_E = emax - emin
                     OpticalCS%state(in)%E_min = emin
                  end if
               end if
            end do
         end if

         num_l_lines = (OpticalCS%numcc+1)/10
         if(10*num_l_lines /= OpticalCS%numcc)num_l_lines = num_l_lines + 1
         do ll = 1, num_l_lines
            minl = (ll-1)*10+1
            maxl = min(OpticalCS%numcc,minl+9)
            read(50,*)
            do in = 1, OpticalCS%nume
               read(50,*)OpticalCS%energy(in),(OpticalCS%optical_cs(in,n), n = minl, maxl)            ! Ecis cross sections for coupled-channels
            end do                                                             ! n=target%istate is elastic
         end do
         do i = 1, OpticalCS%numcc
            read(50,*)
            do in = 1, OpticalCS%nume
               read(50,*)xdummy,(OpticalCS%optical_leg(in,L,i),L = 0, OpticalCS%Max_L)              !  Legendre coefficients for angular distributiuons
            end do
         end do
         close(unit=50)
!--------    Now map the coupled-channels states to the excitation spectrum of the target
!--------    store in OpticalCS%state(n)%istate
         do n = 1, OpticalCS%numcc
            OpticalCS%state(n)%istate = cc_index(n)
            if(cc_index(n) == target%istate .and. OpticalCS%state(n)%state_type == 1)OpticalCS%ielastic = n
         end do
         deallocate(cc_index)


      end if
   end do

!
      particle(k)%lmax = lmax
      particle(k)%nume = nume

!------    Check if we want to scale the elastic cross section -  Engineering fix.
   if(scale_elastic)then
      if(print_me)then
         write(6,*)'***************************************'
         write(6,*)'scale_elastic =',scale_elastic
         write(6,*)'***************************************'
      end if

      mass_target = nucleus(jtarget)%mass + nucleus(jtarget)%state(istate)%energy
      mass_proj = particle(iproj)%mass
      isp = nint(2.0d0*particle(iproj)%spin)
      do ie = 1, particle(iproj)%nume
         E_rel = particle(iproj)%e_grid(ie)
         E_in = E_rel*(mass_target + mass_proj)/mass_target
         sum_direct = 0.0d0
         elastic_cs = 0.0d0
         n_elastic = 1
         do n = 1, OpticalCS%numcc
            if(OpticalCS%state(n)%istate == istate .and.OpticalCS%state(n)%state_type == 1)then
               n_elastic = n
               elastic_cs = OpticalCS%optical_cs(ie,n)
            else
               sum_direct = sum_direct + OpticalCS%optical_cs(ie,n)
            end if
         end do
         absorption_cs = comp_cs(ie,jtarget,istate,iproj)
         total_cs = elastic_cs + sum_direct + absorption_cs
         new_elastic_cs = elastic_cs*elastic_scale*(1.0d0 + elastic_shift*(1.0d0 - exp(-elastic_damp*E_in)))
         new_absorption_cs = total_cs - sum_direct - new_elastic_cs
         T_factor = new_absorption_cs/absorption_cs
         OpticalCS%optical_cs(ie,n_elastic) = new_elastic_cs
!         write(6,*)'***************************************'
!         write(6,*)'T_factor = ',T_factor
!         write(6,*)'***************************************'
         do l = 0, particle(iproj)%lmax
            do j = 0, isp
               particle(iproj)%trans_read(ie,j,l) =               &
                   particle(iproj)%trans_read(ie,j,l)*T_factor
            end do
         end do
!-----    Check if it worked
         new_absorption_cs = comp_cs(ie,jtarget,istate,iproj)
         new_total_cs = new_elastic_cs + sum_direct + new_absorption_cs
 
!         write(6,*)'Incident Energy = ',E_in, 'E_rel = ',particle(iproj)%e_grid(ie)
!         write(6,*)'Before ',total_cs, elastic_cs, sum_direct, absorption_cs
!         write(6,*)'After ',new_total_cs, new_elastic_cs, sum_direct, new_absorption_cs
      end do
   end if



    do k = 1, 6
       nume = particle(k)%nume
       if(.not.particle(k)%in_decay)cycle        !   if not in decay chain cycle over
       particle(k)%nbin = int(particle(k)%max_e/de) + 1
       nbin = particle(k)%nbin
       isp = nint(2.0d0*particle(k)%spin)
       lmax = particle(k)%lmax
       allocate(particle(k)%trans(0:isp,0:lmax,nbin))
       do l = 0, lmax
          do j = 0, isp
             do i= 1, nbin
                energy = real(i,kind=8)*de
                particle(k)%trans(j,l,i) =                              & !   changed order for better access in HF decay
                     tco_interpolate(energy,nume,                       &
                                     particle(k)%e_grid,                &
                                     particle(k)%trans_read(1,j,l))
             end do
         end do
      end do
    end do

   return

end subroutine optical_setup
!
!
!*******************************************************************************
!
real(kind=8) function tco_interpolate(e,nume,e_grid,tco)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine finds the transmission coeeficient at a 
!    specific energy e by interpolating from a list of nume values 
!    in the array tco on a energy grid defined by e_grid 
!    Interpolate using a log-log linear approximation
!    
!    Improved 6 Jan 2021 making use of the constant grid in log(e) to
!    find mid points. Sped up by 6x.
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
   use nodeinfo
   use variable_kinds
   implicit none
   real(kind=8), intent(in) :: e
   integer(kind=4), intent(in) :: nume
   real(kind=8), intent(in) :: e_grid(nume),tco(nume)
!-----------------------------------------------------------------------------
   integer(kind=4) :: i1, i2 
!   real(kind=8) :: tco_1,tco_2
   real(kind=8) :: x1, x2, y1, y2, x, y, a, b 
   real(kind=8), parameter :: tolerance = 1.0d-5

   real(kind=8) :: delta
!-----------------------------------------------------------------------------
   tco_interpolate = 1.0d-9
   if(e <= 0.0d0)return
   tco_interpolate = tco(1)
   if(e <= e_grid(1))return
!   if(abs(log(e) - log(e_grid(i1))) <= tolerance)return

   x = log(e)
   x1 = log(e_grid(1))

   delta = log(e_grid(2)) - x1
   i1 = int((x - x1)/delta) + 1
   i2 = i1 + 1

!   if(i1 >= nume)then
!      if(iproc == 0)write(6,*)'Error in log_tco_interpolate e > egrid(nume)'
!      call MPI_Abort(icomm,101,ierr)
!   end if

   x1 = log(e_grid(i1))
   x2 = log(e_grid(i2))
!------   Put in lower threshold
!   tco_1 = max(tco(i1),1.0d-9)
!   tco_2 = max(tco(i2),1.0d-9)
!   y1 = log(tco_1)  
!   y2 = log(tco_2)
   y1 = log(max(tco(i1),1.0d-9))
   y2 = log(max(tco(i2),1.0d-9))
   a = (y2-y1)/(x2-x1)
   b = y1 - a*x1
   y = a*x + b
   tco_interpolate = exp(y)

!  write(6,*)'LOG'
!  write(6,*)i1,i2
!  write(6,*)x,x1,x2,y1,y2
!  write(6,*)a,b,y

  return
end function tco_interpolate

