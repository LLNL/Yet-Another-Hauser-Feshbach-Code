!
!*****************************************************************************80
!
subroutine print_channels(in, itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par, &
                          num_energies, num_e, max_jx_20, delta_jx_20,                    &
                          de_spec, cs_threshold, reaction_cs, write_error)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write exit channel reaction data the library
!    directory
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        constants
!        nodeinfo
!        nuclei
!        Channel_info
!        particles_def
!
!     Subroutines:
!
!        nucleus_label
!        get_to_eof
!        exit_YAHFC
!
!     External functions:
!
!        poly
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
   use options
   use constants
   use nodeinfo
   use nuclei
   use Channel_info
   use particles_def
   implicit none
   integer(kind=4), intent(in):: in
   integer(kind=4), intent(in):: itarget
   integer(kind=4), intent(in):: istate
   integer(kind=4), intent(in):: ilab
   character(len=100), intent(in) :: file_lab
   integer(kind=4), intent(in):: ilib_dir
   character(len=132), intent(in) :: lib_dir
   character(len=1), intent(in) :: ch_par(0:1)
   integer(kind=4), intent(in) :: num_energies
   integer(kind=4), intent(in) :: num_e
   integer(kind=4), intent(in) :: max_jx_20
   real(kind=8), intent(in) :: delta_Jx_20
   real(kind=8), intent(in) :: de_spec
   real(kind=8), intent(in) :: cs_threshold
   real(kind=8), intent(in) :: reaction_cs(num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: inuc
   integer(kind=4) :: ipi, ipf, i
   integer(kind=4) :: j, jx
   integer(kind=4) :: k
   integer(kind=4) :: m
   integer(kind=4) :: n, nn
   integer(kind=4) :: L, L_max
   integer(kind=4) :: icc, icc_max
   real(kind=8) :: cs
   real(kind=8) :: Q
   real(kind=8) :: x
   real(kind=8) :: e_in, e_out
   real(kind=8) :: xnorm
   real(kind=8) :: Temp
   real(kind=8) :: alf, bet
   real(kind=8) :: sum
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   integer(kind=4) :: ilast
   character(len=132) :: directory
   character(len=132) :: outfile
   integer(kind=4), allocatable, dimension (:) :: ilast_chann
   real(kind=8), allocatable, dimension (:) :: Channel_prob
   character(len=128) :: fstring, temp_string
   real(kind=8) :: part_mult(0:6)

   real(kind=8), allocatable :: Temp_Ang_Dist(:,:,:)
   real(kind=8), allocatable :: Ang_Dis(:,:)
!----------------------------------------------------------------------
   real(kind=8) :: poly
!----------------------------------------------------------------------

   write_error = .false.

   iproj = projectile%particle_type

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Channel Probabilities if a population calculation    ---+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(pop_calc)then
      directory(1:ilib_dir) = lib_dir(1:ilib_dir)
      idir = ilib_dir + 1
      directory(idir:idir) = '/'

      ifile = 12
      outfile(1:ifile) = 'Channel_Prob'
      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')

      write(temp_string,*)num_channels

      if(in == 1)then
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
         if(.not. allocated(ilast_chann))allocate(ilast_chann(num_channels))
         do i = 1, num_channels
            ilast_chann(i) = index(Exit_Channel(i)%Channel_Label,' ') - 1
         end do
         if(pop_calc_prob)then
            write(100,'(''# Probabilities for each exit channel'')')
         else
            write(100,'(''# Cross sections for each exit channel in '',a2)')cs_units
         end if


         fstring = "('#         E_in      ',"//trim(adjustl(temp_string))//'(7x,a10,2x))'
         write(100,fstring)(Exit_Channel(i)%Channel_Label(1:ilast_chann(i)), i = 1, num_channels)

         fstring = "('#   ----------------',"//trim(adjustl(temp_string))//"('   ----------------'))"
         write(100,fstring)
      else
         call get_to_eof(100)
      end if


      if(.not. allocated(Channel_prob))allocate(Channel_prob(num_channels))

      fstring = '(1x,3x,1pe16.7,'//trim(adjustl(temp_string))//'(3x,1pe16.7))'

      e_in = projectile%energy(in)
      do i = 1, num_channels
         Channel_prob(i) = 0.0d0
         do n = 1, Exit_Channel(i)%num_cs
            Channel_prob(i) = Channel_prob(i) + Exit_Channel(i)%Channel_cs(n,in)*       &
                                                reaction_cs(in)*cs_scale
         end do
      end do

!  write(6,*)'fstring ',fstring

      write(100,fstring)e_in, (Channel_prob(i), i = 1, num_channels)

      close(unit=100)

      if(allocated(ilast_chann))deallocate(ilast_chann)
      if(allocated(Channel_prob))deallocate(Channel_prob)

   end if

   if(.not.allocated(Temp_Ang_Dist))allocate(Temp_Ang_Dist(0:max_jx_20,0:num_e,0:6))
   if(.not.allocated(Ang_Dis))allocate(Ang_Dis(0:max_jx_20,0:6))
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
         Q = Exit_Channel(i)%Q_value - nucleus(inuc)%state(j)%energy +                             &
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
            write_error = .true.
            return
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')

         if(in == 1)then
!   write(100,*)i,n
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ') - 1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
               istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
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
            if(inuc == itarget)then
               write(100,'(''# Gammas transitions between discrete states are not included'')')
            end if
            write(100,'(''#'')')
            if(Exit_channel(i)%num_particles > 0)then
               if(explicit_channels)then
                  write(100,'(''# Mass of emitted particles '')')
                  do m = 1, Exit_channel(i)%num_particles
                     k = Exit_channel(i)%decay_particles(m)
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')               &
                        m, particle(k)%mass/mass_u
                  end do
               else
                  write(100,'(''# Mass of emitted particles '')')
                  m = 0
                  do k = 1, 6
                     do nn = 1, Exit_channel(i)%num_part(k)
                        m = m + 1
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                           m,particle(k)%mass/mass_u
                     end do
                  end do
               end if
            end if
            write(100,'(''#'')')
            write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
            write(100,'(''#'')')
            write(100,'(''# Cross section data '')')
            if(pop_calc .and. pop_calc_prob)then
               write(100,1901)(particle(k)%label, k = 0, 6)
               else
                  write(100,1902)cs_units,(particle(k)%label, k = 0, 6)
            end if
 1901       format('#         E_in               Prob',5x,7(8x,'Mult(',a1,')',4x))
 1902       format('#         E_in              xs(',a2,')',5x,7(8x,'Mult(',a1,')',4x))
            write(100,'(''#'',9(''   ----------------''))')
         else
            call get_to_eof(100)
         end if


         e_in = projectile%energy(in)
         cs = Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)*cs_scale
         write(100,'(1x,9(3x,1pe16.7))')e_in, cs, (Exit_Channel(i)%part_mult(k,n),k=0,6)

         close(unit=100)

         if(xs_only)cycle

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Angular distribution and spectrum  ---------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------   Output P(x,Eout|Ein)   -----------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
            call exit_YAHFC(103)
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')

         if(in == 1)then
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
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
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')               &
                        m, particle(k)%mass/mass_u
                  end do
               else
                  write(100,'(''# Mass of emitted particles '')')
                  m = 0
                  do k = 1, 6
                     do nn = 1, Exit_channel(i)%num_part(k)
                        m = m + 1
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                           m,particle(k)%mass/mass_u
                     end do
                   end do
               end if
            end if
            write(100,'(''#'')')
            write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
            write(100,'(''#'')')
            write(100,'(''# Angular Distribution data  P(x,Eout|E_in) '')')
            write(100,'(''# Normalized to unity over integral Int(-1,1) int(0,E_max)  P(x,Eout|Ein)dx,dEout '')')
            write(100,'(''# Normalized using trapaziodial rule'')')
            if(inuc == itarget)then
               write(100,'(''# Gammas transitions between discrete states are not included'')')
            end if
            write(100,'(''# dx = '',1pe16.7,'' dEout = '',1pe16.7,'' MeV'')')delta_jx_20, de_spec
            write(100,'(''# Frame = COM (center-of-momentum)'')')
         else
            call get_to_eof(100)
         end if


         e_in = projectile%energy(in)
!         if(e_in < -Q)cycle
         if(e_in >= -Q)then
            cs = Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)
            if(pop_calc .and. pop_calc_prob)then
               write(100,'(''# E_in = '',1pe16.7,3x,''Decay Probability = '',1pe16.7)')               &
                       e_in,cs
            else
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                       e_in,cs,cs_units
            end if
            if(cs >= cs_threshold)then
               write(100,'(''#'')')
               write(100,1910)(particle(k)%label, k = 0,6)
 1910    format('#         E_in               E_out          x=cos(theta)  ',                         &
             7(3x,('P(x,Eout,',a1,'|E_in)')))
               write(100,'(''#'',10(''   ----------------''))')
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------    Find max value of icc to print out
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----   Normalize the particle spectra and find max icc to print out
               Temp_Ang_Dist(0:max_jx_20,0:num_e,0:6) = 0.0d0
               do k = 0, 6
                  xnorm = 0.0d0
                  do icc = 0, num_e
                     xnorm = xnorm + Exit_Channel(i)%Spect(k,n)%E_spec(icc)*de_spec
                  end do
                  if(xnorm > 1.0d-8)then
                     do icc = 0, num_e
                        Exit_Channel(i)%Spect(k,n)%E_spec(icc) =                                      &
                            Exit_Channel(i)%Spect(k,n)%E_spec(icc)/xnorm
                     end do
                  else
                     do icc = 0, num_e
                        Exit_Channel(i)%Spect(k,n)%E_spec(icc) = 0.0d0
                     end do
                  end if
               end do
               icc_max = num_e
               do icc = num_e - 1, 0, -1
                  sum = 0.0d0
                  do k = 0, 6
                     sum = sum + Exit_Channel(i)%Spect(k,n)%E_spec(icc)
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
                        do L = 0, Exit_Channel(i)%Spect(k,n)%E_Ang_L_Max(icc)
                           sum = sum + Exit_Channel(i)%Spect(k,n)%E_Ang_L(L,icc)*                     &
                                       poly(L,1,alf,bet,x)
                        end do
                        Temp_Ang_Dist(jx,icc,k) = sum
                     end do
                     xnorm = 0.0d0
                     do jx = 0, max_jx_20 - 1
                        xnorm = xnorm + (Temp_Ang_Dist(jx,icc,k) + Temp_Ang_Dist(jx+1,icc,k))*        &
                                delta_jx_20*0.5d0
                     end do
                     if(xnorm > 0.0d0)then
                        do jx = 0, max_jx_20
                           Temp_Ang_Dist(jx,icc,k) = Temp_Ang_Dist(jx,icc,k)/xnorm*                   &
                                               Exit_Channel(i)%Spect(k,n)%E_spec(icc)
                        end do
                     end if
                  end do
                  xnorm = 0.0d0
                  temp = 0.0d0
                  do icc = 0, icc_max
                     sum = 0.0d0
                     do jx = 0, max_jx_20 - 1
                        sum = sum + (Temp_Ang_Dist(jx,icc,k) + Temp_Ang_Dist(jx+1,icc,k))*            &
                                    0.5d0*delta_jx_20
                     end do
                     if(Exit_Channel(i)%Spect(k,n)%E_spec(icc)> 1.0d-8)                               &
                        xnorm = xnorm + sum*de_spec
                     temp = temp + Exit_Channel(i)%Spect(k,n)%E_spec(icc)*de_spec
                  end do
               end do

!----------   Print out data
               do icc = 0, icc_max
                  e_out = real(icc,kind=8)*de_spec
                  if(icc >0)e_out = e_out - de_spec/2.0d0
                  do jx = 0, max_jx_20
                     x = real(jx,kind=8)*delta_jx_20 - 1.0d0
                     write(100,'(1x,10(3x,1pe16.7))')e_in, e_out, x,                                  &
                          (Temp_Ang_Dist(jx,icc,k),k = 0,6)
                  end do
               end do
            end if
         end if

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
            call exit_YAHFC(103)
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')

         if(in == 1)then
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
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
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')             &
                        m, particle(k)%mass/mass_u
                  end do
               else
                  write(100,'(''# Mass of emitted particles '')')
                  m = 0
                  do k = 1, 6
                     do nn = 1, Exit_channel(i)%num_part(k)
                        m = m + 1
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                           m,particle(k)%mass/mass_u
                     end do
                   end do
               end if
            end if
            write(100,'(''#'')')
            write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
            write(100,'(''#'')')
            write(100,'(''# Angular Distribution data  P(x,Eout|Ein) '')')
            write(100,'(''# Normalized to unity over integral Int(-1,1) int(0,E_max)  P(x,Eout|Ein)dx,dEout '')')
            if(inuc == itarget)then
               write(100,'(''# Gammas transitions between discrete states are not included'')')
            end if
            write(100,'(''# dx = '',1pe16.7,'' dEout = '',1pe16.7,'' MeV'')')delta_jx_20, de_spec
            write(100,'(''# Frame = COM (center-of-momentum)'')')
         else
            call get_to_eof(100)
         end if

         e_in = projectile%energy(in)
         icc_max = num_e
         cs = Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)
         if(pop_calc)then
            write(100,'(''# E_in = '',1pe16.7,3x,''Decay Probability = '',1pe16.7)')               &
                    e_in,cs
         else
            write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                    e_in,cs,cs_units
         end if
         if(cs >= cs_threshold)then
            do icc = num_e - 1, 0, -1
               sum = 0.0d0
               do k = 0, 6
                  sum = sum + Exit_Channel(i)%Spect(k,n)%E_spec(icc)
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
               write(100,'(''#'',34x,''P(Eout|E_in)'',7(3x,1pe16.7))')                               &
                     (Exit_Channel(i)%Spect(k,n)%E_spec(icc),k=0,6)
               write(100,'(''#'')')
               write(100,2110)(particle(k)%label, k = 0,6)
 2110    format('#         E_in               E_out           L',7(3x,(5x,'a(L)(',a1,')',4x)))
               write(100,'(''#   ----------------   ----------------     ---'',7(3x,''----------------''))')
               L_max = 0
               do k = 0, 6
                  if(Exit_Channel(i)%Spect(k,n)%E_Ang_L_Max(icc) > L_max)                       &
                     L_max = Exit_Channel(i)%Spect(k,n)%E_Ang_L_Max(icc)
               end do
               do L = 0, L_max
                 write(100,'(1x,3x,1pe16.7,3x,1pe16.7,3x,i5,7(3x,1pe16.7))')e_in, e_out, L,        &
                    (Exit_Channel(i)%Spect(k,n)%E_Ang_L(L,icc),k=0,6)
               end do
            end do
         end if

         close(unit=100)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------------   Output P(x|Ein)   ------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         if(.not.allocated(Ang_Dis))allocate(Ang_Dis(0:max_jx_20,0:6))
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
            call exit_YAHFC(103)
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')

         if(in == 1)then
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
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
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')               &
                        m, particle(k)%mass/mass_u
                  end do
               else
                  write(100,'(''# Mass of emitted particles '')')
                  m = 0
                  do k = 1, 6
                     do nn = 1, Exit_channel(i)%num_part(k)
                        m = m + 1
                        write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                        write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                           m,particle(k)%mass/mass_u
                     end do
                   end do
                  end if
            end if
            write(100,'(''#'')')
            write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
            write(100,'(''#'')')
            write(100,'(''# Angular Distribution data  P(x|E_in) '')')
            write(100,'(''# Normalized to unity over integral Int(-1,1) P(x|E_in)dx '')')
            if(inuc == itarget)then
               write(100,'(''# Gammas transitions between discrete states are not included'')')
            end if
            write(100,'(''# dx = '',1pe16.7)')delta_jx_20
            write(100,'(''# Frame = COM (center-of-momentum)'')')
         else
            call get_to_eof(100)
         end if

         Ang_Dis(0:max_jx_20,0:6) = 0.0d0

          e_in=projectile%energy(in)
          if(e_in < -Q)cycle
          cs = Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)
          write(100,'(''#'')')
          if(pop_calc .and. pop_calc_prob)then
             write(100,'(''# E_in = '',1pe16.7,3x,''Decay Probability = '',1pe16.7)')            &
                    e_in,cs
          else
             write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                     e_in,cs,cs_units
          end if
         if(cs >= cs_threshold)then
            fstring = "('#',9x,'E_in',13x,'cos(theta)',2x,7(3x,(4x,'P(x,',a1,'|E_in)',1x)))"
            write(100,fstring)(particle(k)%label, k = 0,6)
            write(100,'(''#'',9(''   ----------------''))')
!----    Calculate Angular distribution again and then integrate over energy
            alf = 0.0d0
            bet = 0.0d0
            do k = 0, 6
               do jx = 0, max_jx_20
                  Ang_Dis(jx,k) = 0.0d0
                  sum = 0.0d0
                  do icc = 0, num_e
                     x = real(jx,kind=8)*delta_jx_20 - 1.0d0
                     do L = 0, Exit_Channel(i)%Spect(k,n)%E_Ang_L_Max(icc)
                        sum = sum + Exit_Channel(i)%Spect(k,n)%E_Ang_L(L,icc)*                  &
                                   poly(L,1,alf,bet,x)*                                           &
                                   Exit_Channel(i)%Spect(k,n)%E_spec(icc)*de_spec
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
         end if

         close(unit=100)

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
            call exit_YAHFC(103)
         end if
         open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')

         if(in == 1)then
             file_lab2(1:20) = ' '
             ilab2 = ilab
             file_lab2(1:ilab2) = file_lab(1:ilab)
             ilab2 = ilab2 + 1
             ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
             file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
             ilab2 = ilab2 + ilast
             file_lab2(ilab2:ilab2) = ')'
             write(100,'(''# '',a20)')file_lab2
             ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
             write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
                   istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''#'')')
             write(100,'(''# Mass amu = '',1pe23.16,'' MeV'')')mass_u
             write(100,'(''# Mass of target = '',1pe23.16,'' amu'')')nucleus(itarget)%mass/mass_u
             if(pop_calc)then
                write(100,'(''# Population decay, no projectile'')')
             else
                if(iproj >= 0 .and. iproj < 7)write(100,'(''# Mass of projectile = '',1pe23.16,'' amu'')')   &
                   particle(iproj)%mass/mass_u
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
                          write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')         &
                       m, particle(k)%mass/mass_u
                   end do
                else
                   write(100,'(''# Mass of emitted particles '')')
                   m = 0
                   do k = 1, 6
                      do nn = 1, Exit_channel(i)%num_part(k)
                         m = m + 1
                         write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                         write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')         &
                         m,particle(k)%mass/mass_u
                      end do
                   end do
                end if
             end if
             write(100,'(''#'')')
             write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
             write(100,'(''#'')')
             write(100,'(''# Probability density function for emission spectra P(Eout|Ein) '')')
             write(100,'(''# Normalized to unity over integral Int(0,Emax) P(Eout|Ein)dEout '')')
             if(inuc == itarget)then
                write(100,'(''# Gammas transitions between discrete states are not included'')')
             end if
             write(100,'(''# dEout = '',1pe16.7)')de_spec
             write(100,'(''# Frame = COM (center-of-momentum)'')')
         else
            call get_to_eof(100)
         end if

         e_in=projectile%energy(in)
         if(e_in >= -Q)then
            write(100,'(''#'')')
            cs = Exit_Channel(i)%Channel_cs(n,in)*reaction_cs(in)
            if(pop_calc .and. pop_calc_prob)then
               write(100,'(''# E_in = '',1pe16.7,3x,''Decay Probability = '',1pe16.7)')               &
                       e_in,cs
            else
               write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                       e_in,cs,cs_units
            end if
            if(cs >= cs_threshold)then
               write(100,1904)(particle(k)%label, k = 0,6)
  1904   format('#         E_in               E_out     ',7(3x,(2x,'P(Eout,',a1,'|E_in)')))
               write(100,'(''#'',9(''   ----------------''))')
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------    Find max value of icc to print out
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----
               icc_max = num_e
               do icc = num_e - 1, 0, -1
                  sum = 0.0d0
                  do k = 0, 6
                     sum = sum + Exit_Channel(i)%Spect(k,n)%E_spec(icc)
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
                  write(100,'(1x,9(3x,1pe16.7))')e_in, e_out,                                           &
                       (Exit_Channel(i)%Spect(k,n)%E_spec(icc),k = 0,6)
               end do
            end if
         end if

         close(unit=100)

      end do       !   End of loop over channels
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

         if(in == 1)then
            file_lab2(1:20) = ' '
            ilab2 = ilab
            file_lab2(1:ilab2) = file_lab(1:ilab)
            ilab2 = ilab2 + 1
            ilast = index(Exit_Channel(i)%Channel_Label,' ')-1
            file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
            ilab2 = ilab2 + ilast
            file_lab2(ilab2:ilab2) = ')'
            write(100,'(''# '',a20)')file_lab2
            ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
            write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
                  istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
                  nucleus(itarget)%state(istate)%energy
            write(100,'(''# To level = N/A'')')
            if(pop_calc .and. pop_calc_prob)then
               write(100,'(''# Decay Probability data '')')
               write(100,'(''#         E_in               Prob'')')cs_units
            else
               write(100,'(''# Cross section data '')')
               write(100,'(''#         E_in              xs('',a2,'')'')')cs_units
            end if
            write(100,'(''#'',2(''   ----------------''))')
         else
            call get_to_eof(100)
         end if

         e_in = projectile%energy(in)
         cs = Exit_Channel(i)%Channel_cs(0,in)*reaction_cs(in)*cs_scale
         write(100,'(1x,9(3x,1pe16.7))')e_in,cs

         close(unit=100)
      end if

   end do

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

!      j = Exit_Channel(i)%StateLabel(1)
      Q = Exit_Channel(i)%Q_value - nucleus(inuc)%state(1)%energy +                             &
          nucleus(itarget)%state(istate)%energy

      ifile = 8
      outfile(1:ifile) = 'Channel_'
      outfile(ifile+1:ifile+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
      ifile = ifile + ilast
      outfile(ifile+1:ifile+7) = '_cs_Tot'
      ifile = ifile + 7
      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'.dat',status = 'unknown')

      if(in == 1)then
         file_lab2(1:20) = ' '
         ilab2 = ilab
         file_lab2(1:ilab2) = file_lab(1:ilab)
         ilab2 = ilab2 + 1
         ilast = index(Exit_Channel(i)%Channel_Label,' ') - 1
         file_lab2(ilab2:ilab2+ilast) = Exit_Channel(i)%Channel_Label(1:ilast)
         ilab2 = ilab2 + ilast
         file_lab2(ilab2:ilab2) = ')'
         write(100,'(''# '',a20)')file_lab2
         ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
         write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
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
         write(100,'(''# Total Xs in this channel '')')
            if(inuc == itarget)then
               write(100,'(''# Gammas transitions between discrete states are not included'')')
            end if
         write(100,'(''#'')')
         if(Exit_channel(i)%num_particles > 0)then
            if(explicit_channels)then
               write(100,'(''# Mass of emitted particles '')')
               do m = 1, Exit_channel(i)%num_particles
                  k = Exit_channel(i)%decay_particles(m)
                  write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                  write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')               &
                     m, particle(k)%mass/mass_u
               end do
            else
               write(100,'(''# Mass of emitted particles '')')
               m = 0
               do k = 1, 6
                  do nn = 1, Exit_channel(i)%num_part(k)
                     m = m + 1
                     write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                     write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')            &
                        m,particle(k)%mass/mass_u
                  end do
               end do
            end if
         end if
         write(100,'(''#'')')
         write(100,'(''# Q-value = '',f15.8,'' MeV'')')Q
         write(100,'(''#'')')
         write(100,'(''# Cross section data '')')
         if(pop_calc .and. pop_calc_prob)then
            write(100,2901)(particle(k)%label, k = 0, 6)
            else
               write(100,2902)cs_units,(particle(k)%label, k = 0, 6)
         end if
 2901    format('#         E_in               Prob',5x,7(8x,'Mult(',a1,')',4x))
 2902    format('#         E_in              xs(',a2,')',5x,7(8x,'Mult(',a1,')',4x))
         write(100,'(''#'',9(''   ----------------''))')
      else
         call get_to_eof(100)
      end if

      e_in = projectile%energy(in)
      cs = 0.0d0
      part_mult(0:6) = 0.0d0
      do n = 1, Exit_Channel(i)%num_cs
         cs = cs + Exit_Channel(i)%Channel_cs(n,in)
         do k = 0, 6
            part_mult(k) = part_mult(k) + Exit_Channel(i)%part_mult(k,n)*Exit_Channel(i)%Channel_cs(n,in)
         end do
      end do
      if(cs > 1.0d-10)then
         do k = 0,6
            part_mult(k) = part_mult(k)/cs
         end do
      end if
      cs = cs*reaction_cs(in)*cs_scale
      write(100,'(1x,9(3x,1pe16.7))')e_in, cs, (part_mult(k),k=0,6)

      close(unit=100)

   end do

   if(allocated(Temp_Ang_Dist))deallocate(Temp_Ang_Dist)
   if(allocated(Ang_dis))deallocate(Ang_dis)



   return
end subroutine print_channels
