!
!*****************************************************************************80
!
subroutine print_channel_gammas_all_ein(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,      &
                                        num_energies, reaction_cs)
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
   integer(kind=4), intent(in):: itarget
   integer(kind=4), intent(in):: istate
   integer(kind=4), intent(in):: ilab
   character(len=100), intent(in) :: file_lab
   integer(kind=4), intent(in):: ilib_dir
   character(len=132), intent(in) :: lib_dir
   character(len=1), intent(in) :: ch_par(0:1)
   integer(kind=4), intent(in) :: num_energies
   real(kind=8), intent(in) :: reaction_cs(num_energies)
!----------------------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: inuc
   integer(kind=4) :: ipi, ipf, i, in
   integer(kind=4) :: k
   integer(kind=4) :: m
   integer(kind=4) :: nf, nn
   real(kind=8) :: cs
   real(kind=8) :: ratio
   real(kind=8) :: Q
   real(kind=8) :: e_in
   real(kind=8) :: sum
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   integer(kind=4) :: ilast
   character(len=132) :: directory
   character(len=132) :: outfile
   character(len=5) :: resid_label
   integer(kind=4) :: nres
!
!-------------------------------------------------------------------------------------------
!-----   Print out gammas for each channel if track_gammas == .true. for each energy
!-------------------------------------------------------------------------------------------
!
   iproj = projectile%particle_type
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
      file_lab2(ilab2:ilab2) = ')'
      write(100,'(''# '',a20)')file_lab2
      ipi = nint((nucleus(itarget)%state(istate)%parity + 1)/2.0d0)
      write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')') &
            istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
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
            write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')                 &
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
      if(.not. pop_calc)then
         write(100,'(''# Gamma-ray probabilities for each incident energy '')')
      else
         write(100,'(''# Gamma-ray probabilities for each initial excitation energy '')')
      end if


      if(pop_calc .and. pop_calc_prob)then
         write(100,'(''# Fraction of decay for discrete gamma transitions in '',a5)')           &
              resid_label(1:nres)
         write(100,'(''#'')')
      else
         write(100,'(''# Cross section for discrete gamma transitions in '',a5)')resid_label(1:nres)
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
            write(100,1911)k-1, nucleus(inuc)%state(k)%spin,ch_par(ipi),                        &
                              nf-1, nucleus(inuc)%state(nf)%spin, ch_par(ipf),                     &
                              nucleus(inuc)%state(k)%egamma(m)
1911        format('#',i4,1x,f4.1,a1,2x,i4,1x,f4.1,a1,8x,f10.5)

            if(pop_calc .and. pop_calc_prob)then
               write(100,1920)
            else
               write(100,1906)cs_units
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
                  sum = sum + Exit_channel(i)%Channel_cs(nn,in)
               end do
               if(sum > 1.0d-10)ratio = cs/sum
1908   format(1x,f12.6,1x,i4,1x,f4.1,a1,'  ',i4,1x,f4.1,a1,8x,f10.5,3x,e15.7,3x,e15.7,3x,e15.7)
               write(100,1908)e_in, k-1, nucleus(inuc)%state(k)%spin, ch_par(ipi),              &
                              nf-1, nucleus(inuc)%state(nf)%spin, ch_par(ipf),                  &
                              nucleus(inuc)%state(k)%egamma(m),                                 &
                              cs, ratio
            end do
         end do
      end do
      close(unit = 100)
   end do
end subroutine print_channel_gammas_all_ein

