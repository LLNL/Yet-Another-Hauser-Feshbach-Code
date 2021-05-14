!
!*****************************************************************************80
!
subroutine print_channel_gammas(itarget, ilab, file_lab, ilib_dir, lib_dir,            &
                                ch_par, in, e_in, reaction_cs, write_error)
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write exit channel gamma ray data in the library 
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
   use variable_kinds
   use options
   use constants
   use nodeinfo
   use nuclei
   use Channel_info
   use particles_def
   implicit none  
   integer(kind=4), intent(in):: itarget
   integer(kind=4), intent(in):: ilib_dir
   character(len=132), intent(in) :: lib_dir
   character(len=1), intent(in) :: ch_par(0:1)
   integer(kind=4) :: in
   real(kind=8), intent(in) :: e_in
   real(kind=8), intent(in) :: reaction_cs
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: inuc
   integer(kind=4) :: ipi, ipf, i
   integer(kind=4) :: k
   integer(kind=4) :: m
   integer(kind=4) :: nf, nn
   real(kind=8) :: cs
   real(kind=8) :: ratio
   real(kind=8) :: Q
   real(kind=8) :: sum, check_sum
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   integer(kind=4) :: ilast
   character(len=132) :: directory
   character(len=132) :: outfile
   integer(kind=4) :: ilab
   character(len=100) :: file_lab
   character(len=5) :: target_label
   integer(kind=4) :: ntar
   character(len=5) :: resid_label
   integer(kind=4) :: nres
!----------------------------------------------------------------------
   write_error = .false.
   iproj = projectile%particle_type

   call nucleus_label(itarget,ntar,target_label)
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
               write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')                  &
                     m, particle(k)%mass/mass_u
            end do
         else
            write(100,'(''# Mass of emitted particles '')')
            m = 0
            do k = 1, 6
               do nn = 1, Exit_channel(i)%num_part(k)
                  m = m + 1
                  write(100,'(''# Particle # '',i5,'' = '',a1)')m,particle(k)%label
                  write(100,'(''# Mass particle # '',i5,'' = '',1pe23.16,'' amu'')')                &
                        m, particle(k)%mass/mass_u
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
         sum = sum + Exit_channel(i)%Channel_cs(nn,in)
      end do

      if(.not. pop_calc)then
         sum = sum*reaction_cs*cs_scale
         write(100,'(''#'')')
         write(100,'(''# Total cross section to this exit decay channel = '',1pe15.7,1x,a2)')sum,cs_units
         write(100,'(''#'')')
      else
         sum = sum*reaction_cs
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
            cs = Exit_channel(i)%state(k)%cs(m,in)*reaction_cs*cs_scale
            ratio = 0.0d0
            if(sum > 1.0d-10)ratio = Exit_channel(i)%state(k)%cs(m,in)/sum
            check_sum = check_sum + ratio 
            write(100,'(1x,i4,1x,f4.1,a1,''  '',i4,1x,f4.1,a1,8x,f10.5,3x,e15.7,3x,e15.7,3x,e15.7)') &
                                                  k-1, nucleus(inuc)%state(k)%spin, ch_par(ipi),       &
                                                  nf-1, nucleus(inuc)%state(nf)%spin, ch_par(ipf),     &
                                                  nucleus(inuc)%state(k)%egamma(m),                  &
                                                  cs, ratio
         end do
      end do
      close(unit = 100)
   end do                               
   return
end subroutine print_channel_gammas

