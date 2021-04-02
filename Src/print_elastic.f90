!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write all the elastic cross sections to the library 
!    directory
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
subroutine print_elastic(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,     &
                         num_energies, Ang_L_max, max_jx_100, delta_Jx_100,              &
                         SE_cs, SE_Ang, Elastic_cs, Elastic_Ang,                         &
                         nstates, Inelastic_cs, Inelastic_Ang_L, Inelastic_L_max,        &
                         write_error)
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
   integer(kind=4), intent(in) :: Ang_L_max
   integer(kind=4), intent(in) :: max_jx_100
   real(kind=8), intent(in) :: delta_Jx_100
   real(kind=8), intent(in) :: SE_cs(num_energies)
   real(kind=8), intent(in) :: SE_Ang(0:Ang_L_max,num_energies)
   real(kind=8), intent(inout) :: Elastic_cs(num_energies)
   real(kind=8), intent(inout) :: Elastic_Ang(0:Ang_L_max,num_energies)
   integer(kind=4), intent(in) :: nstates
   real(kind=8), intent(in) :: Inelastic_cs(0:nstates,num_energies)
   real(kind=8), intent(in) :: Inelastic_Ang_L(0:Ang_L_max,0:nstates,num_energies)
   integer(kind=4), intent(in) :: Inelastic_L_max(0:nstates,num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: ipi, ipf, in, j, jx, L
   real(kind=8) :: xA
   real(kind=8) :: x
   real(kind=8) :: e_in, cs
   real(kind=8) :: xnorm
   real(kind=8) :: Temp
   real(kind=8) :: alf, bet
   real(kind=8) :: sum
   real(kind=8) :: comp
   real(kind=8) :: P_L
   real(kind=8) :: shape, Coul, Sig_C
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   character(len=132) :: directory
   character(len=132) :: outfile
   real(kind=8), allocatable :: xvalue(:)
   real(kind=8), allocatable :: Ang_Dist(:)
!----------------------------------------------------------------------
   real(kind=8) :: poly
!----------------------------------------------------------------------
   write_error = .false.
   iproj = projectile%particle_type

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   idir = idir + 1
   directory(idir:idir) = particle(iproj)%label
   idir = idir + 1
   directory(idir:idir) = '/'
                                      ! don't have with population calculation
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Shape Elastic scattering cross section  ----------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ifile = 16
   outfile(1:ifile) = 'Shape_Elastic_cs'
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
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                           &
         nucleus(itarget)%state(istate)%energy
   if(iproj == 1)write(100,'(''# Shape Elastic cross section data '')')
   if(iproj > 1)write(100,'(''# <Shape Elastic/Rutherford> '')')
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

   if(xs_only)return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------     Angular distributions     ---------------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
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
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
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
         write(100,'(''# E_in = '',1pe16.7,3x,''OM/Rutherford = '',1pe16.7)')                    &
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
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
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
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
            j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''# Frame = COM'')')

   do in = 1, num_energies
      e_in = projectile%energy(in)
      if(iproj == 1)then
         write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                    e_in,SE_cs(in),cs_units
      elseif(iproj > 1)then
         write(100,'(''# E_in = '',1pe16.7,3x,''OM/Rutherford = '',1pe16.7)')                   &
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
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate - 1, nucleus(itarget)%state(istate)%spin, ch_par(ipi),                          &
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
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
          j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''# Frame = COM'')')

   xA = (nucleus(itarget)%mass+nucleus(itarget)%state(istate)%energy)/particle(iproj)%mass
   Coul = (particle(iproj)%Z*nucleus(itarget)%Z*                                                &
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
            if(L <= Inelastic_L_max(istate,in))comp = comp + Inelastic_Ang_L(L,istate,in)*P_L
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
         call Legendre_expand(max_jx_100+1, xvalue(0), Ang_Dist(0),                             &
                              Ang_L_max, Elastic_Ang(0,in))
         do L = 0, Ang_L_max
            Elastic_Ang(L,in) = Elastic_Ang(L,in)*xnorm
         end do
         Elastic_cs(in) = xnorm*0.5d0
      end if
      if(iproj == 1)then
         write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')') &
                    e_in,(SE_cs(in) + Inelastic_cs(istate,in)),cs_units
      elseif(iproj > 1)then
         write(100,'(''# E_in = '',1pe16.7,3x,'' <EL/Rutherford> = '',1pe16.7)')                &
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
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
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
      e_in = projectile%energy(in)
      write(100,'(1x,4(3x,1pe16.7))')e_in, cs
   end do
   close(unit=100)

!
!---   Also print out Legendre Coefficients
!
   ifile = 15
   outfile(1:ifile) = 'Elastic_Ang_Leg'
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
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         j-1, nucleus(itarget)%state(j)%spin, ch_par(ipi),                                       &
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
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
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
   return
end subroutine print_elastic


!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write all the compound elastic cross sections to 
!    the library directory
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
subroutine print_compound_elastic(itarget, istate, ilab, file_lab, ilib_dir, lib_dir, ch_par,     &
                         num_energies, Ang_L_max, max_jx_100, delta_Jx_100, cs_threshold,         &
                         nstates, Inelastic_cs, Inelastic_Ang_L, Inelastic_L_max, write_error)
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
   integer(kind=4), intent(in) :: Ang_L_max
   integer(kind=4), intent(in) :: max_jx_100
   real(kind=8), intent(in) :: delta_Jx_100
   real(kind=8), intent(in) :: cs_threshold
   integer(kind=4), intent(in) :: nstates
   real(kind=8), intent(in) :: Inelastic_cs(0:nstates,num_energies)
   real(kind=8), intent(in) :: Inelastic_Ang_L(0:Ang_L_max,0:nstates,num_energies)
   integer(kind=4), intent(in) :: Inelastic_L_max(0:nstates,num_energies)
   logical, intent(out) :: write_error
!----------------------------------------------------------------------
   integer(kind=4) :: iproj
   integer(kind=4) :: ipi, ipf, in, j, jx, L
!-rem   real(kind=8) :: xA
   real(kind=8) :: x
   real(kind=8) :: e_in, cs
   real(kind=8) :: xnorm
   real(kind=8) :: Temp
   real(kind=8) :: alf, bet
   real(kind=8) :: sum
!-rem   real(kind=8) :: comp
!-rem   real(kind=8) :: P_L
!-rem   real(kind=8) :: shape, Coul, Sig_C
   integer(kind=4) :: ilab2
   character(len=20) :: file_lab2
   integer(kind=4) :: idir, ifile
   character(len=132) :: directory
   character(len=132) :: outfile
!-rem   real(kind=8), allocatable :: xvalue(:)
   real(kind=8), allocatable :: Ang_Dist(:)
!----------------------------------------------------------------------
   real(kind=8) :: poly
!----------------------------------------------------------------------
   write_error = .false.
   iproj = projectile%particle_type

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   idir = idir + 1
   directory(idir:idir) = particle(iproj)%label
   idir = idir + 1
   directory(idir:idir) = '/'
                                      ! don't have with population calculation
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!--------   Compound Elastic crosss section     -------------------+
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
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
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
      j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''#            E_in              xs'',''('',a2,'')'')')cs_units
   write(100,'(''#'',2(''   ----------------''))')
   do in = 1, num_energies
      cs = Inelastic_cs(istate,in)
      e_in = projectile%energy(in)
      write(100,'(1x,4(3x,1pe16.7))')e_in, cs*cs_scale
   end do
   close(unit=100)

   if(xs_only)return

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
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
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
   write(100,'(''# Final state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')     &
      j-1,nucleus(itarget)%state(j)%spin, ch_par(ipf), nucleus(itarget)%state(j)%energy
   write(100,'(''#'')')
   write(100,'(''# Frame = COM'')')
   alf = 0.0d0
   bet = 0.0d0
   do in = 1, num_energies
      e_in = projectile%energy(in)
      write(100,'(''# E_in = '',1pe16.7,3x,''Cross Section = '',1pe16.7,1x,''('',a2,'')'')')    &
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
   write(100,'(''# Target state = '',i3,3x,''J = '',f4.1,a1,3x,''Ex = '',1pe15.7,'' MeV'')')    &
         istate-1,nucleus(itarget)%state(istate)%spin, ch_par(ipi),                             &
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

   return

end subroutine print_compound_elastic
