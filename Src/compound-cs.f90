!
!*******************************************************************************
!
subroutine compound_xs(e_in, itarget, istate, iproj, sigma,      &
                       num_channel, channel_prob, ichannel)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates reaction cross section
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        particles_def
!        constants 
!        nodeinfo
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!       real(kind=8) :: jhat
!       real(kind=8) :: tco_interpolate
!       real(kind=8) :: KE_com
!       real(kind=8) :: pc_com
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
   use options
   use nuclei
   use particles_def
   use constants 
   use nodeinfo
   implicit none
   real(kind=8), intent(in) :: e_in
   integer(kind=4), intent(in) :: itarget, istate, iproj
   real(kind=8), intent(out) :: sigma
   integer(kind=4), intent(out) :: num_channel
   real(kind=8), pointer, dimension(:) :: channel_prob
   integer(kind=4), pointer, dimension(:,:) :: ichannel
!----------------------------------------------------------
   real(kind=8) cs, cs_fac, channel_xs
   real(kind=8) :: xj
   real(kind=8) :: xI, xI_min, xI_max
   integer(kind=4) :: Ix, Ix_min, Ix_max
   real(kind=8) :: mass_i, mass_t, e_rel
   real(kind=8) :: momentum, wave_number
   real(kind=8) :: spin_proj, spin_target
   integer(kind=4) isp
   integer(kind=4) nume
   real(kind=8) sum,sum1
   integer(kind=4) i, l, nn
   integer(kind=4) cpar, par
   real(kind=8) :: tcoef
!----------   External functions
   real(kind=8) :: jhat
   real(kind=8) :: tco_interpolate
   real(kind=8) :: KE_com
   real(kind=8) :: pc_com
!-----------------------------------------------------------
   spin_target = nucleus(itarget)%state(istate)%spin
   spin_proj = particle(iproj)%spin
   mass_i = particle(iproj)%Mass
   mass_t = nucleus(itarget)%Mass + nucleus(itarget)%state(istate)%energy
!   e_rel = e_in*mass_t/(mass_t + mass_i)
!   mass_rel = mass_i*mass_t/(mass_t + mass_i)
!   momentum = dsqrt(2.0d0*e_rel*mass_rel)
   e_rel = KE_com(mass_i, mass_t, E_in)
   momentum = pc_com(mass_i, mass_t, E_in)
   wave_number = momentum/hbar_c
   cs = pi/wave_number**2*fmsq_eq_barn
   cpar = nint(particle(iproj)%par*nucleus(itarget)%state(istate)%parity)
   isp = nint(2*spin_proj)
   nume = particle(iproj)%nume
   sum1 = 0.0d0

   sum1 = 0.0d0
   do l = 0, particle(iproj)%lmax                                !   loop over angular momentum
      par = nint(cpar*(-1.0d0)**l)
      xj = real(l,kind=8) - spin_proj
      do i = 0, isp                                         !   loop over channel spins
         xj = xj + real(i,kind=8)
         if(xj < 0.0d0)cycle   
         xI_min = abs(xj - spin_target)
         xI_max = xj + spin_target
         Ix_min = max(nint(xI_min-nucleus(1)%jshift),0)
         Ix_max = min(nint(xI_max-nucleus(1)%jshift),nucleus(1)%j_max)
         do Ix = Ix_min, Ix_max
            xI = real(Ix,kind=8) + nucleus(1)%jshift
            cs_fac = jhat(xI)/(jhat(spin_target)*jhat(spin_proj))
            tcoef = tco_interpolate(e_rel,nume,                      &
                                    particle(iproj)%e_grid,          &
                                    particle(iproj)%trans_read(1,i,l))

            channel_xs = cs*cs_fac*tcoef

            sum1 = sum1 + channel_xs
         end do
      end do
   end do

   sigma = sum1

!   if(sigma < 1.0d-7)return
   
   num_channel = 0
   do l = 0, particle(iproj)%lmax                                !   loop over angular momentum
      par = nint(cpar*(-1.0d0)**l)
      xj = real(l,kind=8) - spin_proj
      do i = 0, isp                                         !   loop over channel spins
         xj = xj + real(i,kind=8)
         if(xj < 0.0d0)cycle   
         xI_min = abs(xj - spin_target)
         xI_max = xj + spin_target
         Ix_min = max(nint(xI_min-nucleus(1)%jshift),0)
         Ix_max = min(nint(xI_max-nucleus(1)%jshift),nucleus(1)%j_max)
         do Ix = Ix_min, Ix_max
            xI = real(Ix,kind=8) + nucleus(1)%jshift
            cs_fac = jhat(xI)/(jhat(spin_target)*jhat(spin_proj))
            tcoef = tco_interpolate(e_rel,nume,                      &
                                    particle(iproj)%e_grid,          &
                                    particle(iproj)%trans_read(1,i,l))

            channel_xs = cs*cs_fac*tcoef

            if(channel_xs/sum1 < 1.0d-6)cycle
            num_channel = num_channel + 1
         end do
      end do
   end do

   allocate(channel_prob(num_channel))
   allocate(ichannel(4,num_channel))

   sum = 0.0d0
   nn = 0
   do l = 0, particle(iproj)%lmax                                !   loop over angular momentum
      par = nint(cpar*(-1.0d0)**l)
      xj = real(l,kind=8) - spin_proj
      do i = 0, isp                                             !   loop over channel spins
         xj = xj + real(i,kind=8)
         if(xj < 0.0d0)cycle
         xI_min = abs(xj - spin_target)
         xI_max = xj + spin_target
         Ix_min = max(nint(xI_min-nucleus(1)%jshift),0)
         Ix_max = min(nint(xI_max-nucleus(1)%jshift),nucleus(1)%j_max)
         do Ix = Ix_min, Ix_max
            xI = real(Ix,kind=8) + nucleus(1)%jshift
            cs_fac = jhat(xI)/(jhat(spin_target)*jhat(spin_proj))
            tcoef = tco_interpolate(e_rel,nume,                      &
                                    particle(iproj)%e_grid,          &
                                    particle(iproj)%trans_read(1,i,l))
            channel_xs = cs*cs_fac*tcoef
            if(channel_xs/sum1 < 1.0d-6)cycle
            sum = sum + channel_xs
            nn = nn + 1
            channel_prob(nn) = sum/sigma
            ichannel(1,nn) = l
            ichannel(2,nn) = i
            ichannel(3,nn) = Ix
            ichannel(4,nn) = (par+1)/2
         end do
      end do
   end do

   sigma = sum

   return

end subroutine compound_xs
!
!*******************************************************************************
!
real(kind=8) function comp_cs(ie,itarget,istate,k)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function calculates compound reaction cross section
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        particles_def
!        constants 
!        nodeinfo
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!       real(kind=8) :: jhat
!       real(kind=8) :: pc_com
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
   use options
   use nuclei
   use particles_def
   use constants 
   use nodeinfo
   implicit none
   integer(kind=4), intent(in) :: ie, itarget, istate, k
!----------------------------------------------------------
   real(kind=8) :: sigma
!----------------------------------------------------------
   real(kind=8) :: cs, cs_fac, channel_xs
   real(kind=8) :: xj
   real(kind=8) :: xI, xI_min, xI_max
   integer(kind=4) :: Ix
   real(kind=8) :: mass_i,mass_t
!   real(kind=8) :: mass_rel,e_rel
   real(kind=8) :: momentum,wave_number
   real(kind=8) :: spin_proj,spin_target
   integer(kind=4) isp
   real(kind=8) sum,sum1
   integer(kind=4) i, l
   real(kind=8) :: cpar
   integer(kind=4) par
   integer(kind=4) :: num
   real(kind=8) :: tcoef
   real(kind=8) :: mass_target
!----------   External functions
   real(kind=8) :: jhat
   real(kind=8) :: pc_com
!-----------------------------------------------------------
!   e_rel = particle(k)%e_grid(ie)
   mass_target = nucleus(itarget)%mass + nucleus(itarget)%state(istate)%energy
   spin_target = nucleus(itarget)%state(istate)%spin
   spin_proj = particle(k)%spin
   mass_i = particle(k)%Mass
   mass_t = nucleus(itarget)%Mass + nucleus(itarget)%state(istate)%energy
!   mass_rel = mass_i*mass_t/(mass_t + mass_i)
!   momentum = dsqrt(2.0d0*e_rel*mass_rel)
!   e_rel = KE_com(mass_i, mass_t, E_in)
   momentum = pc_com(mass_i, mass_t, particle(k)%e_grid(ie))
   wave_number = momentum/hbar_c
   cs = pi/wave_number**2*fmsq_eq_barn

   cpar = particle(k)%par*nucleus(itarget)%state(istate)%parity
   isp = nint(2*spin_proj)
   sum1 = 0.0d0

   sum = 0.0d0
   do l = 0, particle(k)%lmax                                !   loop over angular momentum
      par = nint(cpar*(-1.0d0)**l)
      xj = real(l,kind=8) - spin_proj
      do i = 0, isp                                         !   loop over channel spins
         xj = xj + real(i,kind=8)
         if(xj < 0.0d0)cycle   
         xI_min = abs(xj - spin_target)
         xI_max = xj + spin_target
         num = nint(xI_max - xI_min) + 1
         do Ix = 1, num
            xI = real(Ix-1,kind=8) + xI_min
            cs_fac = jhat(xI)/(jhat(spin_target)*jhat(spin_proj))
            tcoef = particle(k)%trans_read(ie,i,l)
            channel_xs = cs*cs_fac*tcoef
            if(channel_xs < 1.0d-6)cycle
            sum = sum + channel_xs
         end do
      end do
   end do
   sigma = sum
   comp_cs = sigma
   return
end function comp_cs
!
!*******************************************************************************
!
real(kind=8) function compound_cs(e_in,ipar,xI,itarget,istate,iproj)
!
!*******************************************************************************
!
!  Discussion:
!
!    This is an alternative function calculates compound reaction cross section
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        particles_def
!        constants 
!        nodeinfo
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!       real(kind=8) :: jhat
!       real(kind=8) :: tco_interpolate
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
   use options
   use nuclei
   use particles_def
   use constants 
   use nodeinfo
   implicit none
   real(kind=8), intent(in) :: e_in
   integer(kind=4), intent(in) :: ipar
   real(kind=8) xI
   integer(kind=4), intent(in) :: itarget,istate,iproj
!----------------------------------------------------------
   real(kind=8) :: cs,cs_fac
   real(kind=8) :: mass_i, mass_t, e_rel
   real(kind=8) :: momentum, wave_number
   real(kind=8) :: spin_proj, spin_target
   integer(kind=4) isp
   integer(kind=4) nume
   real(kind=8) sum
   integer(kind=4) i, l
   integer(kind=4) cpar, par
   real(kind=8) :: tcoef
   real(kind=8) :: xj
   real(kind=8) :: channel_xs
!----------   External functions
   real(kind=8) :: jhat
   real(kind=8) :: tco_interpolate
   real(kind=8) :: KE_com
   real(kind=8) :: pc_com
!-----------------------------------------------------------
   compound_cs = 0.0d0
   par = 2*ipar - 1
   spin_target = nucleus(itarget)%state(istate)%spin
   spin_proj = particle(iproj)%spin
   mass_i = particle(iproj)%Mass
   mass_t = nucleus(itarget)%Mass + nucleus(itarget)%state(istate)%energy
!   e_rel = e_in*mass_t/(mass_t + mass_i)
!   mass_rel = mass_i*mass_t/(mass_t + mass_i)
!   momentum = dsqrt(2.0d0*e_rel*mass_rel)
   e_rel = KE_com(mass_i, mass_t, E_in)
   momentum = pc_com(mass_i, mass_t, E_in)
   wave_number = momentum/hbar_c
   cs = pi/wave_number**2*fmsq_eq_barn
   cpar = nint(particle(iproj)%par*nucleus(itarget)%state(istate)%parity)
   isp = nint(2*spin_proj)
   nume = particle(iproj)%nume
   cs_fac = jhat(xI)/(jhat(spin_target)*jhat(spin_proj))

   sum = 0.0d0
   do l = 0, particle(iproj)%lmax                                !   loop over angular momentum
      if(nint(cpar*(-1.0d0)**l) /= par)cycle
      xj = real(l,kind=8) - spin_proj
      do i = 0, isp                                         !   loop over channel spins
         xj = xj + real(i,kind=8)
         if(xj < 0.0d0)cycle
         if(xI < abs(xj - spin_target) .or. xI > xj + spin_target)cycle    
         tcoef = tco_interpolate(e_rel,nume,                         &
                                 particle(iproj)%e_grid,             &
                                 particle(iproj)%trans_read(1,i,l))

         channel_xs = cs*cs_fac*tcoef

         if(channel_xs < 1.0d-7)cycle
         sum = sum + channel_xs
      end do
   end do

   compound_cs = sum

   return
end function compound_cs
!
!*******************************************************************************
!
subroutine photo_xs(e_in, itarget, istate, sigma,                 &
                    num_channel, channel_prob, ichannel)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates reaction cross section
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        particles_def
!        constants 
!        nodeinfo
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!       real(kind=8) :: jhat
!       real(kind=8) :: EL_absorption
!       real(kind=8) :: ML_absorption
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
   use options
   use nuclei
   use particles_def
   use constants 
   use nodeinfo
   implicit none
   real(kind=8), intent(in) :: e_in
   integer(kind=4), intent(in) :: itarget, istate
   real(kind=8), intent(out) :: sigma
   integer(kind=4), intent(out) :: num_channel
   real(kind=8), pointer, dimension(:) :: channel_prob
   integer(kind=4), pointer, dimension(:,:) :: ichannel
!----------------------------------------------------------
   real(kind=8) :: ex
   real(kind=8) :: cs, channel_xs
   real(kind=8) :: spin_target
   real(kind=8) :: xI, xI_min, xI_max
   integer(kind=4) :: Ix, Ix_min, Ix_max
   integer(kind=4) :: nn
   real(kind=8) :: sum
   real(kind=8) :: xL
   integer(kind=4) :: L, Lmax
   integer(kind=4) :: cpar, par
!----------   External functions
   real(kind=8) :: jhat
   real(kind=8) :: EL_absorption
   real(kind=8) :: ML_absorption
!-----------------------------------------------------------
   spin_target = nucleus(itarget)%state(istate)%spin
   cpar = nint(nucleus(itarget)%state(istate)%parity)
   ex = nucleus(itarget)%sep_e(1)

   cs = 0.0d0
!----   Count channels and allocate arrays
!----   Electric modes
   num_channel = 0
   Lmax = nucleus(itarget)%lmax_E
   do L = 1, Lmax
      xL = real(L,kind=8)
      xI_min = abs(spin_target - L)
      xI_max = spin_target + xL
      Ix_min = nint(xI_min - nucleus(itarget)%jshift)
      Ix_max = nint(xI_max - nucleus(itarget)%jshift)
      Ix_max = min(Ix_max,nucleus(itarget)%j_max)
      num_channel = num_channel + Ix_max - Ix_min + 1
      channel_xs = EL_absorption(itarget, L, e_in, ex)
      cs = cs + channel_xs
   end do
!----   Magnetic modes
   Lmax = nucleus(itarget)%lmax_M
   do L = 1, Lmax
      xL = real(L,kind=8)
      xI_min = abs(spin_target - xL)
      xI_max = spin_target + xL
      Ix_min = nint(xI_min - nucleus(itarget)%jshift)
      Ix_max = nint(xI_max - nucleus(itarget)%jshift)
      Ix_max = min(Ix_max,nucleus(itarget)%j_max)
      num_channel = num_channel + Ix_max - Ix_min + 1
      channel_xs = ML_absorption(itarget, L, e_in)
      cs = cs + channel_xs
   end do

   sigma = cs

   allocate(channel_prob(num_channel))
   allocate(ichannel(4,num_channel))

!----   Calculate the probabilities for each channel
   nn = 0
!----   Electric modes
   sum = 0.0d0
   Lmax = nucleus(itarget)%lmax_E
   do L = 1, Lmax
      xL = real(L,kind=8)
      xI_min = abs(spin_target - xL)
      xI_max = spin_target + xL
      Ix_min = nint(xI_min - nucleus(itarget)%jshift)
      Ix_max = nint(xI_max - nucleus(itarget)%jshift)
      Ix_max = min(Ix_max,nucleus(itarget)%j_max)
      channel_xs = EL_absorption(itarget, L, e_in, ex)
      par = cpar*(-1)**L
      do Ix = Ix_min, Ix_max
         xI = real(Ix,kind=8) + nucleus(itarget)%jshift
         nn = nn + 1
         sum = sum + channel_xs*jhat(xI)/(jhat(spin_target)*jhat(xL))
         channel_prob(nn) = sum/cs
         ichannel(1,nn) = L
         ichannel(2,nn) = 0
         ichannel(3,nn) = Ix
         ichannel(4,nn) = (par+1)/2
      end do
   end do
!----   Magnetic modes
   Lmax = nucleus(itarget)%lmax_M
   do L = 1, Lmax
      xL = real(L,kind=8)
      xI_min = abs(spin_target - xL)
      xI_max = spin_target + xL
      Ix_min = nint(xI_min - nucleus(itarget)%jshift)
      Ix_max = nint(xI_max - nucleus(itarget)%jshift)
      Ix_max = min(Ix_max,nucleus(itarget)%j_max)
      channel_xs = ML_absorption(itarget, L, e_in)
      par = -cpar*(-1)**L
      do Ix = Ix_min, Ix_max
         xI = real(Ix,kind=8) + nucleus(itarget)%jshift
         nn = nn + 1
         sum = sum + channel_xs*jhat(xI)/(jhat(spin_target)*jhat(xL))
         channel_prob(nn) = sum/cs
         ichannel(1,nn) = L
         ichannel(2,nn) = 1
         ichannel(3,nn) = Ix
         ichannel(4,nn) = (par+1)/2
      end do
   end do

   return

end subroutine photo_xs


