!
!*******************************************************************************
!
subroutine Gamma_gamma(icomp, l, Gamma_g, g_error)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine computes Gamma_gamma
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
   use nuclei
   use particles_def
   use constants 
   use nodeinfo
   implicit none
   integer(kind=4), intent(in) :: icomp, l
   real(kind=8), intent(out) :: Gamma_g
   integer(kind=4), intent(out) :: g_error
!-------------------------------------------------------------------------+
   integer(kind=4) k
   integer(kind=4) lmin,le_min,lm_min
   real(kind=8) :: xj, xj_min, xj_max
   real(kind=8) :: j_shift, xj_f, xj_f_min, xj_f_max
   integer(kind=4) n
   integer(kind=4) ip
!   integer(kind=4) i_f,n_f,j_f,j_f_min,j_f_max,ip_f
   integer(kind=4) n_f, j_f, j_f_min, j_f_max, ip_f
   real(kind=8) :: energy, e_f
   real(kind=8) :: trans 
   real(kind=8) :: Ex_max
   real(kind=8) :: S_part
   real(kind=8) :: xspin
   integer(kind=4) ipar
   real(kind=8):: xl, xj_l
   integer(kind=4) :: jjj, num, lem
   real(kind=8) :: delta_e
   real(kind=8) :: temp
!-------------------------------------------------------------------------+
!------     Function declarations
   real(kind=8) :: Gamma_g2
   real(kind=8) :: EL_trans, ML_trans
   integer(kind=4) :: find_ibin
!-------------------------------------------------------------------------+
!------                                                                   +
!--------------------   Start subroutine                                  +
!------                                                                   +
!-------------------------------------------------------------------------+

   g_error = 0
   Gamma_g = 0.0d0
   Gamma_g2 = 0.0d0
   Ex_max = nucleus(icomp)%Ex_max
   S_part = nucleus(icomp)%sep_e(1)                !    set up with too low energy
   if(S_part > Ex_max)then                     !    can't calculate Gamma_g0
      g_error = 1
      Gamma_g = -1.0
      return
   end if
!   n = int((S_part - nucleus(icomp)%e_grid(1))/delta_e) + 1
   energy = S_part
   n = find_ibin(energy, icomp)
   j_shift=nucleus(icomp)%jshift

   xspin = nucleus(icomp)%target_spin

   ipar = nucleus(icomp)%target_ipar
   ip = ((2*ipar - 1)*(-1)**l+1)/2
   xl = real(l, kind=8)
   do jjj = -1, 1, 2
      xj_l = xl + real(jjj,kind=8)/2.0
      if(xj_l < 0.0) cycle
      xj_min = abs(xspin - xj_l)
      xj_max = (xspin + xj_l)
      num = nint(xj_max - xj_min)
      do k = 0, num
         xj = xj_min + real(k, kind=8)
!         i_f = icomp
!---------------------------   gamma decay to continuous level bins
         do n_f = n, 1, -1              !  loop over final excitation energies
            e_f = energy - nucleus(icomp)%e_grid(n_f)
            delta_e = nucleus(icomp)%delta_e(n_f)
            if(e_f + 0.5d0*delta_e <= 0.0d0)cycle     !   bin is misalgined and not allowed
            if(e_f - 0.5d0*delta_e < 0.0d0)then       !   top of bin is above S_n
               temp = e_f + 0.5d0*delta_e
               e_f = 0.5d0*(e_f + 0.5d0*delta_e)
               delta_e = temp
            end if
!---------------------------   Start with Electric decay 
            do lem = 1, nucleus(icomp)%lmax_E                   !  loop over EL decays

               trans = EL_trans(icomp, lem, e_f, energy)

               ip_f = iand((ip+lem),1)                       !  parity of final state
               xj_f_min = abs(xj-dfloat(lem))                !  min final spin
               xj_f_max = xj + dfloat(lem)                     !  max final spin
               j_f_min = nint(xj_f_min - nucleus(icomp)%jshift)      !  min j-index
               j_f_max = min(nint(xj_f_max-nucleus(icomp)%jshift),nucleus(icomp)%j_max)      !  max j-index
               do j_f = j_f_min, j_f_max                    !  loop over final j
                  xj_f = dfloat(j_f) + nucleus(icomp)%jshift
                  if(xj < 1.0d-5.and.xj_f <= 1.0d-5)cycle        !  O -> 0 not allowed
                  Gamma_g = Gamma_g + trans*nucleus(icomp)%bins(j_f,ip_f,n_f)%rho*delta_e
               end do
            end do
!---------------------------   Now Magnetic decay 
            do lem = 1, nucleus(icomp)%lmax_M                      !  loop over ML decays

                trans = ML_trans(icomp, lem, e_f)

                ip_f = iand(ip + lem + 1,1)                   !  parity of final state
                xj_f_min = abs(xj-dfloat(lem))                !  min final spin
                xj_f_max = xj+dfloat(lem)                     !  max final spin
                j_f_min = nint(xj_f_min - nucleus(icomp)%jshift)      !  min j-index
                j_f_max = min(nint(xj_f_max-nucleus(icomp)%jshift),nucleus(icomp)%j_max)      !  max j-index
                do j_f = j_f_min,j_f_max                    !  loop over final j
                   xj_f = dfloat(j_f)+nucleus(icomp)%jshift
                   if(xj < 1.0d-5.and.xj_f <= 1.0d-5)cycle        !  O -> 0 not allowed
                   Gamma_g = Gamma_g+trans*nucleus(icomp)%bins(j_f,ip_f,n_f)%rho*delta_e
               end do
             end do
          end do
      end do
   end do
   do jjj = -1, 1, 2
      xj_l = xl + real(jjj,kind=8)/2
      if(xj_l < 0.0) cycle
      xj_min = abs(xspin - xj_l)
      xj_max = (xspin + xj_l)
      num = nint(xj_max - xj_min)
      do k = 0, num
         xj = xj_min + real(k, kind=8)
!         i_f = icomp
         j_f_max = nucleus(icomp)%j_max + nint(nucleus(icomp)%jshift)
         xj_f_max = real(j_f_max,kind=8)
!---------------------------   gamma decay to discrete states
         do n_f = 1,nucleus(icomp)%num_discrete
            e_f = energy - nucleus(icomp)%state(n_f)%energy
            xj_f = nucleus(icomp)%state(n_f)%spin
            if(xj < 1.0d-5.and.xj_f <= 1.0d-5)cycle        !  O -> 0 not allowed
            lmin = max(1,int(abs(xj_f-xj)))                      !   can't have L=0
            ip_f = iabs(nint((nucleus(icomp)%state(n_f)%parity-1)/2.))
            if(ip  ==  ip_f)then                  !  parity the same even L for E odd L for M
               if(iand(lmin,1) == 0)then
                  le_min=lmin
                  lm_min=lmin+1
               else
                  le_min=lmin+1
                  lm_min=lmin
               end if
            else                                !  parity the same odd L for E even L for M             
               if(iand(lmin,1) == 0)then
                  le_min=lmin+1
                  lm_min=lmin
               else
                  le_min=lmin
                  lm_min=lmin+1
               end if
            end if
!--------------------    Electric
            do lem = le_min,nucleus(icomp)%lmax_E,2
               trans = EL_trans(icomp, lem, e_f, energy)

               Gamma_g = Gamma_g + trans

            end do

!-------------------   Magnetic
            do lem = lm_min, nucleus(icomp)%lmax_M, 2

               trans = ML_trans(icomp, lem, e_f)

               Gamma_g = Gamma_g + trans
            end do
         end do
      end do
   end do
   Gamma_g = Gamma_g*nucleus(icomp)%D0/(2.0d0*pi)*1000.0d0
return
end subroutine Gamma_gamma
