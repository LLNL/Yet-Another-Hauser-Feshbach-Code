!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine fit to experimental values of Gamma_gamma (if known)
!    by adding an additional E1 mode 
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
subroutine fit_nuke_Gamma_gamma
   use variable_kinds
   use options
   use constants
   use nodeinfo
   use nuclei
   implicit none
!-------------------------------------------------------------------------------
   integer(kind=4) :: icomp
   real(kind=8) :: diff, tolerance
   real(kind=8) :: step
   real(kind=8) :: Gamma_g, Gamma_g_old
   real(kind=8) :: Gamma_g_exp, dGamma_g_exp
   integer(kind=4) :: g_error
   integer(kind=4) :: ich
   integer(kind=4) :: num_res, num_res_old
   logical :: converged
!-------------------------------------------------------------------------------
   icomp = 1
   do icomp = 1, num_comp
      converged = .false.
      Gamma_g_exp = nucleus(icomp)%Gamma_g_exp
      call Gamma_gamma(icomp, 0, Gamma_g, g_error)

      if(g_error > 0) converged = .true.                 !  Don't try to fit as Gamma_g can't be computed (Ex_max < Sn)
!----   Don't fit if Gamma_g_exp < 0.0, no data
      if(Gamma_g_exp < 0.0d0)converged = .true.

      if(nucleus(icomp)%dGamma_g_exp > 0.0d0)then
         dGamma_g_exp = nucleus(icomp)%dGamma_g_exp
      else
         dGamma_g_exp = Gamma_g_exp*0.05d0              !  just in case there is data, but error bar is messed up
      end if
      tolerance = min(dGamma_g_exp,0.1d0)             !  tolerance on fit. Make it small so that we can actually get Exp
      if(tolerance < 1.0d-3)tolerance = 0.1d0
      diff = 0.0d0

      if(Gamma_g_exp > 0.0d0)diff = abs(Gamma_g - Gamma_g_exp)                 !  current difference

      num_res_old = nucleus(icomp)%EL_mode(1)%num_gsf

      if(diff > tolerance .and. nucleus(icomp)%fit_gamma_gamma)then
         nucleus(icomp)%EL_mode(1)%num_gsf = nucleus(icomp)%EL_mode(1)%num_gsf + 1
         num_res = nucleus(icomp)%EL_mode(1)%num_gsf
         nucleus(icomp)%EL_mode(1)%gsf(num_res)%er = 5.0d0
         nucleus(icomp)%EL_mode(1)%gsf(num_res)%gr = 5.0d0
         nucleus(icomp)%EL_mode(1)%gsf(num_res)%sr = 0.0d0

         step = nucleus(icomp)%EL_mode(1)%gsf(1)%sr/10000.0

         if(Gamma_g > Gamma_g_exp)step = - step       !   Calculated Gamma_g is too big, make strength negative
         nucleus(icomp)%EL_mode(1)%gsf(num_res)%sr = nucleus(icomp)%EL_mode(1)%gsf(num_res)%sr + step

         Gamma_g_old = Gamma_g
         call Gamma_gamma(icomp, 0, Gamma_g, g_error)

         diff = abs(Gamma_g - Gamma_g_exp)
         if(diff < tolerance)converged = .true.
         if(step > 0.0 .and. Gamma_g > Gamma_g_exp) converged = .true.       !   Shouldn't really get here, but exit a trap
         if(step < 0.0 .and. Gamma_g < Gamma_g_exp) converged = .true.
         do while(.not. converged)
             nucleus(icomp)%EL_mode(1)%gsf(num_res)%sr =                            &
                   nucleus(icomp)%EL_mode(1)%gsf(num_res)%sr + step
            call Gamma_gamma(icomp, 0, Gamma_g, g_error)

            diff = abs(Gamma_g - Gamma_g_exp)
            if(diff < tolerance)then
               converged = .true.
            else
               if(step < 0.0d0 .and. (Gamma_g_old > Gamma_g_exp .and. Gamma_g < Gamma_g_exp))then
                  step = -step/10.0
               elseif(step > 0.0d0 .and. (Gamma_g_old < Gamma_g_exp .and. Gamma_g > Gamma_g_exp))then
                  step = -step/10.0
               end if
            end if
            Gamma_g_old = Gamma_g
         end do
      end if

      nucleus(icomp)%Gamma_g = Gamma_g
      call Gamma_gamma(icomp, 1, Gamma_g, g_error)
      nucleus(icomp)%Gamma_g_1 = Gamma_g

      if(nucleus(icomp)%Gamma_g > -1.0 .and. print_me)then
         num_res = nucleus(icomp)%EL_mode(1)%num_gsf
         write(6,*)
         ich = 1
         if(nucleus(icomp)%atomic_symbol(1:1) == ' ')then
            ich = 2
            write(6,'(5x,i3,a1)')nucleus(icomp)%A,nucleus(icomp)%atomic_symbol(ich:2)
         else
            write(6,'(5x,i3,a2)')nucleus(icomp)%A,nucleus(icomp)%atomic_symbol(ich:2)
         end if
	 
         if(Gamma_g > 0.0d0)then
            if(nucleus(icomp)%Gamma_g_exp > 0.0d0)then
               write(6,'(''Gamma_gamma (l=0)'','' Calc = '',f12.3,'', Exp = '',f12.3,'' +/- '',f12.3)')     &
                  nucleus(icomp)%Gamma_g,nucleus(icomp)%Gamma_g_exp, nucleus(icomp)%dGamma_g_exp
            else
               write(6,'(''Gamma_gamma (l=0)'','' Calc = '',f12.3,'', Exp = UNAVAILABLE'')')                &
                  nucleus(icomp)%Gamma_g
            end if
            if(nucleus(icomp)%Gamma_g_1_exp > 0.0d0)then
               write(6,'(''Gamma_gamma (l=1)'','' Calc = '',f12.3,'', Exp = '',f12.3,'' +/- '',f12.3)')     &
                  nucleus(icomp)%Gamma_g_1,nucleus(icomp)%Gamma_g_1_exp, nucleus(icomp)%dGamma_g_1_exp
            else
               write(6,'(''Gamma_gamma (l=1)'','' Calc = '',f12.3,'', Exp = UNAVAILABLE'')')                &
                  nucleus(icomp)%Gamma_g_1
            end if
            if(num_res > num_res_old)then
               write(6,*)'In order to reproduce Gamma_gamma(l=0), an additional E1 resonance with paramters'
               write(6,*)'Centroid = ',nucleus(icomp)%EL_mode(1)%gsf(num_res)%er
               write(6,*)'Width    = ',nucleus(icomp)%EL_mode(1)%gsf(num_res)%gr
               write(6,*)'Strength = ',nucleus(icomp)%EL_mode(1)%gsf(num_res)%sr
            end if
         elseif(Gamma_g < 0.0d0)then
            if(nucleus(icomp)%Gamma_g_exp > 0.0d0)then
               write(6,'(''Gamma_gamma (l=0)'','' Calc = Not Calculated, Exp = '',f12.3,'' +/- '',f12.3)')  &
                  nucleus(icomp)%Gamma_g,nucleus(icomp)%Gamma_g_exp, nucleus(icomp)%dGamma_g_exp
            else
               write(6,'(''Gamma_gamma (l=0)'','' Calc = Not Calculated, Exp = UNAVAILABLE'')')
            end if
            if(nucleus(icomp)%Gamma_g_1_exp > 0.0d0)then
               write(6,'(''Gamma_gamma (l=1)'','' Calc = Not Calculated, Exp = '',f12.3,'' +/- '',f12.3)')  &
                  nucleus(icomp)%Gamma_g_1,nucleus(icomp)%Gamma_g_1_exp, nucleus(icomp)%dGamma_g_1_exp
            else
               write(6,'(''Gamma_gamma (l=1)'','' Calc = Not Calculated, Exp = UNAVAILABLE'')')
            end if

         end if

         write(6,*)
      end if
   end do
   return
end subroutine fit_nuke_Gamma_gamma
