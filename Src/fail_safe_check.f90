!
!*****************************************************************************
!
subroutine fail_safe_check
!
!*****************************************************************************
!
!  Discussion:
!
!    This subroutine is used to check if certain parameters are unphysical
!    and could lead to dangerous or unreliable results. If found, an error
!    message is given on unit 6 and execution is terminated.
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        constants
!        nodeinfo
!        useful_data
!        nuclei
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        real(kind=8) :: EL_f
!        real(kind=8) :: ML_f
!
!     MPI routines:
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
!*****************************************************************************
!
!----------------------------------------------------------------------
!-------    Use modules
!----------------------------------------------------------------------
   use variable_kinds
   use options
   use constants
   use nodeinfo
   use useful_data
   use nuclei
!----------------------------------------------------------------------
   implicit none
   integer(kind=4) :: inuc
   integer(kind=4) :: i, imax, l, n
   real(kind=8) :: step, e_gamma, gsf
   integer(kind=4) :: ilast
   logical :: terminate_execution, gsf_negative
!---------   External functions
   real(kind=8) :: EL_f
   real(kind=8) :: ML_f
!----------------------------------------------------------------------
!-----   Main loop over all compound nuclei in the calculation
!----------------------------------------------------------------------
   step = 0.01
   terminate_execution = .false.
   do inuc = 1, num_comp
!-----   Check that the matching energy is greater than Delta + 0.2
      ilast = index(nucleus(inuc)%Label,' ')-1
      if(ilast < 0)ilast = 5
      if(nucleus(inuc)%level_param(6) <= nucleus(inuc)%level_param(3) + 0.20d0)then
         if(iproc == 0)then
            write(6,*)
            write(6,'(''*************************************************************'')')
            write(6,'(''For compound nucleus '',a5)')nucleus(inuc)%Label(1:ilast)
            write(6,'(''The matching energy for the level density is less than Delta + 0.2 MeV.'')')
            write(6,'(''This will lead to unstable results. Execution will terminate.'')')
            write(6,'(''*************************************************************'')')
            write(6,*)
         end if
         terminate_execution = .false.
      end if
!-----   Check electromagnetic strength functions
      imax = int(nucleus(inuc)%Ex_max/step) + 1
!-----   Electric transitions
      do L = 1, nucleus(inuc)%lmax_E
         gsf_negative = .false.
         do i = 1, imax
            e_gamma = real(i,kind=8)*step
            gsf = EL_f(inuc, L, e_gamma, nucleus(inuc)%sep_e(1))
            if(gsf < 0.0d0)gsf_negative = .true.
         end do
         if(gsf_negative)then
            if(iproc == 0)then
               write(6,*)
               write(6,'(''*************************************************************'')')
               write(6,'(''For compound nucleus '',a5)')nucleus(inuc)%Label(1:ilast)
               write(6,'(''The E'',i1,'' strength function is negative.'')')L
               write(6,'(''This will lead to unstable results. Execution will terminate.'')')
               write(6,'(''*************************************************************'')')
               write(6,*)
            end if
            terminate_execution = .true.
         end if
      end do
!-----   Magnetic transitions
      do L = 1, nucleus(inuc)%lmax_M
         gsf_negative = .false.
         do i = 1, imax
            e_gamma = real(i,kind=8)*step
            gsf = ML_f(inuc, L, e_gamma)
            if(gsf < 0.0d0)gsf_negative = .true.
         end do
         if(gsf_negative)then
            if(iproc == 0)then
               write(6,*)
               write(6,'(''*************************************************************'')')
               write(6,'(''For compound nucleus '',a5)')nucleus(inuc)%Label(1:ilast)
               write(6,'(''The E'',i1,'' strength function is negative.'')')L
               write(6,'(''This will lead to unstable results. Execution will terminate.'')')
               write(6,'(''*************************************************************'')')
               write(6,*)
            end if
            terminate_execution = .true.
         end if
      end do
!----   Fission barriers
      if(nucleus(inuc)%fission)then
!----    Loop over barriers
         do n = 1, nucleus(inuc)%F_n_barr
!----   Curvature, hbw of barrier < 0
            if(nucleus(inuc)%F_Barrier(n)%hbw < 0.0d0)then
               if(iproc == 0)then
                  write(6,*)
                  write(6,'(''*************************************************************'')')
                  write(6,'(''For compound nucleus '',a5,'' Fision barrier #'',i1)')  &
                     nucleus(inuc)%Label(1:ilast),n
                  write(6,'(''The curvature, hbw, for this barrier is negative.'')')
                  write(6,'(''This will lead to unstable results. Execution will terminate.'')')
                  write(6,'(''*************************************************************'')')
                  write(6,*)
               end if
               terminate_execution = .true.
            end if
!----   E_match < Delta for this barrier
            if(nucleus(inuc)%F_Barrier(n)%level_param(6) <                          &
               nucleus(inuc)%F_Barrier(n)%level_param(3))then
               if(iproc == 0)then
                  write(6,*)
                  write(6,'(''*************************************************************'')')
                  write(6,'(''For compound nucleus '',a5,'' Fision barrier #'',i1)')  &
                     nucleus(inuc)%Label(1:ilast),n
                  write(6,'(''Ematch < Delta + 0.2 MeV for this barrier.'')')
                  write(6,'(''This will lead to unstable results. Execution will terminate.'')')
                  write(6,'(''*************************************************************'')')
                  write(6,*)
               end if
               terminate_execution = .true.
            end if
         end do
      end if     
   end do 

   if(terminate_execution)call MPI_Abort(icomm,1,ierr)

   return

end subroutine fail_safe_check
