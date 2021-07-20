!
!*****************************************************************************
!
subroutine fail_safe_check(delt)
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
!        particles_def
!
!     Subroutines:
!
!        get_lev_den
!        finish_lev_den
!        exit_YAHFC
!
!     External functions:
!
!        real(kind=8) :: EL_f
!        real(kind=8) :: ML_f
!
!     MPI routines:
!
!        MPI_Abort   ----   via exit_YAHFC
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
   use particles_def
!----------------------------------------------------------------------
   implicit none
!----------------------------------------------------------------------
   real(kind=8), intent(in) :: delt
!----------------------------------------------------------------------
   integer(kind=4) :: inuc
   integer(kind=4) :: i, ii, imax, l, n
   real(kind=8) :: step, e_gamma, gsf
   integer(kind=4) :: ilast
   logical :: terminate_execution, gsf_negative
   character(len=132) :: tco_file
   integer(kind=4) :: len_tco
   real(kind=8) Ecut
   integer(kind=4) :: iproj, itarget
   integer(kind=4) :: iZ, iA
   logical :: tco_exist
   integer(kind=4) :: idummy, state_type, nread, ncut
   real(kind=8) :: dummy, state_e, Emin_DWBA
   real(kind=8) :: sum, spin, sg2cut
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

   if(terminate_execution)call exit_YAHFC(1)
!----------------------------------------------------------------------
!----   Check if DWBA calculation has been performed and that the lowest
!----   DWBA state to continuous energy bins is above Ecut for the target
!----   if not, reset Ecut for the target.
   if(do_dwba)then
      iproj = projectile%particle_type
      itarget = target%icomp
      iZ = nucleus(1)%Z - particle(iproj)%Z
      iA = nucleus(1)%A - particle(iproj)%A
!----   check that it is the target
      if(iZ /= nucleus(itarget)%Z .or. iA /= nucleus(itarget)%A)then
          write(6,*)'Not the target nucleus'
          call exit_YAHFC(901)
      end if
      tco_file(1:100) = ' '
      tco_file(1:1) = particle(iproj)%label
      len_tco = 1
      if(nucleus(itarget)%atomic_symbol(1:1) == ' ')then
         tco_file(len_tco+1:len_tco+1) = nucleus(itarget)%atomic_symbol(2:2)
         len_tco = len_tco+1
      else
         tco_file(len_tco+1:len_tco+2) = nucleus(itarget)%atomic_symbol(1:2)
         len_tco = len_tco+2
      end if
      if(nucleus(itarget)%A < 10)then
         write(tco_file(len_tco+1:len_tco+1),'(i1)')nucleus(itarget)%A
         len_tco = len_tco + 1
      elseif(nucleus(itarget)%A < 100)then
         write(tco_file(len_tco+1:len_tco+2),'(i2)')nucleus(itarget)%A
         len_tco = len_tco + 2
      elseif(nucleus(itarget)%A < 1000)then
         write(tco_file(len_tco+1:len_tco+3),'(i3)')nucleus(itarget)%A
         len_tco = len_tco + 3
      end if

      tco_exist = .false.
      inquire(file=tco_file(1:len_tco)//'-CC.data',exist=tco_exist)

!      write(6,*)'tco_exist = ',tco_exist

      if(tco_exist)then
         open(unit=53,file=tco_file(1:len_tco)//'-CC.data',status='unknown')
         read(53,*)idummy, nread
         Emin_DWBA = 1.0d7
         do i = 1, nread
            read(53,*)idummy, idummy, dummy, dummy, dummy, state_e, state_type
            if(state_type < 1 .and. state_e < Emin_DWBA)Emin_DWBA = state_e
         end do
         close(unit=53)
         if(Emin_DWBA < nucleus(itarget)%level_param(7))then
!   write(6,*)'Adjusting Ecut'
            ncut = nucleus(itarget)%ncut
            do ii = nucleus(itarget)%ncut, 1, -1
               if(nucleus(itarget)%state(ii)%energy < Emin_DWBA - 0.5d0*Delt)then
                  ncut = ii
                  exit
               end if
            end do
            nucleus(itarget)%ncut = ncut
            Ecut = nucleus(itarget)%state(ncut)%energy + 0.001d0
            nucleus(itarget)%level_param(7) = Ecut
            nucleus(itarget)%level_ecut = Ecut
            sg2cut = 0.0d0
            sum = 0.0d0
            do i = 1, ncut
               spin = nucleus(itarget)%state(i)%spin
               sg2cut = sg2cut + spin*(spin + 1.0d0)*(2.0d0*spin+1.0d0)
               sum = sum + (2.0d0*spin+1.0d0)
            end do
            sg2cut = sg2cut/(2.0d0*sum)      !   standard model,const. sig**2, integrate over J
            nucleus(itarget)%level_param(8) = sg2cut
!   write(6,*)'new ecut = ',Ecut
!   write(6,*)itarget, iZ,iA
         end if
      end if

   end if

   return

end subroutine fail_safe_check
