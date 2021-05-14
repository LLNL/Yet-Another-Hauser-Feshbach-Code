!
!*******************************************************************************
!
subroutine Fission_data(data_path,len_path,icomp)                
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine to set up fission data for nuclei in the calculation
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        nodeinfo
!        constants
!
!     Subroutines:
!
!       init_barrier_data
!       parse_string
!       lower_case_word
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
!*******************************************************************************
!
   use variable_kinds
   use options
   use nuclei
   use nodeinfo
   use constants
   implicit none
   character(len=200), intent(in) :: data_path
   integer(kind=4), intent(in) :: len_path
   integer(kind=4), intent(in) :: icomp
!-----------------------------------------------------------
   integer(kind=4) :: ia
   real(kind=8) :: Ea, Eb, Delf

   character(len=2) :: symb, symma, symmb
   integer(kind=4) :: izz, iaa
   integer(kind=4) :: i, j

   logical :: finished, reading
   character(len=132) :: command
   integer(kind=4) :: numw
   integer(kind=4) :: startw(66), stopw(66)
   integer(kind=4) :: read_err
   integer(kind=4) :: nchar
!-------------------------------------------------------------------

   ia = nucleus(icomp)%A

   nucleus(icomp)%fission = .false.

!----    First set up fission parameters based on data in 'Fission-barrier.dat'
!----    Then check for data in 'Fission-Parameters.txt', which contains
!----    user specified parameters, likely due to previous fitting
!----
   open(unit=53,file=data_path(1:len_path)//'Fission-barrier.dat',status='old')
   do i=1,4
      read(53,*)
   end do
 1 read(53,'(2(1x,i3),1x,a2,2(f8.3,2x,a2),f8.3)',end=2)    &
              izz,iaa,symb,Ea,symma,Eb,symmb,Delf

   if(izz == nucleus(icomp)%Z .and. iaa == nucleus(icomp)%A)then

      nucleus(icomp)%fission = .true.
      nucleus(icomp)%F_n_barr = 2
      if(Ea <= 0.001)then
         nucleus(icomp)%fission = .true.
         nucleus(icomp)%F_n_barr = 1
         if(.not.allocated(nucleus(icomp)%F_Barrier))allocate(nucleus(icomp)%F_Barrier(1))
         nucleus(icomp)%F_barrier(1)%barrier = Eb



         if(symmb == 'S ')then
             nucleus(icomp)%F_barrier(1)%symmetry = 1
         elseif(symmb == 'GA')then
             nucleus(icomp)%F_barrier(1)%symmetry = 2
         elseif(symmb == 'MA')then
             nucleus(icomp)%F_barrier(1)%symmetry = 3
         end if

         if(iand(iaa,1) == 0)then
            if(iand(izz,1) == 0)then
               nucleus(icomp)%F_Barrier(1)%hbw = 0.9d0
            else
               nucleus(icomp)%F_Barrier(1)%hbw = 0.6d0
            end if
         else
            nucleus(icomp)%F_Barrier(1)%hbw = 0.8d0
         end if
      else
         nucleus(icomp)%fission = .true.
         nucleus(icomp)%F_n_barr = 2
         if(.not.allocated(nucleus(icomp)%F_Barrier))allocate(nucleus(icomp)%F_Barrier(2))
         nucleus(icomp)%F_Barrier(1)%barrier = Ea
         nucleus(icomp)%F_Barrier(2)%barrier = Eb
         if(symma == 'S ')then
             nucleus(icomp)%F_barrier(1)%symmetry = 1
         elseif(symma == 'GA')then
             nucleus(icomp)%F_barrier(1)%symmetry = 3
         elseif(symma == 'MA')then
             nucleus(icomp)%F_barrier(1)%symmetry = 2
         end if
         if(symmb == 'S ')then
             nucleus(icomp)%F_barrier(2)%symmetry = 1
         elseif(symmb == 'GA')then
             nucleus(icomp)%F_barrier(2)%symmetry = 3
         elseif(symmb == 'MA')then
             nucleus(icomp)%F_barrier(2)%symmetry = 2
         end if
         if(iand(iaa,1) == 0)then
            if(iand(izz,1) == 0)then
               nucleus(icomp)%F_Barrier(1)%hbw = 0.9d0
               nucleus(icomp)%F_Barrier(2)%hbw = 0.6d0
            else
               nucleus(icomp)%F_Barrier(1)%hbw = 0.6d0
               nucleus(icomp)%F_Barrier(2)%hbw = 0.4d0
            end if
         else
            nucleus(icomp)%F_Barrier(1)%hbw = 0.8d0
            nucleus(icomp)%F_Barrier(2)%hbw = 0.5d0
         end if
      end if

      call init_barrier_data(icomp)

   end if
   goto 1
 2 continue
   close(unit=53)

!------   Now check the file 'Fission-Parameters.txt'
   open(unit=53,file=data_path(1:len_path)//'Fission-Parameters.txt',status='old')
   finished = .false.
   reading = .false.
   do while(.not. finished)
      read(53,'(a)', iostat = read_err)command
      if(read_err /= 0)finished = .true.                                     !  Hit end of file
      call parse_string(command,numw,startw,stopw)
      if(command(1:1) == '!' .or. command(1:1) == '#')cycle                  !  it is a comment line
      nchar = stopw(numw)

      call lower_case_word(nchar,command(1:stopw(numw)))                     !  convert to lower case
      if(command(startw(1):stopw(1)) == 'end')finished = .true.
      if(.not. reading)then
         if(command(startw(1):stopw(1)) == 'start')then
            read(command(startw(2):stopw(2)),*)izz
            read(command(startw(3):stopw(3)),*)iaa
            if(izz == nucleus(icomp)%Z .and. iaa == nucleus(icomp)%A)then    !  Found nucleus in file
               reading = .true.
               read(command(startw(4):stopw(4)),*)j
               nucleus(icomp)%F_n_barr = j       
               nucleus(icomp)%fission = .true.
               if(allocated(nucleus(icomp)%F_Barrier))deallocate(nucleus(icomp)%F_Barrier)
               allocate(nucleus(icomp)%F_Barrier(nucleus(icomp)%F_n_barr))
               call init_barrier_data(icomp)
            end if
         end if
      else
         if(command(startw(1):stopw(1)) == 'start')exit                      !  hit the next nucleus, so we are done
         call fission_command(command, numw, startw, stopw, icomp)           !  parse data and put in proper place
         nucleus(icomp)%fission_read = .true.
      end if  
   end do
   close(unit=53)
   return
end subroutine Fission_data
!
!*******************************************************************************
!
subroutine Fission_levels(icomp)                
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine to fit matching energy given transition states above the 
!    fission barrier
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        nodeinfo
!
!     Subroutines:
!
!       find_T_E0
!       cumm_rho
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
!*******************************************************************************
!
   use variable_kinds
   use options
   use nuclei
   use nodeinfo
   implicit none
   integer(kind=4), intent(in) :: icomp
   integer(kind=4) :: ia
   real(kind=8) :: spin, sg2cut, sum
   real(kind=8) :: chisq, chi_min
   real(kind=8) :: ematch
   integer(kind=4) :: nfit
   real(kind=8), allocatable :: elv(:)
   real(kind=8), allocatable :: cum_rho(:),dcum_rho(:)
   real(kind=8), allocatable :: cum_fit(:)

   integer(kind=4) i,j, num

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ia = nucleus(icomp)%A


!-------------------------------------------------------------------
   do i = 1, nucleus(icomp)%F_n_barr

!----   Find Temperature and E0 at start -  later fit discrete
      call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,           &
                     nucleus(icomp)%F_Barrier(i)%vib_enh,                  &
                     nucleus(icomp)%F_Barrier(i)%rot_enh)


!-----   Spin cut off parameter at Ecut
!-----   Taken from Reffo based on maximum liklihood arguments
      sg2cut = (0.83*real(nucleus(icomp)%A,kind=8)**0.26)**2
      num = nucleus(icomp)%F_Barrier(i)%num_discrete


      if(num > 0)then
         sg2cut=0.
         sum = 0.0
         do j = 1, num
            spin = nucleus(icomp)%F_Barrier(i)%state_j(j)
            sg2cut = sg2cut + spin*(spin + 1.0)*(2.0*spin+1.0)
            sum = sum + (2.0*spin+1.0)
         end do
         sg2cut = sg2cut/(3.0d0*sum)
      else
       sg2cut = (0.83*real(nucleus(icomp)%A,kind=8)**0.26)**2
      end if
      if(sg2cut < 1.0d-4)sg2cut =                                          &
           (0.83*real(nucleus(icomp)%A,kind=8)**0.26)**2

      nucleus(icomp)%F_Barrier(i)%level_param(9) = nucleus(icomp)%level_param(9)
      nucleus(icomp)%F_Barrier(i)%ecut = 0.0001
      nucleus(icomp)%F_Barrier(i)%level_param(8) = sg2cut
      if(num > 0)nucleus(icomp)%F_Barrier(i)%ecut = nucleus(icomp)%F_Barrier(i)%state_e(num) + 0.001
      nucleus(icomp)%F_Barrier(i)%level_param(7) = nucleus(icomp)%F_Barrier(i)%ecut
      nfit = nucleus(icomp)%F_Barrier(i)%num_discrete
      if(nfit <= 5)cycle                                           !  There are so few discrete states, no point in fitting to the cumulative level density
      if(.not.allocated(elv))allocate(elv(nfit))                   !  allocate cumulative density array
      if(.not.allocated(cum_rho))allocate(cum_rho(nfit))           !  allocate cumulative density array
      if(.not.allocated(dcum_rho))allocate(dcum_rho(nfit))
      if(.not.allocated(cum_fit))allocate(cum_fit(nfit))           !  allocate cumulative density array
      do j = 1, nfit                                               !  Make cumulative rho array
         elv(j) = nucleus(icomp)%F_Barrier(i)%state_e(j)
         cum_rho(j)=real(j,kind=8)
         dcum_rho(j)=1.0d0/sqrt(cum_rho(j))                        !  assume error 1/sqrt
      end do

      ematch = 8.0
      nucleus(icomp)%F_Barrier(i)%level_param(6) = ematch
      call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,            &
                     nucleus(icomp)%F_Barrier(i)%vib_enh,                   &
                     nucleus(icomp)%F_Barrier(i)%rot_enh)

      call cumm_rho(nfit,elv,nucleus(icomp)%A,nucleus(icomp)%level_param,   &
                    nucleus(icomp)%vib_enh,                                 &
                    nucleus(icomp)%rot_enh,cum_fit)
  
      chisq=0.0d0

      do j = 1, nfit
         chisq = chisq + (cum_rho(j)-cum_fit(j))**2/dcum_rho(j)**2
      end do
      ematch = nucleus(icomp)%F_Barrier(i)%level_param(6)
      chi_min = 1.0d10

 10   nucleus(icomp)%F_Barrier(i)%level_param(6) =                          &
              nucleus(icomp)%F_Barrier(i)%level_param(6) - 0.05
      call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,            &
                     nucleus(icomp)%F_Barrier(i)%vib_enh,                   &
                     nucleus(icomp)%F_Barrier(i)%rot_enh)
      if(nucleus(icomp)%F_Barrier(i)%level_param(6) < 3.0)goto 11
      call cumm_rho(nfit,elv,nucleus(icomp)%A,nucleus(icomp)%level_param,   &
                    nucleus(icomp)%vib_enh,                                 &
                    nucleus(icomp)%rot_enh,cum_fit)
      chisq=0.0d0
      do j=1,nfit
         chisq = chisq + (cum_rho(j)-cum_fit(j))**2/dcum_rho(i)**2
      end do


      if(chisq < chi_min)then
         ematch = nucleus(icomp)%F_Barrier(i)%level_param(6)
         chi_min = chisq
      else
         goto 11
      end if
      goto 10
 11   nucleus(icomp)%F_Barrier(i)%level_param(6) = ematch
      call Find_T_E0(ia,nucleus(icomp)%F_Barrier(i)%level_param,           &
                     nucleus(icomp)%F_Barrier(i)%vib_enh,                  &
                     nucleus(icomp)%F_Barrier(i)%rot_enh)
      if(allocated(elv))deallocate(elv)
      if(allocated(cum_rho))deallocate(cum_rho)
      if(allocated(dcum_rho))deallocate(dcum_rho)
      if(allocated(cum_fit))deallocate(cum_fit)
   end do

   return

end subroutine Fission_levels
!
!*******************************************************************************
!
subroutine Fission_transmission(icomp,Ex,xji,ipar,F_trans)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine to compute the fission transmission coefficients at
!    excitation energy ex, angular moment xji, and parity ipar 
!    for nucleus icomp
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        options
!        nuclei
!        nodeinfo
!        constants
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!      real(kind=8) :: HW_trans
!      real(kind=8) :: spin_fac
!      real(kind=8) :: parity_fac
!      logical :: real8_equal
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
   use nodeinfo
   use constants
   implicit none
   integer(kind=4), intent(in) :: icomp
   real(kind=8), intent(in) :: Ex, xji
   integer(kind=4), intent(in) :: ipar
   real(kind=8), intent(out) :: F_trans(4)
!-----------------------------------------------------
   real(kind=8) :: trans(3)
   real(kind=8) :: de
   real(kind=8) :: T_f
   integer(kind=4) :: i, j, j_min
   real(kind=8) :: Ef
   real(kind=8) :: apar,spin_cut,delta,ecut,sg2cut,sig_mod,ematch
   real(kind=8) :: par
   real(kind=8) :: rho,rho_FM,sig2,apu
   real(kind=8) :: tt
   real(kind=8) :: jfac,pfac
   real(kind=8) :: K_vib, K_rot
   integer(kind=4) :: ia
   real(kind=8) :: e_old
   logical converged
   real(kind=8) :: pmode, pe1, pbb
   save e_old

   real(kind=8) :: E0, T, E1
   real(kind=8) :: b, Max_J, F_Barrier, F_Barrier_hbw
   real(kind=8) :: aa, bb, cc
!--------   External Functions  ----------------------
   real(kind=8) :: HW_trans
   real(kind=8) :: spin_fac
   real(kind=8) :: parity_fac
   logical :: real8_equal
!-----------------------------------------------------

   par = 1.0d0
   if(ipar == 0)par = -1.0d0

   de = 0.005
   ia = nucleus(icomp)%A
   do i = 1, nucleus(icomp)%F_n_barr
      apar = nucleus(icomp)%F_barrier(i)%level_param(1)
      spin_cut = nucleus(icomp)%F_barrier(i)%level_param(2)
      delta = nucleus(icomp)%F_barrier(i)%level_param(3)
      ecut = nucleus(icomp)%F_Barrier(i)%ecut
      sg2cut = nucleus(icomp)%F_barrier(i)%level_param(8)
      sig_mod = nucleus(icomp)%F_barrier(i)%level_param(9)
      ematch = nucleus(icomp)%F_barrier(i)%level_param(6)
      E0 = nucleus(icomp)%F_barrier(i)%level_param(15)
      T = nucleus(icomp)%F_barrier(i)%level_param(14)
      if(E0 < 0.0d0)then
         E1 = T*log(1.0d0-exp(E0/T))
      else
         E1 = -5.0d0
      end if
      j_min = nint(E1/de)

      aa = nucleus(icomp)%F_Barrier(i)%barrier_damp(1)
      bb = nucleus(icomp)%F_Barrier(i)%barrier_damp(2)
      cc = nucleus(icomp)%F_Barrier(i)%barrier_damp(3)

      trans(i)=0.0d0


      Max_j = nucleus(icomp)%F_barrier(i)%Max_J
      F_Barrier = nucleus(icomp)%F_Barrier(i)%barrier
      F_Barrier = F_Barrier*aa*exp(-cc**2*(Ex-bb)**2)
      if(Max_J > 0.0d0)then
         b = 1.0d0/(Max_J*(Max_J+1.0d0))
         if(xji <= Max_J)then
            F_Barrier = F_Barrier*(1.0d0 - b*xji*(xji+1.0d0))
         else
            F_Barrier = 0.1d0
         end if
      end if
      F_Barrier = max(F_Barrier,0.1d0)
      F_Barrier_hbw = nucleus(icomp)%F_Barrier(i)%hbw


      do j = 1,nucleus(icomp)%F_Barrier(i)%num_discrete
         if(real8_equal(nucleus(icomp)%F_barrier(i)%state_j(j),xji) .and.   &
            real8_equal(nucleus(icomp)%F_barrier(i)%state_pi(j),par))then
            Ef = nucleus(icomp)%F_barrier(i)%state_e(j)
            tt = HW_trans(Ex,Ef,F_Barrier,F_barrier_hbw)
            trans(i) = trans(i) + tt
         end if
      end do

      pmode = nucleus(icomp)%F_barrier(i)%level_param(16)
      pe1 = nucleus(icomp)%F_barrier(i)%level_param(17)
      pbb = nucleus(icomp)%F_barrier(i)%level_param(18)


      T_f = 0.0d0
      if(ecut > 0.0d0)then
         j = 0
      else
         j = j_min
      end if
      converged = .false.
      do while (.not. converged)
         Ef = real(j,kind=8)*de + ecut  
         call rhoe(Ef,nucleus(icomp)%F_barrier(i)%level_param,                 &
                   nucleus(icomp)%F_barrier(i)%vib_enh,                        &
                   nucleus(icomp)%F_barrier(i)%rot_enh,                        &
                   ia,rho_Fm,apu,sig2,K_vib,K_rot)
!
         jfac = spin_fac(xji,sig2)
         pfac = parity_fac(Ef,xji,ipar,pmode,pe1,pbb)
         rho = rho_FM*jfac*pfac
         tt = HW_trans(Ex,Ef,F_Barrier,F_barrier_hbw)
         tt = tt*rho*de
         if(T_f > 0.0 .and. tt/T_f < 1.0d-6)converged = .true.
         if(tt < 1.0d-7)converged = .true.
         T_f = T_f + tt
         j = j + 1
      end do


      trans(i) = trans(i) + T_f
      if(trans(i) < 1.0d-20)trans(i)=1.0d-20
      F_trans(i) = trans(i)
   end do

   if(nucleus(icomp)%F_n_barr == 1)then
      F_trans(4) = F_trans(1)
   else if(nucleus(icomp)%F_n_barr == 2)then
      F_trans(4) = F_trans(1)*F_trans(2)/(F_trans(1)+F_trans(2))
   else if(nucleus(icomp)%F_n_barr == 3)then
      F_trans(4) = F_trans(1)*(F_trans(2) + F_trans(3))/(F_trans(1) + F_trans(2) + F_trans(3))
   end if

   if(F_trans(4) <= 1.0d-20)F_trans(4) = 0.0d0
   e_old = ex
    return
end subroutine Fission_transmission
!
!*******************************************************************************
!
real(kind=8) function HW_trans(Ex,ei,Barrier,hbw)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the Hill-Wheeler transmission coefficient for 
!    fission barrier
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        constants
!
!     Subroutines:
!
!        None
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
!*******************************************************************************
!
   use variable_kinds
   use constants
   implicit none
   real(kind=8), intent(in) :: Ex, ei, Barrier, hbw
!-------------------------------------------------------------------
   real(kind=8) :: hbwx
   real(kind=8) :: exponent
!-------------------------------------------------------------------
   hbwx = max(hbw,1.0d-2)
   exponent = -2.0d0*pi*(Ex-Barrier-ei)/hbwx
   if(exponent <= -100.0d0)exponent = -100.0d0
   if(exponent >= 100.0d0)exponent = 100.0d0
   HW_trans = 1.0d0/(1.0d0 + exp(exponent))
end function HW_trans
!
!*******************************************************************************
!
subroutine init_barrier_data(i)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine initializes fission parameter info for nucleus i
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nodeinfo
!        options
!        useful_data
!        nuclei
!        particles_def
!        constants
!
!     Subroutines:
!
!        None
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
!*******************************************************************************
!
   use variable_kinds
   use nodeinfo
   use options
   use useful_data
   use nuclei
   use particles_def
   use constants
   implicit none
   integer(kind=4), intent(in) :: i
!-----------------------------------------------------------------------------
   integer(kind=4) :: j
   real(kind=8) :: beta_2
   real(kind=8) :: sig2_perp
   real(kind=8) :: sig2_ax
!-----------------------------------------------------------------------------

   do j = 1, nucleus(i)%F_n_barr
      if(Fiss_Max_J > 0.0d0)then
         nucleus(i)%F_Barrier(j)%Max_J = Fiss_Max_J
         if(iand(nucleus(i)%A,1) == 1)nucleus(i)%F_Barrier(j)%Max_J =        &
                                      nucleus(i)%F_Barrier(j)%Max_J + 0.5d0
      else
         nucleus(i)%F_Barrier(j)%Max_J = 0.0d0
      end if 
!---- Don't reset barrier heights and widths if already set
      if(nucleus(i)%F_Barrier(j)%barrier < 1.0d-1)nucleus(i)%F_Barrier(j)%barrier = 6.0d0
      if(nucleus(i)%F_Barrier(j)%hbw < 1.0d-1)nucleus(i)%F_Barrier(j)%hbw = 0.6d0
      nucleus(i)%F_Barrier(j)%ecut = 0.0d0
!----   If barrier symmetry has not been defined, then make guess for barrier symmetries.
!----   If barrier symmetry is not defined for this call, it is because user has 
!----   reset nucleus(i)%F_n_barr and the data has been reallocated. It is expected that 
!----   the user will reset the barrier symmetry as well. This is initialize based on
!----   the behavior of many actinides.
      if(nucleus(i)%F_Barrier(j)%symmetry == 0)then
         if(nucleus(i)%F_n_barr == 1)then
            nucleus(i)%F_Barrier(j)%symmetry = 1               !  Axially symmetric
         else
            if(j == 1)then
               nucleus(i)%F_Barrier(j)%symmetry = 1                            !  Axially symmetric
               if(nucleus(i)%A >= 144)nucleus(i)%F_Barrier(j)%symmetry = 3     !  Axially asymmetric
            end if
            if(j == 2)nucleus(i)%F_Barrier(j)%symmetry = 2     !  Mass asymmetric
         end if
      end if

      beta_2 = 2.5d0*nucleus(i)%beta(2)
      if(j == 2)then
         beta_2 = 4.5d0*nucleus(i)%beta(2)
      end if
      nucleus(i)%F_Barrier(j)%beta_2 = beta_2

      sig2_perp = (1.0d0 + beta_2/3.0d0)
      sig2_ax = sqrt(pi/2.0d0)*(1.0d0 - 2.0d0*beta_2/3.0d0)

      nucleus(i)%F_Barrier(j)%level_param(1) = nucleus(i)%level_param(1)
      nucleus(i)%F_Barrier(j)%level_param(2) = nucleus(i)%level_param(2)
      nucleus(i)%F_Barrier(j)%level_param(3) = nucleus(i)%level_param(3)
      if(j == 1)nucleus(i)%F_Barrier(j)%level_param(4) = 1.5d0
      if(j == 2)nucleus(i)%F_Barrier(j)%level_param(4) = 0.6d0
      nucleus(i)%F_Barrier(j)%level_param(5) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(6) = nucleus(i)%level_param(6)
      nucleus(i)%F_Barrier(j)%level_param(7) = nucleus(i)%F_Barrier(j)%ecut
      nucleus(i)%F_Barrier(j)%level_param(8) = (0.83*real(nucleus(i)%A,kind=8)**0.26)**2

      nucleus(i)%F_Barrier(j)%level_param(9) = 1.0d0
      nucleus(i)%F_Barrier(j)%level_param(10) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(11) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(12) = nucleus(i)%level_param(12)
      nucleus(i)%F_Barrier(j)%level_param(13) = nucleus(i)%level_param(13)
      nucleus(i)%F_Barrier(j)%level_param(14) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(15) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(16) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(17) = 0.0d0
      nucleus(i)%F_Barrier(j)%level_param(18) = 0.0d0

      nucleus(i)%F_barrier(j)%barrier_damp(1) = 1.0d0
      nucleus(i)%F_barrier(j)%barrier_damp(2) = 0.0d0
      nucleus(i)%F_barrier(j)%barrier_damp(3) = 0.0d0

      if(nucleus(i)%lev_option == 2)then
         nucleus(i)%F_Barrier(j)%level_param(10) = real(nucleus(i)%F_Barrier(j)%symmetry,kind=8)

         nucleus(i)%F_Barrier(j)%level_param(11) = nucleus(i)%level_param(11)
            
         nucleus(i)%F_Barrier(j)%rot_enh(1) = 1.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(2) = 30.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(3) = 5.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(4) = sig2_perp
         nucleus(i)%F_Barrier(j)%rot_enh(5) = sig2_ax
         nucleus(i)%F_Barrier(j)%vib_enh(1) = 1.0d0
         nucleus(i)%F_Barrier(j)%vib_enh(2) = 30.0d0
         nucleus(i)%F_Barrier(j)%vib_enh(3) = 5.0d0
      else
         nucleus(i)%F_Barrier(j)%rot_enh(1) = 0.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(2) = 30.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(3) = 5.0d0
         nucleus(i)%F_Barrier(j)%rot_enh(4) = sig2_perp
         nucleus(i)%F_Barrier(j)%rot_enh(5) = sig2_ax
         nucleus(i)%F_Barrier(j)%vib_enh(1) = 0.0d0
         nucleus(i)%F_Barrier(j)%vib_enh(2) = 30.0d0
         nucleus(i)%F_Barrier(j)%vib_enh(3) = 5.0d0
      end if
      nucleus(i)%F_Barrier(j)%num_discrete = 0
   end do

   return
end subroutine init_barrier_data
!
!*******************************************************************************
!
subroutine fission_command(command, numw, startw, stopw, icomp)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine parses the input command, interprets it, and puts
!    fission parameter data into the proper data slot.
!
!   Dependencies:
!
!     Modules:
!
!        variable_kinds
!        nodeinfo
!        options
!        print_control
!        useful_data
!        nuclei
!        particles_def
!        directory_structure
!        constants
!        pre_equilibrium_no_1
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        None
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
!*******************************************************************************
!
   use variable_kinds
   use nodeinfo
   use options
   use print_control
   use useful_data
   use nuclei
   use particles_def
   use directory_structure
   use constants
   use pre_equilibrium_no_1
   implicit none
!-------------------------------------------------------------------
   character(len=132), intent(in) :: command
   integer(kind=4), intent(in) :: numw
   integer(kind=4), intent(in) :: startw(66), stopw(66)
   integer(kind=4), intent(in) :: icomp
!-------------------------------------------------------------------
   integer(kind=4) :: iZ, iA
   integer(kind=4) :: j, k
   real(kind=8) :: x(66)
   integer(kind=4) :: ilast
   real(kind=8) :: beta_2
   real(kind=8) :: sig2_ax, sig2_perp

   integer(kind=4) :: num
   character(len=50) :: read_file
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_num_barrier')return
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_barrier')then          !   global setting of this parameter
      if(numw < 6)then
         if(iproc == 0)write(6,*)'Error in input for option "f_barrier"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      read(command(startw(6):stopw(6)),*)x(2)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'too many barriers for F_Barrier'
            call MPI_Abort(icomm,101,ierr)
         end if
         nucleus(icomp)%F_Barrier(j)%barrier = x(1)
         nucleus(icomp)%F_Barrier(j)%hbw = x(2)
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_Barrier'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    damping factor for fission barriers
!    F_barrier = F_barrier*x(1)*exp(-x(3)**2*(Ex-x())**2)
!   
   if(command(startw(1):stopw(1)) == 'f_barrier_damp')then          !   global setting of this parameter
      if(numw < 6)then
         if(iproc == 0)write(6,*)'Error in input for option "f_barrier_damp"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(2)
      read(command(startw(6):stopw(6)),*)x(3)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'too many barriers for F_Barrier'
            call MPI_Abort(icomm,101,ierr)
         end if
         nucleus(icomp)%F_Barrier(j)%barrier_damp(2) = x(2)
         nucleus(icomp)%F_Barrier(j)%barrier_damp(3) = x(3)
!-------   Make damping factor = 1.0 at Ex = 0.0
         x(1) = exp((x(3)*x(2))**2)
         nucleus(icomp)%F_Barrier(j)%barrier_damp(1) = x(1)
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_Barrier_damp'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    Set symmetry for barrier
!   
   if(command(startw(1):stopw(1)) == 'f_barrier_symmetry')then          !   global setting of this parameter
      if(numw < 5)then
         if(iproc == 0)write(6,*)'Error in input for option "f_barrier_symmetry"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'Error: requesting too many barriers for F_Barrier for nucleus Z = ',iZ,' A = ',iA
            call MPI_Abort(icomm,101,ierr)
         end if
         if(command(startw(5):stopw(5)) == 's' .or. command(startw(5):stopw(5)) == '1')then
            nucleus(icomp)%F_Barrier(j)%symmetry = 1
            nucleus(icomp)%F_Barrier(j)%level_param(10) = real(nucleus(icomp)%F_Barrier(j)%symmetry,kind=8)
            return
         end if
         if(command(startw(5):stopw(5)) == 'lr-a' .or. command(startw(5):stopw(5)) == '2')then
            nucleus(icomp)%F_Barrier(j)%symmetry = 2
            nucleus(icomp)%F_Barrier(j)%level_param(10) = real(nucleus(icomp)%F_Barrier(j)%symmetry,kind=8)
            if(j == 1 .and. iproc == 0)write(6,*)'WARNING!!!! ----  Setting first barrier to left-right asymmetric'
            return
         end if
         if(command(startw(5):stopw(5)) == 'ta-lr' .or. command(startw(5):stopw(5)) == '3')then
            nucleus(icomp)%F_Barrier(j)%symmetry = 3
            nucleus(icomp)%F_Barrier(j)%level_param(10) = real(nucleus(icomp)%F_Barrier(j)%symmetry,kind=8)
            return
         end if
         if(command(startw(5):stopw(5)) == 'ta-nlr' .or. command(startw(5):stopw(5)) == '4')then
            nucleus(icomp)%F_Barrier(j)%symmetry = 4
            nucleus(icomp)%F_Barrier(j)%level_param(10) = real(nucleus(icomp)%F_Barrier(j)%symmetry,kind=8)
            if(j == 1 .and. iproc == 0)write(6,*)'WARNING!!!! ----  Setting first barrier to triaxial no left-right asymmetry'
            return
         end if
      end if
      if(iproc == 0)then
         write(6,*)'Nucleus not found: f_barrier_symmetry'
         write(6,*)'Keeping the default value for this nucleus and barrier'
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_ecut')then          !   global setting of this parameter
      if(numw < 5)then
         if(iproc == 0)write(6,*)'Error in input for option "f_ecut"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'too many barriers for F_ecut'
            call MPI_Abort(icomm,101,ierr)
         end if
         nucleus(icomp)%F_Barrier(j)%ecut = x(1)
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_ecut'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_aparam')then          !   global setting of this parameter
      if(numw < 5)then
         if(iproc == 0)write(6,*)'Error in input for option "f_lev_aparam"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'too many barriers for F_lev_aparam'
            call MPI_Abort(icomm,101,ierr)
         end if
         nucleus(icomp)%F_Barrier(j)%level_param(1) = x(1)
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_aparam'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_spin')then          !   global setting of this parameter
      if(numw < 5)then
         if(iproc == 0)write(6,*)'Error in input for option "f_lev_spin"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'too many barriers for F_lev_spin'
            call MPI_Abort(icomm,101,ierr)
         end if
         nucleus(icomp)%F_Barrier(j)%level_param(2)=x(1)
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_spin'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(command(startw(1):stopw(1)) == 'f_lev_delta')then          !   global setting of this parameter
      if(numw < 5)then
         if(iproc == 0)write(6,*)'Error in input for option "f_;ev_delta"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'too many barriers for F_lev_delta'
            call MPI_Abort(icomm,101,ierr)
         end if
         nucleus(icomp)%F_Barrier(j)%level_param(3) = x(1)
         nucleus(icomp)%F_Barrier(j)%level_param(6) = 2.5 + 150./real(iA,kind=8) + x(1)
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_delta'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_shell')then
      if(numw < 5)then
         if(iproc == 0)write(6,*)'Error in input for option "f_lev_shell"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'too many barriers for F_lev_shell'
            call MPI_Abort(icomm,101,ierr)
         end if
         nucleus(icomp)%F_Barrier(j)%level_param(4) = x(1)
         return
      end if
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_gamma')then
      if(numw < 5)then
         if(iproc == 0)write(6,*)'Error in input for option "f_lev_gamma"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'too many barriers for F_lev_gamma'
            call MPI_Abort(icomm,101,ierr)
         end if
         nucleus(icomp)%F_Barrier(j)%level_param(5) = x(1)
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_delta'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_rot_enhance')then
      if(numw < 7)then
         if(iproc == 0)write(6,*)'Error in input for option "f_lev_rot_enhance"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      read(command(startw(6):stopw(6)),*)x(2)
      read(command(startw(7):stopw(7)),*)x(3)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'Error in F_lev_rot_enhance index > # of barriers'
            call MPI_Abort(icomm,101,ierr)
         end if
         do k = 1, 3
            if(x(k) >= 0.0d0)nucleus(icomp)%F_Barrier(j)%rot_enh(k) = x(k)
         end do
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_rot_enhance'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_vib_enhance')then
      if(numw < 7)then
         if(iproc == 0)write(6,*)'Error in input for option "f_lev_vib_enhance"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      read(command(startw(6):stopw(6)),*)x(2)
      read(command(startw(7):stopw(7)),*)x(3)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'Error in F_vib_vib_enhance index > # of barriers'
            call MPI_Abort(icomm,101,ierr)
         end if
         do k = 1,3
            if(x(k) >= 0.0d0)nucleus(icomp)%F_Barrier(j)%vib_enh(k) = x(k)
         end do
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_vib_enhance'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_lev_ematch')then
      if(numw < 5)then
         if(iproc == 0)write(6,*)'Error in input for option "f_lev_ematch"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*)'too many barriers for F_lev_ematch'
            call MPI_Abort(icomm,101,ierr)
         end if
         nucleus(icomp)%F_Barrier(j)%level_param(6) = x(1)
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_ematch'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_beta2')then          !   global setting of this parameter
      if(numw < 5)then
         if(iproc == 0)write(6,*)'Error in input for option "f_beta2"'
         return
      end if
      read(command(startw(2):stopw(2)),*)iZ
      read(command(startw(3):stopw(3)),*)iA
      read(command(startw(4):stopw(4)),*)j
      read(command(startw(5):stopw(5)),*)x(1)
      if(iZ == nucleus(icomp)%Z .and. iA == nucleus(icomp)%A)then
         if(j > nucleus(icomp)%F_n_barr)then
            if(iproc == 0)write(6,*) 'WARNING -- too many barriers for F_beta_2 in nucleus ',icomp
            call MPI_Abort(icomm,101,ierr)
         end if
         beta_2 = x(1)
         nucleus(icomp)%F_Barrier(j)%beta_2 = beta_2
         sig2_perp = (1.0d0 + beta_2/3.0d0)
         sig2_ax = sqrt(pi/2.0d0)*(1.0d0 - 2.0d0*beta_2/3.0d0)
         nucleus(icomp)%F_Barrier(j)%rot_enh(4) = sig2_perp
         nucleus(icomp)%F_Barrier(j)%rot_enh(5) = sig2_ax
         return
      end if
      if(iproc == 0)write(6,*)'Nucleus not found: F_lev_aparam'
      return
   end if
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
   if(command(startw(1):stopw(1)) == 'f_barr_levels')then

      read_file(1:50) = ' '
      ilast = stopw(2) - startw(2) + 1

      read_file(1:ilast) = command(startw(2):stopw(2))

      read(command(startw(3):stopw(3)),*)iZ
      read(command(startw(4):stopw(4)),*)iA
      read(command(startw(5):stopw(5)),*)j
      if(nucleus(icomp)%Z == iZ .and. nucleus(icomp)%A == iA)then
         open(unit=8, file = read_file(1:ilast), status='old')
         num = 0
 9       read(8,*,end=10)
         num = num + 1
         goto 9
 10      rewind(8)
         nucleus(icomp)%F_Barrier(j)%num_discrete = num
         if(.not.allocated(nucleus(icomp)%F_barrier(j)%state_e))then
            allocate(nucleus(icomp)%F_barrier(j)%state_e(num))
         else
            deallocate(nucleus(icomp)%F_barrier(j)%state_e)
            allocate(nucleus(icomp)%F_barrier(j)%state_e(num))
         end if
         if(.not.allocated(nucleus(icomp)%F_barrier(j)%state_j))then
            allocate(nucleus(icomp)%F_barrier(j)%state_j(num))
         else
            deallocate(nucleus(icomp)%F_barrier(j)%state_j)
            allocate(nucleus(icomp)%F_barrier(j)%state_j(num))
         end if
         if(.not.allocated(nucleus(icomp)%F_barrier(j)%state_pi))then
            allocate(nucleus(icomp)%F_barrier(j)%state_pi(num))
         else
            deallocate(nucleus(icomp)%F_barrier(j)%state_pi)
            allocate(nucleus(icomp)%F_barrier(j)%state_pi(num))
         end if
         do k = 1, num
            read(8,*)nucleus(icomp)%F_Barrier(j)%state_e(k),                &
                         nucleus(icomp)%F_Barrier(j)%state_j(k),            &
                      nucleus(icomp)%F_Barrier(j)%state_pi(k)
         end do
         close(unit=8)
      end if

      if(iproc == 0)write(6,*)'Nucleus not found: F_Barr_levels'
      return
   end if
!---------------------------------------------------------------
   return
end subroutine fission_command


