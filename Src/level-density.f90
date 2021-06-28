!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_lev_den(data_path,len_path,symb,iz,ia,icomp)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine sets up level density data, parameters, array data
!    based on the selected level density option
!
!   Dependencies:
!
!     Modules:
!
!        options
!        nuclei
!        nodeinfo
!        constants
!        particles_def
!
!     Subroutines:
!
!        Find_T_E0
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
   use options
   use nuclei
   use nodeinfo
   use constants
   use particles_def
   implicit none
   character(len=200), intent(in) :: data_path      ! path where data files are kept
   integer(kind=4), intent(in)  :: len_path             ! length of data_path
   character(len=2), intent(in) :: symb
   integer(kind=4), intent(in) :: iz,ia
   integer(kind=4), intent(in) :: icomp
!----------------------------------------------------------------------
   character(len=5) nuc_name
   integer(kind=4) :: len_nuc
   integer(kind=4) :: max_nuc
   parameter (max_nuc=20000)     ! accomodate up to 100000 nuclei
   integer(kind=4) :: izz,iaa
   real(kind=8) :: x,y,x1,y1, xx, yy
   integer(kind=4) :: ix, ix1
!----------------------------------------------------------------------
   integer(kind=4) :: eof
   real(kind=8) :: aparam
   real(kind=8) :: spin_cut
   real(kind=8) :: delta
   real(kind=8) :: ematch
   real(kind=8) :: ecut
   real(kind=8) :: sg2cut
   real(kind=8) :: shell
   real(kind=8) :: gamma
   real(kind=8) :: afac
   real(kind=8) :: spin
   integer(kind=4) :: maxlev
   real(kind=8) :: b2,b3,b4,b6,beta2,beta3,beta4,beta6
   real(kind=8) :: xZ,xN,xA
   real(kind=8) :: Mass_exp,Mass_LDM
   character(len=132) :: input_line
   logical :: header
   logical :: found
   real(kind=8) :: xxx2, xxx4, xxx6
   integer(kind=4) :: n1, n2
!----------------------------------------------------------------------
   real(kind=8) :: sum
!---------------------------------------------------------------------
   real(kind=8) :: D0exp,dD0exp
   real(kind=8) :: Gamma_g_exp,dGamma_g_exp
   real(kind=8) :: D1exp,dD1exp
   real(kind=8) :: Gamma_g_1_exp,dGamma_g_1_exp
   integer(kind=4) :: i, in
      
!-----   Gibert & cameron parameters 
   real(kind=8) :: deltap(50),deltan(75)
   real(kind=8) :: sp(100),sn(150)
!----   Gilbert & Cameron pairing terms
   data deltap/0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,2.46d0,2.09d0,1.62d0,1.62d0,1.83d0,     &
               1.73d0,1.35d0,1.54d0,1.20d0,1.06d0,1.36d0,1.43d0,1.17d0,1.24d0,1.20d0,     &
               1.28d0,1.28d0,1.35d0,1.36d0,1.19d0,1.14d0,1.12d0,1.58d0,1.17d0,1.18d0,     &
               1.22d0,0.97d0,0.92d0,0.62d0,0.68d0,0.64d0,0.72d0,0.75d0,0.71d0,0.87d0,     &
               0.83d0,0.89d0,0.79d0,0.89d0,0.78d0,0.69d0,0.61d0,0.72d0,0.77d0,0.00d0/
   data deltan/0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,2.67d0,1.80d0,1.67d0,1.86d0,2.04d0,     &
               1.64d0,1.44d0,1.54d0,1.30d0,1.27d0,1.29d0,1.41d0,1.50d0,1.50d0,1.43d0,     &
               1.88d0,1.47d0,1.57d0,1.46d0,0.93d0,0.72d0,1.12d0,1.29d0,0.94d0,1.24d0,     &
               1.25d0,1.14d0,1.32d0,1.15d0,1.24d0,1.43d0,1.09d0,1.20d0,1.04d0,0.70d0,     &
               0.85d0,0.76d0,0.92d0,0.99d0,1.10d0,0.92d0,0.73d0,0.70d0,0.87d0,0.61d0,     &
               0.69d0,0.55d0,0.40d0,0.73d0,0.58d0,0.86d0,1.13d0,0.84d0,0.79d0,0.82d0,     &
               0.71d0,0.41d0,0.38d0,0.67d0,0.61d0,0.78d0,0.67d0,0.67d0,0.79d0,0.60d0,     &
               0.57,0.49,0.43,0.50,0.39/
!---- Gilbert & Cameron shell corrections
    data sp/ 0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,            &
             0.00d0,  0.00d0,  0.00d0,  0.00d0,  0.00d0,            &
            -2.91d0, -4.17d0, -5.72d0, -7.80d0, -8.97d0,            &
            -9.70d0,-10.10d0,-10.70d0,-11.38d0,-12.07d0,            &
           -12.55d0,-13.24d0,-13.93d0,-14.71d0,-15.53d0,            &
           -16.37d0,-17.36d0,-18.52d0,-18.44d0,-18.19d0,            &
           -17.68d0,-17.09d0,-16.65d0,-16.66d0,-16.59d0,            &
           -16.35d0,-16.18d0,-16.41d0,-16.60d0,-16.54d0,            &
           -16.42d0,-16.84d0,-17.22d0,-17.42d0,-17.52d0,            &
           -17.82d0,-18.19d0,-18.58d0,-19.11d0,-19.83d0,            &
           -19.14d0,-18.35d0,-17.40d0,-16.54d0,-15.68d0,            &
           -14.75d0,-13.71d0,-12.87d0,-12.18d0,-11.61d0,            &
           -11.09d0,-10.78d0,-10.53d0,-10.41d0,-10.21d0,            &
            -9.85d0, -9.36d0, -8.97d0, -8.56d0, -8.13d0,            &
            -7.68d0, -7.33d0, -7.11d0, -7.16d0, -7.05d0,            &
            -6.81d0, -6.56d0, -6.95d0, -7.52d0, -8.03d0,            &
            -8.41d0, -8.86d0, -7.71d0, -6.38d0, -5.47d0,            &
            -4.78d0, -4.37d0, -4.17d0, -4.12d0, -4.29d0,            &
            -4.61d0, -5.04d0, -5.48d0, -5.96d0, -6.40d0,            &
            -6.87d0, -7.20d0, -7.74d0,  0.00d0,  0.00d0/
   data sn/ 0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,            &
            0.00d0, 0.00d0, 0.00d0, 0.00d0, 0.00d0,            &
            6.80d0, 7.53d0, 7.55d0, 7.21d0, 7.44d0,            &
            8.07d0, 8.94d0, 9.81d0,10.60d0,11.39d0,            &
           12.54d0,13.68d0,14.34d0,14.19d0,13.83d0,            &
           13.50d0,13.00d0,12.13d0,12.60d0,13.26d0,            &
           14.13d0,14.92d0,15.60d0,16.38d0,17.08d0,            &
           17.55d0,17.98d0,18.33d0,18.56d0,18.71d0,            &
           18.65d0,18.55d0,18.52d0,18.34d0,18.01d0,            &
           17.38d0,16.56d0,15.62d0,14.38d0,12.88d0,            &
           13.24d0,13.71d0,14.40d0,15.16d0,15.89d0,            &
           16.43d0,16.97d0,17.59d0,18.08d0,18.72d0,            &
           19.22d0,19.51d0,19.73d0,19.91d0,20.06d0,            &
           20.16d0,20.09d0,19.83d0,19.41d0,19.06d0,            &
           18.66d0,17.73d0,17.03d0,16.44d0,16.00d0,            &
           15.33d0,14.49d0,13.42d0,12.28d0,11.14d0,            &
           10.10d0, 9.09d0,10.00d0,10.64d0,11.18d0,            &
           11.70d0,12.22d0,12.71d0,13.05d0,12.99d0,            &
           12.62d0,12.11d0,11.66d0,11.21d0,10.81d0,            &
           10.38d0,10.03d0, 9.65d0, 9.38d0, 8.99d0,            &
            8.62d0, 8.33d0, 8.10d0, 7.82d0, 7.56d0,            &
            7.33d0, 7.15d0, 6.83d0, 6.69d0, 6.55d0,            &
            6.53d0, 6.49d0, 6.39d0, 5.82d0, 5.26d0,            &
            4.53d0, 3.83d0, 3.08d0, 2.37d0, 1.72d0,            &
            1.05d0, 0.27d0,-0.69d0,-1.69d0,-2.58d0,            &
           -3.16d0,-1.72d0,-0.41d0, 0.71d0, 1.66d0,            &
            2.62d0, 3.22d0, 3.76d0, 4.10d0, 4.46d0,            &
            4.83d0, 5.09d0, 5.18d0, 5.17d0, 5.10d0,            &
            5.05d0, 5.04d0, 5.03d0, 4.99d0, 4.98d0,            &
            5.11d0, 5.27d0, 5.39d0, 5.37d0, 5.30d0/
!-----   Start calculation

!------------   Number of neutrons
   in = ia - iz
!------   set experimental D0 values to zero - remain zero if no exp.
!-----------   Find data file 
   if(ia < 9)then
      write(nuc_name(1:1),'(i1)')ia
      len_nuc = 1
   elseif(ia < 100)then
      write(nuc_name(1:2),'(i2)')ia
      len_nuc = 2
   else
      write(nuc_name(1:3),'(i3)')ia
      len_nuc = 3
   end if
   if(symb(1:1) == ' ')then
      write(nuc_name(len_nuc+1:len_nuc+1),'(a1)')symb(2:2)
      len_nuc = len_nuc+1
   else
      write(nuc_name(len_nuc+1:len_nuc+2),'(a2)')symb(1:2)
      len_nuc = len_nuc+2
   end if


!-----   First get experimental D0 and Gamma_gamma
   if(nucleus(icomp)%D0exp < 0.0d0)then
      open(unit=53,file=data_path(1:len_path)//'Level-Density-Data-RIPL3-L0.dat',   &
           status='old')
      D0exp = -1.0d0
      dD0exp = -1.0d0
      Gamma_g_exp = -1.0
      dGamma_g_exp = -1.0
      header = .true.
      do while(header)
         read(53,'(a)')input_line
         if(input_line(1:1) /= '#')then
            backspace(53)
            header = .false.
         end if
      end do
      eof = 0
      do i = 1,max_nuc
         read(53,'(1x,i2,4x,i3,1x,f4.1,1x,f7.3,2(1x,e9.2),1x,2(1x,f5.2),1x,i5,1x,i5)', &
                   iostat = eof)izz,iaa,xx,yy,x,x1,y,y1,ix,ix1
         if(eof == 1)exit
         if(izz == iz .and. iaa == ia-1)then
            D0exp = x
            dD0exp = x1
            Gamma_g_exp = real(ix)
            dGamma_g_exp = real(ix1)
            exit
         end if
      end do
      close(unit=53)
      D0exp = D0exp*1000.0d0
      dD0exp = dD0exp*1000.0d0
      nucleus(icomp)%D0exp = D0exp
      nucleus(icomp)%dD0exp = dD0exp
      nucleus(icomp)%Gamma_g_exp = Gamma_g_exp
      nucleus(icomp)%dGamma_g_exp = dGamma_g_exp
      if(D0exp <= 1.0d-6 )nucleus(icomp)%fit_D0 = .false.
   end if

!-----   Now, get experimental D1 and Gamma_gamma

   if(nucleus(icomp)%D1exp < 0.0d0)then
      open(unit=53,file=data_path(1:len_path)//'Level-Density-Data-RIPL3-L1.dat',   &
          status='old')
      D1exp = -1.0
      dD1exp = -1.0
      Gamma_g_1_exp =  -1.0
      dGamma_g_1_exp = -1.0
      header = .true.
      do while(header)
         read(53,'(a)')input_line
         if(input_line(1:1) /= '#')then
            backspace(53)
            header = .false.
         end if
      end do
      eof = 0
      do i = 1, max_nuc
         read(53,'(1x,i2,4x,i3,1x,f4.1,1x,f7.3,2(1x,e9.2),1x,2(1x,f5.2),1x,i5,1x,i5)',   &
                   iostat = eof)izz,iaa,xx,yy,x,x1,y,y1,ix,ix1
         if(eof == 1)exit
         if(izz == iz .and. iaa ==  ia-1)then
            D1exp = x
            dD1exp = x1
            Gamma_g_1_exp = real(ix)
            dGamma_g_1_exp = real(ix1)
            exit
         end if
      end do
      close(unit=53)
      D1exp = D1exp*1000.0d0
      dD1exp = dD1exp*1000.0d0
      nucleus(icomp)%D1exp = D1exp
      nucleus(icomp)%dD1exp = dD1exp
      nucleus(icomp)%Gamma_g_1_exp = Gamma_g_1_exp
      nucleus(icomp)%dGamma_g_1_exp = dGamma_g_1_exp
   end if
!
!-----   Next get GS deformation from Moller-Nix neutron separation energy
!-----   if larger than .25 call nucleus deformed
!-----   Check if they have already been set, if so then we can skip
!-----   First, check moller-nix file, then check cc_file for whatever is there.
   sum = 0.0d0
   do i = 2,6
      sum = abs(nucleus(icomp)%beta(i))
   end do

   if(sum <= 1.0d-5)then
      open(unit=53,file=data_path(1:len_path)//'mass-frdm95.dat',status='old')
      beta2 = 0.0d0
      beta3 = 0.0d0
      beta4 = 0.0d0
      beta6 = 0.0d0
      header = .true.
      do while(header)
         read(53,'(a)')input_line
         if(input_line(1:1) /= '#')then
            backspace(53)
            header = .false.
         end if
      end do
      do i = 1, max_nuc
         read(53,'(2(1x,i3),6x,f9.3,1x,f9.3,10x,4(1x,f7.3))',end=100)       &
               izz,iaa,x1,y1,b2,b3,b4,b6
         if(izz == iz .and. iaa == ia)then

            beta2 = b2
            beta3 = b3
            beta4 = b4
            beta6 = b6
            goto 100
         end if
      end do
 100  close(unit=53)
!----   Check the cc_file
      open(unit=53,file = data_path(1:len_path)//'Coupled-Channels.txt',status='old')
      found = .false.
      do while(.not. found)
         read(53,'(a)')input_line
         if(input_line(1:3) == 'END' .or. input_line(1:3) == 'End' &  
            .or. input_line(1:3) =='end')exit
         read(53,*)izz, iaa
         if(izz == iz .and. iaa == ia)then
            read(53,*)xxx2, xxx4, xxx6
            if(xxx2 > -1.0d0)then
               beta2 = xxx2
               beta4 = xxx4
               beta6 = xxx6
               beta3 = 0.0d0
            end if
            found = .true.
            exit
         else
            read(53,*)
            read(53,*)n1, n2
            do i = 1, n1 + n2
               read(53,*)
            end do
         end if
      end do
      close(unit=53)

      nucleus(icomp)%beta(2) = beta2
      nucleus(icomp)%beta(3) = beta3
      nucleus(icomp)%beta(4) = beta4
      nucleus(icomp)%beta(5) = 0.0d0
      nucleus(icomp)%beta(6) = beta6
   end if

   beta2 = nucleus(icomp)%beta(2)
   beta3 = nucleus(icomp)%beta(3)
   beta4 = nucleus(icomp)%beta(4)
   beta6 = nucleus(icomp)%beta(6)


!      write(6,*)(nucleus(icomp)%coll_enh(i), i = 1, 3)
!-------    Deformation parameter in Gilbert & Cameron
   xZ = real(iz,kind=8)
   xN = real(in,kind=8)
   xA = real(ia,kind=8)
   if(nucleus(icomp)%lev_option == 0)then
      nucleus(icomp)%level_model = 'G & C     '
      afac = 0.142
      if(beta2 > 0.2)afac = 0.120
!-------  Set up Gilbert & cameron parameters
      aparam = (0.00917d0*(sp(iz)+sn(in))+afac)*real(ia,kind=8)
      spin_cut = 0.0888d0
      delta = 0.0d0
      if(nucleus(icomp)%pair_model == 0)then
         if(iand(iz,1) == 0)delta = delta + deltap(iz/2)
         if(iand(in,1) == 0)delta = delta + deltan(in/2)
      end if

      if(nucleus(icomp)%pair_model == 1)then
         if(iand(iz,1) == 0 .and.  iand(in,1) == 0)delta = 24.0d0/sqrt(xA)
         if(iand(iz,1) == 0 .and. iand(in,1) == 1)delta = 12.0d0/sqrt(xA)
         if(iand(iz,1) == 1 .and. iand(in,1) == 0)delta = 12.0d0/sqrt(xA)
      end if

      if(nucleus(icomp)%pair_model == 2)delta = nucleus(icomp)%Delta_exp

      ematch = 2.5d0 + 150.0d0/real(ia,kind=8) + delta               !  at this point use the G&C starting value
 
      shell = 0.0d0
      gamma = 0.054d0
      ematch = 2.5d0 + 150.0d0/real(ia,kind=8) + delta
      maxlev = nucleus(icomp)%ncut
      ecut = nucleus(icomp)%level_ecut

!-----   Spin cut off parameter at Ecut
!-----   Taken from Reffo based on maximum liklihood arguments
      if(maxlev > 0)then
         sg2cut = 0.0d0
         sum = 0.0d0
         do i = 1, maxlev
            spin = nucleus(icomp)%state(i)%spin
            sg2cut = sg2cut + spin*(spin + 1.0d0)*(2.0d0*spin+1.0d0)
            sum = sum + (2.0d0*spin+1.0d0)
         end do
         sg2cut = sg2cut/(2.0d0*sum)      !   standard model,const. sig**2, integrate over J
!         sg2cut = sg2cut/(3.0d0*sum)     !   TALYS, and other arguments
      else
         sg2cut = (0.83d0*real(nucleus(icomp)%A,kind=8)**0.26d0)**2
      end if
      if(sg2cut < 1.0d-4)sg2cut = (0.83d0*real(nucleus(icomp)%A,kind=8)**0.26d0)**2
      if(ecut >= nucleus(icomp)%sep_e(1))nucleus(icomp)%fit_D0 = .false.

!-------   Deformation factors for spin-cutoff parameter for rotational enhancement
      nucleus(icomp)%sig2_perp = 0.0d0
      nucleus(icomp)%sig2_ax = 0.0d0
!----------------------------------------------------------------
      nucleus(icomp)%level_param(1) = aparam
      nucleus(icomp)%level_param(2) = spin_cut
      nucleus(icomp)%level_param(3) = delta
      nucleus(icomp)%level_param(4) = shell
      nucleus(icomp)%level_param(5) = gamma
      nucleus(icomp)%level_param(6) = ematch
      nucleus(icomp)%level_param(7) = ecut
      nucleus(icomp)%level_param(8) = sg2cut
      nucleus(icomp)%level_param(9) = 1.0d0
      nucleus(icomp)%level_param(10) = 0.0d0
      nucleus(icomp)%level_param(11) = 0.0d0
      nucleus(icomp)%level_param(12) = spin_cut*xA**(2.d0/3.d0)
      nucleus(icomp)%level_param(13) = 0.0d0
      nucleus(icomp)%level_param(16) = 0.0d0
      nucleus(icomp)%level_param(17) = 0.0d0
      nucleus(icomp)%level_param(18) = 0.0d0
      nucleus(icomp)%level_param(19) = real(nucleus(icomp)%lev_option,kind=8)   !  option is needed for later use

      nucleus(icomp)%vib_enh(1) = 0.0d0
      nucleus(icomp)%vib_enh(2) = 30.0d0
      nucleus(icomp)%vib_enh(3) = 5.0d0
      nucleus(icomp)%rot_enh(1) = 0.0d0
      nucleus(icomp)%rot_enh(2) = 30.0d0
      nucleus(icomp)%rot_enh(3) = 5.0d0
      call Find_T_E0(ia,nucleus(icomp)%level_param,                  &
                     nucleus(icomp)%vib_enh,                         &
                     nucleus(icomp)%rot_enh)
   end if
!------------------------------------------------------------------------
   if(nucleus(icomp)%lev_option == 1)then
      nucleus(icomp)%level_model = 'Talys-Def '
      aparam = 0.0692559d0*xA+0.282769d0*xA**(2.0d0/3.0d0)
      spin_cut = 0.01389d0
!  ----  first use of delta is for pairing in the LDM
      delta = 0.0d0
      if(iand(iz,1) == 0 .and. iand(in,1) == 0)delta = -11.0d0/sqrt(xA)
      if(iand(iz,1) == 1 .and. iand(in,1) == 1)delta = 11.0d0/sqrt(xA)
      gamma = 0.433090d0/real(ia,kind=8)**(1.0d0/3.0d0)
      shell = 0.0d0
! ------   Use difference between experimental Mass and Liquid-drop mass for shell
      Mass_LDM = -15.677d0*(1.0d0-1.79d0*((xN-xZ)/xA)**2)*xA +                 &
                  18.56d0*(1.0d0-1.79d0*((xN-xZ)/xA)**2)*xA**(2.0d0/3.0d0) +   &
                   0.717d0*xZ**2/xA**(1.0d0/3.0d0) -                           &
                   1.21129d0*xZ**2/xA +                                        &
                   delta +                                                     &
                   xZ*particle(2)%ME + xN*particle(1)%ME
      Mass_exp = nucleus(icomp)%ME
      shell = Mass_exp - Mass_LDM
!------    Now set delta for the Fermi-gas model
      delta = 0.0d0
      if(nucleus(icomp)%pair_model == 0)then
         if(iand(iz,1) == 0)delta = delta + deltap(iz/2)
         if(iand(in,1) == 0)delta = delta + deltan(in/2)
      end if

      if(nucleus(icomp)%pair_model == 1)then
         if(iand(iz,1) == 0 .and. iand(in,1) == 0)delta = 24.0d0/sqrt(xA)
         if(iand(iz,1) == 0 .and. iand(in,1) == 1)delta = 12.0d0/sqrt(xA)
         if(iand(iz,1) == 1 .and. iand(in,1) == 0)delta = 12.0d0/sqrt(xA)
      end if

      if(nucleus(icomp)%pair_model == 2)delta = nucleus(icomp)%Delta_exp

      ematch = 2.5d0 + 150.0d0/real(ia,kind=8) + delta               !  at this point use the G&C starting value
                                                      !  we fit it to the cummaltive level density
      maxlev = nucleus(icomp)%ncut
      ecut = nucleus(icomp)%level_ecut

!-------   Deformation factors for spin-cutoff parameter for rotational enhancement
      nucleus(icomp)%sig2_perp = 0.0d0
      nucleus(icomp)%sig2_ax = 0.0d0


      if(maxlev > 0)then
         sg2cut = 0.
         sum = 0.0
         do i = 1, maxlev
            spin = nucleus(icomp)%state(i)%spin
            sg2cut = sg2cut + spin*(spin + 1.0)*(2.0*spin+1.0)
            sum = sum + (2.0d0*spin+1.0d0)
         end do
         sg2cut = sg2cut/(3.0d0*sum)
      else
         sg2cut = (0.83d0*real(nucleus(icomp)%A,kind=8)**0.26d0)**2
      end if
      if(sg2cut < 1.0d-4)sg2cut =                                 &
         (0.83d0*real(nucleus(icomp)%A**(5.0d0/3.d0),kind=8)**0.26d0)**2

      if(ecut >= nucleus(icomp)%sep_e(1))nucleus(icomp)%fit_D0 = .false.
      nucleus(icomp)%level_param(1) = aparam
      nucleus(icomp)%level_param(2) = spin_cut
      nucleus(icomp)%level_param(3) = delta
      nucleus(icomp)%level_param(4) = shell
      nucleus(icomp)%level_param(5) = gamma
      nucleus(icomp)%level_param(6) = ematch
      nucleus(icomp)%level_param(7) = ecut
      nucleus(icomp)%level_param(8) = sg2cut
      nucleus(icomp)%level_param(9) = 1.0d0
      nucleus(icomp)%level_param(10) = 0.0d0
      nucleus(icomp)%level_param(11) = 0.0d0
      nucleus(icomp)%level_param(12) = spin_cut*xA**(5.0d0/3.d0)
      nucleus(icomp)%level_param(13) = 1.0d0
      nucleus(icomp)%level_param(16) = 0.0d0
      nucleus(icomp)%level_param(17) = 0.0d0
      nucleus(icomp)%level_param(18) = 0.0d0
      nucleus(icomp)%level_param(19) = real(nucleus(icomp)%lev_option,kind=8)   !  option is needed for later use

      nucleus(icomp)%vib_enh(1) = 0.0d0
      nucleus(icomp)%vib_enh(2) = 30.0d0
      nucleus(icomp)%vib_enh(3) = 5.0d0
      nucleus(icomp)%rot_enh(1) = 0.0d0
      nucleus(icomp)%rot_enh(2) = 30.0d0
      nucleus(icomp)%rot_enh(3) = 5.0d0
      nucleus(icomp)%rot_enh(4) = nucleus(icomp)%sig2_perp
      nucleus(icomp)%rot_enh(5) = nucleus(icomp)%sig2_ax
      call Find_T_E0(ia,nucleus(icomp)%level_param,                     &
                     nucleus(icomp)%vib_enh,                            &
                     nucleus(icomp)%rot_enh)
   end if


   if(nucleus(icomp)%lev_option >= 2 .and. nucleus(icomp)%lev_option <= 5)then
      nucleus(icomp)%level_model = 'Collective'
      aparam = 0.0207305d0*xA + 0.229537d0*xA**(2.0d0/3.0d0)
      spin_cut = 0.01389d0
!  ----  first use of delta is for pairing in the LDM
      delta = 0.0d0
      if(iand(iz,1) == 0 .and. iand(in,1) == 0)delta = -11.0d0/sqrt(xA)
      if(iand(iz,1) == 1 .and. iand(in,1) == 1)delta = 11.0d0/sqrt(xA)
      gamma = 0.473625d0/real(ia,kind=8)**(1.0d0/3.0d0)
      shell = 0.0d0
! ------   Use difference between experimental Mass and Liquid-drop mass for shell
      Mass_LDM = -15.677d0*(1.0d0-1.79d0*((xN-xZ)/xA)**2)*xA +               &
                  18.56d0*(1.0d0-1.79d0*((xN-xZ)/xA)**2)*xA**(2.0d0/3.0d0) + &
                  0.717d0*xZ**2/xA**(1.0d0/3.0d0) -                          &
                  1.21129d0*xZ**2/xA +                                       &
                  delta +                                                    &
                  xZ*particle(2)%ME + xN*particle(1)%ME
      Mass_exp = nucleus(icomp)%ME
      shell = Mass_exp - Mass_LDM
!------    Now set delta for the Fermi-gas model
      delta = 0.0d0
      if(nucleus(icomp)%pair_model == 0)then
         if(iand(iz,1) == 0)delta = delta + deltap(iz/2)
         if(iand(in,1) == 0)delta = delta + deltan(in/2)
      end if

      if(nucleus(icomp)%pair_model == 1)then
         if(iand(iz,1) == 0 .and. iand(in,1) == 0)delta = 24.0d0/sqrt(xA)
         if(iand(iz,1) == 0 .and. iand(in,1) == 1)delta = 12.0d0/sqrt(xA)
         if(iand(iz,1) == 1 .and. iand(in,1) == 0)delta = 12.0d0/sqrt(xA)
      end if

      if(nucleus(icomp)%pair_model == 2)delta = nucleus(icomp)%Delta_exp

      ematch = 2.5d0 + 150.0d0/real(ia,kind=8) + delta               !  at this point use the G&C starting value
                                                                     !  we fit it to the cummaltive level density
      maxlev = nucleus(icomp)%ncut
      ecut = nucleus(icomp)%level_ecut

!-------   Deformation factors for spin-cutoff parameter for rotational enhancement
      if(nucleus(icomp)%lev_option == 2 .or. nucleus(icomp)%lev_option == 3)then
         nucleus(icomp)%sig2_perp = (1.0d0+beta2/3.0d0)
         nucleus(icomp)%sig2_ax = (1.0d0-2.0d0*beta2/3.0d0)
      elseif(nucleus(icomp)%lev_option == 4 .or. nucleus(icomp)%lev_option == 5)then 
         nucleus(icomp)%sig2_perp = (1.0d0+sqrt(5.0d0/(4.0d0*pi))*beta2)
         nucleus(icomp)%sig2_ax = (1.0d0-0.5d0*sqrt(5.0d0/(4.0d0*pi))*beta2)
      end if

      if(maxlev > 0)then
         sg2cut = 0.0d0
         sum = 0.0d0
         do i = 1, maxlev
            spin = nucleus(icomp)%state(i)%spin
            sg2cut = sg2cut + spin*(spin + 1.0d0)*(2.0d0*spin+1.0d0)
            sum = sum + (2.0d0*spin+1.0d0)
         end do
         sg2cut = sg2cut/(3.0d0*sum)
      else
        sg2cut = (0.83d0*real(nucleus(icomp)%A,kind=8)**0.26d0)**2
      end if
      if(sg2cut < 1.0d-4)sg2cut = (0.83*real(nucleus(icomp)%A,kind=8)**0.26)**2
      if(ecut >= nucleus(icomp)%sep_e(1))nucleus(icomp)%fit_D0 = .false.
      nucleus(icomp)%level_param(1) = aparam
      nucleus(icomp)%level_param(2) = spin_cut
      nucleus(icomp)%level_param(3) = delta
      nucleus(icomp)%level_param(4) = shell
      nucleus(icomp)%level_param(5) = gamma
      nucleus(icomp)%level_param(6) = ematch
      nucleus(icomp)%level_param(7) = ecut
      nucleus(icomp)%level_param(8) = sg2cut
      nucleus(icomp)%level_param(9) = 1.0d0
      nucleus(icomp)%level_param(10) = 1.0d0
      nucleus(icomp)%level_param(11) = 2.0d0
      nucleus(icomp)%level_param(12) = spin_cut*xA**(5.0d0/3.d0)
      nucleus(icomp)%level_param(13) = 1.0d0
      nucleus(icomp)%level_param(16) = 0.0d0
      nucleus(icomp)%level_param(17) = 0.0d0
      nucleus(icomp)%level_param(18) = 0.0d0

      nucleus(icomp)%level_param(19) = real(nucleus(icomp)%lev_option,kind=8)   !  option is needed for later use

      nucleus(icomp)%vib_enh(1) = 1.0d0
      nucleus(icomp)%vib_enh(2) = 30.0d0
      nucleus(icomp)%vib_enh(3) = 5.0d0
      nucleus(icomp)%rot_enh(1) = 1.0d0
      nucleus(icomp)%rot_enh(2) = 30.0d0
      nucleus(icomp)%rot_enh(3) = 5.0d0
      nucleus(icomp)%rot_enh(4) = nucleus(icomp)%sig2_perp
      nucleus(icomp)%rot_enh(5) = nucleus(icomp)%sig2_ax
      call Find_T_E0(ia,nucleus(icomp)%level_param,nucleus(icomp)%vib_enh,   &
                     nucleus(icomp)%rot_enh)

   end if

   spin_cut = nucleus(icomp)%level_param(2)

   return
end subroutine get_lev_den
!
!*******************************************************************************
!
subroutine finish_lev_den(icomp)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine finishes up the level densities. It is called after
!    all options are read in so that it updates all level density information,
!    refits D0 if needed, and fits Ematch if needed
!
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!        options
!        nuclei
!        constants
!        particles_def
!
!     Subroutines:
!
!        find_prob
!        unpack_data
!
!     External functions:
!
!        real(kind=8) :: random_64
!        real(kind=8) :: random_32
!
!    MPI routines:
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
   use nodeinfo
   use options
   use nuclei
   use constants
   use particles_def
   implicit none
!-----------   Input Data   --------------------------------------------
   integer(kind=4), intent(in) :: icomp
!-----------   Internal Data -------------------------------------------
   real(kind=8) :: e,s1,s2
   real(kind=8) :: xI
   integer(kind=4) :: ipar
   real(kind=8) :: U
   real(kind=8) :: D0,D0exp,dD0exp, D1
   real(kind=8) :: aparam,a_Sn
   real(kind=8) :: spin_cut,delta,shell,gamma
   real(kind=8) :: xA
   integer(kind=4) :: ia
!-----------   External functions   ------------------------------------
   real(kind=8) :: lev_space
   real(kind=8) :: aparam_u
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Calcuate D0 with these level density parameters              +
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   nucleus(icomp)%sig2_perp = 0.0d0
   nucleus(icomp)%sig2_ax = 0.0d0

   if(nucleus(icomp)%lev_option == 2 .or. nucleus(icomp)%lev_option == 3)then
      nucleus(icomp)%sig2_perp = (1.0d0+nucleus(icomp)%beta(2)/3.0d0)
      nucleus(icomp)%sig2_ax = (1.0d0-2.0d0*nucleus(icomp)%beta(2)/3.0d0)
   elseif(nucleus(icomp)%lev_option == 4 .or. nucleus(icomp)%lev_option == 5)then
      nucleus(icomp)%sig2_perp = (1.0d0+sqrt(5.0d0/(4.0d0*pi))*nucleus(icomp)%beta(2))
      nucleus(icomp)%sig2_perp = (1.0d0-0.5d0*sqrt(5.0d0/(4.0d0*pi))*nucleus(icomp)%beta(2))
   end if

   nucleus(icomp)%rot_enh(4) = nucleus(icomp)%sig2_perp
   nucleus(icomp)%rot_enh(5) = nucleus(icomp)%sig2_ax

!---------------------------------------------------------------------
   e = nucleus(icomp)%sep_e(1)          !   neutron separation energy
   s1 = nucleus(icomp)%s1
   s2 = nucleus(icomp)%s2
   ipar = nucleus(icomp)%ipar
   U = e - nucleus(icomp)%level_param(3)
   D0exp = nucleus(icomp)%D0exp
   dD0exp = nucleus(icomp)%dD0exp
   if(D0exp < 1.0d-6)nucleus(icomp)%fit_D0 = .false.
   xA = nucleus(icomp)%A
   ia = int(xA)
   xI = nucleus(icomp)%target_spin
   ipar = nucleus(icomp)%target_ipar

!-----   Make sure that the spin cutoff parameter gets updated
   if(nucleus(icomp)%lev_option == 0)then
      nucleus(icomp)%level_param(12) = nucleus(icomp)%level_param(2)*xA**(2.d0/3.d0)
   elseif(nucleus(icomp)%lev_option == 1)then
      nucleus(icomp)%level_param(12) = nucleus(icomp)%level_param(2)*xA**(5.d0/3.d0)
   elseif(nucleus(icomp)%lev_option == 2)then
      nucleus(icomp)%level_param(12) = nucleus(icomp)%level_param(2)*xA**(5.d0/3.d0)
   end if



   if(nucleus(icomp)%level_param(7) < nucleus(icomp)%sep_e(1))then


      nucleus(icomp)%D0 = lev_space(e,xI,ipar,0,                   &
                                    nucleus(icomp)%level_param,    &
                                    nucleus(icomp)%vib_enh,        &
                                    nucleus(icomp)%rot_enh,ia)
      D0 = nucleus(icomp)%D0
      if(nucleus(icomp)%fit_D0)then
         if(D0 < D0exp-dD0exp/2 .or. D0 > D0exp+dD0exp/2)then     !  D0 value is outside expt. range
            call fit_D0(D0exp,e,xI,ipar,                    &
                        nucleus(icomp)%level_param,         &
                        nucleus(icomp)%fit_aparam,          &
                        nucleus(icomp)%vib_enh,             &
                        nucleus(icomp)%rot_enh,ia)
         end if
      end if
      D0 = lev_space(e,xI,ipar,0,nucleus(icomp)%level_param,                &
                     nucleus(icomp)%vib_enh,nucleus(icomp)%rot_enh,ia)    !   if ematch is abvoe sep-e, then check D0
      nucleus(icomp)%D0 = D0
      D1 = lev_space(e,xI,ipar,1,nucleus(icomp)%level_param,                &
                   nucleus(icomp)%vib_enh,nucleus(icomp)%rot_enh,ia)    !   if ematch is abvoe sep-e, then check D0
      nucleus(icomp)%D1 = D1
   else
      nucleus(icomp)%D0 = -1.0d0
      nucleus(icomp)%D1 = -1.0d0
   end if

   if(nucleus(icomp)%fit_ematch)then
      call fit_lev_den(icomp)

      call Find_T_E0(ia,nucleus(icomp)%level_param,                   &
                     nucleus(icomp)%vib_enh,                          &
                     nucleus(icomp)%rot_enh)
   end if


   aparam = nucleus(icomp)%level_param(1)
   spin_cut = nucleus(icomp)%level_param(2)
   delta = nucleus(icomp)%level_param(3)
   shell = nucleus(icomp)%level_param(4)
   gamma = nucleus(icomp)%level_param(5)
   a_Sn = aparam_u(u,aparam,shell,gamma)
   nucleus(icomp)%a_Sn = a_Sn

   if(nint(nucleus(icomp)%level_param(9)) == 0)then
      nucleus(icomp)%sig2 = spin_cut*xA**(2.0d0/3.0d0)
      nucleus(icomp)%sig2_Sn = nucleus(icomp)%sig2*sqrt(a_Sn*U)
   else
      nucleus(icomp)%sig2 = spin_cut*xA**(5.0d0/3.0d0)
      nucleus(icomp)%sig2_Sn = spin_cut*xA**(5.0d0/3.0d0)*          &
                               sqrt(a_Sn*U)/aparam
   end if

   if(print_me)write(6,*)'finished level_param for nucleus',icomp
   return
end subroutine finish_lev_den
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine fit_lev_den(icomp)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine fits the matching energy to the experimental cummlative
!    level density
!
!   Dependencies:
!
!     Modules:
!
!        nuclei
!        lev_fit
!        nodeinfo
!
!     Subroutines:
!
!        cum_rho
!        Find_T_E0
!
!     External functions:
!
!        lev_space
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
   use nuclei
   use lev_fit
   use nodeinfo
   implicit none
   integer(kind=4), intent(in) :: icomp
!----------------------------------------------------------------------
!---------   Parameters to fit cumulative level density
   real(kind=8), allocatable :: elv(:)
   real(kind=8), allocatable :: cum_rho(:),dcum_rho(:)
   real(kind=8), allocatable :: cum_fit(:)
   real(kind=8) :: chisq
   integer(kind=4) :: nca, nfit
   real(kind=8) :: D0,D0old
   real(kind=8) :: D0_exp,dD0_exp
!-----------------------------
   integer(kind=4) :: i,k
   integer(kind=4) :: ia
   real(kind=8) :: chi_min,ematch,emin,emax
   real(kind=8) :: xI
   integer(kind=4) :: ipar
   integer(kind=4) :: numf
   real(kind=8) :: ratio
!-------------  External functions
   real(kind=8) :: lev_space

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------   Set up fit to reproduce cumulative level density up to Ecut +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   D0_exp = nucleus(icomp)%D0exp
   dD0_exp = nucleus(icomp)%dD0exp
   nfit = nucleus(icomp)%ncut
   if(nfit <= 5)return                            !  Their are so few discretes, no point in fitting to the cumulative level density
   if(.not.allocated(elv))allocate(elv(nfit))     !  allocate cumulative density array
   if(.not.allocated(cum_rho))allocate(cum_rho(nfit))     !  allocate cumulative density array
   if(.not.allocated(dcum_rho))allocate(dcum_rho(nfit))
   if(.not.allocated(cum_fit))allocate(cum_fit(nfit))     !  allocate cumulative density array
   do i = 1,nfit                                    !  Make cumulative rho array
      elv(i) = nucleus(icomp)%state(i)%energy
      cum_rho(i) = real(i,kind=8)
      dcum_rho(i) = 1.0d0/dsqrt(cum_rho(i))           !  assume error 1/sqrt
   end do
   D0old = nucleus(icomp)%D0                         !  original D0
   sep_e_l = nucleus(icomp)%sep_e(1)                 !  variables stored in module and used by subroutine func1
   s1_lev = nucleus(icomp)%s1
   s2_lev = nucleus(icomp)%s2
   ipar_lev = nucleus(icomp)%ipar
   AA_lev = nucleus(icomp)%A
   ia = AA_lev
   nca = 9                                          !  dimension of covariance matrix
   xI = nucleus(icomp)%target_spin
   ipar = nucleus(icomp)%target_ipar
   D0 = lev_space(sep_e_l,xI,ipar,0,                  &
                  nucleus(icomp)%level_param,         &
                  nucleus(icomp)%vib_enh,             &
                  nucleus(icomp)%rot_enh,AA_lev)    !   if ematch is above sep-e, then check D0

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------    set up params to fit - depending on number of discrete    +
!-------    if # discrete states is >= 10 fit both Delta and Ematch   +
!-------    if # discrete states is < 10 fit only Ematch              +
!-------    Note if we fit Delta, this will likely change D0          +
!-------    Hence, after we fit to the cumulative level density,      +
!-------    we redetermine a [dc(1)] to get old D0                    +
!-------    Then refit Delta and Ematch.                              +
!-------    Rinse and repeat as many times as needed                  +
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   call cumm_rho(nfit,elv,ia,nucleus(icomp)%level_param,      &
                 nucleus(icomp)%vib_enh,                      &
                 nucleus(icomp)%rot_enh,cum_fit)

   chisq = 0.0d0

   do i = 2, nfit
      chisq = chisq+(cum_rho(i)-cum_fit(i))**2/dcum_rho(i)**2
   end do
   chisq = chisq/real(nfit-2,kind=8)
   if(.not.nucleus(icomp)%fit_ematch)goto 119                       !  they have changed, so keep trying to fit
!------   Check to see if fit makes sense, and if known spectrum does. It is 
!------   possible that the known spectrum is incomplete, and we should actually
!------   use many fewer states. Check cumulative fit and if fitted value exceeds
!------   actual by 50% then reduce the number of states we keep.
!      ematch=dc(6)
   ematch = nucleus(icomp)%level_param(6)
   chi_min = 1.0d10
!
!   Min and max range of E_match. Must be greater than Delta = level_param(3)
!   and E_cut = level_param(7)
!   Search for E_match in 10 keV steps

   emin  = max(nucleus(icomp)%level_param(7),                   &
               nucleus(icomp)%level_param(3)) + 0.2
   emax = nucleus(icomp)%sep_e(1) - 0.2d0          !   neutron separation energy
   numf = int((emax-emin)/0.01d0) + 1

!************************************************************************
!                                                                       *
!   Brute force loop over possible values of E_match. Find value with   *
!   best match to cummlative level density for known discrete states    *
!                                                                       *
!************************************************************************
   do k = 1, numf
      nucleus(icomp)%level_param(6) = emin + real(k-1,kind=8)*0.01d0
      call Find_T_E0(ia,nucleus(icomp)%level_param,             &
                     nucleus(icomp)%vib_enh,                    &
                     nucleus(icomp)%rot_enh)
      call cumm_rho(nfit,elv,ia,nucleus(icomp)%level_param,     &
                    nucleus(icomp)%vib_enh,                     &
                    nucleus(icomp)%rot_enh,cum_fit)
      chisq = 0.0d0
      do i = 1,nfit
         chisq = chisq + (cum_rho(i)-cum_fit(i))**2/dcum_rho(i)**2
      end do
      chisq = chisq/real(nfit-2,kind=8)
      if(chisq < chi_min)then
         ematch = nucleus(icomp)%level_param(6)
         chi_min = chisq
      end if
   end do
   nucleus(icomp)%level_param(6) = ematch
   call Find_T_E0(ia,nucleus(icomp)%level_param,                &
                  nucleus(icomp)%vib_enh,                       &
                  nucleus(icomp)%rot_enh)
   call cumm_rho(nfit,elv,ia,nucleus(icomp)%level_param,        &
                 nucleus(icomp)%vib_enh,                        &
                 nucleus(icomp)%rot_enh,cum_fit)
   ratio = cum_fit(nfit)/cum_rho(nfit)

   nucleus(icomp)%cum_rho_ratio = ratio

   if(ratio > 1.25d0 .or. ratio < 0.75d0)nucleus(icomp)%ematch_warn = .true.


 119  if(allocated(elv))deallocate(elv)
   if(allocated(cum_rho))deallocate(cum_rho)
   if(allocated(dcum_rho))deallocate(dcum_rho)
   if(allocated(cum_fit))deallocate(cum_fit)
   return
end subroutine fit_lev_den
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine rho_J_par_e(E, xJ, ipar, level_param, vib_enhance, rot_enhance, ia, &
                       rho, apu, sig2, K_vib, K_rot)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine returns level density rho(E)
!    and spin-cutoff parameter,sig at excitation energy E
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        rho_FT
!        rho_BFM
!
!     External functions:
!
!        real(kind=8) :: aparam_u
!        real(kind=8) :: sig2_param
!        real(kind=8) :: enhance_rot
!        real(kind=8) :: enhance_vib
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
   implicit none
   real(kind=8), intent(in) :: e
   real(kind=8), intent(in) :: xJ
   integer(kind=4), intent(in) :: ipar
   real(kind=8), intent(in) :: level_param(20)
   real(kind=8), intent(in) :: vib_enhance(3), rot_enhance(5)
   integer(kind=4), intent(in) :: ia
!-----
   real(kind=8), intent(out) :: rho
   real(kind=8), intent(out) :: apu
   real(kind=8), intent(out) :: sig2
   real(kind=8), intent(out) :: K_vib, K_rot 
!-----
   real(kind=8) :: sig2_perp
   real(kind=8) :: aparam, spin_cut, delta, shell, gamma
   real(kind=8) :: ematch, ecut, sg2cut
   real(kind=8) :: A
   real(kind=8) :: Em, U, E0
   real(kind=8) :: T
   real(kind=8) :: de
   real(kind=8) :: low_e_mod
   real(kind=8) :: enhance
   integer(kind=4) :: vib_mode, rot_mode
   real(kind=8) :: pmode, pe1, pbb
   real(kind=8) :: jfac, pfac
!-----    External Functions-----------------------------------------------------
   real(kind=8) :: aparam_u
   real(kind=8) :: sig2_param
   real(kind=8) :: enhance_rot
   real(kind=8) :: enhance_vib
   real(kind=8) :: spin_fac
   real(kind=8) :: parity_fac  
!--------------------------------------------------------------------------------
!-----------   Atomic mass -> real value
   A = real(ia,kind=8)
!-----------   level-density parameters --------------------

   aparam = level_param(1)
   spin_cut = level_param(2)
   delta = level_param(3)
   shell = level_param(4)
   gamma = level_param(5)
   ematch = level_param(6)
   ecut = level_param(7)
   sg2cut = level_param(8)
   low_e_mod = level_param(9)

   pmode = level_param(16)
   pe1 = level_param(17)
   pbb = level_param(18)

   rot_mode = nint(level_param(10))
   vib_mode = nint(level_param(11))

!------------    start of calculation   --------------------
   de = 0.01d0
   Em = ematch
   if(E < ematch)then                     !  finite temperature
      T = level_param(14)
      E0 = level_param(15)              !  constant temperature

      sig2 = sig2_param(E,level_param,A)

      if(nint(level_param(19)) == 3 .or. nint(level_param(19)) == 5)then
         sig2 = sig2*rot_enhance(4)
      end if
 
      U = E - delta
      apu = aparam_u(U,aparam,shell,gamma)
      K_rot = 1.0d0
      K_vib = 1.0d0

      call rho_FT(E, E0, T, rho)

      enhance = K_rot*K_vib
      rho = rho*enhance
   else                                   !  fermi Gas
      U = E - delta
      sig2 = sig2_param(E,level_param,A)
      sig2_perp = rot_enhance(4)*sig2
      apu = aparam_u(U,aparam,shell,gamma)

      K_rot = enhance_rot(rot_mode, sig2, sig2_perp, rot_enhance, E)
      K_vib = enhance_vib(vib_mode, A, U, E, apu, Shell, vib_enhance)

      if(nint(level_param(19)) == 3 .or. nint(level_param(19)) == 5)then
         sig2 = sig2*rot_enhance(4)
      end if

      call rho_BFM(U, apu, sig2, rho)

      enhance = K_rot*K_vib
      rho = rho*enhance
   end if

   jfac = spin_fac(xJ, sig2)
   pfac = parity_fac(E, xJ, ipar, pmode, pe1, pbb)

   rho = rho*jfac*pfac
   

   return
end subroutine rho_J_par_E
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine rhoe(E, level_param, vib_enhance, rot_enhance, ia, rho,     &
                apu, sig2, K_vib, K_rot)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine returns level density rho(E)
!    and spin-cutoff parameter,sig at excitation energy E
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        rho_FT
!        rho_BFM
!
!     External functions:
!
!        real(kind=8) :: aparam_u
!        real(kind=8) :: sig2_param
!        real(kind=8) :: enhance_rot
!        real(kind=8) :: enhance_vib
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
   implicit none
   real(kind=8), intent(in) :: e
   real(kind=8), intent(in) :: level_param(20)
   real(kind=8), intent(in) :: vib_enhance(3), rot_enhance(5)
   integer(kind=4), intent(in) :: ia
!-----
   real(kind=8), intent(out) :: rho
   real(kind=8), intent(out) :: apu
   real(kind=8), intent(out) :: sig2
   real(kind=8), intent(out) :: K_vib, K_rot 
!-----
   real(kind=8) :: sig2_perp
   real(kind=8) :: aparam, spin_cut, delta, shell, gamma
   real(kind=8) :: ematch, ecut, sg2cut
   real(kind=8) :: A
   real(kind=8) :: Em, U, E0
   real(kind=8) :: T
   real(kind=8) :: de
   real(kind=8) :: low_e_mod
   real(kind=8) :: enhance
   integer(kind=4) :: vib_mode, rot_mode
!-----    External Functions-----------------------------------------------------
   real(kind=8) :: aparam_u
   real(kind=8) :: sig2_param
   real(kind=8) :: enhance_rot
   real(kind=8) :: enhance_vib
!--------------------------------------------------------------------------------
!-----------   Atomic mass -> real value
   A = real(ia,kind=8)
!-----------   level-density parameters --------------------
   aparam = level_param(1)
   spin_cut = level_param(2)
   delta = level_param(3)
   shell = level_param(4)
   gamma = level_param(5)
   ematch = level_param(6)
   ecut = level_param(7)
   sg2cut = level_param(8)
   low_e_mod = level_param(9)

   rot_mode = nint(level_param(10))
   vib_mode = nint(level_param(11))

!------------    start of calculation   --------------------
   de = 0.01d0
   Em = ematch
   if(E < ematch)then                     !  finite temperature
      T = level_param(14)
      E0 = level_param(15)              !  constant temperature

      sig2 = sig2_param(E,level_param,A)
 
      U = E - delta
      apu = aparam_u(U,aparam,shell,gamma)
      K_rot = 1.0d0
      K_vib = 1.0d0

      if(nint(level_param(19)) == 3 .or. nint(level_param(19)) == 5)then
         sig2 = sig2*rot_enhance(4)
      end if

      call rho_FT(E, E0, T, rho)

      enhance = K_rot*K_vib
      rho = rho*enhance

   else                                   !  fermi Gas
      U = E - delta
      sig2 = sig2_param(E,level_param,A)
      sig2_perp = rot_enhance(4)*sig2
      apu = aparam_u(U,aparam,shell,gamma)

      K_rot = enhance_rot(rot_mode, sig2, sig2_perp, rot_enhance, E)
      K_vib = enhance_vib(vib_mode, A, U, E, apu, Shell, vib_enhance)

      if(nint(level_param(19)) == 3 .or. nint(level_param(19)) == 5)then
         sig2 = sig2*rot_enhance(4)
      end if

      call rho_BFM(U, apu, sig2, rho)

      enhance = K_rot*K_vib
      rho = rho*enhance

   end if

   return
end subroutine rhoe
!
!--------------------------------------------------------------------
!
real(kind=8) function sig2_param(E,level_param,A)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns the spin cutoff parameter
!
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!        constants
!
!     Subroutines:
!
!        exit_YAHFC
!
!     External functions:
!
!        aparam_u
!
!     MPI routines:
!
!        MPI_Abort   ----    via exit_YAHFC
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
   use nodeinfo
   use constants
   implicit none
!------  Passed in
   real(kind=8), intent(in) :: E, level_param(20), A
!------  Passed out
!------  Internal
   real(kind=8) :: aparam,apu,delta,shell,gamma,spin_cut,sg2cut
   real(kind=8) :: ecut,ematch,low_e_mod
   real(kind=8) :: U, Um, sig2, sig2_em, deriv, sig, sig2_min
   integer(kind=4) :: i, sig_model
!------  External functions  ---------------------------------------------------
   real(kind=8) :: aparam_u
!------  Start Calculation -----------------------------------------------------
   aparam = level_param(1)
   spin_cut = level_param(2)
   delta = level_param(3)
   shell = level_param(4)
   gamma = level_param(5)
   ematch = level_param(6)
   ecut = level_param(7)
   sg2cut = level_param(8)
   low_e_mod = level_param(9)
   U = E - delta
   Um = Ematch - delta
   sig = level_param(12)
   sig_model = nint(level_param(13))
   sig2_min = min((0.83d0*A**0.26d0)**2.0d0,sg2cut)

!   if( E < ecut)then
!      sig2_param = sg2cut
!      return
!   end if

   if(E < ematch )then
      apu = aparam_u(Um,aparam,shell,gamma)
      sig2_em = sig*sqrt(max(0.2d0,Um*apu))/aparam
      if(sig_model == 0)sig2_em = sig*sqrt(max(0.2d0,Um*apu))/aparam
      if(sig_model == 1)sig2_em = sig*sqrt(max(0.2d0,Um)/apu)
      sig2_em = max(sig2_em,sig2_min)
      deriv = (sig2_em - sg2cut)/(ematch - ecut)
      sig2 = sig2_em - deriv*(ematch - e)
      if(sig2 < sig2_min) sig2 = sig2_min
      if(sig2 < 0.0)then
         write(6,*)'A =',A, E
         do i = 1, 9
            write(6,*)i,level_param(i)
         end do
         write(6,*)'sig_mod',sig_model
         write(6,*)spin_cut, apu, aparam
         write(6,'(''sig2'',7(1x,f10.5))')sg2cut,sig2_em,deriv,    &
            e,ecut,ematch,sig2
         if(iproc == 0)write(6,*)'sig2 < 0'
         call exit_YAHFC(401)
      end if
   else
      apu = aparam_u(U,aparam,shell,gamma)
      sig2 = (0.83*A**0.26)**2
      if(U > 0.0d0)then
         if(sig_model == 0)sig2 = sig*sqrt(apu*U)/aparam
         if(sig_model == 1)sig2 = sig*sqrt(U/apu)
         sig2 = max(sig2,(0.83*A**0.26)**2)
      end if
   end if
   sig2_param = sig2

   return
end function sig2_param
!
!--------------------------------------------------------------------
!
subroutine rho_FT(e,E0,T,rho)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine returns the level density using the 
!    finite-temperature model
!
!   Dependencies:
!
!     Modules:
!
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
   use constants
   implicit none
   real(kind=8), intent(in) :: e
   real(kind=8), intent(in) :: E0, T
   real(kind=8), intent(out) :: rho
!-------------------------------------------------------------------
   rho = exp((e-E0)/T)/T
   return
end subroutine rho_FT

!
!--------------------------------------------------------------------
!
subroutine rho_BFM(U, apu, sig2, rho)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine returns the level density using the 
!    back-shifted fermi gas model
!
!   Dependencies:
!
!     Modules:
!
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
   use constants
   implicit none
   real(kind=8), intent(in) :: U
   real(kind=8), intent(in) :: apu
   real(kind=8), intent(in) :: sig2
   real(kind=8), intent(out) :: rho
!-------------------------------------------------------------------
   real(kind=8) :: rho_F, rho_0
   real(kind=8) :: T,an,ap
   real(kind=8) :: U1,exponent1,exponent2
!-------------------------------------------------------------------
   exponent1 = 2.0d0*sqrt(apu*U)
   U1 = U
!   if(U <= 0.0d0)U1 = 1.0d-6
   if(U <= 0.0d-6)U1 = 1.0d-6
   rho_F = exp(exponent1)/                                            &
         (sqrt(2.0d0*sig2)*12.0d0*apu**(0.25)*U1**(1.25))
   T = sqrt(U/apu)
   an = apu/2.0d0
   ap = apu/2.0d0
   exponent2 = 4.0d0*ap*an*T**2
   rho_0 = exp(exponent1+1.0d0)*(an+ap)**2/(24.0d0*sqrt(sig2)*sqrt(an*ap))
   rho = rho_F*rho_0/(rho_F+rho_0)

!   write(40,*)U,rho_f,rho_0,1.0d0/(rho_f+

!   if(nint(low_e_mod) == 0)then
!      rho = rho_F
!   elseif(nint(low_e_mod) == 1)then
!      rho = 1.0d0/(rho_F + rho_0)
!   end if
   return
end subroutine rho_BFM
!
!*******************************************************************************
!
real(kind=8) function spin_fac(xj,sg2)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns the spin-dependence factor for the level density
!
!   Dependencies:
!
!     Modules:
!
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
   use constants
   implicit none
   real(kind=8), intent(in) :: xj,sg2
!------------------------------------------------------------------------
   spin_fac = (d_two*xj+d_one)*exp(-(xJ+d_half)**2/(d_two*sg2))/   &
            (d_two*sg2)
   return
end function spin_fac
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function parity_fac(E,xJ,ipar,mode,e1,bb)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns the parity-dependence factor for the level density
!
!   Dependencies:
!
!     Modules:
!
!        None
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
   implicit none
   real(kind=8), intent(in) :: E, xj, mode, e1, bb
   integer(kind=4) :: ipar
   integer(kind=4) :: ippar
   real(kind=8) :: xj_fac
!-------------------------------------------------------------------------------
   xj_fac = 1.0d0
   if(xj < 0.0d0)xj_fac = 0.0d0
   if(nint(mode) == 0)then
      parity_fac=0.5d0
      return
   else
      ippar = nint(mode) - 1
      if(ippar == ipar)then
         parity_fac = 0.0d0
         if(E >= e1)then
            parity_fac = 0.5*tanh(bb*(E - e1))
         end if
         return
      else
         parity_fac = 1.0d0
         if(E >= e1)then
            parity_fac = 1.0d0 - 0.5*tanh(bb*(E - e1))
         end if
         return
     end if
   end if
   return
end function parity_fac
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function lev_space(e, xj, ipar, l, level_param,       &
                                vib_enh, rot_enh, ia)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns the level spacing
!    at the separation energy for neutrons with orbital angular momentum l 
!    on A-1 target nucleus
!
!   Dependencies:
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        rhoe
!
!     External functions:
!
!        real(kind=8) :: spin_fac
!        real(kind=8) :: parity_fac
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
   implicit none
   real(kind=8), intent(in) :: e, xj
   integer(kind=4), intent(in) :: l
   integer(kind=4), intent(in) :: ipar
   real(kind=8), intent(in) :: level_param(20)
   real(kind=8), intent(in) :: vib_enh(3), rot_enh(5)
   integer(kind=4), intent(in) :: ia
   real(kind=8) :: mode, e1, bb
! -------------------------------------------------------
   integer(kind=4) :: i, k
   real(kind=8) :: rrho,rho0,apu,sg2
   real(kind=8) :: xl
   integer(kind=4) :: iparn
   integer(kind=4) :: j
   integer(kind=4) :: iJ, il
   real(kind=8) :: K_vib, K_rot, xI
!------------    External functions     -----------------------------
   real(kind=8) :: spin_fac
   real(kind=8) :: parity_fac
!--------------------------------------------------------------------  
   mode = level_param(16)
   e1 = level_param(17)
   bb = level_param(18)
          
   call rhoe(e,level_param,vib_enh,rot_enh,ia,rrho,apu,sg2,k_vib,K_rot)
!---------   Factor of two is to account for parity (only one here)
   xl = real(l, kind=8)
   iparn = ((2*ipar - 1)*(-1)**l+1)/2


   rho0 = 0.0d0
   iJ = nint(2.0d0*xj)
   il = 2*l
   do i = iabs(iJ - il), iJ + il, 2
      do j = -1, 1, 2
         k = i + j
         if(k < 0)cycle
         xI = real(k,kind=8)/2.0d0
         rho0 = rho0 + rrho*spin_fac(xI,sg2)*                    &
                       parity_fac(e,xI,iparn,mode,e1,bb)
     end do     
   end do

!---------  invert rho multiply by 1000 for spacing in keV per level
   lev_space = 1.0d6/rho0
   return
end function lev_space
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function aparam_u(u,aparam,shell,gamma)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns the excitation-energy dependent a-parameter
!
!   Dependencies:
!
!     Modules:
!
!        None
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
   implicit none
   real(kind=8), intent(in) :: u,aparam,shell,gamma
!------------------------------------------------------------------
   if(u >= 1.0d-6)then
      aparam_u = aparam*(1.0d0 + shell*(1.0d0 - exp(-gamma*u))/u)
   else
      aparam_u = aparam*(1.0d0 + shell*gamma)
   end if
   return
end function aparam_u
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
real(kind=8) function Ux(energy,delta)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns pairing shifted excitation energy
!
!   Dependencies:
!
!     Modules:
!
!        None
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
   implicit none
   real(kind=8), intent(in) :: energy,delta
!--------------------------------------------------------------------
   Ux = 0.0d0
   if(energy-delta > 0.0d0)Ux = energy - delta
   return
end function Ux
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine cumm_rho(nfit,elv,ia,level_param,vib_enh,rot_enh,cumrho)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function returns the cummulative level density
!
!   Dependencies:
!
!   Dependencies:
!
!     Modules:
!
!        rhoe
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
   use constants
   implicit none
   integer(kind=4), intent(in) :: nfit
   real(kind=8), intent(in) :: elv(nfit)
   integer(kind=4), intent(in) :: ia
   real(kind=8), intent(in) :: level_param(20)
   real(kind=8), intent(in) :: vib_enh(3), rot_enh(5)
   real(kind=8), intent(out) :: cumrho(nfit)
!----------------------------------------------------------------------
   real(kind=8) :: rrho, apu, sg2
   real(kind=8) :: e,de
   integer(kind=4), parameter :: nmax = 10000
   real(kind=8),allocatable :: cumr(:)
   integer(kind=4) :: i, j
   real(kind=8) :: deriv
   real(kind=8) :: K_vib, K_rot
   real(kind=8) :: E0, T
!----------------------------------------------------------------------
   rrho = 0.0d0
   cumrho(1:nfit) = 0.0d0
   if(.not. allocated(cumr))allocate(cumr(0:nmax))
   cumr(0:nmax) = 0.0d0
   de = (elv(nfit)+0.02d0)/real(nmax,kind=8)
   T = level_param(14)
   E0 = level_param(15)
   if(E0 > 0.0d0)then
      cumr(0) = exp(-E0/T)
   else
      cumr(0) = 1.0d0
   end if
   do i = 1, nmax
      e = real(i,kind=8)*de
      call rhoe(e, level_param, vib_enh, rot_enh, ia, rrho,     &
                apu, sg2, K_vib, K_rot)
      cumr(i) = cumr(i-1) + rrho*de
 
   end do
!-----  Now interpolate points
      cumrho(1) = cumr(0)
   do i = 2, nfit
      j = int(elv(i)/de)
      e = real(j,kind=8)*de
      if(j > 0)then
         deriv = (cumr(j+1) - cumr(j-1))/(2.*de)
         cumrho(i) = cumr(j) + deriv*(elv(i)-e)
      else
         cumrho(i) = cumrho(i-1)
      end if
   end do
   return
end subroutine cumm_rho
!
!*******************************************************************************
!
subroutine fit_D0(D0exp,sep_e,xj,ipar,lev_param,fit_aparam,vib_enh,rot_enh,iA)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine fits D0 by adjusting the level-density a-parameter 
!    or the shell corrections
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        Find_T_E0
!
!     External functions:
!
!        lev_space
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
   implicit none
   real(kind=8), intent(in) :: D0exp
   real(kind=8), intent(in) :: sep_e
   real(kind=8), intent(in) :: xj
   integer(kind=4), intent(in) :: ipar
   real(kind=8), intent(inout) :: lev_param(20)
   logical, intent(in) :: fit_aparam
   real(kind=8), intent(inout) :: vib_enh(3), rot_enh(5)
   integer(kind=4), intent(in) :: iA
!------------------------------------------------------------------------
   real(kind=8) :: deriv
   real(kind=8) :: ap,am,da,anew
   real(kind=8) :: Wp,Wm,dW,Wnew
   real(kind=8) :: D0,D0p,D0m
   integer(kind=4) :: icount
   real(kind=8) :: a_old
!------   External Functions         ----------------------------------------------
   real(kind=8) :: lev_space
!----------------------------------------------------------------------------------
   if(D0exp <= 0.0d0)return
   a_old = lev_param(1)
   icount = 0
   da = 0.05d0
   if(fit_aparam)then
      D0 = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
      if(D0 < D0exp)then
         ap = lev_param(1)
         D0p = D0
 1       am = ap-da
         lev_param(1) = am
         call Find_T_E0(ia,lev_param,vib_enh,rot_enh)

         D0m = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
         if(abs(D0m-D0exp)/D0exp < 1.0e-4)then       ! it agrees
            return
         end if
         if(D0m < D0exp)then           !   still less than so keep going
            ap = am
            D0p = D0m
            goto 1
         end if
      else
         am = lev_param(1)
         D0m = D0
 2       ap = am + da
         lev_param(1) = ap
         call Find_T_E0(ia,lev_param,vib_enh,rot_enh)
         D0p = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
         if(abs(D0p-D0exp)/D0exp < 1.0e-4)then       ! it agrees
            return
         end if
         if(D0p > D0exp)then           !   still greater than so keep going
            am = ap
            D0m = D0p
            goto 2
         end if
      end if
      lev_param(1) = ap
      call Find_T_E0(ia,lev_param,vib_enh,rot_enh)
      D0p = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
      lev_param(1) = am
      call Find_T_E0(ia,lev_param,vib_enh,rot_enh)
      D0m = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
      deriv = (D0p-D0m)/(ap-am)
      anew = (D0exp-D0m)/deriv + am
      lev_param(1) = anew
      call Find_T_E0(ia,lev_param,vib_enh,rot_enh)
      D0 = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
      return
   end if
   dW = 0.01d0
   if(.not. fit_aparam)then
      D0 = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
      if(D0 < D0exp)then
         Wp = lev_param(4)
         D0p = D0
 11      Wm = Wp-dW
         lev_param(4) = Wm
         call Find_T_E0(ia,lev_param,vib_enh,rot_enh)
         D0m = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
         if(abs(D0m-D0exp)/D0exp < 1.0e-4)then       ! it agrees
            return
         end if
         if(D0m < D0exp)then           !   still less than so keep going
            Wp = Wm
            D0p = D0m
            goto 11
         end if
      else
         Wm = lev_param(4)
         D0m = D0
 12      Wp = Wm + dW
         lev_param(4) = Wp
         call Find_T_E0(ia,lev_param,vib_enh,rot_enh)
         D0p = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
         if(abs(D0p-D0exp)/D0exp < 1.0e-4)then       ! it agrees
            return
         end if
         if(D0p > D0exp)then           !   still greater than so keep going
            Wm = Wp
            D0m = D0p
            goto 12
         end if
      end if
      lev_param(4) = Wp
      call Find_T_E0(ia,lev_param,vib_enh,rot_enh)
      D0p = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
      lev_param(4) = Wm
      call Find_T_E0(ia,lev_param,vib_enh,rot_enh)
      D0m = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
      deriv = (D0p-D0m)/(Wp-Wm)
      Wnew = (D0exp - D0m)/deriv + Wm
      lev_param(4) = Wnew
      call Find_T_E0(ia,lev_param,vib_enh,rot_enh)
      D0 = lev_space(sep_e,xj,ipar,0,lev_param,vib_enh,rot_enh,iA)
      return
   end if
   return
end subroutine fit_D0
!
!*******************************************************************************
!
real(kind=8) function enhance_vib(mode, A, U, E, apu, Shell, enh)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the vibrational enhancement factor
!
!     Modules:
!
!        None
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
   implicit none
   integer(kind=4), intent(in) :: mode
   real(kind=8), intent(in) :: A, U, E, apu, Shell, enh(3)
!-----------------------------------------------------
   integer(kind=4) :: L
   real(kind=8) :: T, K_v, ap
   real(kind=8) :: C, w(6), xn(6), g(6), dS, dU, pi, L_hat
   real(kind=8) :: damp
   integer(kind=4) :: vib_mode
!----------   Calculation   --------------------------
   pi = 2.0d0*asin(1.0d0)
   T = 0.0d0
   vib_mode = mode
   if(mode < 3)then
      if( U >= 0.0d0)then
         T = sqrt(U/apu)
      end if
   else
      vib_mode = vib_mode - 2
      if( U >= 0.0d0)then
         ap = A/13.
         T = sqrt(U/ap)
      end if
   end if
   enhance_vib = 1.0d0
   K_v = 1.0d0
   if(vib_mode == 1)then
      K_v = 1.0d0
      K_v = exp(0.05555*A**(2./3.)*T**(4./3.))
      K_v = max(K_v,1.0d0)
      K_v = 1.0d0 + (K_v - 1.0d0)*damp(E, enh(2), enh(3))*enh(1)
      enhance_vib = K_v 
      return
   end if
   if(vib_mode == 2) then
      K_v = 1.0d0
      C = 0.0075d0*A**(1.0d0/3.0d0)
      w(2) =  65.0d0*A**(-5.0d0/6.0d0)/(1.0d0 + 0.05d0*Shell)  
      w(3) = 100.0d0*A**(-5.0d0/6.0d0)/(1.0d0 + 0.05d0*Shell)  
      dS = 0.0d0
      dU = 0.0d0
      do L = 2,3
         L_hat = (2.0d0*real(L,kind=8) + 1.0d0)
         g(L) = C*(w(L)**2 + 4.0d0*(pi*T)**2)
         xn(L) = exp(-0.5d0*g(L)/w(L))/(exp(w(L)/T) - 1.0d0)
         dS = dS + L_hat*((1.0d0+xn(L))*log(1.0d0 + xn(L)) -    &
                   xn(L)*log(xn(L)))
         dU = dU + L_hat*w(L)*xn(L)
      end do
      K_v = exp(dS - dU/T)
      enhance_vib = K_v
   end if     

   return
end function enhance_vib
!
!*******************************************************************************
!
real(kind=8) function enhance_rot(mode, sig2, sig2_perp, enh, Ex)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the rotational enhancement factor
!
!     Modules:
!
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
   use constants
   integer(kind=4), intent(in) :: mode
   real(kind=8), intent(in) :: sig2
   real(kind=8), intent(in) :: sig2_perp
   real(kind=8), intent(in) :: enh(5)
   real(kind=8), intent(in) :: Ex
   real(kind=8) :: K_r
   real(kind=8) :: damp
!----------------------------------------------------------------------------------
   K_r = 0.0d0
!   if(mode > 0)K_r = enh(4)*sig2
   if(mode > 0)K_r = sig2_perp
   if(mode == 2 )then
      K_r = 2.0d0*K_r
   end if
   if(mode == 3)then
      K_r = K_r*sqrt(pi/2.0d0)*enh(5)*sqrt(sig2)
   end if
   if(mode == 4)then
      K_r = 2.0d0*K_r*sqrt(pi/2.0d0)*enh(5)*sqrt(sig2)
   end if
   K_r = max(K_r,1.0d0)
   K_r = 1.0d0 + (K_r - 1.0d0)*damp(Ex, enh(2), enh(3))*enh(1)
   enhance_rot = K_r
   return
end function enhance_rot 
!
!*******************************************************************************
!        
real(kind=8) function damp(E, B, C)
!
!*******************************************************************************
!
!  Discussion:
!
!    This function computes the damping factor for the rotational and vibrational
!    enhancement factors
!
!     Modules:
!
!        None
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
   real(kind=8), intent(in) :: E, B, C
!------------------------------------------------------------------------
   real(kind=8) :: exponent
!----------------------------------------------------------------------------------
   damp = 0.0d0      
   exponent = (E-B)/C
   if(exponent <= -100.0d0)exponent = -100.0d0
   if(exponent >= 100.0d0)exponent = 100.0d0
   damp = 1.0d0/(1.0d0 + exp(exponent))
   return
end function damp
!
!*******************************************************************************
!
subroutine Find_T_E0(ia, level_param, vib_enhance, rot_enhance)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates T and E0 for Gilbert and Cameron from
!    continuity of rho and d rho/dE at ematch
!
!     Modules:
!
!        None
!
!     Subroutines:
!
!        rho_BFM
!
!     External functions:
!
!        real(kind=8) :: aparam_u
!        real(kind=8) :: sig2_param
!        real(kind=8) :: enhance_rot
!        real(kind=8) :: enhance_vib    
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
   implicit none
!---------------  data passed in 
   integer(kind=4), intent(in) :: ia    
   real(kind=8), intent(inout) :: level_param(20)
   real(kind=8), intent(inout) :: vib_enhance(3), rot_enhance(5)
!---------------   Internal data      
   real(kind=8) :: A
   real(kind=8) :: aparam, spin_cut, delta, shell, gamma
   real(kind=8) :: ematch, ecut, sg2cut, low_e_mod
   real(kind=8) :: sig2_perp
   integer(kind=4) :: rot_mode, vib_mode
   real(kind=8) :: enhance, K_rot, K_vib
      
   real(kind=8) :: E, U, de
      
   real(kind=8) :: deriv
      
   real(kind=8) :: apu, sig2
   real(kind=8) :: T, E0
      
   real(kind=8) rho_F0, rho_Fm, rho_Fp
      
!------     External functions     ------------------------------------------------
   real(kind=8) :: aparam_u
   real(kind=8) :: sig2_param
   real(kind=8) :: enhance_rot
   real(kind=8) :: enhance_vib    
!----------------------------------------------------------------------------------
   A = real(ia, kind=8)
      
   rho_F0 = 0.0d0
   rho_Fm = 0.0d0
   rho_Fp = 0.0d0
      
   aparam = level_param(1)
   spin_cut = level_param(2)
   delta = level_param(3)
   shell = level_param(4)
   gamma = level_param(5)
   ematch = level_param(6)
   ecut = level_param(7)
   sg2cut = level_param(8)
   low_e_mod = level_param(9)

   rot_mode = nint(level_param(10))
   vib_mode = nint(level_param(11))

!------------    start of calculation   --------------------
   de = 0.01d0

   level_param(6) = 0.0d0         !   set ematch to zero to get sig2
                                     !   without interpolation
   E = ematch
   U = E - delta

   apu = aparam_u(U,aparam,shell,gamma)
   sig2 = sig2_param(E,level_param,A)
   sig2_perp = rot_enhance(4)*sig2

   K_rot = enhance_rot(rot_mode, sig2, sig2_perp, rot_enhance, E)
   K_vib = enhance_vib(vib_mode, A, U, E, apu, Shell, vib_enhance)

   if(nint(level_param(19)) == 3 .or. nint(level_param(19)) == 5) sig2 = sig2_perp

   call rho_BFM(U, apu, sig2, rho_F0)

   enhance = K_rot*K_vib
   rho_F0 = rho_F0*enhance

   E = ematch - de
   U = E - delta

   apu = aparam_u(U,aparam,shell,gamma)
   sig2 = sig2_param(E,level_param,A)
   sig2_perp = rot_enhance(4)*sig2

   K_rot = enhance_rot(rot_mode, sig2, sig2_perp, rot_enhance, E)
   K_vib = enhance_vib(vib_mode, A, U, E, apu, Shell, vib_enhance)

   if(nint(level_param(19)) == 3 .or. nint(level_param(19)) == 5) sig2 = sig2_perp

   call rho_BFM(U, apu, sig2, rho_Fm)

   enhance = K_rot*K_vib
   rho_Fm = rho_Fm*enhance


   E = ematch + de
   U = E - delta

   apu = aparam_u(U,aparam,shell,gamma)
   sig2 = sig2_param(E,level_param,A)
   sig2_perp = rot_enhance(4)*sig2

   K_rot = enhance_rot(rot_mode, sig2, sig2_perp, rot_enhance, E)
   K_vib = enhance_vib(vib_mode, A, U, E, apu, Shell, vib_enhance)

   if(nint(level_param(19)) == 3 .or. nint(level_param(19)) == 5) sig2 = sig2_perp

   call rho_BFM(U, apu, sig2, rho_Fp)

   enhance = K_rot*K_vib
   rho_Fp = rho_Fp*enhance


   deriv = (rho_Fp - rho_Fm)/(2.0d0*de)
   T = rho_F0/deriv
   E0 = ematch - T*log(T*rho_F0)              !  constant temperature

   level_param(6) = ematch                    !   reset to correct value

   level_param(14) = T
   level_param(15) = E0

   return
end subroutine Find_T_E0
