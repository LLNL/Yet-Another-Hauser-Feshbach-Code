!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine prints out pre-equilibrium emission spectra to appropriate
!    file in the directory 
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
!*****************************************************************************80
!
subroutine print_preeq_spectra(ilib_dir, lib_dir, ifile, file_name,                   &
                               e_in, de_spec2, num_e, Preeq_spect,                    &
                               Preeq_spect_full, smear)
   use variable_kinds
   use options
   use constants
   use nodeinfo
   use nuclei
   use Channel_info
   use particles_def
   implicit none  
   integer(kind=4), intent(in):: ilib_dir
   character(len=132), intent(in) :: lib_dir
   integer(kind=4), intent(in):: ifile
   character(len=40), intent(in) :: file_name
   real(kind=8), intent(in) :: e_in
   real(kind=8), intent(in) :: de_spec2
   integer(kind=4), intent(in):: num_e
   real(kind=8), intent(in) :: Preeq_spect(0:6,0:num_e)
   real(kind=8), intent(in) :: Preeq_spect_full(0:6,0:num_e)
   real(kind=8), intent(in) :: smear(0:6)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer(kind=4) :: i, j, k
   real(kind=8) :: e1, e2
   real(kind=8) :: energy
   real(kind=8) :: xnorm1, xnorm2
   integer(kind=4) :: idir
   character(len=132) :: directory
   real(kind=8), allocatable :: resp(:,:)
   real(kind=8), allocatable :: resp2(:,:)   
!-------------------------------------------------------------------------------------------

   if(.not. allocated(resp))allocate(resp(0:6,0:num_e))
   if(.not. allocated(resp2))allocate(resp2(0:6,0:num_e))
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Fold in a resolution smear
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do k = 0,6
      resp(k,0:num_e) = 0.0d0 
      resp2(k,0:num_e) = 0.0d0 
      do i = 0, num_e
         e1 = real(i,kind=8)*de_spec2
         do j = 0, num_e
            e2 = real(j,kind=8)*de_spec2
            resp(k,j) = resp(k,j) +                                                             &
              Preeq_spect(k,i)*exp(-((e1-e2)/smear(1))**2/2.0d0)/sqrt(2.0d0*pi*smear(1)**2)
           resp2(k,j) = resp2(k,j) +                                                            &
              Preeq_spect_full(k,i)*exp(-((e1-e2)/smear(1))**2/2.0d0)/sqrt(2.0d0*pi*smear(1)**2)
         end do
      end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Normalize for distribution
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      xnorm1 = 0.0d0
      xnorm2 = 0.0d0
      do i = 0, num_e
         e1 = real(i,kind=8)*de_spec2
         xnorm1 = xnorm1 + resp(k,i)*de_spec2
         xnorm2 = xnorm2 + Preeq_Spect(k,i)*de_spec2
      end do
      if(xnorm1 > 0.0d0)then
         do i = 0, num_e
            e1 = real(i,kind=8)*de_spec2
            resp(k,i) = resp(k,i)*xnorm2/xnorm1
            if(resp(k,i) < 1.0d-20)resp(k,i) = 1.0d-20
         end do
      end if
      xnorm1 = 0.0d0
      xnorm2 = 0.0d0
      do i = 0, num_e
         e1 = real(i,kind=8)*de_spec2
         xnorm1 = xnorm1 + resp2(k,i)*de_spec2
         xnorm2 = xnorm2 + Preeq_Spect_full(k,i)*de_spec2
      end do
      if(xnorm1 > 0.0d0)then
         do i = 0, num_e
            e1 = real(i,kind=8)*de_spec2
            resp2(k,i) = resp2(k,i)*xnorm2/xnorm1
            if(resp2(k,i) < 1.0d-20)resp2(k,i) = 1.0d-20
         end do
      end if
   end do

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   idir = idir + 1
   directory(idir:idir+6) = 'Spectra'
   idir = idir + 7
   directory(idir:idir) = '/'

   open(unit=24,file=                                                                           &
       directory(1:idir)//file_name(1:ifile)//'-Preeq-Spectrum.dat',status='unknown')
   write(24,'(''#Pre-equilibirum spectrum for each particle type '')')
   write(24,'(''#Calculated as barns/MeV '')')
   write(24,'(''#Incident energy = '',f10.4)')e_in
   if(fission)then
      write(24,'(''#First set includes fission events, second set excludes fission events'')')
      write(24,'(''#'',10x,14(1x,13x,i2))')(k, k = 0, nucleus(1)%num_decay-1),                  &
                                           (k, k = 0, nucleus(1)%num_decay-1)
   else
      write(24,'(''#'',10x,14(1x,13x,i2))')(k, k = 0, nucleus(1)%num_decay-1)
   end if
   do i = 0, num_e
      energy = real(i,kind=8)*de_spec2
      if(fission)then
         write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,                                           &
                   (resp2(k,i),k = 0, nucleus(1)%num_decay-1),                                  &
                   (resp(k,i),k = 0, nucleus(1)%num_decay-1)
         energy = energy + de_spec2
         write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,                                           &
                    (resp2(k,i),k = 0, nucleus(1)%num_decay-1),                                 &
                    (resp(k,i),k = 0, nucleus(1)%num_decay-1)
      else
         write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,(resp2(k,i),k = 0, nucleus(1)%num_decay-1)
         energy = energy + de_spec2
         write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,(resp2(k,i),k = 0, nucleus(1)%num_decay-1)
      end if
   end do
   close(unit=24)
   if(allocated(resp))deallocate(resp)
   if(allocated(resp2))deallocate(resp2)
   return
end subroutine print_preeq_spectra
