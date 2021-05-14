!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine prints out (proj,Xk) spectra (k = g,n,p,d,t,h,a) to appropriate
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
subroutine print_x_particle_spectra(ilib_dir, lib_dir, ifile, file_name,                        &
                                    e_in, de_spec2, max_num, num_e,                             &
                                    x_particle_spectrum, smear)
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
   integer(kind=4), intent(in):: max_num
   integer(kind=4), intent(in):: num_e
   real(kind=8), intent(in) :: x_particle_spectrum(0:num_e,0:6)
   real(kind=8), intent(in) :: smear(0:6)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer(kind=4) :: i, j, k
   real(kind=8) :: e1, e2
   real(kind=8) :: energy
   real(kind=8) :: xnorm1, xnorm2
   integer(kind=4) :: idir
   character(len=132) :: directory
   real(kind=8), allocatable :: x_particle_spectrum_res(:,:)
   real(kind=8), allocatable :: x_particle_spectrum_res_dist(:,:)   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   if(.not. allocated(x_particle_Spectrum_res))allocate(x_particle_Spectrum_res(0:num_e,0:6))
   if(.not. allocated(x_particle_Spectrum_res_dist))allocate(x_particle_Spectrum_res_dist(0:num_e,0:6))
   x_particle_Spectrum_res(0:num_e,0:6) = 0.0d0 
   x_particle_Spectrum_res_dist(0:num_e,0:6) = 0.0d0 

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'
   idir = idir + 1
   directory(idir:idir+6) = 'Spectra'
   idir = idir + 7
   directory(idir:idir) = '/'

   open(unit=24,file=                                                                           &
        directory(1:idir)//file_name(1:ifile)//'-xParticle-Spectrum.dat',status='unknown')
   write(24,'(''#Outgoing spectra for each particle type '')')
   write(24,'(''#Incident energy = '',f10.4)')e_in
   write(24,'(''#Calculated as barns/MeV '')')
   write(24,'(''#k='',8x,7(10x,i2,4x))')(k, k = 0, 6)
   write(24,'(''#'',10x,7(3x,''----------------''))')
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Fold in a resolution smear
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do k = 0, 6
      do i = 0, num_e
         e1 = real(i,kind=8)*de_spec2
         do j = 0, num_e
            e2 = real(j,kind=8)*de_spec2
            x_particle_spectrum_res(j,k) = x_particle_spectrum_res(j,k) +                       &
                 x_particle_Spectrum(i,k)*exp(-((e1-e2)/smear(k))**2/2.0d0)/                    &
            sqrt(2.0d0*pi*smear(k)**2)
         end do
      end do
   end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Normalize the distribution
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do k = 0, 6
      xnorm1 = 0.0d0
      xnorm2 = 0.0d0
      do i = 0, num_e
         e1 = real(i,kind=8)*de_spec2
         xnorm1 = xnorm1 + x_particle_Spectrum_res(i,k)*de_spec2
         xnorm2 = xnorm2 + x_particle_Spectrum(i,k)*de_spec2
      end do
      if(xnorm1 > 0.0d0)then
         do i = 0, num_e
            e1 = real(i,kind=8)*de_spec2
            x_particle_Spectrum_res_dist(i,k) = x_particle_Spectrum_res(i,k)/xnorm1
            if(x_particle_Spectrum_res_dist(i,k) < 1.0d-20)                                     &
               x_particle_Spectrum_res_dist(i,k) = 1.0d-20
         end do
      end if
      if(xnorm2 > 0.0d0)then
         do i = 0, num_e
            e1 = real(i,kind=8)*de_spec2
            x_particle_Spectrum_res(i,k) = x_particle_Spectrum_res(i,k)*xnorm2/xnorm1
            if(x_particle_Spectrum_res(i,k) < 1.0d-20)x_particle_Spectrum_res(i,k) = 1.0d-20
         end do
      end if
   end do

   do i = 0, max_num
      energy = real(i,kind=8)*de_spec2
      write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,                                              &
         (x_particle_Spectrum(i,k)*cs_scale,k=0, 6)   
      energy = energy + de_spec2
      write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,                                              &
         (x_particle_Spectrum(i,k)*cs_scale,k=0, 6)   
   end do
   close(unit=24)
   if(allocated(x_particle_Spectrum_res))deallocate(x_particle_Spectrum_res)
   if(allocated(x_particle_Spectrum_res_dist))deallocate(x_particle_Spectrum_res_dist)
   return
end subroutine print_x_particle_spectra
