!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine prints out (proj,Xk) spectra (k = g,n,p,d,t,h,a) to appropriate
!    file in the directory 
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
!   do i = 0, max_num
!      energy = real(i,kind=8)*de_spec2
!      write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,                                              &
!         (x_particle_Spectrum_res(i,k)*cs_scale,k=0, 6),                                        &
!         (x_particle_Spectrum_res_dist(i,k),k=0, 6)    
!      energy = energy + de_spec2
!      write(24,'(1x,f10.4,14(1x,1pe15.7))')energy,                                              &
!         (x_particle_Spectrum_res(i,k)*cs_scale,k=0, 6),                                        &
!         (x_particle_Spectrum_res_dist(i,k),k=0, 6)    
!   end do
   close(unit=24)
   if(allocated(x_particle_Spectrum_res))deallocate(x_particle_Spectrum_res)
   if(allocated(x_particle_Spectrum_res_dist))deallocate(x_particle_Spectrum_res_dist)
   return
end subroutine print_x_particle_spectra
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine prints out pre-equilibrium emission spectra to appropriate
!    file in the directory 
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
subroutine print_preeq_spectra(ilib_dir, lib_dir, ifile, file_name,                             &
                               e_in, de_spec2, num_e, Preeq_spect, Preeq_spect_full, smear)
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
!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine prints out direct (coupled-channels) and dwba emission spectra 
!    to appropriate file in the directory 
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
subroutine print_direct_spectra(ilib_dir, lib_dir, ifile, file_name,                            &
                                e_in, de_spec2, max_num, num_e,                                 &
                                direct_Spectrum, dwba_Spectrum, smear)
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
   real(kind=8), intent(in) :: direct_Spectrum(0:num_e)
   real(kind=8), intent(in) :: dwba_spectrum(0:num_e)
   real(kind=8), intent(in) :: smear(0:6)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   integer(kind=4) :: i, j
   real(kind=8) :: e1, e2
   real(kind=8) :: energy
   real(kind=8) :: xnorm1, xnorm2, xnorm3, xnorm4
   integer(kind=4) :: idir
   character(len=132) :: directory
   real(kind=8), allocatable :: res1(:)
   real(kind=8), allocatable :: res2(:)   


   if(.not. allocated(res1))allocate(res1(0:num_e))
   if(.not. allocated(res2))allocate(res2(0:num_e))

   res1(0:num_e) = 0.0d0 
   res2(0:num_e) = 0.0d0 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Fold in a resolution smear
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   do i = 0, num_e
      e1 = real(i,kind=8)*de_spec2
      do j = 0, num_e
         e2 = real(j,kind=8)*de_spec2
         res1(j) = res1(j) +                                                                    &
            direct_Spectrum(i)*exp(-((e1-e2)/smear(1))**2/2.0d0)/sqrt(2.0d0*pi*smear(1)**2)
         res2(j) = res2(j) +                                                                    &
            dwba_Spectrum(i)*exp(-((e1-e2)/smear(1))**2/2.0d0)/sqrt(2.0d0*pi*smear(1)**2)
      end do
   end do
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Normalize for distribution
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   xnorm1 = 0.0d0
   xnorm2 = 0.0d0
   xnorm3 = 0.0d0
   xnorm4 = 0.0d0
   do i = 0, num_e
      e1 = real(i,kind=8)*de_spec2
      xnorm1 = xnorm1 + res1(i)*de_spec2
      xnorm2 = xnorm2 + res2(i)*de_spec2
      xnorm3 = xnorm3 + direct_Spectrum(i)*de_spec2
      xnorm4 = xnorm4 + dwba_Spectrum(i)*de_spec2
   end do
   if(xnorm1 > 0.0d0)then
      do i = 0, num_e
         e1 = real(i,kind=8)*de_spec2
         res1(i) = res1(i)*xnorm3/xnorm1
         if(res1(i) < 1.0d-20)res1(i) = 1.0d-20
      end do
   end if
   if(xnorm2 > 0.0d0)then
      do i = 0, num_e
         e1 = real(i,kind=8)*de_spec2
         res2(i) = res2(i)*xnorm4/xnorm2
         if(res2(i) < 1.0d-20)res2(i) = 1.0d-20
      end do
   end if

   if(iproc == 0)then
      directory(1:ilib_dir) = lib_dir(1:ilib_dir)
      idir = ilib_dir + 1
      directory(idir:idir) = '/'
      idir = idir + 1
      directory(idir:idir+6) = 'Spectra'
      idir = idir + 7
      directory(idir:idir) = '/'

      open(unit=24,file=                                                                        &
          directory(1:idir)//file_name(1:ifile)//'-Direct-Spectrum.dat',status='unknown')
      write(24,'(''#Outgoing spectra for direct and DWBA reactions '')')
      write(24,'(''#Calculated as barns/MeV '')')
      write(24,'(''#Incident energy = '',f10.4)')e_in
      write(24,'(''#'',10x,10x,''Direct'',12x,''DWBA'')')
      do i = 0, max_num
         energy = real(i,kind=8)*de_spec2
         write(24,'(1x,f10.4,14(1x,1pe15.7))')energy, direct_Spectrum(i), dwba_Spectrum(i)
         energy = energy + de_spec2
         write(24,'(1x,f10.4,14(1x,1pe15.7))')energy, direct_Spectrum(i), dwba_Spectrum(i)
      end do
!      do i = 0, max_num
!         energy = real(i,kind=8)*de_spec2
!          write(24,'(1x,f10.4,14(1x,1pe15.7))')energy, res1(i), res2(i)
!         energy = energy + de_spec2
!         write(24,'(1x,f10.4,14(1x,1pe15.7))')energy, res1(i), res2(i)
!      end do
      close(unit=24)
   end if
   if(allocated(res1))deallocate(res1)
   if(allocated(res2))deallocate(res2)
   return
end subroutine print_direct_spectra
