!
!*******************************************************************************
!
subroutine get_binding_energy(data_path, len_path,         &
                              iz1, ia1, me, be, sep_e)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine reads in binding energies and stores in arrays
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
  use nodeinfo
  use useful_data
  use particles_def
  use nodeinfo
  implicit none
  character(len=200), intent(in) :: data_path
  integer(kind=4), intent(in) :: len_path
  integer(kind=4), intent(in) :: iz1,ia1
  real(kind=8), intent(out) :: me, be, sep_e(6)
!-------------------------------------------------------------------
  integer(kind=4) :: in1
  integer(kind=4) :: istart(0:400)        !  starting position for A
  integer(kind=4) :: num_nuc
  integer(kind=4), parameter :: max_nuc=20000        ! accomodate up to 100000 nuclei
  real(kind=8) :: mass_ex(max_nuc)
  real(kind=8) binding_e(max_nuc)
  integer(kind=4) :: iz(max_nuc),ia(max_nuc)
  integer(kind=4) :: inn,izz,iaa,izb,iz2,in2,ia2
  integer(kind=4) :: i,j,k
  real(kind=8) :: xtemp
  real(kind=8) :: mass_ex_aw,mass_ex_mn
  real(kind=8) :: u,mh,mn
  integer(kind=4) :: eof
  save mass_ex, binding_e, istart, iz, ia
!-------------------------------------------------------------------
  u = 931.494013d0
  mh = 938.782982d0
  mn = 939.565336d0
  istart(0) = 1
  if(i_bind == 0)then
     open(unit=50,file=data_path(1:len_path)//'mass-frdm95.dat',status='old')
     izb=0
     num_nuc = 0
!---   Read in file and store in array
     do i = 1, 4
        read(50,*)
     end do
     do i = 1, max_nuc
        read(50,'(2(1x,i3),6x,f9.3,1x,f9.3)', iostat = eof)    &
             izz,iaa,mass_ex_aw,mass_ex_mn
        if(eof == 1)exit
        if(mass_file == 'aw')then
           if(abs(mass_ex_aw) > 1.0d-4)then   !  fail safe in case exp is not known
              mass_ex(i) = mass_ex_aw
           else
              mass_ex(i) = mass_ex_mn
           end if
        elseif(mass_file == 'mn')then
           if(abs(mass_ex_mn) > 1.0-4)then         !  use exp when MN=0.0
              mass_ex(i) = mass_ex_mn
           else
              mass_ex(i) = mass_ex_aw
           end if
        end if
        inn = iaa-izz
        if(izz > izb)then
           istart(izz) = i
           izb = izz
        end if
        iz(i) = izz
        ia(i) = iaa
!-----------------     Convert to binding energy
        binding_e(i) = real(izz,kind=8)*particle(2)%ME +       &
            real(iaa-izz,kind=8)*particle(1)%ME - mass_ex(i)
        num_nuc = num_nuc+1
     end do
     close(unit=50)
     istart(iz(num_nuc)+1) = num_nuc+1     
  end if
  i_bind = 1
!----  With binding energies, now compute separation energies
  in1 = ia1-iz1
  do j = istart(iz1), istart(iz1+1) - 1
     if(ia1 == ia(j))then
        me = mass_ex(j)
        be = binding_e(j)
        goto 101
     else
        me = -1.01d6
     end if
  end do
  return
 101  continue
!----   now for the separation energies
  do j = 1, 6
     iz2 = iz1 - particle(j)%Z
     ia2 = ia1 - particle(j)%A
     in2 = ia2 - iz2
     xtemp = 0.0
     do k = istart(iz2), istart(iz2+1) - 1
        if(ia2 == ia(k))xtemp = mass_ex(k)
     end do
     sep_e(j) = xtemp + particle(j)%ME - me
  end do
  return
end
