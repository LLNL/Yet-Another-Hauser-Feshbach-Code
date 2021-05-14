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
!   Dependencies:
!
!     Modules:
!
!        nodeinfo
!        useful_data
!        particles_def
!        nodeinfo
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
  use nodeinfo
  use useful_data
  use particles_def
  use nodeinfo
  implicit none
  character(len=200), intent(in) :: data_path
  integer(kind=4), intent(in) :: len_path
  integer(kind=4), intent(in) :: iz1,ia1
  real(kind=8), intent(out) :: me, be, sep_e(0:6)
!-------------------------------------------------------------------
  integer(kind=4) :: in1
  integer(kind=4), allocatable, dimension(:) :: istart    !  starting position for Z
  integer(kind=4) :: num_nuc, max_z
  integer(kind=4) :: max_nuc                              ! Number of nuclei in file
  real(kind=8), allocatable :: mass_ex(:)
  real(kind=8), allocatable :: binding_e(:)
  integer(kind=4),allocatable :: iz(:), ia(:)
  integer(kind=4) :: inn, izz, iaa, izb, iz2, in2, ia2
  integer(kind=4) :: i, j, k
  real(kind=8) :: xtemp
  real(kind=8) :: mass_ex_aw, mass_ex_mn
  integer(kind=4) :: eof
  integer :: read_error
  character(len=132) :: line

  save mass_ex, binding_e, istart, iz, ia
!-------------------------------------------------------------------

!-------------------------------------------------------------------
!---   If arrays are not allocated, then subroutine hasn't been called before
!---   Therefore, read in mass data
!-------------------------------------------------------------------
!

  if(.not. allocated(mass_ex))then
     open(unit=50,file=data_path(1:len_path)//'mass-frdm95.dat',status='old')
!---   Count how many entries and allocate arrays
     read_error = 0
     max_nuc = 0
     max_z = 0
     do while(read_error == 0)
        read(50,'(a)',iostat = read_error)line
        if(line(1:1) == '#')cycle
        max_nuc = max_nuc + 1
        read(line,'(1x,i3)')izz
        if(izz > max_z)max_z = izz
     end do

     if(.not. allocated(mass_ex))allocate(mass_ex(max_nuc))
     if(.not. allocated(binding_e))allocate(binding_e(max_nuc))
     if(.not. allocated(iz))allocate(iz(max_nuc))
     if(.not. allocated(ia))allocate(ia(max_nuc))
     if(.not. allocated(istart))allocate(istart(0:max_z+1))

     rewind(50)
     izb=0
     num_nuc = 0
!---   Read in file and store in array
     read_error = 0
     istart(0) = 0
     i = 0
     do while(read_error == 0)
        read(50,'(a)',iostat = read_error)line
        if(line(1:1) == '#')cycle
        read(line,'(2(1x,i3),6x,f9.3,1x,f9.3)', iostat = eof)    &
             izz,iaa,mass_ex_aw,mass_ex_mn
        i = i + 1
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
        num_nuc = num_nuc + 1
     end do
     close(unit=50)
     istart(iz(num_nuc)+1) = num_nuc + 1     
  end if
  i_bind = 1
!----  With binding energies, now compute separation energies
  in1 = ia1-iz1
  me = -1.01d6
  if(iz1 < 1)return
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
  sep_e(0) = 0.0d0
  do j = 1, 6
     sep_e(j) = 1.0d3          !  default it to a large value if not found
     iz2 = iz1 - particle(j)%Z
     if(iz2 < 1)cycle
     ia2 = ia1 - particle(j)%A
     if(ia2 < 2)cycle
     in2 = ia2 - iz2
     xtemp = 0.0
     if(istart(iz2) > 0)then
        do k = istart(iz2), istart(iz2+1) - 1
           if(ia2 == ia(k))xtemp = mass_ex(k)
        end do
        sep_e(j) = xtemp + particle(j)%ME - me
     else
     end if
  end do
  return
end subroutine get_binding_energy
