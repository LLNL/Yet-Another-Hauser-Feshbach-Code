!
!*****************************************************************************80
!
!  Discussion:
!
!    This Subroutine to write data on nuclear states: properties and how they 
!    decay to the library directory
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
subroutine print_nuke_data(ilib_dir, lib_dir)
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
!----------------------------------------------------------------------
   integer(kind=4) :: i, iA, k, m
   integer(kind=4) :: num_states
   integer(kind=4) :: idir, ifile
   character(len=132) :: directory
   character(len=132) :: outfile
   character(len=14) iso_label
   character(len=16) mod_label

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Core Directory     -------------------------------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   directory(1:ilib_dir) = lib_dir(1:ilib_dir)
   idir = ilib_dir + 1
   directory(idir:idir) = '/'

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------    Decay pattern for each for each nucleus   --------------+
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   do i = 1, num_comp
      iA = nucleus(i)%A
      ifile = 0
      if(iA < 10)then
         ifile = 1
         write(outfile(1:ifile),'(i1)')iA
      elseif(iA < 99)then
         ifile = 2
         write(outfile(1:ifile),'(i2)')iA
      elseif(iA < 1000)then
         ifile = 3
         write(outfile(1:ifile),'(i3)')iA
      end if
      if(nucleus(i)%atomic_symbol(1:1) == ' ')then
         outfile(ifile+1:ifile+1) = nucleus(i)%atomic_symbol(2:2)
         ifile = ifile + 1
      else
         outfile(ifile+1:ifile+2) = nucleus(i)%atomic_symbol(1:2)
         ifile = ifile + 2
      end if
      open(unit=100, file = directory(1:idir)//outfile(1:ifile)//'-decay-properties.dat',status = 'unknown')            
      write(100,'(''#Decay properties of discrete states for'',a5)')outfile(1:ifile)
      write(100,'(''#Mass = '',1pe23.16,'' MeV'')')nucleus(i)%mass
      write(100,'(''#Mass = '',1pe23.16,'' amu'')')nucleus(i)%mass/mass_u
      write(100,'(''#AMU = '',1pe23.16,'' MeV'')')mass_u
      write(100,'(''#'')')
      num_states = nucleus(i)%ncut
      if(all_discrete_states)num_states = nucleus(i)%num_discrete
      do k = 1, num_states
         write(100,'(''    #       Ex     J   par #decay     Isomer'')')
         write(100,'(''  ---     ------  --- ----  -----  --------------'')')
         iso_label = 'Isomer = False'
         if(nucleus(i)%state(k)%isomer)iso_label = 'Isomer = True '
         mod_label = 'Modified = False'
         if(nucleus(i)%state(k)%state_modified)mod_label = 'Modified = True '
         write(100,'(1x,i4,1x,f10.4,2(1x,f4.1),1x,i5,3x,a14,3x,a14)')                     &
                     k-1,nucleus(i)%state(k)%energy,                                      &
                     nucleus(i)%state(k)%spin,nucleus(i)%state(k)%parity,                 &
                     nucleus(i)%state(k)%nbranch, iso_label, mod_label
         write(100,'(''#      i --->   f'',12x,''branch'',6x,''prob_gamma'',9x,''prob_ic'',9x,''Modified'')')
         write(100,'(''#    ---      ---'',3(5x,''-------------''),5x,''--------'')')
         do m = 1, nucleus(i)%state(k)%nbranch
            write(100,'(4x,i4,'' --->'',i4,3(3x,1pe15.7),5x,l4)')                         &
                  k-1,nucleus(i)%state(k)%ibranch(m)-1,                                   &
                      nucleus(i)%state(k)%branch(m),                                      &
                      nucleus(i)%state(k)%p_gamma(m),                                     &
                      nucleus(i)%state(k)%p_ic(m),                                        &
                      nucleus(i)%state(k)%branch_modified(m)
         end do
         write(100,'(''#'')')
      end do
   end do
  return
end subroutine print_nuke_data
