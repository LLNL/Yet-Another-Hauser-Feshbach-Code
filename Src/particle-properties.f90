!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
subroutine particle_data
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine puts particle data into appropriate arrays for later use
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
   use variable_kinds
   use constants
   use particles_def
   use directory_structure
   implicit none
!--------------   data for emitted particle types
   allocate(particle(-2:7))
   particle(-2:7)%name = '        '
!--------------   Blank - used for population calculations
   particle(-2)%Z = 0
   particle(-2)%A = 0
   particle(-2)%spin = 0.0
   particle(-2)%par = 0.0
   particle(-2)%ME = 0.0d0
   particle(-2)%mass = 0.0d0
   particle(-2)%label = 'X'
   particle(-2)%name = 'PopDecay'
   particle(-2)%opt_pot_set = .false.
   particle(-2)%max_opt_pot = 1
   particle(-2)%om_option = 0  
   particle(-2)%nume = 0
!--------------   electron
   particle(-1)%Z = 0
   particle(-1)%A = 0
   particle(-1)%spin = 0.5
   particle(-1)%par = 1.0
   particle(-1)%ME = 0.00000
   particle(-1)%mass = 5.485799090d-4*mass_u
   particle(-1)%label = 'e'
   particle(-1)%name = 'electron'
   particle(-1)%opt_pot_set = .false.
   particle(-1)%max_opt_pot = -1
   particle(-1)%om_option = -1
   particle(-1)%nume = 0
!--------------   photon
   particle(0)%Z = 0
   particle(0)%A = 0
   particle(0)%spin = 1.0
   particle(0)%par = 1.0
   particle(0)%ME = 0.00000
   particle(0)%mass = 0.00000
   particle(0)%label = 'g'
   particle(0)%name = 'photon'
   particle(0)%opt_pot_set = .false.
   particle(0)%max_opt_pot = -1
   particle(0)%om_option = -1  
   particle(0)%nume = 0
!--------------   neutron
   particle(1)%Z = 0
   particle(1)%A = 1
   particle(1)%spin = 0.5
   particle(1)%par = 1.0
   particle(1)%ME = 8.071323
   particle(1)%mass = real(1,kind=8)*mass_u + particle(1)%ME
   particle(1)%label = 'n'
   particle(1)%name = 'neutron'
   particle(1)%opt_pot_set = .false.
   particle(1)%max_opt_pot = 3
   particle(1)%om_option = 0  
   particle(1)%nume = 0
!--------------   proton
   particle(2)%Z = 1
   particle(2)%A = 1
   particle(2)%spin = 0.5
   particle(2)%par = 1.0
   particle(2)%ME = 7.288969
   particle(2)%mass = real(1,kind=8)*mass_u + particle(2)%ME
   particle(2)%label = 'p'
   particle(2)%name = 'proton'
   particle(2)%opt_pot_set = .false.
   particle(2)%max_opt_pot = 3
   particle(2)%om_option = 0  
   particle(2)%nume = 0
!--------------   deuteron
   particle(3)%Z = 1
   particle(3)%A = 2
   particle(3)%spin = 1.0
   particle(3)%par = 1.0
   particle(3)%ME = 13.135720
   particle(3)%mass = real(2,kind=8)*mass_u + particle(3)%ME
   particle(3)%label = 'd'
   particle(3)%name = 'deuteron'
   particle(3)%opt_pot_set = .false.
   particle(3)%max_opt_pot = 1
   particle(3)%om_option = 0 
   particle(3)%nume = 0
!--------------   triton
   particle(4)%Z = 1
   particle(4)%A = 3
   particle(4)%spin = 0.5
   particle(4)%par = 1.0
   particle(4)%ME = 14.949794
   particle(4)%mass = real(3,kind=8)*mass_u + particle(4)%ME
   particle(4)%label = 't'
   particle(4)%name = 'triton'
   particle(4)%opt_pot_set = .false.
   particle(4)%max_opt_pot = 1
   particle(4)%om_option = 0 
   particle(4)%nume = 0
!--------------   Helium-3
   particle(5)%Z = 2
   particle(5)%A = 3
   particle(5)%spin = 0.5
   particle(5)%par = 1.0
   particle(5)%ME = 14.931204
   particle(5)%mass = real(3,kind=8)*mass_u + particle(5)%ME
   particle(5)%label = 'h'
   particle(5)%name = 'He3'
   particle(5)%opt_pot_set = .false.
   particle(5)%max_opt_pot = 1
   particle(5)%om_option = 0  
   particle(5)%nume = 0
!--------------   Alpha
   particle(6)%Z = 2
   particle(6)%A = 4
   particle(6)%spin = 0.0
   particle(6)%par = 1.0
   particle(6)%ME = 2.424911
   particle(6)%mass = real(4,kind=8)*mass_u + particle(6)%ME
   particle(6)%label = 'a'
   particle(6)%name = 'alpha'
   particle(6)%opt_pot_set = .false.
   particle(6)%max_opt_pot = 1
   particle(6)%om_option = 0  
   particle(6)%nume = 0
!--------------   Blank - used to define fission
   particle(7)%Z = 0
   particle(7)%A = 0
   particle(7)%spin = 0.0
   particle(7)%par = 0.0
   particle(7)%ME = 0.0d0
   particle(7)%mass = 0.0d0
   particle(7)%label = 'f'
   particle(7)%name = 'fission'
   particle(7)%opt_pot_set = .false.
   particle(7)%max_opt_pot = 1
   particle(7)%om_option = 0  
   particle(7)%nume = 0
   return
end subroutine particle_data
