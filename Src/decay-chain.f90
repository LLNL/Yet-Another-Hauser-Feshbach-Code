subroutine set_up_decay_chain(Z_p, A_p, Z_t, A_t, num_comp)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine to sets up the entire decay chain for the calculation
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
   use nodeinfo
   use options
   use print_control
   use useful_data
   use nuclei
   use Channel_info
   use particles_def
   use directory_structure
   use constants
   implicit none
   integer(kind=4), intent(in) :: Z_p, A_p, Z_t, A_t
   integer(kind=4), intent(out) :: num_comp
!-------------------------------------------------------
   integer(kind=4) :: inuc, inucp
   integer(kind=4) :: Z_i,A_i
   integer(kind=4) :: i, k
   integer(kind=4) :: ia,ih,it,id,ip,in
   integer(kind=4) :: Z_f,N_f, A_f
   integer(kind=4) :: Z_pp, N_pp, A_pp
   real(kind=8) :: D_p, D_n
   real(kind=8) :: me,be,sep(6),me_f,be_f,sep_f(6)
   real(kind=8) :: sep_tot
   real(kind=8) :: emax
   integer(kind=4) ::  N_i,d_i,alpha_i
   integer(kind=4) ::  NZ
   real(kind=8) :: e_rel
   integer, parameter :: max_nuc = 200
   integer(kind=4) :: storeZA(2,max_nuc)                 !  Temporary storage array for compound nuclei
   real(kind=8) ::  store_exmax(max_nuc)                 !  Temporary storage array for compound nuclei
   real(kind=8) :: exmax
   real(kind=8) :: em_proj
   real(kind=8) :: Coul
   logical found
   integer(kind=4) :: npart
   integer(kind=4) :: num
   integer(kind=4) :: num_temp
   integer(kind=4) :: num_paths
   integer(kind=4) :: max_num
   integer(kind=4) :: num_decay
   integer(kind=4) :: nump(6)
   integer(kind=4) :: ichannel
   real(kind=8) :: max_e
!-------------------------------------------------------
   Z_i=Z_p+Z_t
   A_i=A_p+A_t
 
   N_i = A_i - Z_i
   NZ = min(Z_i,N_i)
   d_i = NZ
   alpha_i = NZ/2

   if(max_particle(1) < 0)max_particle(1) = N_i - 2
   if(max_particle(2) < 0)max_particle(2) = Z_i - 2
   if(max_particle(3) < 0)max_particle(3) = d_i - 2
   if(max_particle(6) < 0)max_particle(6) = alpha_i - 2

   Coulomb_Barrier(1:6) = 0.0d0

   do k = 1, 6
      Coulomb_Barrier(k) = 0.5*1.44*(Z_i-particle(k)%Z)*particle(k)%Z/               &
                           (1.3*((A_i-particle(k)%A)**(1./3.) + particle(k)%A**(1./3.)))
   end do
      

   call get_binding_energy(data_path,len_path,              &
                           Z_i,A_i,me,be, sep)

   if(projectile%particle_type >= 0 .and. projectile%particle_type <= 6)then
      em_proj = projectile%e_max
      if(em_proj < de) em_proj = de/2.
      e_rel = em_proj*dfloat(A_t)/dfloat(A_t+A_p)
      emax = e_rel + sep(projectile%particle_type)
   else
      emax = 0.0
      do i=1, num_pop_e
         if(Pop_data(i)%Ex_pop + 3.5*Pop_data(i)%dEx_pop > emax) then
            emax = Pop_data(i)%Ex_pop + 3.0*Pop_data(i)%dEx_pop
         end if
      end do
   end if

!----------------------------------------------------------------!
!--------------   Find out how many compound nuclei there are    !
!--------------   And number of possible exit channels           !
!----------------------------------------------------------------!
   num_paths=0
   num_temp=0
   max_num = 0

   particle(0:6)%in_decay = .false.

   num_comp = 0
   
!-------   Loop over particle types
   do ia = 0, max_particle(6)
      do ih = 0, max_particle(5)
         do it = 0, max_particle(4)
            do id = 0, max_particle(3)
               do ip = 0, max_particle(2)
                  do in = 0, max_particle(1)
                     Z_f = Z_i - ip - id - it - ih*2 - ia*2
                     A_f = A_i - in - ip - id*2 - it*3 - ih*3 - ia*4
                     npart = in + ip + id + it + ih + ia
                     if(Z_f < 2 .or. A_f < 4)exit
!-----   Get Binding energy and compute total mass excess
!-----   See if channel is open
                     call get_binding_energy(data_path,len_path,      &
                                             Z_f,A_f,me_f,be_f,sep_f)
                     if(me_f <= -1.01d6)exit
                     sep_tot= me_f - me + ia*particle(6)%ME +    &
                                          ih*particle(5)%ME +    &
                                          it*particle(4)%ME +    &
                                          id*particle(3)%ME +    &
                                          ip*particle(2)%ME +    &
                                          in*particle(1)%ME
                     Coul = 0.0d0
                     if(Apply_Coulomb_Barrier)                   &
                        Coul = in*Coulomb_Barrier(1) +           &
                               ip*Coulomb_Barrier(2) +           &
                               id*Coulomb_Barrier(3) +           &
                               it*Coulomb_Barrier(4) +           &
                               ih*Coulomb_Barrier(5) +           &
                               ia*Coulomb_Barrier(6)
                     if(emax < sep_tot + Coul)cycle

                     exmax = emax - sep_tot
                     found = .false.
!------   Check temporary storage   -------
                     do i = 1, num_comp
                        if(Z_f == storeZA(1,i) .and. A_f == storeZA(2,i))then
                           if(exmax > store_exmax(i))store_exmax(i) = exmax
                           found=.true.
                           exit
                        end if
                     end do
!------   Not found in temporary storage, so put into list
                     if(.not. found)then
                        num_comp = num_comp + 1
                        if(num_comp > max_nuc) stop "Error too many compound nuclei in set_up_decay_chain"
                        storeZA(1,num_comp) = Z_f
                        storeZA(2,num_comp) = A_f
                        store_exmax(num_comp) = exmax
                     end if
                     num_paths = num_paths + nint(factorial(npart)/(factorial(in)*factorial(ip)*factorial(id)*    &
                                 factorial(it)*factorial(ih)*factorial(ia)))
                     num = in + ip + id + it + ih + ia
                     if(num > max_num) max_num = num
                  end do
               end do
            end do
         end do
      end do
   end do

   if(iproc == 0)write(6,*)'Total number of compound nuclei',num_comp
   if(iproc == 0)write(6,*)'Total number of decay paths',num_paths

!--------------------------------------!
!------    Allocate compound nuclei    !
!--------------------------------------!
   allocate(nucleus(num_comp))
   num_channels = num_paths
!--------------------------------------!
!------    Allocate Exit channles      !
!--------------------------------------!
   if(.not.allocated(Exit_channel))then
      allocate(Exit_channel(num_channels))
   else
      stop 'channels already allocated'
   end if

!--------------------------------------!
!-----    Set up compound nuclei       !
!-----    use same algorithm as above  !
!--------------------------------------!
   inuc = 0
   max_num = 0
   ichannel = 0
   do ia = 0, max_particle(6)
      nump(6) = ia
      do ih = 0, max_particle(5)
         nump(5) = ih
         do it = 0, max_particle(4)
            nump(4) = it
            do id = 0, max_particle(3)
               nump(3) = id
               do ip = 0, max_particle(2)
                  nump(2) = ip
                  do in = 0, max_particle(1)
                     nump(1) = in
                     Z_f = Z_i - ip - id - it - ih*2 - ia*2
                     A_f = A_i - in - ip - id*2 - it*3 - ih*3 - ia*4
                     npart = in + ip + id + it + ih + ia
                     if(Z_f < 2 .or. A_f < 4)exit
                     call get_binding_energy(data_path,len_path,       &
                                             Z_f,A_f,me_f,be_f,sep_f)
                     if(me_f <= -1.01d6)exit
                     sep_tot = me_f - me + ia*particle(6)%ME +    &
                                           ih*particle(5)%ME +    &
                                           it*particle(4)%ME +    &
                                           id*particle(3)%ME +    &
                                           ip*particle(2)%ME +    &
                                           in*particle(1)%ME
                     Coul = 0.0d0
                     if(Apply_Coulomb_Barrier)                   &
                        Coul = in*Coulomb_Barrier(1) +           &
                               ip*Coulomb_Barrier(2) +           &
                               id*Coulomb_Barrier(3) +           &
                               it*Coulomb_Barrier(4) +           &
                               ih*Coulomb_Barrier(5) +           &
                               ia*Coulomb_Barrier(6)
                     if(emax < sep_tot + Coul)cycle

                     exmax = emax - sep_tot
                     found = .false.
!----   Check if already made
                     do i = 1, inuc
                        if(Z_f == nucleus(i)%Z .and. A_f == nucleus(i)%A)then
                           if(exmax > store_exmax(i))store_exmax(i) = exmax
                           found=.true.
!----   Found, but check if exmax is greater than current, if so, use this instead
                           if(exmax > nucleus(i)%Ex_max)nucleus(i)%Ex_max = exmax
                           exit
                        end if
                     end do
!----   Not made previously, so add to list and fill data arrays
                     if(.not.found)then
                        inuc = inuc + 1
                        if(inuc > num_comp) stop "inuc > num_comp"
                        nucleus(inuc)%Z = Z_f
                        nucleus(inuc)%A = A_f
                        nucleus(inuc)%D0exp=-1.0d0
                        nucleus(inuc)%dD0exp=-1.0d0
                        nucleus(inuc)%D1exp=-1.0d0
                        nucleus(inuc)%dD1exp=-1.0d0
                        nucleus(inuc)%Gamma_g = -1.0d0
                        nucleus(inuc)%Gamma_g_exp = -1.0d0
                        nucleus(inuc)%dGamma_g_exp = -1.0d0
                        nucleus(inuc)%Gamma_g_1 = -1.0d0
                        nucleus(inuc)%Gamma_g_1_exp = -1.0d0
                        nucleus(inuc)%dGamma_g_1_exp = -1.0d0
                        nucleus(inuc)%beta(2) = 0.0d0
                        nucleus(inuc)%beta(3) = 0.0d0
                        nucleus(inuc)%beta(4) = 0.0d0
                        nucleus(inuc)%beta(5) = 0.0d0
                        nucleus(inuc)%beta(6) = 0.0d0
			

                        if(A_f > 20)then
                           nucleus(inuc)%lev_option = 1
                           if(A_f > 130) nucleus(inuc)%lev_option = 2
                           nucleus(inuc)%fit_aparam=.false.   !   When fitting to D0 do we adjust a-parameter
                                                              !   or shell correction 
                        else
                           nucleus(inuc)%lev_option = 0
                           nucleus(inuc)%fit_aparam=.true.    !   When fitting to D0 do we adjust a-parameter
                                                              !   or shell correction 
                        end if
                        nucleus(inuc)%pair_model = 1
                        nucleus(inuc)%level_param(1:11)=0.0d0
                        nucleus(inuc)%fit_D0=.true.       !   Fit to D0 (if known)  
                        nucleus(inuc)%fit_ematch=.true.   !   Fit ematch to cummlative level density
                                                       !----------   Now set up connections in the primary array nucleus so that
                                                       !----------   the HF denominators can be calculated
                                                       !             to overide set to false with 
                                                       !             option lev_fit_ematch for each nucleus
                                                       !             or globally with fit_ematch; 0 for false, 1 for true
                        nucleus(inuc)%fission_read=.false.       !  Set true once fission parameters from default are read  
                        nucleus(inuc)%atomic_symbol=symb(Z_f)
                        nucleus(inuc)%BE = be_f
                        nucleus(inuc)%ME = me_f
                        nucleus(inuc)%Mass = real(A_f, kind=8)*mass_u + me_f
                        nucleus(inuc)%Ex_max = exmax
                        nucleus(inuc)%Kinetic_energy = 0.0
                        nucleus(inuc)%dKinetic_energy = 0.0
                        nucleus(inuc)%nump(1:6) = 0
                        nucleus(inuc)%num_part_decay = 0
                        if(iand(A_f,1) == 0)then
                           nucleus(inuc)%jshift=0.0                       !   j-shift for level densities
                        else                                              !   =0.5 for odd-A nuclei
                           nucleus(inuc)%jshift=0.5                       !   =0 for even-A nuclei
                        end if
                        nucleus(inuc)%sep_e(0:6)=0.0d0
                        do i=1,6
                           nucleus(inuc)%sep_e(i)=sep_f(i)            !  separation energy neutron to alpha
                        end do
                     end if

                     num = in + ip + id + it + ih + ia
                     if(num > max_num) max_num = num
!----  Set up all the channels for this final compound nucleus - note ichannel is incremented in make_channels
                     call make_channels(num, nump, inuc, ichannel)
                  end do
               end do
            end do
         end do
      end do
   end do



   do inuc = 1, num_comp
!-----   First allow gammas - always allowed
      nucleus(inuc)%num_decay = 1
!-----   Start by counting decays from this nucleus
      
      do k = 1, 6
         if(max_particle(k) == 0)cycle                    ! If max_particle == 0 then decay of this particl is not allowed
         Coul = Coulomb_Barrier(k)
         if(nucleus(inuc)%Ex_max > nucleus(inuc)%sep_e(k) + Coul)then        !   Check separation energy
           Z_f = nucleus(inuc)%Z - particle(k)%Z
           A_f = nucleus(inuc)%A - particle(k)%A
           N_f = A_f - Z_f
           found = .false.                                            !  search if decay is possible to listed compound nuclei
           do inucp = 1, num_comp
              if((Z_f == nucleus(inucp)%Z) .and. (A_f == nucleus(inucp)%A))then
                 found = .true.
                 exit
              end if 
           end do
           if(found)nucleus(inuc)%num_decay = nucleus(inuc)%num_decay + 1
         end if
      end do
!-----    Now allocate arrays
      if(nucleus(inuc)%num_decay > 0)then
         allocate(nucleus(inuc)%decay_to(nucleus(inuc)%num_decay))
         allocate(nucleus(inuc)%decay_particle(nucleus(inuc)%num_decay))
      end if
!------    Fill arrays for decays: particle type and to which nucleus for each decay
!-----   First allow gammas - always allowed
      num_decay = 1
      k = 0
      nucleus(inuc)%decay_particle(num_decay) = k           !   particle type
      nucleus(inuc)%decay_to(num_decay) = inuc             !  decay to which nucleus
      particle(k)%in_decay = .true.
      do k = 1, 6
         if(max_particle(k) == 0)cycle                    ! If max_particle == 0 then decay of this particl is not allowed
         Coul = Coulomb_Barrier(k)
         if(nucleus(inuc)%Ex_max > nucleus(inuc)%sep_e(k) + Coul)then        !   Check separation energy
            Z_f = nucleus(inuc)%Z - particle(k)%Z
            A_f = nucleus(inuc)%A - particle(k)%A
             found = .false.                                          !  search if decay is possible to listed compound nuclei
            do inucp = 1, num_comp
               if((Z_f == nucleus(inucp)%Z) .and. (A_f == nucleus(inucp)%A))then
                  found = .true.
                  exit
               end if 
            end do
            if(found)then
               num_decay = num_decay + 1
               nucleus(inuc)%decay_particle(num_decay) = k           !   particle type
               nucleus(inuc)%decay_to(num_decay) = inucp             !  decay to which nucleus
               particle(k)%in_decay = .true.
            end if
         end if
      end do
   end do

!------------   Gammas are always in the decay chain
   particle(0:6)%max_e = 0.0
   do i=1,num_comp
      do k=0,6
         max_e=nucleus(i)%Ex_max - nucleus(i)%sep_e(k)
         if(max_e > particle(k)%max_e)                         &
            particle(k)%max_e = max_e
      end do
   end do

!-------   Calculate proton and neutron pairing energies from binding energies
!---------   Neutron pairing energy
   do inuc = 1, num_comp
      Z_f = nucleus(inuc)%Z
      A_f = nucleus(inuc)%A
      N_f = A_f - Z_f

      D_n = 0.0d0
      Z_pp = Z_f
      N_pp = N_f - 2
      A_pp = Z_pp + N_pp
      call get_binding_energy(data_path,len_path,   &
                              Z_pp,A_pp,me_f,be_f,sep_f)
      D_n = be_f
      N_pp = N_f - 1
      A_pp = Z_pp + N_pp
      call get_binding_energy(data_path,len_path,   &
                              Z_pp,A_pp,me_f,be_f,sep_f)
      D_n = D_n -3.0d0*be_f
      N_pp = N_f
      A_pp = Z_pp + N_pp
      call get_binding_energy(data_path,len_path,   &
                              Z_pp,A_pp,me_f,be_f,sep_f)
      D_n = D_n +3.0d0*be_f
      N_pp = N_f + 1
      A_pp = Z_pp + N_pp
      call get_binding_energy(data_path,len_path,   &
                              Z_pp,A_pp,me_f,be_f,sep_f)
      D_n = D_n - be_f
      D_n = 0.25d0*D_n
      if(iand(N_f,1) == 1)D_n = -1.0d0*D_n
!-------   Proton pairing energy
      D_p = 0.0d0
      N_pp = N_f
      Z_pp = Z_f - 2
      A_pp = Z_pp + N_pp
      call get_binding_energy(data_path,len_path,   &
                              Z_pp,A_pp,me_f,be_f,sep_f)
      D_p = be_f
      Z_pp = Z_f - 1
      A_pp = Z_pp + N_pp
      call get_binding_energy(data_path,len_path,   &
                              Z_pp,A_pp,me_f,be_f,sep_f)
      D_p = D_p - 3.0d0*be_f
      Z_pp = Z_f
      A_pp = Z_pp + N_pp
      call get_binding_energy(data_path,len_path,   &
                              Z_pp,A_pp,me_f,be_f,sep_f)
      D_p = D_p + 3.0d0*be_f
      Z_pp = Z_f + 1
      A_pp = Z_pp + N_pp
      call get_binding_energy(data_path,len_path,   &
                              Z_pp,A_pp,me_f,be_f,sep_f)
      D_p = D_p - be_f
      D_p = 0.25d0*D_p
      if(iand(Z_f,1) == 1)D_p = -1.0d0*D_p
      nucleus(inuc)%Delta_exp_p = D_p
      nucleus(inuc)%Delta_exp_n = D_n
      
      if(iand(Z_f,1) == 0 .and. iand(N_f,1) == 0)nucleus(inuc)%Delta_exp = (D_p + D_n)
      if(iand(Z_f,1) == 1 .and. iand(N_f,1) == 0 .or.             &
         iand(Z_f,1) == 0 .and. iand(N_f,1) == 1)nucleus(inuc)%Delta_exp = 0.5d0*(D_p + D_n)
     if(iand(Z_f,1) == 1 .and. iand(N_f,1) == 1)nucleus(inuc)%Delta_exp = 0.0d0


   end do




    return
!----------   Now set up connections in the primary array nucleus so that
!----------   the HF denominators can be calculated
end subroutine set_up_decay_chain
!
!*******************************************************************************
!
subroutine make_channels(num_part_tot, num_part, inuc, ichannel)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine defines the channels and gives the appropriate label
!    based on emitted particles
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
   use Channel_info
   use particles_def
   implicit none

   integer(kind=4), intent(in) :: num_part_tot
   integer(kind=4), intent(in) :: num_part(6)
   integer(kind=4), intent(in) :: inuc
   integer(kind=4), intent(out) :: ichannel
!---------------------------------------------------
   integer(kind=4) :: sum(6)
   logical too_many
   integer(kind=4) :: k, n, nloops, low(16), hi(16), idex(16)

   integer(kind=4) :: channel_code

   nloops = num_part_tot 

   low(1:nloops) = 1
   hi(1:nloops) = 6

   if(nloops == 0)then
      ichannel = ichannel + 1
      channel_code = 0
      Exit_channel(ichannel)%Channel_code = channel_code
      Exit_channel(ichannel)%num_particles = num_part_tot
      Exit_channel(ichannel)%Channel_Label(1:20) = ' '
      Exit_channel(ichannel)%Channel_Label(1:1) = 'g'   
      do n = 1, nloops
         Exit_channel(ichannel)%decay_particles(n) = idex(n)
         Exit_Channel(ichannel)%Channel_Label(n:n) = particle(idex(n))%label            
      end do
      Exit_channel(ichannel)%Final_nucleus = inuc
      return
   end if


   do n = 1, nloops
      idex(n) = min(low(n), hi(n))
      if(idex(n) < low(n)) return       !      didn't do anything
   end do


   do while(idex(1) <= hi(1))

      sum(1:6) = 0
      do n = 1, nloops
         sum(idex(n)) = sum(idex(n)) + 1
      end do
      too_many = .false.
      do k = 1, 6
         if(sum(k) > num_part(k))too_many = .true.
      end do

      if(.not. too_many)then
         ichannel = ichannel + 1
         channel_code = 0
         do n = 1, nloops
            channel_code = ior(channel_code,ishft(idex(n),(n-1)*3))
         end do
         Exit_channel(ichannel)%Channel_code = channel_code
         Exit_channel(ichannel)%num_particles = num_part_tot
         Exit_channel(ichannel)%Channel_Label(1:20) = ' '
         allocate(Exit_channel(ichannel)%decay_particles(nloops))
         do n = 1, nloops
            Exit_channel(ichannel)%decay_particles(n) = idex(n)
            Exit_Channel(ichannel)%Channel_Label(n:n) = particle(idex(n))%label
         end do
         Exit_channel(ichannel)%Final_nucleus = inuc
      end if        
      idex(nloops) = idex(nloops) + 1
      do n = nloops, 2, -1
         if(idex(n) > hi(n))then
            idex(n-1) = idex(n-1) + 1
            idex(n) = low(n)
         end if
      end do
   end do

   return
end subroutine make_channels
