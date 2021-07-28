!
!*******************************************************************************
!
subroutine run_fresco(ener, fresco_dir, len_fresco, fresco_name, iendf, fname, err_name, symb,      &
                      pindex, mass_proj, iZ, iA, namet, mass_target, beta, deformed, J_gs, K_band,  &
                      V_pot, R_pot, a_pot, RC, iradius, rela, ncc, nex, if_state,                         &
                      cc_state_par, cc_state_type, cc_state_k, cc_state_kpp,                        &
                      cc_state_j, cc_state_e, cc_state_str)
!
!*******************************************************************************
!
!  Discussion:
!
!    This routine is called from make_fresco_tco to run fresco for the
!    case requested, i.e., a particular optical potential at a given
!    incident energy
!
!   Dependencies:
!
!     Modules:
!
!        constants
!        particles_def
!        options
!
!     Subroutines:
!
!        None
!
!     External functions:
!
!        clebr
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
!  Revised:
!      19 July 2021
!  Author:
!      Ian Thompson, LLNL
!  Changes:
!    1. Relativistic specification 'rela' from each optical potential routine
!    2. Now 5 (not 3) columns of potential forms available.
!       Array 'ipotk' specifies Fresco kind of each column.
!    3. Add opt_pot = 4 for the Capote/Soukhovitskii dispersive optical model potential
!
!*******************************************************************************
   use constants
   use particles_def
   use options
   implicit none
   real(kind=8), intent(in) :: ener
   character(len=132), intent(in) :: fresco_dir
   integer(kind=4), intent(in) :: len_fresco
   character(len=100), intent(in) :: fresco_name
   integer(kind=4), intent(in) :: iendf
   character(len=6), intent(in) :: fname
   character(len=100), intent(in) :: err_name
   character(len=2), intent(in) :: symb
   integer(kind=4), intent(in) :: pindex
   real(kind=8), intent(in) :: mass_proj
   integer(kind=4), intent(in) :: iZ, iA
   character(len=8), intent(in) :: namet
   real(kind=8), intent(in) :: mass_target
   real(kind=8), intent(in) :: beta(6)
   logical, intent(in) :: deformed
   real(kind=8), intent(in) :: J_gs
   real(kind=8), intent(in) :: K_band
   real(kind=8), intent(in) :: V_pot(2,3)
   real(kind=8), intent(in) :: R_pot(2,3)
   real(kind=8), intent(in) :: A_pot(2,3)
   real(kind=8), intent(in) :: RC
   integer(kind=4), intent(in) :: iradius
   integer(kind=4), intent(in) :: ncc, nex, if_state
   integer(kind=4), intent(in) :: cc_state_par(nex), cc_state_type(nex)
   integer(kind=4), intent(in) :: cc_state_k(nex), cc_state_kpp(nex)
   real(kind=8), intent(in) :: cc_state_j(nex), cc_state_e(nex)
   real(kind=8), intent(in) :: cc_state_str(nex)
!-------    End of in and out variables
!----------------------------------------------------------------------------------------
   integer(kind=4) :: ipot,ipk,ipotk(5)
   real(kind=8) :: th_min, th_max, th_inc
   integer(kind=4) :: zpart, apart
   real(kind=8) :: A, Z, Ap, AAp
   real(kind=8) :: ac
   real(kind=8) :: absend
   integer(kind=4) :: kp
   integer(kind=8) :: ifirst, ilast, ipot_end
   integer(kind=4) :: i, k
   character(len=2) :: opt_label,rela
   real(kind=8) :: K_state, J_state
   real(kind=8) :: xk_factor
   character(len=15) char_energy
   character(len=132) command, line
   integer(kind=4) :: istart, istop, iend_err
   real(kind=8) :: jtmin, jtmax
   real(kind=8) :: rmatch
   integer(kind=4) :: iter
   character(len=1) :: label
   integer(kind=4) :: len_dir
   real(kind=8) :: zzero
   real(kind=8) :: spin
   real(kind=8) :: mass_rel, e_rel, momentum, wave_number
   real(kind=8) :: hcm, hcm_check
!----------    External functions   -------------------
   real(kind=8) :: clebr

   write(char_energy,'(e15.7)')ener

   len_dir = index(fresco_dir,' ') - 1
   if(len_dir == 0)len_dir = 132

   zzero = 0.0d0
   jtmin = 0
   jtmax = particle(pindex)%lmax
   absend = 0.001d0
   rmatch = 20.0d0

   A = real(iA,kind=8)
   Z = real(iZ,kind=8)
   ac = A**(1.0d0/3.0d0)

   mass_rel = mass_proj*mass_target*mass_u/(mass_proj + mass_target)
   e_rel = ener*mass_target/(mass_target + mass_proj)
   momentum = dsqrt(2.0d0*e_rel*mass_rel)
   wave_number = momentum/hbar_c
   hcm_check = 0.2/wave_number
   hcm = 0.1d0
   if(hcm_check < 0.1d0)then
      hcm = 0.05d0
   elseif(hcm_check < 0.05d0)then
      hcm = 0.01d0
   elseif(hcm_check < 0.01d0)then
      hcm = 0.005d0
   end if

   th_min = 0.0d0
   th_max = 180.0d0
   th_inc = 1.0d0

   label = particle(pindex)%label
   zpart = particle(pindex)%Z
   apart = particle(pindex)%A
   Ap = real(apart,kind=8)
   AAp = 0.0d0
   if(iradius == 1) AAp = Ap
   spin = particle(pindex)%spin
!---------    Coupled channels control
!---------    ncc = total number of coupled channles
!---------    nex = total number of excited states in calculaiton including CCBA states
!---------    if ncc /= nex run one iteration after coupled channels for DWBA states.
  iter = 0
  if(ncc /= nex)iter = 1


!----    Open and create fresco input file for this energy
   open(unit=20,file = fresco_dir(1:len_dir)//'/'//fresco_name(1:len_fresco)//'.in', status='unknown')
!----   write fresco input
   if(particle(pindex)%opt_pot < 10)then
      write(opt_label(1:1),'(i1)')particle(pindex)%opt_pot
      ipot_end = 1
   else
      write(opt_label(1:2),'(i2)')particle(pindex)%opt_pot
      ipot_end = 2
   end if

     write(20,'(a)') fname(1:iendf)//' with potential #'//opt_label(1:ipot_end)//', at E_lab ='//char_energy
     write(20,'(a)') 'NAMELIST'
     write(20,13) hcm, rmatch, rela
13   format(' &Fresco  hcm= ',f6.4, ' rmatch= ',f7.3,' rela="',a2,'"')
     absend = 1.0d-4
     if(ener <= 0.5d0)absend = 1.0d-6
     write(20,'(''    jtmin= '',f5.2,'' jtmax = '',f5.2,'' absend = '',f10.7)')jtmin,jtmax,absend
     write(20,14) th_min, th_inc, th_max, ncc
14      format('    thmin= ',f3.1,' thinc= ',f3.1,' thmax= ',f5.1,' iblock= ',i3)
     write(20,'(''    chans= 1 smats= 2 xstabl= 1 tcfile= 3 iter= '',i2)')iter
     write(20,15) ener
15      format('    elab=',e15.7,'  pel= 1 exl= 1 lab= 1 lin= 1 lex= 1 /')
     write(20,*)
     if(ncc == nex)then
        write(20,16) label, mass_proj, zpart, nex
16         format('&Partition  namep=''',a1,'       '' massp= ',f12.8,' zp= ',i3,' nex=',i3)
     else
        write(20,166) label, mass_proj, zpart, nex
166        format('&Partition  namep=''',a1,'       '' massp= ',f12.8,' zp= ',i3,' nex=',i3,4x,'mixpot=2')
     end if
     write(20,17) namet, mass_target, iZ
17      format('            namet=''',a8,''' masst= ',f12.8,' zt= ',i3,' qval=  0.000/')
     write(20,18)spin,1, cc_state_kpp(1), cc_state_j(1), cc_state_par(1), cc_state_e(1), K_band
18      format('&States jp= ',f3.1,' ptyp=',i2,' ep=  0.000000  cpot=',i3,' jt=',f4.1,' ptyt=',i2,' et=',f8.4,' kkt = ',f4.1,'/')

   do i = 2, ncc
      write(20,21)cc_state_kpp(i), cc_state_j(i), cc_state_par(i), cc_state_e(i),K_band
21    format('&States copyp= 1                       cpot=',i3,' jt=',f4.1,' ptyt=',i2,' et=',f8.4,' kkt = ',f4.1,'/')
   end do
   write(20,*)
   do i = ncc + 1, nex
      write(20,21)cc_state_kpp(i), cc_state_j(i), cc_state_par(i), cc_state_e(i),K_band
   end do

   write(20,'(''&Partition /'')')
   write(20,*)

30 format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:3)=',3f10.5,'/')
31 format('&POT kp=',i3,' type =',i2,' shape=',i2,' p(1:6)=',6f10.5,'/')
32 format('&POT /'/)
!---------    Write data to generate potential
   kp = 1                           !   potential for coupled-channels states
!---------    First print Coulomb
   write(20,30) kp, 0, 0, A, AAp, RC
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!---------   Now loop over indvidual components of nuclear potentials
!---------   ipot = 1, 6,  giving ipk=ipotk(ipot) for Fresco kind
!---------   ipk = 1 : Volume
!---------   ipk = 2 : derivative - surface
!---------   ipk = 3 : Spin-orbit
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ipotk(:) = (/ 1, 2, 3, 1, 3/)
     do ipot = 1, 5
        ipk = ipotk(ipot)
        if(abs(V_pot(1,ipot)) > 1.0d-6)then
           write(20,31) kp,ipk,0,V_pot(1,ipot),R_pot(1,ipot),a_pot(1,ipot),0.,0.,0.
           if(ipk /= 3 .and. deformed )write(20,31)kp,11,ifresco_shape,(beta(k)*ac*R_pot(1,ipot),k=1,6)
        end if
        if(abs(V_pot(2,ipot)) > 1.0d-6)then
           write(20,31) kp,ipk,0,0.,0.,0.,V_pot(2,ipot),R_pot(2,ipot),a_pot(2,ipot)
           if(ipk /= 3 .and. deformed) write(20,31)kp,11,ifresco_shape,(beta(k)*ac*R_pot(2,ipot),k=1,6)
        end if
     end do
     if(ncc /= nex)then
        write(20,*)
        kp = 2
        write(20,30) kp,0,0,A,AAp,RC
        ipot = 1
        ipk = ipotk(ipot)
        write(20,31) kp,ipk,0,V_pot(1,ipot),R_pot(1,ipot),a_pot(1,ipot),V_pot(2,ipot),R_pot(2,ipot),a_pot(2,ipot)
        ipot = 2
        ipk = ipotk(ipot)
        if(abs(V_pot(1,ipot)) + abs(V_pot(2,ipot)) > 1.0d-6)then
         write(20,31) kp,-ipk,0,V_pot(1,ipot),R_pot(1,ipot),a_pot(1,ipot),V_pot(2,ipot),R_pot(2,ipot),a_pot(2,ipot)
         endif
        ipot = 4
        ipk = ipotk(ipot)
        if(abs(V_pot(1,ipot)) + abs(V_pot(2,ipot)) > 1.0d-6)then
         write(20,31) kp,-ipk,0,V_pot(1,ipot),R_pot(1,ipot),a_pot(1,ipot),V_pot(2,ipot),R_pot(2,ipot),a_pot(2,ipot)
         endif
        write(20,31)kp,13,10,0.0,0.5,0.5,0.5,0.0,0.0

 33  format('  &step ib=',i3,1x,'ia= ',i2,1x,'k=',i2,' str=',f10.6,1x,'/')
 34  format('  &step /')
        do i = ncc + 1, nex
           xk_factor = sqrt(2.0d0*J_gs + 1.0d0)
           if(cc_state_type(i) == 0)then
              K_state = real(cc_state_k(i),kind=8)
              J_state = real(cc_state_j(i),kind=8)
              xk_factor = sqrt(2.0d0*J_gs + 1.0d0)*clebr(J_gs,K_band,K_state,zzero,J_state,K_band)
           end if
           write(20,33)i,if_state,cc_state_k(i),cc_state_str(i)*ac*R_pot(1,1)*xk_factor*cc_scale
        end do
        write(20,34)
        ipot = 3
        ipk = ipotk(ipot)
        write(20,31) kp,ipk,0,V_pot(1,ipot),R_pot(1,ipot),a_pot(1,ipot),V_pot(2,ipot),R_pot(2,ipot),a_pot(2,ipot)
     end if
     write(20,*)
     write(20,32)
     write(20,'(''&Overlap /'')')
     write(20,'(''&Coupling /'')')

   close(unit=20)

!-------     Now run fresco
   write(6,*)'**********************************************'
   write(6,*)'Calling unix system command to execute fresco'

   line(1:132) = ' '
   ifirst = 1
   ilast = 12
   line(ifirst:ilast) = 'Calculating '
   ifirst = 13
   ilast = 13
   line(ifirst:ilast) = particle(pindex)%label
   ifirst = 14
   ilast = 16
   line(ifirst:ilast) = ' + '
   ifirst = 17
   if(iA < 10)then
      ilast = 18
      write(line(ifirst:ilast),'(i1)')iA
   elseif(iA < 100)then
      ilast = 19
      write(line(ifirst:ilast),'(i2)')iA
   elseif(iA < 1000)then
      ilast = 20
      write(line(ifirst:ilast),'(i3)')iA
   end if
   ifirst = ilast + 1
   if(symb(1:1).ne.' ')then
      ilast = ifirst + 1
      line(ifirst:ilast) = symb
   else
      ilast = ifirst
      line(ifirst:ilast) = symb(2:2)
   end if


   write(6,'(a)')line(1:ilast)//', at E_lab ='//char_energy

   iend_err = index(err_name,' ') - 1

   command(1:132) = ' '
   istart = 1
   istop = 10
   command(istart:istop) = 'frescox < '
   istart = istop + 1
   istop = istop + len_dir + 1 + len_fresco + 3
   command(istart:istop) = fresco_dir(1:len_dir)//'/'//fresco_name(1:len_fresco)//'.in'
   istart = istop + 1
   istop = istart + 2
   command(istart:istop) = ' > '
   istart = istop + 1
   istop = istart + len_dir + 1 + len_fresco + 4
   command(istart:istop) = fresco_dir(1:len_dir)//'/'//fresco_name(1:len_fresco)//'.out'
   istart = istop + 1
   istop = istart + 4
   command(istart:istop) = ' 2>> '
   istart = istop + 1
   istop = istart + len_dir + 1 + iend_err - 1
   command(istart:istop) = fresco_dir(1:len_dir)//'/'//err_name(1:iend_err)

   write(6,'(a)')command(1:istop)
   call system(command(1:istop))
   return

end subroutine run_fresco
