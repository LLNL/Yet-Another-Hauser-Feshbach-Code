!
!*******************************************************************************
!
subroutine KD_potential(part_type, iA, iZ, energy, V_pot, R_pot, a_pot, RC, D3, iradius)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates the parameters of the 
!    Koning & Delaroche optical potential
!
!  Reference: 
!  
!     A.J. Koning and J.P. Delaroche  
!     Local and global optical models from 1 keV to 200 MeV
!     Nuclear Physics A 713 (2003) 310
!          
!  Licensing:
!
!    This code is distributed under the GNU LGPL version 2 license. 
!
!  Date:
!     Written: 14 October 2004 by Jutta Escher 
!     Modified: 18 Dec 2008 by Ian Thompson
!     Modified: 25 September 2019 by Erich Ormand
!----    Modified 29 November 2017 by WEO to work with YAHFC front end                                +
!----    converted to FORTRAN90 format                                     +
!----    changed to implicit none with explict statement for all variables +
!----    converted arithmetic to double                                    +
!----    converted NA, NN, and NZ to real variables, used in formulae      +
!----    especially previous use of NA**XX                                 +
!
!*******************************************************************************
!
   use nodeinfo
   use variable_kinds
   implicit none
   integer(kind=4), intent(in) :: part_type, iA, iZ
   real(kind=8), intent(in) :: energy
   real(kind=8), intent(out) :: V_pot(2,3)
   real(kind=8), intent(out) :: R_pot(2,3)
   real(kind=8), intent(out) :: a_pot(2,3)
   real(kind=8), intent(out) :: rc, d3
   integer(kind=4), intent(out) :: iradius
!
!    (1,1)   Real Volume
!    (2,1)   Imaginary Volume
!    (1,2)   Real Surface
!    (2,2)   Imaginary Surface
!    (1,3)   Real Spin-orbit
!    (2,3)   Imaginary Spin-orbit
!
!------------------------------------------------------------------
!   real(kind=8) :: vv, rv, av, wv, rw, aw, wd, rd, ad
!   real(kind=8) :: vso, rso, aso, wso, wrso, waso
!------------------------------------------------------------------
   integer(kind=4) :: iN
   real(kind=8) :: diff
   real(kind=8) :: xA, xZ, xN
   real(kind=8) :: v1, v2, v3, v4, w1, w2, d1, d2
   real(kind=8) :: vso_1, vso_2, wso_1, wso_2
   real(kind=8) :: eFermi, vc
   real(kind=8) :: del_vc

!     iradius = 0    R = r0*A_T^1/3
!     iradius = 1    R = r0*(A_Target^1/3 + A_proj^1/3)
   iradius = 0

   iN = iA - iZ
   xA = real(iA,kind=8)
   xZ = real(iZ,kind=8)
   xN = real(iN,kind=8)

   diff = (xN - xZ)/xA

! part_type = 1 (neutrons)
   if(part_type > 2)then
      if(iproc == 0)write(6,*)'Error in KDParam, use only for protons and neutrons, not k = '
      call MPI_Abort(icomm,301,ierr)
   end if
   v1 = 59.30d0 - 21.0d0*diff - 2.4d-2*xA
   v2 = 7.228d-3 - 1.48d-6*xA
   v3 = 1.994d-5 - 2.0d-8*xA
   v4 = 7.1d-9
   w1 = 12.195d0 + 1.67d-2*xA
   W2 = 73.55d0 + 7.95d-2*xA
   d1 = 16.0d0 - 16.0d0*diff
   d2 = 1.80d-2 + 3.802d-3/(1.0d0 + EXP((xA - 156.0d0)/8.0d0))
   d3 = 11.5
   vso_1 = 5.922d0 + 3.0d-3*xA
   vso_2 = 4.0d-3
   wso_1 = -3.1d0
   wso_2 = 160.0d0
   eFermi = -11.2814d0 + 2.646d-2*xA
   rc = 0.0d0
   rc = 1.198d0 + 0.697d0/(xA**(0.6666666666d0)) + 12.994d0/(xA**(1.6666666666d0))
   vc = 0.0d0
! part_type = 2 (protons)
   if(part_type == 2) then
      v1 = 59.30d0 + 21.0d0*diff - 2.4d-2*xA
      v2 = 7.067d-3 + 4.23d-6*xA
      v3 = 1.729d-5 - 1.136d-8*xA
      v4 = 7.1d-9
      w1 = 14.667d0 + 9.629d-3*xA
      W2 = 73.55d0 + 7.95d-2*xA
      d1 = 16.0d0 + 16.0d0*diff
      d2 = 1.80d-2 + 3.802d-3/(1.0d0 + exp((xA - 156.0d0)/8.0d0))
      d3 = 11.5d0
      vso_1 = 5.922d0 + 3.0d-3*xA
      vso_2 = 4.0d-3
      wso_1 = -3.1d0
      wso_2 = 160.0d0
      eFermi = -8.4075d0 + 1.378d-2*xA
      rc = 1.198d0 + 0.697d0/(xA**(0.6666666666d0)) + 12.994d0/(xA**(1.6666666666d0))
      vc = 1.73d0*xZ/(rc*(xA**0.333333333d0))
   end if
   del_vc = vc*v1*(v2 - 2.0d0*V3*(energy - eFermi) + 3.0d0*V4*(energy - eFermi)**2)

   V_pot(1:2,1:3) = 0.0d0

   V_pot(1,1) = v1*(1.0d0 - v2*(energy - eFermi) + v3*(energy - eFermi)**2 - &
                v4*(energy - eFermi)**3) + del_vc
   R_pot(1,1) = 1.3039d0 - 0.4054d0/(xA**0.333333333d0)
   a_pot(1,1) = 0.6778d0 - 1.487d-4*xA

   V_pot(2,1) = w1*(energy - eFermi)**2/((energy - eFermi)**2 + W2**2)
   R_pot(2,1) = R_pot(1,1)
   a_pot(2,1) = a_pot(1,1)

   V_pot(2,2) = d1*(energy - eFermi)**2*exp(-d2*(energy - eFermi))/((energy - eFermi)**2 + d3**2)
   R_pot(2,2) = 1.3424d0 - 1.585d-2*(xA**0.333333333d0)
   if(part_type == 1) a_pot(2,2) = 0.5446d0 - 1.656d-4*xA     
   if(part_type == 2) a_pot(2,2) = 0.5187d0 - 5.206d-4*xA

   V_pot(1,2) = 0.0d0
   R_pot(1,2) = 1.3424d0 - 1.585d-2*(xA**0.333333333d0)
   if(part_type == 1) a_pot(1,2) = 0.5446d0 - 1.656d-4*xA     
   if(part_type == 2) a_pot(1,2) = 0.5187d0 - 5.206d-4*xA

   V_pot(1,3) = vso_1*exp(-vso_2*(energy - eFermi))
   R_pot(1,3) = 1.1854d0 - 0.647d0/(xA**0.333333333d0)      
   a_pot(1,3) = 0.59d0

   V_pot(2,3) = wso_1*(energy - eFermi)**2/((energy - eFermi)**2 + wso_2**2)
   R_pot(2,3) = R_pot(1,3)
   a_pot(2,3) = a_pot(1,3)

   return
end subroutine KD_potential
!
!*******************************************************************************
!
subroutine maslov_03_potential(E, V_pot, R_pot, a_pot, RC, iradius)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates the parameters for the maslov potential
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
   use variable_kinds
   implicit none
   real(kind=8), intent(in) :: E
   real(kind=8), intent(out) :: V_pot(2,3)
   real(kind=8), intent(out) :: R_pot(2,3)
   real(kind=8), intent(out) :: a_pot(2,3)
   real(kind=8), intent(out) :: RC
   integer(kind=4), intent(out) :: iradius
!----------------------------------------------------------
!   real(kind=8) :: VR, RR, AR
!   real(kind=8) :: W, RW, AW
!   real(kind=8) :: VD, RVD, AVD
!   real(kind=8) :: WD, RD, AD
!   real(kind=8) :: VSO, RSO, ASO
!   real(kind=8) :: WSO, WRSO, WASO

!     iradius = 0    R = r0*A_T^1/3
!     iradius = 1    R = r0*(A_Target^1/3 + A_proj^1/3)
   iradius = 0

   V_pot(1:2,1:3) = 0.0d0

   V_pot(1,1) = 45.93d0 - 0.28d0*E + 5.73d-4*E*E
   R_pot(1,1) = 1.26d0
   a_pot(1,1) = 0.63d0

   if(E < 8.0d0) then
      V_pot(2,2) = 3.14d0 + 0.436d0*E
   else
      V_pot(2,2) = 6.628d0
   end if
   R_pot(2,2) = 1.26d0
   a_pot(2,2) = 0.52d0

   V_pot(2,1) = 0.0d0
   R_pot(2,1) = R_pot(1,1)
   a_pot(2,1) = a_pot(1,1)

   RC = 0.0d0

   V_pot(1,3) = 6.2d0
   R_pot(1,3) = 1.120d0
   a_pot(1,3) = 0.47d0

   V_pot(2,3) = 0.0d0
   R_pot(2,3) = R_pot(1,3)
   a_pot(2,3) = a_pot(1,3)

   V_pot(1,2) = 0.0d0
   R_pot(1,2) = R_pot(1,1)
   a_pot(1,2) = a_pot(1,1)

   return
end subroutine maslov_03_potential
!
!*******************************************************************************
!
subroutine soukhovitskii_potential(part_type, iA, iZ, OM_option,        &
                                   energy, V_pot, R_pot, a_pot, RC, iradius)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates the parameters for the Soukhovitskii optical
!    model potential
!
!  Reference:
!    Soukhovitskii et al, J. Phys. G30, p. 905 (2004)
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
   use variable_kinds
   use directory_structure
   use useful_data
   use options
   implicit none
!--------------------------------------------------------------------
   integer(kind=4) :: part_type, iA, iZ, OM_option
   real(kind=8), intent(in) :: energy
   real(kind=8), intent(out) :: V_pot(2,3)
   real(kind=8), intent(out) :: R_pot(2,3)
   real(kind=8), intent(out) :: a_pot(2,3)
   real(kind=8), intent(out) :: RC
   integer(kind=4), intent(out) :: iradius
!--------------------------------------------------------------------
   real(kind=8) :: v, rv, av, vd, rvd, avd
   real(kind=8) :: w, rw, aw, wd, rwd, awd
   real(kind=8) :: vso, rvso, avso, wso, rwso, awso
!--------------------------------------------------------------------
   integer(kind=4) :: iN
   real(kind=8) :: xA, xZ, xN
   real(kind=8) :: asym, efermi, ef
   real(kind=8) :: Cviso, V0r, Var, vrdisp, v1r, v2r, lambdaR
   real(kind=8) :: viso, Ccoul, phicoul, rr, Crr, widr
   real(kind=8) :: w1loc, w2loc, wddisp
   real(kind=8) :: Cwiso, Wad, d1loc, d2loc, d3loc
   real(kind=8) :: vso1loc, vso2loc, wso1loc, wso2loc
   real(kind=8) :: onethird
   real(kind=8) :: me, be, sep(0:6)
   real(kind=8) :: phase

!     iradius = 0    R = r0*A_T^1/3
!     iradius = 1    R = r0*(A_Target^1/3 + A_proj^1/3)
   iradius = 0

   if(part_type > 2)then
      write(6,*)'Error in soukhovitskii_potential: part_type > 2'
      call MPI_Abort(icomm,301,ierr)
   end if
   onethird = 1.0d0/3.0d0
   iN = iA - iZ
   xA = real(iA,kind=8)
   xZ = real(iZ,kind=8)
   xN = real(iN,kind=8)

   asym = (xA - 2.0d0*xZ)/xA
   phase = 1.0d0
   eFermi = 0.0d0
   if(part_type == 1)then
      if(OM_option == 0)then
         call get_binding_energy(data_path, len_path,     &
                                 iZ, IA, me, be, sep)
         eFermi = -0.5d0*sep(1)
         call get_binding_energy(data_path, len_path,     &
                                 iZ, IA+1, me, be, sep)
         eFermi = eFermi - 0.5d0*sep(1)
      elseif(OM_option > 0)then
         eFermi = -11.2814d0 + 0.02646d0*xA
      end if
      phase = -1.0d0
   elseif(part_type == 2)then
      if(OM_option == 0)then
         call get_binding_energy(data_path, len_path,     &
                                 iZ, IA, me, be, sep)
         eFermi = -0.5d0*sep(2)
         call get_binding_energy(data_path, len_path,     &
                                 iZ+1, IA+1, me, be, sep)
         eFermi = eFermi - 0.5d0*sep(2)
      elseif(OM_option > 0)then
         eFermi = -8.4075d0 + 0.01378d0*xA
      end if
      phase = 1.0d0
   end if
   ef = energy - eFermi
   Cviso = 10.5d0
   V0r = -41.45d0
   Var = -0.06667d0
   vrdisp = 92.44d0
   v1r = 0.03d0
   v2r = 2.05d-4
   lambdaR = 3.9075d-3
   viso = 1.0d0 + phase*Cviso*asym/(V0r + Var*(xA-232.0d0)+vrdisp)
   v = (V0r + Var*(xA-232.0d0) + v1r*ef + v2r*(ef**2) + vrdisp*exp(-lambdaR*ef))*viso

   if(part_type == 2) then
     Ccoul = 0.9d0
     phicoul = (lambdaR*vrdisp*exp(-lambdaR*ef) - v1r - 2.0d0*v2r*ef)*viso
     v = v + Ccoul*xZ/(xA**onethird)*phicoul
   end if

   rr = 1.245d0
   Crr = 0.05d0
   widr = 100.0d0
   rv = rr*(1.0d0 - Crr*ef**2/(ef**2 + widr**2))
   av = (0.660d0 + 2.53d-4*energy)
   w1loc = 14.74d0
   w2loc = 81.63d0
   w = w1loc*ef**2/(ef**2 + w2loc**2)
   rw = 1.2476d0
   aw = 0.594d0
   vd = 0.0d0
   rvd = 1.2080d0
   avd = 0.614d0
   wddisp = 17.38d0
   Cwiso = 24.d0
   Wad = 0.03833d0
   d1loc = (wddisp + Wad*(xA-232.0d0) + phase*Cwiso*asym)
   d2loc = 0.01759d0
   d3loc = 11.79d0
   wd = d1loc*ef**2*exp(-d2loc*ef)/(ef**2 + d3loc**2)
   rwd = rvd
   awd = avd
   vso1loc = 5.86d0
   vso2loc = 0.0050d0
   vso = vso1loc*exp(-vso2loc*ef)
   rvso = 1.1213d0
   avso = 0.59d0
   wso1loc = -3.1d0
   wso2loc = 160.0d0
   wso = wso1loc*ef**2/(ef**2 + wso2loc**2)
   rwso = rvso
   awso = avso
   if(part_type == 1) then
     rc = 1.2643d0
   elseif(part_type == 2) then
     rc = 1.2643d0
   end if

   V_pot(1:2,1:3) = 0.0d0

   V_pot(1,1) = v
   R_pot(1,1) = rv
   a_pot(1,1) = av
   V_pot(2,1) = w
   R_pot(2,1) = rw
   a_pot(2,1) = aw
   V_pot(1,2) = vd
   R_pot(1,2) = rvd
   a_pot(1,2) = avd
   V_pot(2,2) = wd
   R_pot(2,2) = rwd
   a_pot(2,2) = awd
   V_pot(1,3) = vso
   R_pot(1,3) = rvso
   a_pot(1,3) = avso
   V_pot(2,3) = wso
   R_pot(2,3) = rwso
   a_pot(2,3) = awso

   return
end subroutine soukhovitskii_potential
!
!*******************************************************************************
!
subroutine perey_d_potential(part_type, iA, iZ, E,                   &
                             V_pot, R_pot, a_pot, RC, iradius)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates the parameters for the Perey
!    deuteron potential
!    Use with caution for energies below 12 MeV!
!
!    Reference:
!
!     C. M. Perey and F. G. Perey, At. Data and Nuc. Data Tables, Vol. 17,
!     No. 1, 1976, p. 1.
!          
!  Licensing:
!
!    This code is distrid under the GNU LGPL version 2 license. 
!
!  Date:
!
!    28 April 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   use nodeinfo
   use variable_kinds
   implicit none
   integer(kind=4), intent(in) :: part_type, iA, iZ
   real(kind=8), intent(in) :: E
   real(kind=8), intent(out) :: V_pot(2,3)
   real(kind=8), intent(out) :: R_pot(2,3)
   real(kind=8), intent(out) :: a_pot(2,3)
   real(kind=8), intent(out) :: RC
   integer(kind=4), intent(out) :: iradius
!----------------------------------------------------------
!   real(kind=8) :: VR, RR, AR
!   real(kind=8) :: W, RW, AW
!   real(kind=8) :: VD, RVD, AVD
!   real(kind=8) :: WD, RD, AD
!   real(kind=8) :: VSO, RSO, ASO
!   real(kind=8) :: WSO, WRSO, WASO
   integer(kind=4) :: iN
   real(kind=8) :: xZ, xN, xA, xA13
   real(kind=8) :: diff

!     iradius = 0    R = r0*A_T^1/3
!     iradius = 1    R = r0*(A_Target^1/3 + A_proj^1/3)
   iradius = 0

   if(part_type /= 3)then
      if(iproc == 0)then
         write(6,*)'Error in perey_d_potential. Attempting to call for '
         write(6,*)'an incident particle other than deutrons'
      end if
      call MPI_Abort(icomm,301,ierr)
   end if

   iN = iA - iZ
   xZ = real(iZ,kind=8)
   xN = real(iN,kind=8)
   xA = real(iA,kind=8)
   xA13 = xA**(1./3.)

   diff = (xN - xZ)/xA

   V_pot(1:2,1:3) = 0.0d0
!
!-------    Real volume potential
!
   V_pot(1,1) = 81.0d0 - 0.22d0*E + xZ/xA13
   R_pot(1,1) = 1.15d0
   a_pot(1,1) = 0.81d0
!
!-------   Imaginary surface potential
!
   V_pot(2,2) = 14.4d0 + 0.24d0*E
   R_pot(2,2) = 1.34d0
   a_pot(2,2) = 0.68d0
!
!-------   Coulomb radius
!
   RC = 1.15d0

   return
end subroutine perey_d_potential
!
!*******************************************************************************
!
subroutine becchetti_t_potential(part_type, iA, iZ, E,           &
                                 V_pot, R_pot, a_pot, RC, iradius)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates the parameters for the Perey
!    deuteron potential
!
!    Reference:
!
!    Quoted in C. M. Perey and F. G. Perey, At. Data and Nuc. Data Tables,
!    Vol. 17, No. 1, 1976, p. 1.
!          
!  Licensing:
!
!    This code is distributed under the GNU LGPL version 2 license. 
!
!  Date:
!
!    28 April 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   use nodeinfo
   use variable_kinds
   implicit none
   integer(kind=4), intent(in) :: part_type, iA, iZ
   real(kind=8), intent(in) :: E
   real(kind=8), intent(out) :: V_pot(2,3)
   real(kind=8), intent(out) :: R_pot(2,3)
   real(kind=8), intent(out) :: a_pot(2,3)
   real(kind=8), intent(out) :: RC
   integer(kind=4), intent(out) :: iradius
!----------------------------------------------------------
!   real(kind=8) :: VR, RR, AR
!   real(kind=8) :: W, RW, AW
!   real(kind=8) :: VD, RVD, AVD
!   real(kind=8) :: WD, RD, AD
!   real(kind=8) :: VSO, RSO, ASO
!   real(kind=8) :: WSO, WRSO, WASO
   integer(kind=4) :: iN
   real(kind=8) :: xZ, xN, xA
   real(kind=8) :: diff

!     iradius = 0    R = r0*A_T^1/3
!     iradius = 1    R = r0*(A_Target^1/3 + A_proj^1/3)
   iradius = 0

   if(part_type /= 4)then
      if(iproc == 0)then
         write(6,*)'Error in becchetti_t_potential. Attempting to call for '
         write(6,*)'an incident particle other than tritons'
      end if
      call MPI_Abort(icomm,301,ierr)
   end if

   iN = iA - iZ
   xZ = real(iZ,kind=8)
   xN = real(iN,kind=8)
   xA = real(iA,kind=8)

   diff = (xN - xZ)/xA

   V_pot(1:2,1:3) = 0.0d0
!
!-------    Real volume potential
!
   V_pot(1,1) = 165.0d0 - 0.17d0*E - 6.4d0*diff
   R_pot(1,1) = 1.20d0
   a_pot(1,1) = 0.72d0
!
!-------   Imaginary volume potential
!
   V_pot(2,1) = 46.0d0 - 0.33d0*E - 110.0d0*diff
   R_pot(2,1) = 1.40d0
   a_pot(2,1) = 0.84d0
!
!-------    Real Spin-orbit potential
!
   V_pot(1,3) = 2.5
   R_pot(1,3) = 1.20d0
   a_pot(1,3) = 0.72d0
!
!-------   Coulomb radius
!
   RC = 1.30d0

   return
end subroutine becchetti_t_potential
!
!*******************************************************************************
!
subroutine becchetti_h_potential(part_type, iA, iZ, E,           &
                                 V_pot, R_pot, a_pot, RC, iradius)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates the parameters for the Perey
!    He-3 potential
!
!    Reference:
!
!    Quoted in C. M. Perey and F. G. Perey, At. Data and Nuc. Data Tables,
!    Vol. 17, No. 1, 1976, p. 1.
!          
!  Licensing:
!
!    This code is distributed under the GNU LGPL version 2 license. 
!
!  Date:
!
!    28 April 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   use nodeinfo
   use variable_kinds
   implicit none
   integer(kind=4), intent(in) :: part_type, iA, iZ
   real(kind=8), intent(in) :: E
   real(kind=8), intent(out) :: V_pot(2,3)
   real(kind=8), intent(out) :: R_pot(2,3)
   real(kind=8), intent(out) :: a_pot(2,3)
   real(kind=8), intent(out) :: RC
   integer(kind=4), intent(out) :: iradius
!----------------------------------------------------------
!   real(kind=8) :: VR, RR, AR
!   real(kind=8) :: W, RW, AW
!   real(kind=8) :: VD, RVD, AVD
!   real(kind=8) :: WD, RD, AD
!   real(kind=8) :: VSO, RSO, ASO
!   real(kind=8) :: WSO, WRSO, WASO
   integer(kind=4) :: iN
   real(kind=8) :: xZ, xN, xA
   real(kind=8) :: diff

!     iradius = 0    R = r0*A_T^1/3
!     iradius = 1    R = r0*(A_Target^1/3 + A_proj^1/3)
   iradius = 0

   if(part_type /= 5)then
      if(iproc == 0)then
         write(6,*)'Error in becchetti_h_potential. Attempting to call for '
         write(6,*)'an incident particle other than He-3'
      end if
      call MPI_Abort(icomm,301,ierr)
   end if

   iN = iA - iZ
   xZ = real(iZ,kind=8)
   xN = real(iN,kind=8)
   xA = real(iA,kind=8)

   diff = (xN - xZ)/xA

   V_pot(1:2,1:3) = 0.0d0
!
!-------    Real volume potential
!
   V_pot(1,1) = 151.9d0 - 0.17d0*E + 50d0*diff
   R_pot(1,1) = 1.20d0
   a_pot(1,1) = 0.72d0
!
!-------   Imaginary volume potential
!
   V_pot(2,1) = 41.7d0 - 0.33d0*E + 44.0d0*diff
   R_pot(2,1) = 1.40d0
   a_pot(2,1) = 0.88d0
!
!-------    Real Spin-orbit potential
!
   V_pot(1,3) = 2.5
   R_pot(1,3) = 1.20d0
   a_pot(1,3) = 0.72d0
!
!-------   Coulomb radius
!
   RC = 1.30d0

   return
end subroutine becchetti_h_potential
!
!*******************************************************************************
!
subroutine avrigeanu_a_potential(part_type, iA, iZ, E,                  &
                                 V_pot, R_pot, a_pot, RC, iradius)
!
!*******************************************************************************
!
!  Discussion:
!
!    This subroutine calculates the parameters for the Perey
!    deuteron potential
!
!    Reference:
!
!    V. Avrigeanu, P. E. Hodgson, and M. Avrigeanu, Phys. Rev. C49, 2136 (1994).
!          
!  Licensing:
!
!    This code is distributed under the GNU LGPL version 2 license. 
!
!  Date:
!
!    28 April 2020
!
!  Author:
!
!      Erich Ormand, LLNL
!
!*******************************************************************************
!
   use nodeinfo
   use variable_kinds
   implicit none
   integer(kind=4), intent(in) :: part_type, iA, iZ
   real(kind=8), intent(in) :: E
   real(kind=8), intent(out) :: V_pot(2,3)
   real(kind=8), intent(out) :: R_pot(2,3)
   real(kind=8), intent(out) :: a_pot(2,3)
   real(kind=8), intent(out) :: RC
   integer(kind = 4), intent(out) :: iradius
!----------------------------------------------------------
!   real(kind=8) :: VR, RR, AR
!   real(kind=8) :: W, RW, AW
!   real(kind=8) :: VD, RVD, AVD
!   real(kind=8) :: WD, RD, AD
!   real(kind=8) :: VSO, RSO, ASO
!   real(kind=8) :: WSO, WRSO, WASO
   integer(kind=4) :: iN
   real(kind=8) :: xZ, xN, xA, xA13
   real(kind=8) :: diff

!     iradius = 0    R = r0*A_T^1/3
!     iradius = 1    R = r0*(A_Target^1/3 + A_proj^1/3)
   iradius = 0

   if(part_type /= 6)then
      if(iproc == 0)then
         write(6,*)'Error in avrigeanu_a_potential. Attempting to call for '
         write(6,*)'an incident particle other than alphas'
      end if
      call MPI_Abort(icomm,301,ierr)
   end if

   iN = iA - iZ
   xZ = real(iZ,kind=8)
   xN = real(iN,kind=8)
   xA = real(iA,kind=8)
   xA13 = xA**(1./3.)

   diff = (xN - xZ)/xA

   V_pot(1:2,1:3) = 0.0d0
!
!-------    Real volume potential
!
   V_pot(1,1) = 101.1 - 0.248d0*E + 6.051*xZ/xA13
   R_pot(1,1) = 1.20d0
   a_pot(1,1) = 0.72d0
!
!-------   Imaginary volume potential
!
   if( E < 73.0d0)then
      V_pot(2,1) = 12.64 + 0.2d0*E - 1.706*xA13
   else
      V_pot(2,1) = 26.82 + 0.006*E - 1.706*xA13
   end if
   R_pot(2,1) = 1.57d0
   a_pot(2,1) = 0.692 - 0.02*xA13
!
!-------   Coulomb radius
!
   RC = 1.245d0

   return
end subroutine avrigeanu_a_potential
