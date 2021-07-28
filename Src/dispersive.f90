  !     *******************************************************
  !     START of dispersive PACK
  !     *******************************************************
  !==========================================================================
  !     AUTHOR: Dr. Roberto Capote Noy
  !
  !     e-mail: r.capotenoy@iaea.org ; rcapotenoy@yahoo.com;
  !
  !     DISPERSIVE OPTICAL MODEL POTENTIAL PACKAGE
  !
  !     Analytical dispersive integrals are included
  !     see Quesada JM, Capote R et al,
  !             Computer Physics Communications 153(2003) 97
  !             Phys. Rev. C67(2003) 067601
  !
  !     Dispersive integral's derivatives calculated by Dr.J.M.Quesada
  !     
  !     Further coding July 2021 by Ian Thompson.
  !         Conversion to double precision everywhere.
  !
  !==========================================================================

  REAL*8 FUNCTION DOM_INT_Wv (Ef,Ep,Av,Bv,n,Einc,DerivIntWv)

    IMPLICIT NONE
    REAL*8 Ef,Ep,Av,Bv,E,pi,Einc
    REAL*8 E0,Ex,Eplus,Emin,Rs,ResEmin,ResEplus
    REAL*8 DerEmin, DerEplus, Rds, DerivIntWv
    DOUBLE COMPLEX Pj,I,Zj,Ztmp
    DOUBLE COMPLEX Fs,Ds
    INTEGER N,j,IS

    DATA I/(0.d0,1.d0)/

    pi=4.0*atan(1.0)

    IS = 1
    E = Einc
    IF(Einc.LE.Ef) THEN
      E=2.0*Ef-Einc
      IS = -1
    ENDIF

    E0 = Ep - Ef
    Ex = E  - Ef
    Eplus = Ex + E0
    Emin  = Ex - E0
    DOM_INT_Wv = 0.0
    DerivIntWv = 0.0

    ResEmin  =  Emin**n / (Emin**n + Bv**n)

    DerEmin  =  Emin**(n-1)* &
                ( Emin**n + Bv**n*(1.0 + n*log(ABS(Emin)) ) ) &
                / (Emin**n + Bv**n)**2

    ResEplus = -Eplus**n / (Eplus**n + Bv**n)

    DerEplus = -Eplus**(n-1) * &
                   ( Eplus**n + Bv**n*(1.0+n*log(Eplus)) ) &
                   / (Eplus**n + Bv**n)**2

    Fs = (0.d0,0.d0)
    Ds = (0.d0,0.d0)
    DO j=1,n
      Ztmp = I*(2*j-1)/dble(n)*pi
      Pj = Bv*exp(Ztmp)
      Zj = Pj * (2*Pj +Eplus -Emin) * Ex
      Zj = Zj / ( (Pj+E0) * (Pj+Eplus) * (Pj-Emin) )
      Fs = Fs + Zj*log(-Pj)
      Ds = Ds + 2*Pj*(Ex*Ex + (Pj+E0)**2)*log(-Pj) &
                  /( (Pj+Eplus)**2 * (Pj-Emin)**2 )
    ENDDO

    IF(ABS(IMAG(Fs)).gt.1.d-4) STOP 'Too big imag part in Wv'
    IF(ABS(IMAG(Ds)).gt.1.d-4) STOP 'Too big imag deriv in Wv'
    Rs  = DBLE(Fs)
    Rds = DBLE(Ds)
    DOM_INT_Wv = -Av/pi*IS* &
        (Rs/n  + ResEplus*log(Eplus) + ResEmin*log(ABS(Emin)))
    DerivIntWv =  Av/pi*IS*( Rds/n + DerEplus + DerEmin)

    RETURN
  END FUNCTION DOM_INT_Wv

  REAL*8 FUNCTION DOM_INT_Ws (Ef,Ep,As,Bs,Cs,m,Einc,DerivIntWs)

    IMPLICIT NONE
    REAL*8 Ef,Ep,As,Bs,Cs,E,Einc
    DOUBLE COMPLEX I,Pj,Zj,Ztmp
    REAL*8 Emin,Eplus,E0,Ex,pi
    REAL*8 Rs,ResEmin,ResEplus
    REAL*8 DerivIntWs,DerEmin,DerEplus,Rds
    INTEGER m,j,IS
    DOUBLE COMPLEX Fs,Ds,zfi
    REAL*8 EIn
    
    DATA I/(0.d0,1.d0)/

    pi=4.0*atan(1.0)

    IS = 1
    E = Einc
    IF(Einc.LE.Ef) THEN
      E=2.0*Ef-Einc
      IS = -1
    ENDIF

    E0 = Ep - Ef
    Ex = E  - Ef
    Eplus = Ex + E0
    Emin  = Ex - E0
    DOM_INT_Ws = 0.0
    DerivIntWs = 0.0
    ResEmin  =  Emin**m / (Emin**m + Bs**m)

    DerEmin  = -Emin**(m-1) * &
                 ( Emin**m + Bs**m + ( -Cs*Emin**(m+1) + &
                Bs**m *(-Cs*Emin+m) ) * exp(-Cs*Emin)*EIn(Cs*Emin) ) &
                / (Emin**m + Bs**m)**2

    ResEplus = -Eplus**m / (Eplus**m + Bs**m)

    DerEplus =  Eplus**(m-1) * &
                  ( Eplus**m + Bs**m + ( Cs*Eplus**(m+1) + &
                   Bs**m *(Cs*Eplus+m) ) * exp(Cs*Eplus)*EIn(-Cs*Eplus) ) &
                   / (Eplus**m + Bs**m)**2

    Fs = (0.d0,0.d0)
    Ds = (0.d0,0.d0)
    DO j=1,m
         Ztmp = I*(2*j-1)/dble(m)*pi
         Pj = Bs*exp(Ztmp)
         Zj = Pj * (2*Pj +Eplus -Emin) * Ex
         Zj = Zj / (Pj+E0) / (Pj+Eplus) / (Pj-Emin)
         Fs = Fs + Zj* zfi(-Pj*Cs)
         Ds = Ds + 2*Pj*(Ex*Ex + (Pj+E0)**2)*zfi(-Pj*Cs) &
                  /( (Pj+Eplus)**2 * (Pj-Emin)**2 )
    ENDDO

    IF(ABS(IMAG(Fs)).GT.1.d-4) STOP 'Too big imag part in Ws'
    IF(ABS(IMAG(Ds)).GT.1.d-4) STOP 'Too big imag deriv in Ws'
    Rs = DBLE(Fs)
    Rds = DBLE(Ds)

    DOM_INT_Ws = As/pi*IS*(Rs/m &
                         - ResEplus*exp(Cs*Eplus)*EIn(-Cs*Eplus) &
                         - ResEmin*exp(-Cs*Emin)*EIn(Cs*Emin) )
    RETURN
  END FUNCTION DOM_INT_Ws

  REAL*8 FUNCTION WV(A,B,Ep,Ef,E,n)

    IMPLICIT NONE
    REAL*8  A,B,Ep,Ef,E,ee
    INTEGER n
    WV=0.d0
    IF(E.LE.Ef) E=2.d0*Ef-E
    IF(E.LT.Ep) RETURN
    ee=(E-Ep)**n
    WV=A*ee/(ee+B**n)
    RETURN
  END FUNCTION WV

  REAL*8 FUNCTION WDD(A,B,C,Ep,Ef,E,m)
    IMPLICIT NONE
    REAL*8 A,B,C,Ep,Ef,E,ee,arg
    INTEGER m
    WDD=0.d0
    IF(E.LE.Ef) E=2.d0*Ef-E
    IF(E.LT.Ep) RETURN
    arg=C*(E-Ep)
    IF(arg.GT.15) RETURN
    ee=(E-Ep)**m
    WDD=A*ee/(ee+B**m)*EXP(-arg)
    RETURN
  END FUNCTION WDD

  REAL*8 FUNCTION DOM_int_T1(Ef,Ea,E)

    IMPLICIT NONE
    REAL*8 E,Ea,Ef,Ex,Ea2,Eax,Pi,T11,T12,T13
    Pi=4.0*ATAN(1.0)
    Ex=E-Ef
    Ea2=Ea**2
    Eax=Ex+Ea
    T11 = 0.5d0*log(Ea)/Ex
    T12 =  ( (2*Ea+Ex)*log(Ea)+0.5d0*pi*Ex ) &
             /(2.*(Eax**2 + Ea2))
    T13 = -Eax**2*log(Eax)/(Ex*(Eax**2+Ea2))
    DOM_int_T1 = Ex/Pi*(T11+T12+T13)
    RETURN
  END FUNCTION DOM_int_T1

  REAL*8 FUNCTION DOM_int_T2(Ef,Ea,E)

    IMPLICIT NONE
    REAL*8 E,Ea,Ef,EL,Pi
    Pi=4.0*ATAN(1.0)
    EL=Ef+Ea
    DOM_int_T2= 1.0 / Pi * ( &
             sqrt(abs(Ef)) * atan( (2*sqrt(EL*abs(Ef)))/(EL-abs(Ef)) ) &
        +    EL**1.50/(2*Ef)*log(Ea/EL) )
    IF(E.GT.EL) THEN
      DOM_int_T2 = DOM_int_T2 + 1.0/Pi* ( &
         sqrt(E) * log( (sqrt(E)+sqrt(EL)) / (sqrt(E)-sqrt(EL)) ) + &
         1.50*sqrt(EL)*log((E-EL)/Ea) + EL**1.50/(2*E)*log(EL/(E-EL)) )
    ELSEIF(E.EQ.EL) THEN
      DOM_int_T2 = DOM_int_T2 + 1.0/Pi*1.50*sqrt(EL) &
        *log((2**(4.0/3.0)*EL)/Ea)
    ELSEIF(E.GT.0.d0 .AND. E.LE.EL) THEN
      DOM_int_T2 = DOM_int_T2 + 1.0/Pi * ( &
        sqrt(e) * log( (sqrt(E)+sqrt(EL)) / (sqrt(EL)-sqrt(E)) ) + &
        1.50*sqrt(EL)*log((EL-E)/Ea)+EL**1.50/(2.0*E)*log(EL/(EL-E)) )
    ELSEIF(abs(E)<1e-10) then
      DOM_int_T2 = DOM_int_T2 + 1.0/Pi*( 1.5*sqrt(EL) &
        * log(EL/Ea) + 0.50*sqrt(EL) )
    ELSE
      DOM_int_T2 = DOM_int_T2 + 1.0/Pi * ( &
        -sqrt(abs(E))*atan( 2*(sqrt(EL*abs(E))) / (EL-abs(E)) ) + &
        1.50*sqrt(EL)*log((EL-E)/Ea)+EL**1.50/(2.0*E)*log(EL/(EL-E)) )
    ENDIF
  	!write(101,*) E,DOM_int_T2,Ef,Ea
    RETURN
  END FUNCTION DOM_int_T2

  DOUBLE COMPLEX FUNCTION zfi(za)

    IMPLICIT NONE
    REAL*8 aj
    DOUBLE COMPLEX za,y
    INTEGER m,i
    zfi=0.d0
    IF (za.EQ.0.d0) RETURN
    IF (abs(real(za)+18.5d0).GE.25.d0) GO TO 3
    IF (SQRT(625.d0-(DBLE(za)+18.5d0)**2)/1.665d0.LT.ABS(imag(za))) GO TO 3
    zfi=-.57721566490153d0-log(za)
    y=1.d0
    DO 1 m=1,2000
      aj=m
      y=-y*za/aj
      IF (ABS(y).lt.1.d-15*ABS(zfi)) GO TO 2
        zfi=zfi-y/aj
      1 continue
      2 zfi=EXP(za)*zfi
        RETURN
      3 DO 4 i=1,20
          aj=21-i
          zfi=aj/(za+zfi)
            zfi=aj/(1.d0+zfi)
          4 continue
          zfi=1.d0/(zfi+za)
    RETURN
  END FUNCTION zfi

  REAL*8 FUNCTION EIn(X)

    IMPLICIT NONE
    REAL*8 FAC, H, X
    INTEGER N
    EIn = 0.57721566490153d0+LOG(ABS(X))
    FAC = 1.0
    DO N = 1,100
      H = FLOAT(N)
      FAC = FAC*H
      EIn = EIn + X**N/(H*FAC)
    ENDDO
    RETURN
  END FUNCTION EIn


  SUBROUTINE dispers2(A,Z,k,eopt, &
  v,rvv,avv, dv,drv,dav, dvs,drs,das, w,rw,aw, wd,rwd,awd, &
  vso,rvso,avso,dvso, wso,rwso,awso, &
  Vlin,Vdep,lambdaHF,Cviso,Vso0,lambdaso,Ccoul, &
  AAv,BBv,W0l,W0dep,BBs,CCs,Cwiso,Wso0,BBso, &
  Ea,alpha,eferm,Ades, &
  rHFl,rHFdep,aHFl,aHFdep,rv,avpot, &
  rsl,rsdep,as, &
  rso,aso)
  
    implicit REAL*8(A-H,O-Z)
    REAL*8 eopt,asym,eferm,f,Cviso,Ccoul,Cwiso
    REAL*8 lambdaHF,lambdaso,Ades
    pi = 4.d0*atan(1.d0)
    !
    ! *** Parameters definitions of Soukhovitskii, Capote, Quesada, Chiba and Martyanov (Nov 25, 2015) ***
    ! *** with Asymmetrical W energy-dependence
    ! *** dispers2:  T1 integral with correct coefficient
    ! *** done by Ian Thompson
    ! k         : designator for particle
    ! Z         : charge number of residual nucleus
    ! A         : mass number of residual nucleus
    ! eopt      : incident energy
    ! asym      : asymmetry parameter
    ! eferm     : Fermi energy
    ! f         : eopt-eferm

    Au = A-Ades
    asym=(A-2.*Z)/A
    V0 = Vlin + Vdep*Au
    rHF = rHFL + rHFdep * Au
    aHF = aHFl + aHFdep * Au
    W0 = W0l + W0dep*Au
    av = avpot
    rs = rsl + rsdep * Au
    AAHF = V0 * (1 + (-1)**k * Cviso*asym/V0)
    AAs  = W0 * (1 + (-1)**k * Cwiso*asym/W0)
    eoffset = 0.
    IF (k==2) eoffset =  Ccoul * Z/A**(1./3.)
    Eeff = eopt - eoffset
    f = Eeff  - eferm
    v = AAHF * EXP(-lambdaHF*f)
    vso=Vso0*EXP(-lambdaso*f)

    ! sources of dispersive terms
    w = AAv * f*f/(f*f + BBv**2)
    IF(f < -Ea) THEN
      fe  = f + Ea
      w = w * (1 - fe*fe/(fe*fe + Ea*Ea))
    ELSE IF(f>Ea) THEN
      w = w + alpha * (SQRT(Eeff) + (eferm+Ea)**1.50/(2*Eeff) &
           - 1.50*SQRT(eferm+Ea))
    endif
    wd = AAs * f*f/(f*f + BBs**2) * EXP( -CCs * ABS(f))
    wso=Wso0* f*f/(f*f + BBso**2)
    ! dispersive terms to add for real volume and real surface forms
    drv= rv;  dav = av
    drs = rs; das = as
    dvs = DOM_INT_Ws (eferm,eferm,AAs,BBs,CCs,2,Eeff,DerivIntWs)
    DWv = DOM_INT_Wv (eferm,eferm,AAv,BBv,2,Eeff,DerivIntWv)
    T1 = DOM_int_T1(eferm,Ea,Eeff) * AAv * Ea*Ea/(Ea*Ea + BBv**2)
    T2 = DOM_int_T2(eferm,Ea,Eeff) * alpha
    dv = DWv + T1 + T2
    dvso = DOM_INT_Wv (eferm,eferm,Wso0,BBso,2,Eeff,DerivIntWv)


    ! name translations
    rvv = rHF ; avv = aHF
    rvso = rso; avso = aso
    rwso = rso; awso = aso
    rw  = rv  ; aw   = av
    rwd = rs  ; awd  = as
    RETURN
  END SUBROUTINE dispers2

   
  ! ***********************************************************
  ! *                    END of dispersive                    *
  ! ***********************************************************

