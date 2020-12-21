!----   Distributed as part of the RIPL-3 level library to read
!----   RIPL level libraries

      PROGRAM READ_LEVELS
!
! PARAMETER definitions
!
      INTEGER, PARAMETER :: LEVMAX = 600, GAMMAX = 1200, BRMAX = 30
!
! COMMON variables
!
      INTEGER :: A, I, J, JJ, K, M, N, NC, NGAm, NLEv, NMAx, W, Z
      LOGICAL :: BOTh, EVEn, HASq
      REAL(KIND=4), DIMENSION(BRMAX) :: BR
      REAL(KIND=4), DIMENSION(GAMMAX) :: CC_gam, E_Gam, RI1_gam, RI_gam
      REAL(KIND=4) :: DIFf, EG, IT_percent, L1, LL1, PF, PI, RISum,     &
                    & RITot, SF, SI, SN, SP, SUM, TOTsum
      REAL(KIND=8) :: DMIx, EGD, ICC
      CHARACTER(7), DIMENSION(LEVMAX, BRMAX) :: DMOdes_lev
      REAL(KIND=4), DIMENSION(LEVMAX, BRMAX) :: DPErcent_lev
      REAL(KIND=4), DIMENSION(LEVMAX) :: E_Lev, JPMy_lev, JSMy_lev,     &
                & T_Lev
      INTEGER, DIMENSION(GAMMAX) :: FINal_gam, INItial_gam
      CHARACTER(100) :: FNAme
      CHARACTER(1) :: IMPos
      CHARACTER(1), DIMENSION(LEVMAX) :: JPEstimate_lev
      CHARACTER(18), DIMENSION(LEVMAX) :: JPText_lev
      REAL(KIND=8), DIMENSION(GAMMAX) :: MIXr_gam, PROb_gam
      CHARACTER(10), DIMENSION(GAMMAX) :: MULt_gam
      INTEGER, DIMENSION(LEVMAX) :: NOBr_lev, NOG_lev
      CHARACTER(2), DIMENSION(LEVMAX, BRMAX) :: PFIx_lev
      CHARACTER(10) :: STRmul, TEXt
      CHARACTER(5) :: SYMb
      CHARACTER(4), DIMENSION(LEVMAX) :: UNCertaint_lev
      COMMON /INP   / SYMb, JPEstimate_lev, UNCertaint_lev, JPText_lev, &
                    & PFIx_lev, DMOdes_lev, IMPos, FNAme, MULt_gam, A,  &
                    & Z, NLEv, NGAm, NMAx, NC, M, I, K, J, N, NOG_lev,  &
                    & NOBr_lev, INItial_gam, FINal_gam, W, JJ, BOTh,    &
                    & HASq, EVEn, SN, SP, E_Gam, RI_gam, RI1_gam,       &
                    & CC_gam, SUM, DIFf, EG, SI, SF, PI, PF, RITot, BR, &
                    & E_Lev, JSMy_lev, JPMy_lev, T_Lev, DPErcent_lev,   &
                    & LL1, L1, IT_percent, RISum, TOTsum, MIXr_gam,     &
                    & PROb_gam, EGD, DMIx, ICC, STRmul, TEXt
!
! Local variables
!
      LOGICAL :: data_exists
!
!
!Generate all z & a & call getza subroutine
!The  getza routine provide a logical flag. If data exists it is true else it is
!false
      WRITE(*,*) 'YOU MAY INPUT ATOMIC(>0) OR MASS(<0) NUMBER TO RETRIEVE LEVELS:'
	  READ(*,*) Z
	  WRITE(*,*) 
	  WRITE(*,*) 'Isot    A    Z  Nlev Ngam Ncut  NC        Sn        Sp'

	  IF(Z.GT.0) THEN
        DO A = 1, 294
          CALL GETZA(data_exists)
        ENDDO
      ENDIF
	  IF(Z.LT.0) THEN
 	    A=-Z 
        DO Z = 0, 118
          CALL GETZA(data_exists)
        ENDDO
      ENDIF

	  STOP 'READ_LEVELS OK'
      END

      SUBROUTINE GETZA(Data_exists)
!
! PARAMETER definitions
!
      INTEGER, PARAMETER :: LEVMAX = 600, GAMMAX = 1200, BRMAX = 30
!
! COMMON variables
!
      INTEGER :: A, I, J, JJ, K, M, N, NC, NGAm, NLEv, NMAx, W, Z
      LOGICAL :: BOTh, EVEn, HASq
      REAL(KIND=4), DIMENSION(BRMAX) :: BR
      REAL(KIND=4), DIMENSION(GAMMAX) :: CC_gam, E_Gam, RI1_gam, RI_gam
      REAL(KIND=4) :: DIFf, EG, IT_percent, L1, LL1, PF, PI, RISum,     &
                    & RITot, SF, SI, SN, SP, SUM, TOTsum
      REAL(KIND=8) :: DMIx, EGD, ICC
      CHARACTER(7), DIMENSION(LEVMAX, BRMAX) :: DMOdes_lev
      REAL(KIND=4), DIMENSION(LEVMAX, BRMAX) :: DPErcent_lev
      REAL(KIND=4), DIMENSION(LEVMAX) :: E_Lev, JPMy_lev, JSMy_lev,     &
                & T_Lev
      INTEGER, DIMENSION(GAMMAX) :: FINal_gam, INItial_gam
      CHARACTER(100) :: FNAme
      CHARACTER(1) :: IMPos
      CHARACTER(1), DIMENSION(LEVMAX) :: JPEstimate_lev
      CHARACTER(18), DIMENSION(LEVMAX) :: JPText_lev
      REAL(KIND=8), DIMENSION(GAMMAX) :: MIXr_gam, PROb_gam
      CHARACTER(10), DIMENSION(GAMMAX) :: MULt_gam
      INTEGER, DIMENSION(LEVMAX) :: NOBr_lev, NOG_lev
      CHARACTER(2), DIMENSION(LEVMAX, BRMAX) :: PFIx_lev
      CHARACTER(10) :: STRmul, TEXt
      CHARACTER(5) :: SYMb
      CHARACTER(4), DIMENSION(LEVMAX) :: UNCertaint_lev
      COMMON /INP   / SYMb, JPEstimate_lev, UNCertaint_lev, JPText_lev, &
                    & PFIx_lev, DMOdes_lev, IMPos, FNAme, MULt_gam, A,  &
                    & Z, NLEv, NGAm, NMAx, NC, M, I, K, J, N, NOG_lev,  &
                    & NOBr_lev, INItial_gam, FINal_gam, W, JJ, BOTh,    &
                    & HASq, EVEn, SN, SP, E_Gam, RI_gam, RI1_gam,       &
                    & CC_gam, SUM, DIFf, EG, SI, SF, PI, PF, RITot, BR, &
                    & E_Lev, JSMy_lev, JPMy_lev, T_Lev, DPErcent_lev,   &
                    & LL1, L1, IT_percent, RISum, TOTsum, MIXr_gam,     &
                    & PROb_gam, EGD, DMIx, ICC, STRmul, TEXt
!
! Dummy arguments
!
      LOGICAL :: Data_exists
!
! Local variables
!
      INTEGER :: aa, zz
!
      Data_exists = .FALSE.
!Check if data request is in range
      IF(Z.LT.0 .OR. Z.GT.118)THEN
         WRITE(*, *)'Data with this charge number Z does not exists Z=',&
                  & Z
         RETURN
      ENDIF
      IF(A.LT.1 .OR. A.GT.294)THEN
         WRITE(*, *)'Data with this mass number A does not exists A=', A
         RETURN
      ENDIF
!open file with the give Z
      IF(Z.LT.10)WRITE(FNAme, '(a4,i1,a4)')'z00', Z, '.dat'
      IF(Z.LT.100 .AND. Z.GE.10)WRITE(FNAme, '(a3,i2,a4)')'z0', Z,     &
                                     &'.dat'
      IF(Z.GT.99)WRITE(FNAme, '(a2,i3,a4)')'z', Z, '.dat'
!write(fname,'(a2,i3.0,a4)') 'z',Z,'.dat'
      OPEN(2, FILE = FNAme, STATUS = 'old',ERR=200)
 100  READ(2, '(a5,6i5,2f12.6)', END = 200)SYMb, aa, zz, NLEv, NGAm,    &
         & NMAx, NC, SN, SP

      J = 0
      JJ = 0
      IF(A.EQ.aa .AND. Z.EQ.zz)THEN
         Data_exists = .TRUE.
  	     WRITE(*,'(a5,6i5,2f12.6)')SYMb,aa,zz,NLEv,NGAm,NMAx,NC,SN,SP
		 CLOSE(2)
		 RETURN
      ENDIF
      DO JJ = 1, NLEv
         READ(2,                                                        &
     &'(i3,1x,f10.6,1x,f5.1,i3,1x,(1pe10.2),i3,1x,a1,1x,a4,1x,a18,i3,10(&
     &1x,a2,1x,0pf10.4,1x,a7))')I, E_Lev(I), JSMy_lev(I), W, T_Lev(I),  &
                              & NOG_lev(I), JPEstimate_lev(I),          &
                              & UNCertaint_lev(I), JPText_lev(I),       &
                              & NOBr_lev(I),                            &
                              & (PFIx_lev(I, M), DPErcent_lev(I, M),    &
                              & DMOdes_lev(I, M), M = 1, NOBr_lev(I))
         JPMy_lev(I) = W
         DO K = 1, NOG_lev(I)
            J = J + 1
!           initial_gam(j), &
            READ(2, '(39x,i4,1x,f10.3,3(1x,e10.3))')FINal_gam(J),       &
               & E_Gam(J), RI_gam(J), RI1_gam(J), CC_gam(J)
         ENDDO
      ENDDO
      GOTO 100
 200  CLOSE(2)
      RETURN
      END


      FUNCTION EVEN(I)
!
! Dummy arguments
!
      LOGICAL :: EVEN
      INTEGER :: I
!
! Local variables
!
      REAL :: FLOAT
      INTEGER :: INT
      REAL(KIND=4) :: x
!
      x = I
      EVEN = .FALSE.
      IF(x/2. - FLOAT(INT(x/2.+.1)).LT.0.1)EVEN = .TRUE.
      END
