MODULE VISCO
  USE ISO_C_BINDING

  ! ----------------------------------------------------------------------- !
  ! PROCEDURES TO APPLY VISCOELASTIC RELAXATION
  !
  ! NOTE
  ! ----
  ! SYMMETRIC SECOND ORDER TENSOR ORDERING CORRESPONDS TO ABAQUS EXPLICIT,
  ! THAT IS, A SYMMETRIC SECOND ORDERING TENSOR A IS STORED AS
  !
  !                          | A(1)  A(4)  A(6) |
  !                          |       A(2)  A(5) |
  !                          |             A(3) |
  !
  ! THIS ORDERING DIFFERS FROM ABAQUS STANDARD THAT USES
  !
  !                          | A(1)  A(4)  A(5) |
  !                          |       A(2)  A(6) |
  !                          |             A(3) |
  !
  ! WHY CHOOSE THE FIRST? HABIT, AND IT IS THE ORDERING USED AT SANDIA, ANSYS,
  ! LSDYNA, ...
  !
  ! FOR USE IN MATMODLAB, THIS FILE MUST BE LINKED WITH mmlfio.f90
  ! ----------------------------------------------------------------------- !
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: PROPCHECK, VISCOINI, VISCORELAX
  INTEGER, PARAMETER :: DP=SELECTED_REAL_KIND(14)
  INTEGER, PARAMETER :: NTENS=6
  INTEGER, PARAMETER ::  KXX=1, KYY=2, KZZ=4, KXY=4, KYZ=5, KZX=6
  REAL(KIND=DP), PARAMETER :: TOLER=3.0001E+00_DP
  REAL(KIND=DP), PARAMETER :: ZERO=0._DP, ONE=1._DP, TWO=2._DP, THREE=3._DP
  REAL(KIND=DP), PARAMETER :: SIX=6._DP, TEN=10._DP, THIRTY=30._DP
  REAL(KIND=DP), PARAMETER :: HALF=.5_DP, P3RD=ONE/THREE, P23RD=TWO/THREE
  REAL(KIND=DP), PARAMETER :: P6TH=ONE/SIX
  REAL(KIND=DP), PARAMETER :: P1M3=1.E-3_DP, P1M10=1.E-10_DP, P1M30=1.E-30_DP
  REAL(KIND=DP), PARAMETER :: POINT03=.03_DP, BIGNUM=1.E+30_DP
  !REAL(KIND=DP), PARAMETER :: I3x3(3,3)=RESHAPE((/ ONE, ZERO, ZERO,&
  !                                                 ZERO, ONE, ZERO,&
  !                                                 ZERO, ZERO, ONE/), (/3,3/))
  REAL(KIND=DP), PARAMETER :: I6(6)=(/ONE,ONE,ONE,ZERO,ZERO,ZERO/)
  REAL(KIND=DP), PARAMETER :: MACHINE_EPSILON=EPSILON(ONE)

  ! PROPERTIES
  ! ----------
  INTEGER, PARAMETER :: NUI=24
  INTEGER, PARAMETER :: IPWC1=1  ! WLF C1
  INTEGER, PARAMETER :: IPWC2=2  ! WLF C2
  INTEGER, PARAMETER :: IPWTR=3  ! WLF TREF
  INTEGER, PARAMETER :: IPGOO=4  ! PRONY SHEAR INFINITY
  INTEGER, PARAMETER :: IPG01=IPGOO+1  ! PRONY SHEAR COEFFICIENTS (10)
  INTEGER, PARAMETER :: IPG02=IPGOO+2
  INTEGER, PARAMETER :: IPG03=IPGOO+3
  INTEGER, PARAMETER :: IPG04=IPGOO+4
  INTEGER, PARAMETER :: IPG05=IPGOO+5
  INTEGER, PARAMETER :: IPG06=IPGOO+6
  INTEGER, PARAMETER :: IPG07=IPGOO+7
  INTEGER, PARAMETER :: IPG08=IPGOO+8
  INTEGER, PARAMETER :: IPG09=IPGOO+9
  INTEGER, PARAMETER :: IPG10=IPGOO+10
  INTEGER, PARAMETER :: IPT=IPG10 ! SHEAR RELAX TIME (10)
  INTEGER, PARAMETER :: IPT01=IPT+1
  INTEGER, PARAMETER :: IPT02=IPT+2
  INTEGER, PARAMETER :: IPT03=IPT+3
  INTEGER, PARAMETER :: IPT04=IPT+4
  INTEGER, PARAMETER :: IPT05=IPT+5
  INTEGER, PARAMETER :: IPT06=IPT+6
  INTEGER, PARAMETER :: IPT07=IPT+7
  INTEGER, PARAMETER :: IPT08=IPT+8
  INTEGER, PARAMETER :: IPT09=IPT+9
  INTEGER, PARAMETER :: IPT10=IPT+10

  ! STATE DEPENDENT VARIABLES
  ! ----- --------- ---------
  INTEGER, PARAMETER :: NSDV=68
  !     (1) : WLF AVERAGE SHIFT FACTOR
  !     (2) : AVERAGE NUMERICAL SHIFT FACTOR
  !   (3:8) : INSTANTANEOUS DEVIATORIC PK2 STRESS COMPONENTS AT START OF
  !           CURRENT TIME STEP WRITTEN USING THE INITIAL CONFIGURATION AS
  !           THE REFERENCE STATE
  !  (9:14) : XX,YY,ZZ,XY,YZ,ZX COMPONENTS OF VISCOELASTIC DEVIATORIC 2ND
  !           PIOLA KIRCHHOFF (PK2) STRESS FOR 1ST PRONY TERM USING THE
  !           INITIAL CONFIGURATION AS THE REFERENCE STATE
  ! (15:20) : VISCO DEV PK2 STRESS FOR 2ND PRONY TERM
  ! (21:26) : VISCO DEV PK2 STRESS FOR 3RD PRONY TERM
  ! (27:32) : VISCO DEV PK2 STRESS FOR 4TH PRONY TERM
  ! (33:38) : VISCO DEV PK2 STRESS FOR 5TH PRONY TERM
  ! (39:44) : VISCO DEV PK2 STRESS FOR 6TH PRONY TERM
  ! (45:50) : VISCO DEV PK2 STRESS FOR 7TH PRONY TERM
  ! (51:56) : VISCO DEV PK2 STRESS FOR 8TH PRONY TERM
  ! (57:62) : VISCO DEV PK2 STRESS FOR 9TH PRONY TERM
  ! (63:68) : VISCO DEV PK2 STRESS FOR 10TH PRONY TERM

CONTAINS

  ! ************************************************************************* !

  SUBROUTINE SHIFTFAC(DTIME, TIME, NPROP, PROPS, TEMPOLD, DTEMP, F, &
       NSTATEV, STATEV)
    ! ----------------------------------------------------------------------- !
    ! INITIALIZE THERMAL VISCOELASTIC VARIABLES AT BEGINNING OF STEP
    ! ----------------------------------------------------------------------- !
    INTEGER :: NPROP, NSTATEV
    REAL(KIND=DP), INTENT(IN) :: DTIME, TIME, PROPS(NPROP), TEMPOLD, DTEMP, F(3,3)
    REAL(KIND=DP), INTENT(INOUT) :: STATEV(NSTATEV)
    REAL(KIND=DP) :: WLF_C1, WLF_C2, WLF_TREF, AT, LOGAT, ANAVG, TIMEAVG
    REAL(KIND=DP) :: TEMPDIFF, TEMPAVG, TEMP
    LOGICAL :: TS_FLAG
    ! ----------------------------------------------------------------------- !

    ! RETRIEVE THE WLF PARAMETERS - THERMAL ANALYSIS
    WLF_C1   = PROPS(IPWC1)
    WLF_C2   = PROPS(IPWC2)
    WLF_TREF = PROPS(IPWTR)

    TEMP = TEMPOLD + DTEMP

    ! COMPUTE THE NUMERICAL SHIFT FUNCTION AT THE TIMESTEP MIDPOINT
    TIMEAVG = TIME - HALF * DTIME
    CALL NUMERICAL_SHIFT(TIMEAVG, ANAVG)

    ! EVALUATE THE WLF SHIFT FACTOR AT THE AVERAGE TEMP OF THE STEP
    AT = ONE ! DEFAULT SHIFT FACTOR
    TS_FLAG = .TRUE.
    IF (TS_FLAG) THEN
       TEMPAVG = HALF * (TEMPOLD + TEMP)
       TEMPDIFF = TEMPAVG - WLF_TREF
       IF (ABS(WLF_C1) < P1M10) THEN
          AT = ONE
       ELSE IF (WLF_C2 + TEMPDIFF <= P1M30) then
          AT = P1M30
       ELSE
          LOGAT = WLF_C1 * TEMPDIFF / (WLF_C2 + TEMPDIFF)
          IF (LOGAT > THIRTY) THEN
             AT = BIGNUM
          ELSE IF (LOGAT < -THIRTY) THEN
             AT = P1M30
          ELSE
             AT = TEN ** LOGAT
          END IF
       END IF
    END IF

    ! STORE THE NUMERICAL SHIFT FACTOR FOR WLF
    STATEV(1) = ONE / AT
    STATEV(2) = ANAVG

    RETURN
  END SUBROUTINE SHIFTFAC

  ! ************************************************************************* !

  SUBROUTINE NUMERICAL_SHIFT(T, A)
    REAL(KIND=DP), INTENT(IN) :: T
    REAL(KIND=DP), INTENT(OUT) :: A
    A = ONE
  END SUBROUTINE NUMERICAL_SHIFT

  ! ************************************************************************* !

  SUBROUTINE VISCORELAX(DTIME, TIME, TEMPOLD, DTEMP, NPROP, PROPS, F, &
       NSTATEV, STATEV, SIGO, SIG, CFAC) BIND(C, NAME='viscorelax_')
    ! ----------------------------------------------------------------------- !
    ! VISCOELASTIC RELAXATION. THIS ROUTINE COMPUTES THE DEVIATORIC PART OF
    ! THE CAUCHY STRESS AT THE END OF THE CURRENT TIME STEP.
    !
    ! PARAMETERS
    ! ----------
    ! SIGO : INSTANTANEOUS CAUCHY STRESS AT END OF TIME STEP (PK2 STRESS AT
    !        END OF TIME STEP WITH REFERENCE STATE EQUAL TO CONFIGURATION AT
    !        END OF TIME STEP)
    !    F : TOTAL DEFORMATION GRADIENT TENSOR TO CURRENT CONFIGURATION
    !        (INCLUDES BOTH MECHANICAL & THERMAL PARTS)
    !  SIG : CAUCHY STRESS AT THE END OF TIME STEP (PK2 STRESS AT END OF TIME
    !        STEP WITH REFERENCE STATE EQUAL TO CONFIGURATION AT END OF TIME
    !        STEP) INCLUDING ALL VISCOELASTIC EFFECTS
    !
    ! NOTES
    ! -----
    ! 1) THIS PROCEDURE IS MEANT TO BE COMPATIBLE WITH ARBITRARY ELASTIC
    !    MODELS. ACOORDINGLY, THE PROPERTIES NEEDED IN THIS PROCEDURE ARE
    !    LIKELY A SUBSET OF THE MATERIAL'S PROPS ARRAY. SO, THIS PROCEDURE
    !    SHOULD BE CALLED AS
    !
    !    CALL VISCORELAX(DTIME, TEMPOLD, DTEMP, 24, PROPS(I), F, ...
    !
    !    WHERE I IS THE INTEGER LOCATION OF WLF C1 IN THE MATERIAL'S PROPS
    !    ARRAY PRESUMING, OF COURSE, THAT THE REST OF THE PROPS ARE ORDERED
    !    ACCORDINGLY
    !
    ! 2) LIKEWISE, THE STATEV ARRAY IS LIKELY JUST A SUBSET OF THE MATERIAL'S
    !    STATEV ARRAY. SO, THE REST OF THE CALL SHOULD LOOK LIKE
    !
    !    ..., 68, STATEV(L), STRESSO, STRESS)
    !
    !    WHERE L IS THE INTEGER LOCATION OF THE FIRST STATE VARIABLE OF THIS
    !    PROCEDURE (THE TEMPERATURE AVERAGED WLF SHIFT)
    !
    ! 3) PROCEDURE CURRENTLY SUPPORTS ONLY SHEAR RELAXATION
    ! ----------------------------------------------------------------------- !
    INTEGER(c_int), INTENT(IN) :: NPROP, NSTATEV
    REAL(c_double), INTENT(IN) :: DTIME, TIME, TEMPOLD, DTEMP, PROPS(NPROP)
    REAL(c_double), INTENT(IN) :: SIGO(NTENS), F(3,3)
    REAL(c_double), INTENT(OUT) :: SIG(NTENS)
    REAL(c_double), INTENT(INOUT) :: STATEV(NSTATEV), CFAC(2)
    INTEGER :: J, K, L
    REAL(KIND=DP) :: PK2ODEV(NTENS), PK2O(NTENS), C(NTENS)
    REAL(KIND=DP) :: RATIO, E, S, DTRED, PRESSURE
    REAL(KIND=DP) :: SODEV(NTENS), SDEV(NTENS), PK2DEV(NTENS), TR
    CHARACTER(LEN=120) :: JNKSTR
    INTEGER :: INTV(2)
    REAL(KIND=DP) :: REALV(1)
    CHARACTER(LEN=8) :: CHARV(1)
    ! ----------------------------------------------------------------------- !

    CFAC = ZERO

    IF (NSTATEV /= NSDV) THEN
       JNKSTR = "VISCORELAX: EXPECTED NSTATEV=%I GOT %I"
       INTV(1) = NSDV
       INTV(2) = NSTATEV
       print*,nsdv
       print*,nstatev
       print*,statev
       CALL VISCOERR(-3, JNKSTR, INTV, REALV, CHARV)
    END IF
    IF (NPROP /= NUI) THEN
       JNKSTR = "VISCORELAX: EXPECTED NPROP=%I GOT %I"
       INTV(1) = NUI
       INTV(2) = NPROP
       CALL VISCOERR(-3, JNKSTR, INTV, REALV, CHARV)
    END IF

    ! GET THE SHIFT FACTORS (STORED IN STATEV)
    CALL SHIFTFAC(DTIME, TIME, NPROP, PROPS, TEMPOLD, DTEMP, F, NSTATEV, STATEV)

    ! CHANGE REFERENCE STATE ON SODEV FROM CONFIGURATION AT END OF CURRENT
    ! TIME STEP TO INITIAL CONFIGURATION
    C = ASARRAY(MATMUL(TRANSPOSE(F), F))
    PK2O = PULL(F, SIGO)
    PK2ODEV = DEV(PK2O, C)

    ! REDUCED TIME STEP
    DTRED = DTIME / STATEV(2) / STATEV(1)

    ! LOOP OVER THE PRONY TERMS
    DO K = 1,10
       J = (K - 1) * 6
       ! COMPUTE NEEDED VISCOELASTIC FACTORS
       RATIO = DTRED / PROPS(IPT+K)
       E = EXP(-RATIO)

       IF (RATIO > P1M3) THEN
          ! EXPLICIT CALCULATION OF (1-EXP(-RATIO))/RATIO
          S = (ONE - E) / RATIO
       ELSE
          ! TAYLOR SERIES CALCULATION OF (1 - EXP(-RATIO))/RATIO
          S = ONE - HALF * RATIO + P6TH * RATIO ** 2
       END IF

       ! UPDATE THE VISCOELASTIC STATE VARIABLE HISTORY FOR KTH PRONY TERM
       DO L = 1,6
          STATEV(8+J+L) = E * STATEV(8+J+L) &
                        + PROPS(IPGOO+K) * (S - E) * STATEV(2+L) &
                        + PROPS(IPGOO+K) * (ONE - S) * PK2ODEV(L)
       END DO
       CFAC(1) = CFAC(1) + (ONE - S) * PROPS(IPGOO+K)
    END DO

    ! COMPUTE DECAYING DEVIATORIC STRESS
    PK2DEV = ZERO
    DO L=1,6
       DO K=1,10
          J = (K - 1) * 6
          PK2DEV(L) = PK2DEV(L) + STATEV(8+J+L)
       END DO
    END DO

    ! CHANGE REFERENCE STATE ON DECAYING PORTION OF DEVIATORIC STRESS FROM
    ! INITIAL CONFIGURATION TO CONFIGURATION AT END OF CURRENT TIME STEP
    SDEV = PUSH(F, PK2DEV)

    ! ELIMINATE THE PRESSURE ARISING FROM THE REFERENCE STATE CHANGES USED IN
    ! COMPUTING THE DECAYING PORTION OF THE DEVIATORIC STRESS
    TR = SUM(SDEV(1:3))
    IF (ABS(TR) > MACHINE_EPSILON) THEN
       SDEV(1) = SDEV(1) - TR / THREE
       SDEV(2) = SDEV(2) - TR / THREE
       SDEV(3) = -(SDEV(1) + SDEV(2))
    END IF

    ! COMPUTE TOTAL DEVIATORIC STRESS
    SODEV = DEV(SIGO)
    SDEV = SODEV - SDEV

    ! TOTAL STRESS
    SIG = SDEV
    PRESSURE = -P3RD * SUM(SIGO(1:3))
    SIG(1:3) = SIG(1:3) - PRESSURE

    ! INSTANTANEOUS DEVIATORIC STRESS WITH ORIGINAL CONFIGURATION AS REFERENCE
    ! STATE
    STATEV(3:8) = PK2ODEV

    RETURN
  END SUBROUTINE VISCORELAX

  ! ************************************************************************* !

  SUBROUTINE VISCOINI(NPROP, PROPS, NSTATEV, STATEV) BIND(C, NAME='viscoini_')
    ! ----------------------------------------------------------------------- !
    ! INITIALIZE RELAXATION STATE DEPENDENT VARIABLES
    ! ----------------------------------------------------------------------- !
    INTEGER(c_int), INTENT(IN) :: NPROP, NSTATEV
    REAL(c_double), INTENT(IN) :: PROPS(NPROP)
    REAL(c_double), INTENT(INOUT) :: STATEV(NSTATEV)
    ! ----------------------------------------------------------------------- !
    STATEV = ZERO
    ! STATE DEPDENDENT VARIABLES
    ! ----- ---------- ---------
    ! (1) : WLF AVERAGE SHIFT FACTOR
    ! (2) : AVERAGE NUMERICAL SHIFT FACTOR
    STATEV(1) = ONE
    STATEV(2) = ONE
    RETURN
  END SUBROUTINE VISCOINI

  ! ************************************************************************* !

  SUBROUTINE PROPCHECK(NPROP, PROPS) BIND(C, NAME='propcheck_')
    ! ----------------------------------------------------------------------- !
    ! CHECK PROPERTY ARRAY FOR VISCOELASTIC MODEL
    ! ----------------------------------------------------------------------- !
    INTEGER(c_int), INTENT(IN) :: NPROP
    REAL(c_double), INTENT(INOUT) :: PROPS(NPROP)
    CHARACTER(LEN=120) :: STR1, STR2, STR3
    REAL(KIND=DP) :: PSUM
    INTEGER :: I, J, IFLG
    INTEGER :: INTV(1)
    REAL(KIND=DP) :: REALV(1)
    CHARACTER(LEN=8) :: CHARV(1)
    ! ----------------------------------------------------------------------- !

    ! CHECK SUM OF PRONY SERIES COEFFICIENTS
    PSUM  = SUM(PROPS(IPGOO:IPG10))
    IF (ANY(PROPS(IPGOO:IPG10) < 0.)) THEN
       WRITE(STR1,'(A)' ) 'EXPECTED ALL SHEAR PRONY SERIES COEFFICIENTS > 0'
       CALL VISCOERR(-3, STR1, INTV, REALV, CHARV)
    END IF

    IF (ABS(PSUM) < P1M10) THEN
       IFLG = -1
       WRITE(STR1,'(A)' ) 'SUM OF NORMALIZED SHEAR PRONY SERIES COEFFICIENTS'
       WRITE(STR2,'(A)' ) 'INCLUDING INFINITY TERM IS ZERO. NORMALIZED INFINITY'
       WRITE(STR3,'(A)' ) 'COEFFICIENT SET TO ONE FOR ELASTIC RESPONSE.'
       CALL VISCOERR(-1, STR1, INTV, REALV, CHARV)
       CALL VISCOERR(-1, STR2, INTV, REALV, CHARV)
       CALL VISCOERR(-1, STR3, INTV, REALV, CHARV)
       PROPS(IPGOO) = ONE

    ELSE IF (ABS(PSUM - ONE) > P1M3) THEN
       WRITE(STR1, '(A)') 'EXPECTED SUM OF NORMALIZED SHEAR PRONY SERIES'
       WRITE(STR2, '(A,1PE15.7)') &
            'COEFFICIENTS INCLUDING INF TERM TO BE 1, GOT ', PSUM
       IF (ABS(PSUM - ONE) < POINT03) THEN
          IFLG = -1
       ELSE
          IFLG = -3
       END IF
       CALL VISCOERR(-1, STR1, INTV, REALV, CHARV)
       CALL VISCOERR(IFLG, STR2, INTV, REALV, CHARV)
    END IF

    ! VERIFY THAT ALL RELAXATION TIMES ARE POSITIVE
    J = 1
    DO I = IPT01, IPT10
       IF (PROPS(I) <= ZERO) THEN
          WRITE(STR1,'(A,I3,A)' ) &
               'SHEAR RELAXATION TIME TERM ',J,' <=0, SETTING TO 1'
          CALL VISCOERR(-1, STR1, INTV, REALV, CHARV)
          PROPS(I) = ONE
       END IF
       J = J + 1
    END DO

    RETURN
  END SUBROUTINE PROPCHECK

  ! ************************************************************************* !

  FUNCTION PUSH(F, A)
    ! ----------------------------------------------------------------------- !
    ! PUSH TRANSFORMATION [A'] = 1/DET[F] [F] [A] [F]^T
    ! ----------------------------------------------------------------------- !
    REAL(KIND=DP), INTENT(IN) :: F(3,3), A(6)
    REAL(KIND=DP) :: PUSH(6)
    REAL(KIND=DP) :: JAC

    ! DET[F]
    JAC = F(1,1) * (F(2,2) * F(3,3) - F(2,3) * F(3,2)) &
        + F(1,2) * (F(2,3) * F(3,1) - F(2,1) * F(3,3)) &
        + F(1,3) * (F(2,1) * F(3,2) - F(2,2) * F(3,1))

    ! Compute [C]=( 1/det[F] ) [F] [A] [F]^T
    PUSH(1) = (A(1)*F(1,1)*F(1,1) +     A(2)*F(1,2)*F(1,2) &
            + A(3)*F(1,3)*F(1,3) + TWO*A(4)*F(1,1)*F(1,2) &
            + TWO*A(5)*F(1,2)*F(1,3) + TWO*A(6)*F(1,1)*F(1,3)) / JAC

    PUSH(2) = (A(1)*F(2,1)*F(2,1) + A(2)*F(2,2)*F(2,2) &
            + A(3)*F(2,3)*F(2,3) + TWO*A(4)*F(2,2)*F(2,1) &
            + TWO*A(5)*F(2,2)*F(2,3) + TWO*A(6)*F(2,3)*F(2,1)) / JAC

    PUSH(3) = (A(1)*F(3,1)*F(3,1) + A(2)*F(3,2)*F(3,2) &
            + A(3)*F(3,3)*F(3,3) + TWO*A(4)*F(3,1)*F(3,2) &
            + TWO*A(5)*F(3,3)*F(3,2) + TWO*A(6)*F(3,3)*F(3,1)) / JAC

    PUSH(4) = (A(1)*F(1,1)*F(2,1) + A(2)*F(2,2)*F(1,2) &
            + A(3)*F(2,3)*F(1,3) + A(4)*(F(1,1)*F(2,2) + F(1,2)*F(2,1)) &
            + A(5)*( F(1,2)*F(2,3) + F(2,2)*F(1,3)) &
            + A(6)*( F(1,1)*F(2,3) + F(2,1)*F(1,3))) / JAC

    PUSH(5) = (A(1)*F(3,1)*F(2,1) + A(2)*F(2,2)*F(3,2) &
            + A(3)*F(3,3)*F(2,3) + A(4)*(F(2,2)*F(3,1) + F(2,1)*F(3,2)) &
            + A(5)*( F(2,2)*F(3,3) + F(2,3)*F(3,2)) &
            + A(6)*( F(2,3)*F(3,1) + F(3,3)*F(2,1))) / JAC

    PUSH(6) = (A(1)*F(1,1)*F(3,1) + A(2)*F(1,2)*F(3,2) &
            + A(3)*F(3,3)*F(1,3) + A(4)*( F(1,2)*F(3,1) + F(1,1)*F(3,2)) &
            + A(5)*( F(3,3)*F(1,2) + F(3,2)*F(1,3)) &
            + A(6)*( F(1,1)*F(3,3) + F(3,1)*F(1,3))) / JAC

    RETURN
  END FUNCTION PUSH

  ! ************************************************************************* !

  FUNCTION INV_6X1(A)
    ! ----------------------------------------------------------------------- !
    ! 3x3 MATRIX INVERSE (STORED AS 6x1 ARRAY)
    ! ----------------------------------------------------------------------- !
    REAL(KIND=DP) :: INV_6X1(6)
    REAL(KIND=DP), INTENT(IN) :: A(6)
    REAL(KIND=DP) :: DETA
    ! ----------------------------------------------------------------------- !
    ! DET[A]
    DETA = A(1)*( A(2)*A(3) - A(5)*A(5) ) &
         + A(4)*( A(5)*A(6) - A(4)*A(3) ) &
         + A(6)*( A(4)*A(5) - A(2)*A(6) )

    ! [A]^-1
    INV_6X1(1) = (A(2)*A(3) - A(5)*A(5)) / DETA
    INV_6X1(2) = (A(1)*A(3) - A(6)*A(6)) / DETA
    INV_6X1(3) = (A(1)*A(2) - A(4)*A(4)) / DETA
    INV_6X1(4) = (A(5)*A(6) - A(3)*A(4)) / DETA
    INV_6X1(5) = (A(4)*A(6) - A(1)*A(5)) / DETA
    INV_6X1(6) = (A(4)*A(5) - A(2)*A(6)) / DETA
  END FUNCTION INV_6X1

  ! ************************************************************************* !

  FUNCTION INV_3x3(A)
    ! ----------------------------------------------------------------------- !
    ! 3x3 MATRIX INVERSE
    ! ----------------------------------------------------------------------- !
    REAL(KIND=DP) :: INV_3x3(3,3)
    REAL(KIND=DP), INTENT(IN) :: A(3,3)
    REAL(KIND=DP) :: DETA
    ! ----------------------------------------------------------------------- !
    ! DET[A]
    DETA = A(1,1) * (A(2,2) * A(3,3) - A(2,3) * A(3,2)) &
         + A(1,2) * (A(2,3) * A(3,1) - A(2,1) * A(3,3)) &
         + A(1,3) * (A(2,1) * A(3,2) - A(2,2) * A(3,1))
    ! [A]^-1
    INV_3x3(1,1) = (A(2,2)*A(3,3) - A(2,3)*A(3,2)) / DETA
    INV_3x3(2,2) = (A(1,1)*A(3,3) - A(3,1)*A(1,3)) / DETA
    INV_3x3(3,3) = (A(1,1)*A(2,2) - A(1,2)*A(2,1)) / DETA
    INV_3x3(1,2) = (A(3,2)*A(1,3) - A(3,3)*A(1,2)) / DETA
    INV_3x3(2,3) = (A(2,1)*A(1,3) - A(1,1)*A(2,3)) / DETA
    INV_3x3(3,1) = (A(2,1)*A(3,2) - A(2,2)*A(3,1)) / DETA
    INV_3x3(2,1) = (A(2,3)*A(3,1) - A(3,3)*A(2,1)) / DETA
    INV_3x3(3,2) = (A(1,2)*A(3,1) - A(1,1)*A(3,2)) / DETA
    INV_3x3(1,3) = (A(1,2)*A(2,3) - A(2,2)*A(1,3)) / DETA
    RETURN
  END FUNCTION INV_3x3

  ! ************************************************************************* !

  FUNCTION PULL(F, A)
    ! ----------------------------------------------------------------------- !
    ! PULL TRANSFORMATION [A] = DET[F] [F]^-1 [A'] [F]^-T
    ! ----------------------------------------------------------------------- !
    REAL(KIND=DP) :: PULL(6)
    REAL(KIND=DP), INTENT(IN) :: F(3,3), A(6)
    REAL(KIND=DP) :: FINV(3,3)
    ! ----------------------------------------------------------------------- !
    FINV = INV_3x3(F)
    PULL = PUSH(FINV, A)
  END FUNCTION PULL

  ! ************************************************************************* !

  FUNCTION ASARRAY(AIN)
    ! ----------------------------------------------------------------------- !
    ! RETURN SECOND ORDER SYMMETRIC TENSOR STORED AS 6x1
    ! ----------------------------------------------------------------------- !
    REAL(KIND=DP) :: ASARRAY(6)
    REAL(KIND=DP), INTENT(IN) :: AIN(3,3)
    REAL(KIND=DP) :: A(3,3)
    ! ----------------------------------------------------------------------- !
    A = HALF * (AIN + TRANSPOSE(AIN))
    ASARRAY(1) = A(1,1); ASARRAY(4) = A(1,2); ASARRAY(6) = A(1,3)
                         ASARRAY(2) = A(2,2); ASARRAY(5) = A(2,3)
                                              ASARRAY(3) = A(3,3)
  END FUNCTION ASARRAY

  ! ************************************************************************* !

  FUNCTION DEV(A, METRIC)
    ! ----------------------------------------------------------------------- !
    ! DEVIATORIC PART OF SECOND ORDER SYMMETRIC TENSOR
    ! ----------------------------------------------------------------------- !
    REAL(KIND=DP) :: DEV(6)
    REAL(KIND=DP), INTENT(IN) :: A(6)
    REAL(KIND=DP), INTENT(IN), OPTIONAL :: METRIC(6)
    REAL(KIND=DP) :: ISOA(6)
    ! ----------------------------------------------------------------------- !
    IF (PRESENT(METRIC)) THEN
       ISOA = ISO(A, METRIC)
    ELSE
       ISOA = ISO(A)
    END IF
    DEV = A - ISOA
  END FUNCTION DEV

  ! ************************************************************************* !

  FUNCTION ISO(A, METRIC)
    ! ----------------------------------------------------------------------- !
    ! ISOTROPIC PART OF SECOND ORDER SYMMETRIC TENSOR
    ! ----------------------------------------------------------------------- !
    REAL(KIND=DP) :: ISO(6)
    REAL(KIND=DP), INTENT(IN) :: A(6)
    REAL(KIND=DP), INTENT(IN), OPTIONAL :: METRIC(6)
    REAL(KIND=DP) :: X(6), I(6)
    REAL(KIND=DP), PARAMETER :: W(6)=(/ONE,ONE,ONE,TWO,TWO,TWO/)
    ! ----------------------------------------------------------------------- !
    IF (PRESENT(METRIC)) THEN
       X = INV_6X1(METRIC)
       I = METRIC
    ELSE
       X = I6
       I = I6
    END IF
    ISO = P3RD * SUM(W * A * I) * X
  END FUNCTION ISO

END MODULE VISCO

! *************************************************************************** !

SUBROUTINE VISCOERR(I, MSG, INTV, REALV, CHARV)
  ! BORROWS FROM ABAQUS' STDB_ABQERR. FOR USE IN ABAQUS, JUST REPLACE CODE
  ! BELOW WITH CALL TO STDB_ABQERR WITH SAME ARGUMENTS AS THIS PROCEDURE
  IMPLICIT NONE
  INTEGER :: I
  CHARACTER(LEN=120) :: MSG
  INTEGER :: INTV(1), IDUM
  REAL(8) :: REALV(1), RDUM
  CHARACTER(LEN=8) :: CHARV(1), CDUM(1)
  CHARACTER(LEN=200) :: JNKSTR
  EXTERNAL LOG_MESSAGE
  EXTERNAL LOG_ERROR
  EXTERNAL LOG_WARNING
  RDUM=REALV(1)
  IDUM=INTV(1)
  CDUM=CHARV(1)
  IF (I == -3) THEN
     JNKSTR = "VISCO: BOMBED: " // MSG
     CALL LOG_ERROR(JNKSTR)
  ELSE IF (I == -1) THEN
     CALL LOG_WARNING(MSG)
  ELSE
     CALL LOG_MESSAGE(MSG)
  END IF
  RETURN
END SUBROUTINE VISCOERR


