MODULE utils
#ifdef USE_INTEL
#ifdef FMAD
  USE, INTRINSIC :: ieee_arithmetic, ONLY: IEEE_FMA !, IEEE_MAX_NUM, IEEE_MIN_NUM
#endif
#endif
  USE params
#ifndef FMAD
#define FMAD(a,b,c) ((a)*(b)+(c))
#endif
  IMPLICIT NONE

#ifndef USE_INTEL
  INTERFACE
     PURE FUNCTION C_FMA(A, B, C) BIND(C,NAME='fma')
       USE, INTRINSIC :: iso_c_binding
       IMPLICIT NONE
       REAL(KIND=c_double), INTENT(IN), VALUE :: A, B, C
       REAL(KIND=c_double) :: C_FMA
     END FUNCTION C_FMA
  END INTERFACE

#ifdef MAXD
  INTERFACE
     PURE FUNCTION C_FMAX(A, B) BIND(C,NAME='fmax')
       USE, INTRINSIC :: iso_c_binding
       IMPLICIT NONE
       REAL(KIND=c_double), INTENT(IN), VALUE :: A, B
       REAL(KIND=c_double) :: C_FMAX
     END FUNCTION C_FMAX
  END INTERFACE
#endif

#ifdef MIND
  INTERFACE
     PURE FUNCTION C_FMIN(A, B) BIND(C,NAME='fmin')
       USE, INTRINSIC :: iso_c_binding
       IMPLICIT NONE
       REAL(KIND=c_double), INTENT(IN), VALUE :: A, B
       REAL(KIND=c_double) :: C_FMIN
     END FUNCTION C_FMIN
  END INTERFACE
#endif
#endif

#ifdef HYPOT
  INTERFACE
     PURE FUNCTION HYPOTwX87(A, B) BIND(C,NAME='HYPOTwX87')
       USE, INTRINSIC :: iso_c_binding
       IMPLICIT NONE
       REAL(KIND=c_double), INTENT(IN), VALUE :: A, B
       REAL(KIND=c_double) :: HYPOTwX87
     END FUNCTION HYPOTwX87
  END INTERFACE
#endif

  ! INTERFACE
  !    PURE FUNCTION DASUM2(A, B) BIND(C,NAME='DASUM2')
  !      USE, INTRINSIC :: iso_c_binding
  !      IMPLICIT NONE
  !      REAL(KIND=c_double), INTENT(IN), VALUE :: A, B
  !      REAL(KIND=c_double) :: DASUM2
  !    END FUNCTION DASUM2
  ! END INTERFACE

  ! INTERFACE
  !    PURE FUNCTION DASUM4(A, B, C, D) BIND(C,NAME='DASUM4')
  !      USE, INTRINSIC :: iso_c_binding
  !      IMPLICIT NONE
  !      REAL(KIND=c_double), INTENT(IN), VALUE :: A, B, C, D
  !      REAL(KIND=c_double) :: DASUM4
  !    END FUNCTION DASUM4
  ! END INTERFACE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef NDEBUG
  INTEGER FUNCTION OPEN_LOG(FN, S)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: FN
    INTEGER, INTENT(IN) :: S

    CHARACTER(LEN=LEN_TRIM(FN)+12) :: F
    INTEGER :: U, I

    WRITE (F,'(A,A,I11.11)') TRIM(FN), '.', S
    OPEN(NEWUNIT=U, IOSTAT=I, FILE=F, STATUS='REPLACE', ACCESS='SEQUENTIAL', ACTION='WRITE')
    IF (I .EQ. 0) THEN
       OPEN_LOG = U
    ELSE ! error
       OPEN_LOG = -1
    END IF
  END FUNCTION OPEN_LOG
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE REAL(KIND=DWP) FUNCTION QUIET_NAN(PAYLOAD)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: PAYLOAD

    QUIET_NAN = TRANSFER(IOR(INT(PAYLOAD,DWP), D_QNAN_MASK), 0.0_DWP)
  END FUNCTION QUIET_NAN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE LOGICAL FUNCTION VERIFY_MIN(STRONG)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: STRONG
#ifdef MIND
    VERIFY_MIN = .TRUE.
#else
    IF (STRONG) THEN
       VERIFY_MIN = (&
            (MIN(QUIET_NAN(-1), 0.0_DWP) .EQ. 0.0_DWP) .AND. &
            (MIN(0.0_DWP, QUIET_NAN(-1)) .EQ. 0.0_DWP))
    ELSE ! weak
       VERIFY_MIN = (MIN(QUIET_NAN(-1), 0.0_DWP) .EQ. 0.0_DWP)
    END IF
#endif
  END FUNCTION VERIFY_MIN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE LOGICAL FUNCTION VERIFY_MAX(STRONG)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: STRONG
#ifdef MAXD
    VERIFY_MAX = .TRUE.
#else
    IF (STRONG) THEN
       VERIFY_MAX = (&
            (MAX(QUIET_NAN(0), -1.0_DWP) .EQ. -1.0_DWP) .AND. &
            (MAX(-1.0_DWP, QUIET_NAN(0)) .EQ. -1.0_DWP))
    ELSE ! weak
       VERIFY_MAX = (MAX(QUIET_NAN(0), -1.0_DWP) .EQ. -1.0_DWP)
    END IF
#endif
  END FUNCTION VERIFY_MAX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE LOGICAL FUNCTION VERIFY_MIN_MAX(STRONG)
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: STRONG

    VERIFY_MIN_MAX = (VERIFY_MIN(STRONG) .AND. VERIFY_MAX(STRONG))
  END FUNCTION VERIFY_MIN_MAX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef HYPOT
  PURE REAL(KIND=DWP) FUNCTION HYPOTwFMA(A, B)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A, B

    REAL(KIND=DWP) :: AA, AB, X, Y, Y_X

    AA = ABS(A)
    AB = ABS(B)
    X = MAX(AA, AB)
    Y = MIN(AA, AB)
#ifdef MAXD
    Y_X = MAXD(Y / X, D_ZERO)
#else
    Y_X = MAX(Y / X, D_ZERO)
#endif
    HYPOTwFMA = X * SQRT(FMAD(Y_X, Y_X, D_ONE))
  END FUNCTION HYPOTwFMA

  ELEMENTAL REAL(KIND=DWP) FUNCTION HYPOTw128(A, B)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A, B

    REAL(KIND=QWP) :: X, Y

    X = REAL(A,QWP)
    Y = REAL(B,QWP)
    X = X * X
    Y = Y * Y
    HYPOTw128 = REAL(SQRT(X + Y),DWP)
  END FUNCTION HYPOTw128
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef ABSZ
  PURE REAL(KIND=DWP) FUNCTION ABSwFMA(A)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A

    ABSwFMA = HYPOTwFMA(REAL(A), AIMAG(A))
  END FUNCTION ABSwFMA

  PURE REAL(KIND=DWP) FUNCTION ABSwX87(A)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A

    ABSwX87 = HYPOTwX87(REAL(A), AIMAG(A))
  END FUNCTION ABSwX87

  ELEMENTAL REAL(KIND=DWP) FUNCTION ABSw128(A)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A

    ABSw128 = HYPOTw128(REAL(A), AIMAG(A))
  END FUNCTION ABSw128
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE utils
