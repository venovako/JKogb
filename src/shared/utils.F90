MODULE utils
  USE params
  IMPLICIT NONE

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

#ifdef HYPOT
#ifndef USE_GNU
  INTERFACE
     PURE FUNCTION HYPOTwX87(A, B) BIND(C,NAME='HYPOTwX87')
       USE, INTRINSIC :: iso_c_binding
       IMPLICIT NONE
       REAL(KIND=c_double), INTENT(IN), VALUE :: A, B
       REAL(KIND=c_double) :: HYPOTwX87
     END FUNCTION HYPOTwX87
  END INTERFACE
#endif
#endif

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
  ELEMENTAL REAL(KIND=DWP) FUNCTION HYPOTwFMA(A, B)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A, B

    REAL(KIND=DWP) :: AA, AB, X, Y, Y_X

    AA = ABS(A)
    AB = ABS(B)
    X = MAX(AA, AB)
    Y = MIN(AA, AB)
    Y_X = MAXD(Y / X, D_ZERO)

    HYPOTwFMA = X * SQRT(Y_X * Y_X + D_ONE)
  END FUNCTION HYPOTwFMA

#ifdef USE_GNU
  ELEMENTAL REAL(KIND=DWP) FUNCTION HYPOTwX87(A, B)
    IMPLICIT NONE
    REAL(KIND=DWP), INTENT(IN) :: A, B

    REAL(KIND=XWP) :: X, Y

    X = REAL(A,XWP)
    Y = REAL(B,XWP)
    X = X * X
    Y = Y * Y
    HYPOTwX87 = REAL(SQRT(X + Y),DWP)
  END FUNCTION HYPOTwX87
#endif
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef ABSZ
  ELEMENTAL REAL(KIND=DWP) FUNCTION ABSwFMA(A)
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A

    ABSwFMA = HYPOTwFMA(REAL(A), AIMAG(A))
  END FUNCTION ABSwFMA

#ifdef USE_GNU
  ELEMENTAL REAL(KIND=DWP) FUNCTION ABSwX87(A)
#else
  PURE REAL(KIND=DWP) FUNCTION ABSwX87(A)
#endif
    IMPLICIT NONE
    COMPLEX(KIND=DWP), INTENT(IN) :: A

    ABSwX87 = HYPOTwX87(REAL(A), AIMAG(A))
  END FUNCTION ABSwX87
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE utils
