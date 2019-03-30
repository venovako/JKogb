MODULE BLAS_UTILS
#ifdef USE_MKL
#ifdef USE_INTEL
#ifdef USE_X200
  USE IFCORE
#endif
  USE MKL_SERVICE
#endif
#endif
  USE PARAMS
  USE VN_BLAS_F
  USE OMP_LIB
  IMPLICIT NONE

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION BLAS_PREPARE()
    IMPLICIT NONE
#ifdef USE_MKL
#ifdef USE_INTEL
#ifdef USE_X200
    IF (FOR_GET_HBW_AVAILABILITY() .NE. FOR_K_HBW_AVAILABLE) THEN
       BLAS_PREPARE = FOR_SET_FASTMEM_POLICY(FOR_K_FASTMEM_RETRY_WARN)
       BLAS_PREPARE = MKL_SET_MEMORY_LIMIT(MKL_MEM_MCDRAM, 0)
    END IF
#endif
#endif
#endif
    BLAS_PREPARE = VN_BLAS_PREPARE()
  END FUNCTION BLAS_PREPARE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER FUNCTION BLAS_SET_NUM_THREADS(NT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NT
    BLAS_SET_NUM_THREADS = VN_BLAS_SET_NUM_THREADS(NT)
  END FUNCTION BLAS_SET_NUM_THREADS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  PURE INTEGER FUNCTION BLAS_IZAMAX(N, ZX, INCX)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: INCX, N
    COMPLEX(KIND=DWP), INTENT(IN) :: ZX(*)

    REAL(KIND=DWP) :: DMAX, DTMP
    INTEGER :: I, IX

    BLAS_IZAMAX = 0
    IF ((N .LT. 1) .OR. (INCX .LE. 0)) RETURN
    BLAS_IZAMAX = 1
    IF (N .EQ. 1) RETURN

    IF (INCX .EQ. 1) THEN
       DMAX = ABS(ZX(1))
       DO I = 2, N
          DTMP = ABS(ZX(I))
          IF (DTMP .GT. DMAX) THEN
             BLAS_IZAMAX = I
             DMAX = DTMP
          END IF
       END DO
    ELSE
       IX = 1
       DMAX = ABS(ZX(1))
       DO I = 2, N
          IX = IX + INCX
          DTMP = ABS(ZX(IX))
          IF (DTMP .GT. DMAX) THEN
             BLAS_IZAMAX = I
             DMAX = DTMP
          END IF
       END DO
    END IF
  END FUNCTION BLAS_IZAMAX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE BLAS_UTILS
