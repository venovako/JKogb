MODULE blas_utils
#ifdef USE_MKL
#ifdef USE_INTEL
#ifdef USE_X200
  USE ifcore
#endif
  USE mkl_service
#endif
#endif
  USE params
  USE vn_blas_f
  !$ USE omp_lib
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

END MODULE blas_utils
