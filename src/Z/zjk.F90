PROGRAM ZJK
  USE BINIO
  USE ZSTEP
  USE OMP_LIB
  IMPLICIT NONE

  CHARACTER(LEN=FNL,KIND=c_char) :: FN
  INTEGER :: N, N_2, ID_MAG, ID_CMP, ID_TRU, INFO, NN, NT, SZ(3), FD(3)
  TYPE(ZPROC) :: R

  COMPLEX(KIND=DWP), ALLOCATABLE, TARGET :: ZARR(:,:)
  COMPLEX(KIND=DWP), POINTER, CONTIGUOUS :: A(:,:), U(:,:), Z(:,:)
  REAL(KIND=DWP), ALLOCATABLE, TARGET :: DARR(:)
  REAL(KIND=DWP), POINTER, CONTIGUOUS :: S(:)
  INTEGER, ALLOCATABLE, TARGET :: IARR(:)
  INTEGER, POINTER, CONTIGUOUS :: J(:), P(:), Q(:), STEP(:)
  TYPE(AW), ALLOCATABLE, TARGET :: DZ(:)

  IF (.NOT. VERIFY_MIN_MAX(.FALSE.)) STOP 'MIN and/or MAX do NOT handle NaNs properly!'
  CALL SetCtrlC

  CALL READCL(FN, N, N_2, ID_MAG, ID_CMP, ID_TRU, INFO)
  IF (INFO .NE. 0) STOP 'zjk.exe FN N N_2 [ID_MAG [ID_CMP [ID_TRU]]]'
#ifndef NDEBUG
  WRITE (ULOG,'(A,A)')   '    FN=', TRIM(FN)
  WRITE (ULOG,'(A,I11)') '     N=', N
#endif
  IF (N .LE. 1) STOP 'N < 2'

  INFO = JSTEP_LEN(N, N_2)
  N_2 = INFO
#ifndef NDEBUG
  WRITE (ULOG,'(A,I11)') '   N_2=', N_2
#endif
  IF (INFO .LE. 0) STOP 'JSTEP_LEN'

  ! number of threads
  NT = MIN(MAX(1, INT(OMP_GET_MAX_THREADS())), N_2)
#ifndef NDEBUG
  WRITE (ULOG,'(A,I11)') '    NT=', NT
#endif

  CALL ZPROC_INIT(NT, ID_MAG, ID_CMP, ID_TRU, R, INFO)
#ifndef NDEBUG
  WRITE (ULOG,'(A,I11)') 'ID_MAG=', ID_MAG
  WRITE (ULOG,'(A,I11)') 'ID_CMP=', ID_CMP
  WRITE (ULOG,'(A,I11)') 'ID_TRU=', ID_TRU
#endif
  IF (INFO .NE. 0) STOP 'ZPROC_INIT'

  CALL ZOPEN_YJ_RO(FN, N, N, SZ, FD, INFO)
  IF (INFO .NE. 0) THEN
     WRITE (ULOG,'(A,I11)') 'INFO=', INFO
     STOP 'ZOPEN_YJ_RO'
  END IF

  ALLOCATE(ZARR(N, 3 * N))
#ifndef NDEBUG
  !DIR$ VECTOR ALWAYS
  ZARR = Z_ZERO
#endif
  A => ZARR(:,1:N)
  U => ZARR(:,N+1:2*N)
  Z => ZARR(:,2*N+1:3*N)

  ALLOCATE(DARR(N))
#ifndef NDEBUG
  !DIR$ VECTOR ALWAYS
  DARR = D_ZERO
#endif
  S => DARR(1:N)

  ALLOCATE(IARR(N * N + N_2))
#ifndef NDEBUG
  !DIR$ VECTOR ALWAYS
  IARR = 0
#endif
  NN = (N * (N - 1)) / 2
  J => IARR(1:N)
  P => IARR(N+1:NN+N)
  Q => IARR(NN+N+1:2*NN+N)
  STEP => IARR(2*NN+N+1:2*NN+N+N_2)

  ALLOCATE(DZ(NN))
#ifndef NDEBUG
  !DIR$ VECTOR ALWAYS
  DO INFO = 1, NN
     DZ(INFO)%W = D_ZERO
     DZ(INFO)%P = 0
     DZ(INFO)%Q = 0
     DZ(INFO)%B = 0
  END DO
#endif

  CALL ZREAD_YJ(FD, A, J, N, N, SZ, INFO)
  IF (INFO .NE. 0) THEN
     WRITE (ULOG,'(A,I11)') 'INFO=', INFO
     STOP 'ZREAD_YJ'
  END IF
  CALL BCLOSEN(FD, 3)

  CALL R%TRU(N, P, Q, NN, INFO)
  IF (INFO .NE. NN) THEN
     WRITE (ULOG,'(A,I11)') 'INFO=', INFO
     STOP 'R%TRU'
  END IF
#ifndef NDEBUG
  WRITE (ULOG,'(A,I11)') '    NN=', NN
#endif

  FD(1) = GET_THREAD_NS()
  CALL ZSTEP_LOOP(NT, N, U, N, A, N, Z, N, J, NN, P, Q, R, DZ, N_2, STEP, INFO)
  FD(1) = GET_THREAD_NS() - FD(1)
  IF (INFO .GE. 0) THEN
     WRITE (UOUT,'(A,I10,A,F12.6,A)') 'Executed ', INFO, ' steps with transformations in ', (FD(1) * DNS2s), ' s'
     FLUSH(UOUT)
  ELSE ! error
     WRITE (ULOG,'(A,I10,A,F12.6,A)') 'ERROR ', INFO, ' after ', (FD(1) * DNS2s), ' s'
     FLUSH(ULOG)
  END IF

  IF (ALLOCATED(DZ)) DEALLOCATE(DZ)
  STEP => NULL()
  Q => NULL()
  P => NULL()
  J => NULL()
  IF (ALLOCATED(IARR)) DEALLOCATE(IARR)

  CALL ZOPEN_UZS_RW(FN, N, N, SZ, FD, INFO)
  IF (INFO .NE. 0) THEN
     WRITE (ULOG,'(A,I11)') 'INFO=', INFO
     STOP 'ZOPEN_UZS_RW'
  END IF

  CALL ZWRITE_UZS(FD, U, Z, S, N, N, INFO)
  IF (INFO .NE. 0) THEN
     WRITE (ULOG,'(A,I11)') 'INFO=', INFO
     STOP 'ZWRITE_UZS'
  END IF

  CALL BCLOSEN(FD, 3)

  S => NULL()
  IF (ALLOCATED(DARR)) DEALLOCATE(DARR)
  Z => NULL()
  U => NULL()
  A => NULL()
  IF (ALLOCATED(ZARR)) DEALLOCATE(ZARR)

  STOP !'zjk.exe successfully terminated'
CONTAINS
#include "readcl.F90"
#include "bio.F90"
END PROGRAM ZJK
