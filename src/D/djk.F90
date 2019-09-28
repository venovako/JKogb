PROGRAM DJK
  USE BINIO
  USE DSTEP
  USE OMP_LIB
  IMPLICIT NONE

  CHARACTER(LEN=FNL,KIND=c_char) :: FN
  INTEGER :: N, N_2, ID_MAG, ID_CMP, ID_TRU, INFO, NN, NM, NT, SZ(3), FD(3)
  TYPE(DPROC) :: R

  REAL(KIND=DWP), ALLOCATABLE, TARGET :: DARR(:,:)
  REAL(KIND=DWP), POINTER, CONTIGUOUS :: A(:,:), U(:,:), Z(:,:), S(:)
  INTEGER, ALLOCATABLE, TARGET :: IARR(:)
  INTEGER, POINTER, CONTIGUOUS :: J(:), P(:), Q(:), STEP(:)
  TYPE(AW), ALLOCATABLE, TARGET :: DZ(:)

  ! IF (.NOT. VERIFY_MIN_MAX(.FALSE.)) STOP 'MIN and/or MAX do NOT handle NaNs properly!'
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

  CALL DPROC_INIT(NT, ID_MAG, ID_CMP, ID_TRU, R, INFO)
#ifndef NDEBUG
  WRITE (ULOG,'(A,I11)') 'ID_MAG=', ID_MAG
  WRITE (ULOG,'(A,I11)') 'ID_CMP=', ID_CMP
  WRITE (ULOG,'(A,I11)') 'ID_TRU=', ID_TRU
#endif
  IF (INFO .NE. 0) STOP 'DPROC_INIT'

  CALL DOPEN_YJ_RO(FN, N, N, SZ, FD, INFO)
  IF (INFO .NE. 0) THEN
     WRITE (ULOG,'(A,I11)') 'INFO=', INFO
     STOP 'DOPEN_YJ_RO'
  END IF

  ALLOCATE(DARR(N, 3 * N + 1))
  A => DARR(:,1:N)
  U => DARR(:,N+1:2*N)
  Z => DARR(:,2*N+1:3*N)
  S => DARR(:,3*N+1)

  ALLOCATE(IARR(N * N + N_2))
  NN = (N * (N - 1)) / 2
  J => IARR(1:N)
  P => IARR(N+1:NN+N)
  Q => IARR(NN+N+1:2*NN+N)
  STEP => IARR(2*NN+N+1:2*NN+N+N_2)

  INFO = MOD(NN, NT)
  IF (INFO .GT. 0) INFO = NT - INFO
  NM = 2 * (NN + INFO)
  ALLOCATE(DZ(NM))
#ifndef NDEBUG
  WRITE (ULOG,'(A,I11)') '    NM=', NM
#endif

  CALL DREAD_YJ(FD, A, J, N, N, SZ, INFO)
  IF (INFO .NE. 0) THEN
     WRITE (ULOG,'(A,I11)') 'INFO=', INFO
     STOP 'DREAD_YJ'
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
  CALL DSTEP_LOOP(NT, N, U, N, A, N, Z, N, J, S, NN, P, Q, R, NM, DZ, N_2, STEP, INFO)
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

  CALL DOPEN_UZS_RW(FN, N, N, SZ, FD, INFO)
  IF (INFO .NE. 0) THEN
     WRITE (ULOG,'(A,I11)') 'INFO=', INFO
     STOP 'DOPEN_UZS_RW'
  END IF

  CALL DWRITE_UZS(FD, U, Z, S, N, N, INFO)
  IF (INFO .NE. 0) THEN
     WRITE (ULOG,'(A,I11)') 'INFO=', INFO
     STOP 'DWRITE_UZS'
  END IF

  CALL BCLOSEN(FD, 3)

  S => NULL()
  Z => NULL()
  U => NULL()
  A => NULL()
  IF (ALLOCATED(DARR)) DEALLOCATE(DARR)

  STOP !'djk.exe successfully terminated'
CONTAINS
#include "readcl.F90"
#include "bio.F90"
END PROGRAM DJK
