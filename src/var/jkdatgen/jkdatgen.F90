PROGRAM JKDATGEN
  USE HDF5
  USE H5LT
  IMPLICIT NONE

  ! command-line parameters
  INTEGER :: LGIX, SEED, N, LDG
  CHARACTER(LEN=256) :: FIL, GRP

  LOGICAL :: FEX
  INTEGER(HID_T) :: FID, GID

  INTEGER :: ISEED(4), NRANK, NPLUS, INFO(2)
  REAL(KIND=8), ALLOCATABLE :: LAM(:), A(:,:), G(:,:)
  INTEGER, ALLOCATABLE :: IPIV(:), JVEC(:), IPL(:), INVP(:)

  CALL READCL(LGIX, SEED, N, LDG, FIL, GRP, INFO)
  IF (INFO(1) .NE. 0) THEN
    PRINT *, INFO(1), INFO(2)
    STOP 'readcl'
  END IF

  ALLOCATE(LAM(N))
  ALLOCATE(A(LDG,N))
  ALLOCATE(G(LDG,N))
  ALLOCATE(IPIV(N))
  ALLOCATE(JVEC(N))
  ALLOCATE(IPL(N))
  ALLOCATE(INVP(N))

  CALL GENDAT(LGIX, SEED, N, LAM, ISEED, NRANK, A, G, LDG, IPIV, JVEC, NPLUS, IPL, INVP, INFO)
  IF (INFO(1) .NE. 0) THEN
    PRINT *, INFO(1), INFO(2)
    STOP 'gendat'
  END IF

  CALL H5OPEN_F(INFO(1))
  IF (INFO(1) .NE. 0) THEN
    PRINT *, INFO(1)
    STOP 'h5open_f'
  END IF

  INQUIRE(FILE=FIL, EXIST=FEX)
  IF (FEX) THEN
    CALL H5FOPEN_F(FIL, H5F_ACC_RDWR_F, FID, INFO(1))
  ELSE
    CALL H5FCREATE_F(FIL, H5F_ACC_TRUNC_F, FID, INFO(1))
  END IF
  IF (INFO(1) .NE. 0) THEN
    PRINT *, INFO(1)
    IF (FEX) THEN
      STOP 'h5open_f'
    ELSE
      STOP 'h5create_f'
    END IF
  END IF

  CALL H5GCREATE_F(FID, GRP, GID, INFO(1))
  IF (INFO(1) .NE. 0) THEN
    PRINT *, INFO(1)
    STOP 'h5gcreate_f'
  END IF

  CALL H5WRDS(GID, N, LAM, ISEED, A, G, LDG, NRANK, IPIV, JVEC, NPLUS, IPL, INVP, INFO)

  DEALLOCATE(INVP)
  DEALLOCATE(IPL)
  DEALLOCATE(JVEC)
  DEALLOCATE(IPIV)
  DEALLOCATE(G)
  DEALLOCATE(A)
  DEALLOCATE(LAM)

  IF (INFO(1) .NE. 0) THEN
    PRINT *, INFO(1), INFO(2)
    STOP 'h5wrds'
  END IF

  CALL H5GCLOSE_F(GID, INFO(1))
  IF (INFO(1) .NE. 0) THEN
    PRINT *, INFO(1)
    STOP 'h5gclose_f'
  END IF

  CALL H5FCLOSE_F(FID, INFO(1))
  IF (INFO(1) .NE. 0) THEN
    PRINT *, INFO(1)
    STOP 'h5fclose_f'
  END IF

  CALL H5CLOSE_F(INFO(1))
  IF (INFO(1) .NE. 0) THEN
    PRINT *, INFO(1)
    STOP 'h5close_f'
  END IF

CONTAINS

  SUBROUTINE READCL(LGIX, SEED, N, LDG, FIL, GRP, INFO)

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: LGIX, SEED, N, LDG
    CHARACTER(LEN=*), INTENT(OUT) :: FIL, GRP
    INTEGER, INTENT(OUT) :: INFO(2)

    CHARACTER(LEN=256) :: CAS

    INFO(1) = 0
    INFO(2) = 0

    IF (COMMAND_ARGUMENT_COUNT() .LT. 6) THEN
      PRINT *, 'jkdatgen.exe LAMGEN SEEDIX N LDG FILE.h5 GROUP [ LAMGEN_PARAMS ]'
      PRINT *, '>> COMMAND LINE (INPUT) ARGUMENTS <<'
      PRINT *, 'LAMGEN : spectrum [\Lambda] generator to use (see below): 1 or 2'
      PRINT *, 'SEEDIX : index of hard-coded pRNG seed (see seedix.f90): 1 or 2'
      PRINT *, 'N      : order of the output matrix: > 0'
      PRINT *, 'LDG    : leading dimension of the output matrix: >= N'
      PRINT *, 'FILE.h5: output HDF5 file (may exist): max 256 chars'
      PRINT *, 'GROUP  : output HDF5 group (must NOT exist): max 256 chars'
      PRINT *, 'LAMGEN | LAMGEN_PARAMS'
      PRINT *, '1      : LAMFILE.txt: max 256 chars, >= N lines [each line = one eigenvalue]'
      PRINT *, '2      : IDIST EPS SCALE'
      PRINT *, ' IDIST : 1 [uniform (0,1)], 2 [uniform(-1,1)], or 3 [normal(\mu=0,\sigma=1)]'
      PRINT *, ' EPS   : pseudorandom \lambda''_i survives iff |\lambda''_i| > EPS'
      PRINT *, ' SCALE : final \lambda_i = \lambda''_i * SCALE'
      PRINT *, '<< OUTPUT DATASETS IN FILE.h5/GROUP >>'
      PRINT *, 'IDADIM : integer(4) { LDG, N, NRANK, NPLUS }'
      PRINT *, ' NRANK : numerical column-rank: <= N [but should fail when rank-deficient]'
      PRINT *, ' NPLUS : number of positive eigenvalues: <= N [+/-0 counts as non-positive]'
      PRINT *, 'LAMBDA : double precision(N): the eigenvalues [prescribed or pseudorandom]'
      PRINT *, 'ISEED  : integer(4): initial seed for (d|x)laran pRNG (see LaPACK dlaran.f)'
      PRINT *, 'A      : double precision(LDG,N): the generated symmetric indefinite matrix'
      PRINT *, 'G      : double precision(LDG,N): the `lower'' Bunch-Parlett factor of'
      PRINT *, '       : P A P^T == G'' Q Q^T J'' Q Q^T G''^T == G J G^T'
      PRINT *, 'IPIV   : integer(N): complete pivoting permutation vector P (see xsybpc.f90)'
      PRINT *, 'JVEC   : integer(N): diag(J) == diag(I_{NPLUS},-I_{NRANK-NPLUS},0_{N-NRANK})'
      PRINT *, 'IPL    : integer(N): Q [J-partitioning permutation vector (see jpart.f90)]'
      PRINT *, 'INVP   : integer(N): Q^T [inverse of IPL (see jpart.f90)]'
      INFO(1) = -1
      INFO(2) = COMMAND_ARGUMENT_COUNT() - 6
      RETURN
    END IF

    CALL GET_COMMAND_ARGUMENT(1, CAS)
    READ (CAS,*) LGIX
    IF (LGIX .LT. 1) THEN
      INFO(1) = 1
      INFO(2) = LGIX
      RETURN
    END IF

    CALL GET_COMMAND_ARGUMENT(2, CAS)
    READ (CAS,*) SEED
    IF (SEED .LT. 1) THEN
      INFO(1) = 2
      INFO(2) = SEED
      RETURN
    END IF

    CALL GET_COMMAND_ARGUMENT(3, CAS)
    READ (CAS,*) N
    IF (N .LT. 1) THEN
      INFO(1) = 3
      INFO(2) = N
      RETURN
    END IF

    CALL GET_COMMAND_ARGUMENT(4, CAS)
    READ (CAS,*) LDG
    IF (LDG .LT. N) THEN
      INFO(1) = 4
      INFO(2) = LDG
      RETURN
    END IF

    CALL GET_COMMAND_ARGUMENT(5, FIL)
    IF (LEN_TRIM(FIL) .EQ. 0) THEN
      INFO(1) = 5
      RETURN
    END IF

    CALL GET_COMMAND_ARGUMENT(6, GRP)
    IF (LEN_TRIM(GRP) .EQ. 0) THEN
      INFO(1) = 6
      RETURN
    END IF

  END SUBROUTINE READCL

  SUBROUTINE H5WRDS(GID, N, LAM, ISEED, A, G, LDG, NRANK, IPIV, JVEC, NPLUS, IPL, INVP, INFO)

    IMPLICIT NONE

    INTEGER(HID_T), INTENT(IN) :: GID
    INTEGER, INTENT(IN) :: N, ISEED(4), LDG, NRANK, IPIV(N), JVEC(N), NPLUS, IPL(N), INVP(N)
    REAL(KIND=8), INTENT(IN) :: LAM(N), A(LDG,N), G(LDG,N)
    INTEGER, INTENT(OUT) :: INFO(2)

    INTEGER :: IDADIM(4)
    INTEGER(HSIZE_T) :: DIMS(2)

    INFO(1) = 0
    INFO(2) = 0

    IDADIM = (/ LDG, N, NRANK, NPLUS /)
    DIMS(1) = 4
1   CALL H5LTMAKE_DATASET_INT_F(GID, 'IDADIM', 1, DIMS, IDADIM, INFO(2))
    IF (INFO(2) .NE. 0) THEN
      INFO(1) = 1
      RETURN
    END IF

    DIMS(1) = N
2   CALL H5LTMAKE_DATASET_DOUBLE_F(GID, 'LAMBDA', 1, DIMS, LAM, INFO(2))
    IF (INFO(2) .NE. 0) THEN
      INFO(1) = 2
      RETURN
    END IF

    DIMS(1) = 4
3   CALL H5LTMAKE_DATASET_INT_F(GID, 'ISEED', 1, DIMS, ISEED, INFO(2))
    IF (INFO(2) .NE. 0) THEN
      INFO(1) = 3
      RETURN
    END IF

    DIMS(1) = LDG
    DIMS(2) = N
4   CALL H5LTMAKE_DATASET_DOUBLE_F(GID, 'A', 2, DIMS, A, INFO(2))
    IF (INFO(2) .NE. 0) THEN
      INFO(1) = 4
      RETURN
    END IF

    DIMS(1) = LDG
    DIMS(2) = N
5   CALL H5LTMAKE_DATASET_DOUBLE_F(GID, 'G', 2, DIMS, G, INFO(2))
    IF (INFO(2) .NE. 0) THEN
      INFO(1) = 5
      RETURN
    END IF

    DIMS(1) = N
6   CALL H5LTMAKE_DATASET_INT_F(GID, 'IPIV', 1, DIMS, IPIV, INFO(2))
    IF (INFO(2) .NE. 0) THEN
      INFO(1) = 6
      RETURN
    END IF

    DIMS(1) = N
7   CALL H5LTMAKE_DATASET_INT_F(GID, 'JVEC', 1, DIMS, JVEC, INFO(2))
    IF (INFO(2) .NE. 0) THEN
      INFO(1) = 7
      RETURN
    END IF

    DIMS(1) = N
8   CALL H5LTMAKE_DATASET_INT_F(GID, 'IPL', 1, DIMS, IPL, INFO(2))
    IF (INFO(2) .NE. 0) THEN
      INFO(1) = 8
      RETURN
    END IF

    DIMS(1) = N
9   CALL H5LTMAKE_DATASET_INT_F(GID, 'INVP', 1, DIMS, INVP, INFO(2))
    IF (INFO(2) .NE. 0) THEN
      INFO(1) = 9
      RETURN
    END IF

  END SUBROUTINE H5WRDS

END PROGRAM JKDATGEN
