    ! in shared/params.F90

    REAL(KIND=DWP), PARAMETER :: CHFIX = 1.25_DWP ! 5/4
    REAL(KIND=DWP), PARAMETER :: SHFIX = 0.75_DWP ! 3/4

    ! in D/dtransf.F90

    ! two hyperbolic cosines should be the same, but...
    ! can also carry different signs if a J-unitary transformation of the form
    ! | -1 0 | and | 1  0 |
    ! |  0 1 | /or | 0 -1 |
    ! has been applied to the right (not at the present)
    IF (H .AND. ((ABS(Z(1,1)) .GT. CHFIX) .OR. (ABS(Z(2,2)) .GT. CHFIX))) THEN
       INFO = INFO + 32
       Z(1,1) = SIGN(CHFIX, Z(1,1))
       Z(2,1) = SIGN(SHFIX, Z(2,1))
       Z(1,2) = SIGN(SHFIX, Z(1,2))
       Z(2,2) = SIGN(CHFIX, Z(2,2))
    ELSE ! explicitly make A diagonal
       A(2,1) = D_ZERO
       A(1,2) = D_ZERO
    END IF

    ! in D/dstep.F90

    IF (IAND(IT(I), 32) .EQ. 0) THEN
       A(Q,P) = D_ZERO
       A(P,Q) = D_ZERO
    END IF
