subroutine xlaev2(a, b, c, rt1, rt2, cs1, sn1)

  implicit none

!  -- LAPACK auxiliary routine (version 3.2) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

!     .. Scalar Arguments ..
  real(kind=16), intent(in) :: a, b, c
  real(kind=16), intent(out) :: cs1, rt1, rt2, sn1

!  Purpose
!  =======
!
!  XLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!  eigenvector for RT1, giving the decomposition
!
!     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
!     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
!
!  Arguments
!  =========
!
!  A       (input) REAL(KIND=16)
!          The (1,1) element of the 2-by-2 matrix.
!
!  B       (input) REAL(KIND=16)
!          The (1,2) element and the conjugate of the (2,1) element of
!          the 2-by-2 matrix.
!
!  C       (input) REAL(KIND=16)
!          The (2,2) element of the 2-by-2 matrix.
!
!  RT1     (output) REAL(KIND=16)
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) REAL(KIND=16)
!          The eigenvalue of smaller absolute value.
!
!  CS1     (output) REAL(KIND=16)
!  SN1     (output) REAL(KIND=16)
!          The vector (CS1, SN1) is a unit right eigenvector for RT1.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  CS1 and SN1 are accurate to a few ulps barring over/underflow.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================

!     .. Parameters ..
  real(kind=16), parameter :: one = +1.0e+0_16, two = +2.0e+0_16, zero = +0.0e+0_16, half=+0.5e+0_16

!     .. Local Scalars ..
  integer :: sgn1, sgn2
  real(kind=16) :: ab, acmn, acmx, acs, adf, cs, ct, df, rt, sm, tb, tn

!     .. Executable Statements ..

!     Compute the eigenvalues

  sm = a + c
  df = a - c
  adf = abs(df)
  tb = b + b
  ab = abs(tb)
  if (abs(a) .gt. abs(c)) then
     acmx = a
     acmn = c
  else
     acmx = c
     acmn = a
  end if
  if (adf .gt. ab) then
     rt = adf * sqrt(one + ( ab / adf )**2)
  else if (adf .lt. ab) then
     rt = ab * sqrt(one + ( adf / ab )**2)
  else

!        Includes case AB=ADF=0

     rt = ab*sqrt(two)
  end if
  if (sm .lt. zero) then
     rt1 = half * (sm - rt)
     sgn1 = -1

!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.

     rt2 = (acmx / rt1) * acmn - (b / rt1) * b
  else if (sm .gt. zero) then
     rt1 = half * (sm + rt)
     sgn1 = 1

!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.

     rt2 = (acmx / rt1) * acmn - (b / rt1) * b
  else

!        Includes case RT1 = RT2 = 0

     rt1 = half * rt
     rt2 = -half * rt
     sgn1 = 1
  end if

!     Compute the eigenvector

  if (df .ge. zero) then
     cs = df + rt
     sgn2 = 1
  else
     cs = df - rt
     sgn2 = -1
  end if
  acs = abs(cs)
  if (acs .gt. ab) then
     ct = -tb / cs
     sn1 = one / sqrt(one + ct * ct)
     cs1 = ct * sn1
  else
     if (ab .eq. zero) then
        cs1 = one
        sn1 = zero
     else
        tn = -cs / tb
        cs1 = one / sqrt(one + tn * tn)
        sn1 = tn * cs1
     end if
  end if
  if (sgn1 .eq. sgn2) then
     tn = cs1
     cs1 = -sn1
     sn1 = tn
  end if

end subroutine xlaev2
