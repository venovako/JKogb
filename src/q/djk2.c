#include "djk2.h"

static void qhsvd2d(extended A[static 2][2], extended U[static 2][2], fint info[static 1])
{
  *info = FINT_C(0);
  A[1][0] = A[0][1] = 0.0L;

  if (copysignl(1.0L, A[0][0]) == -1.0L) {
    U[0][0] = -U[0][0];
    U[1][0] = -U[1][0];
    A[0][0] = -A[0][0];
  }

  if (copysignl(1.0L, A[1][1]) == -1.0L) {
    U[0][1] = -U[0][1];
    U[1][1] = -U[1][1];
    A[1][1] = -A[1][1];
  }
}

static void CA1(const extended B[static 2][2], extended A[static 2])
{
  const extended C[2] =
    {
     (fmal(B[1][0], A[1], A[0]) * B[0][0]),
     (fmal(B[0][1], A[0], A[1]) * B[1][1])
    };
  A[0] = C[0];
  A[1] = C[1];
}

static void CA2(const extended B[static 2][2], extended A[static 2][2])
{
  const extended C[2][2] =
    {
     {
      (fmal(B[1][0], A[0][1], A[0][0]) * B[0][0]),
      (fmal(B[0][1], A[0][0], A[0][1]) * B[1][1])
     },
     {
      (fmal(B[1][0], A[1][1], A[1][0]) * B[0][0]),
      (fmal(B[0][1], A[1][0], A[1][1]) * B[1][1])
     }
    };
  A[0][0] = C[0][0];
  A[0][1] = C[0][1];
  A[1][0] = C[1][0];
  A[1][1] = C[1][1];
}

static void AC(extended A[static 2][2], const extended B[static 2][2])
{
  const extended C[2][2] =
    {
     {
      (fmal(A[1][0], B[0][1], A[0][0]) * B[0][0]),
      (fmal(A[1][1], B[0][1], A[0][1]) * B[0][0])
     },
     {
      (fmal(A[0][0], B[1][0], A[1][0]) * B[1][1]),
      (fmal(A[0][1], B[1][0], A[1][1]) * B[1][1])
     }
    };
  A[0][0] = C[0][0];
  A[0][1] = C[0][1];
  A[1][0] = C[1][0];
  A[1][1] = C[1][1];
}

static void qhsvd2u(const bool h, extended A[static 2][2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
{
  *info = FINT_C(0);

  if (h) {
  }
  else {
  }

  if (*info)
    *info = FINT_C(0);
  else {
  }

  qhsvd2d(A, U, info);
  *info = FINT_C(1);
}

static void qhsvd2l(extended A[static 2][2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
{
  *info = FINT_C(0);

  if (*info)
    *info = FINT_C(0);
  else {
  }

  qhsvd2d(A, U, info);
  *info = FINT_C(1);
}

static void qhsvd2g(const bool h, extended A[static 2][2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
{
}

static void qhsvd2s(const fint h, extended A[static 2][2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
{
}

static void qhsvd2_(extended A[static 2][2], const fint J[static 2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
{
  if ((A[0][1] == 0.0L) && (A[1][0] == 0.0L))
    qhsvd2d(A, U, info);
  else /* A general */
    qhsvd2g((J[0] != J[1]), A, U, Z, info);

  if (*info >= FINT_C(0))
    qhsvd2s((J[0] + J[1]), A, U, Z, info);
}

void dhsvd2_(double A[static 2][2], const fint J[static 2], double U[static 2][2], double Z[static 2][2], fint info[static 1])
{
  if (!(fabs(A[0][0]) <= DBL_MAX))
    *info = FINT_C(-1);
  else if (!(fabs(A[0][1]) <= DBL_MAX))
    *info = FINT_C(-2);
  else if (!(fabs(A[1][0]) <= DBL_MAX))
    *info = FINT_C(-3);
  else if (!(fabs(A[1][1]) <= DBL_MAX))
    *info = FINT_C(-4);
  else /* A OK */
    *info = FINT_C(0);

  if (*info)
    return;

  switch (J[0]) {
  case FINT_C(-1):
  case FINT_C(1):
    break;
  default:
    *info = FINT_C(-5);
    return;
  }

  switch (J[1]) {
  case FINT_C(-1):
  case FINT_C(1):
    break;
  default:
    *info = FINT_C(-6);
    return;
  }

  extended A_[2][2] = { { A[0][0], A[0][1] }, { A[1][0], A[1][1] } };
  extended U_[2][2] = { { 1.0L, 0.0L }, { 0.0L, 1.0L } };
  extended Z_[2][2] = { { 1.0L, 0.0L }, { 0.0L, 1.0L } };

  qhsvd2_(A_, J, U_, Z_, info);

  if (*info >= FINT_C(0)) {
    A[0][0] = (double)(A_[0][0]);
    A[0][1] = (double)(A_[0][1]);
    A[1][0] = (double)(A_[1][0]);
    A[1][1] = (double)(A_[1][1]);

    U[0][0] = (double)(U_[0][0]);
    U[0][1] = (double)(U_[0][1]);
    U[1][0] = (double)(U_[1][0]);
    U[1][1] = (double)(U_[1][1]);

    Z[0][0] = (double)(Z_[0][0]);
    Z[0][1] = (double)(Z_[0][1]);
    Z[1][0] = (double)(Z_[1][0]);
    Z[1][1] = (double)(Z_[1][1]);
  }
}
