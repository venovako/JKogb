#include "djk2.h"

static void qhsvd2_(extended A[static 2][2], const fint J[static 2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
{
  if ((A[0][1] == 0.0L) && (A[1][0] == 0.0L))
    /* A diagonal */;
  else
    /* A general */;
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
