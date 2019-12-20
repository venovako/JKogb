#include "djk2.h"

static inline void qswp(extended a[static 1], extended b[static 1])
{
  const extended t = *a;
  *a = *b;
  *b = t;
}

static inline void qassgn1(extended A[static 2], const extended B[static 2])
{
  A[0] = B[0];
  A[1] = B[1];
}

static inline void qassgn2(extended A[static 2][2], const extended B[static 2][2])
{
  A[0][0] = B[0][0];
  A[0][1] = B[0][1];
  A[1][0] = B[1][0];
  A[1][1] = B[1][1];
}

static inline void dassgn2(double A[static 2][2], const extended B[static 2][2])
{
  A[0][0] = (double)(B[0][0]);
  A[0][1] = (double)(B[0][1]);
  A[1][0] = (double)(B[1][0]);
  A[1][1] = (double)(B[1][1]);
}

static inline void CA1(const extended B[static 2][2], extended A[static 2])
{
  const extended C[2] =
    {
     (fmal(B[1][0], A[1], A[0]) / B[0][0]),
     (fmal(B[0][1], A[0], A[1]) / B[1][1])
    };
  qassgn1(A, C);
}

static inline void CA2(const extended B[static 2][2], extended A[static 2][2])
{
  const extended C[2][2] =
    {
     {
      (fmal(B[1][0], A[0][1], A[0][0]) / B[0][0]),
      (fmal(B[0][1], A[0][0], A[0][1]) / B[1][1])
     },
     {
      (fmal(B[1][0], A[1][1], A[1][0]) / B[0][0]),
      (fmal(B[0][1], A[1][0], A[1][1]) / B[1][1])
     }
    };
  qassgn2(A, C);
}

static inline void AC(extended A[static 2][2], const extended B[static 2][2])
{
  const extended C[2][2] =
    {
     {
      (fmal(A[1][0], B[0][1], A[0][0]) / B[0][0]),
      (fmal(A[1][1], B[0][1], A[0][1]) / B[0][0])
     },
     {
      (fmal(A[0][0], B[1][0], A[1][0]) / B[1][1]),
      (fmal(A[0][1], B[1][0], A[1][1]) / B[1][1])
     }
    };
  qassgn2(A, C);
}

static inline void qhsvd2d(extended A[static 2][2], extended U[static 2][2], fint info[static 1])
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

static inline void qhsvd2u(const bool h, extended A[static 2][2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
{
  extended tu = 0.0L, cu = 1.0L, tz = 0.0L, cz = 1.0L;
  *info = FINT_C(0);

  if (h) {
    if (A[1][1] != 0.0L) {
      const extended
        x = (A[1][0] / A[0][0]),
        y = (A[1][1] / A[0][0]);
      extended t2 = ((fabsl(x) <= y) ? (scalbnl(x, 1) * y) : (scalbnl(y, 1) * x));
      t2 /= fmal((y - x), (y + x), 1.0L);
      tu = (!(fabsl(t2) <= LDBL_MAX) ? copysignl(1.0L, t2) : (t2 / (1.0L + sqrtl(fmal(t2, t2, 1.0L)))));
      cu = sqrtl(fmal(tu, tu, 1.0L)); /* 1.0L / */
      tz = fmal(y, tu, -x);
      if (fabsl(tz) >= 1.0L) {
        *info = FINT_C(-8);
        return;
      }
      cz = sqrtl(fmal(-tz, tz, 1.0L)); /* 1.0L / */
    }
    else if (fabsl(A[1][0]) < fabsl(A[0][0])) {
      *info = FINT_C(1);
      tz = -(A[1][0] / A[0][0]);
      if (fabsl(tz) >= 1.0L) {
        *info = FINT_C(-10);
        return;
      }
      cz = sqrtl(fmal(-tz, tz, 1.0L)); /* 1.0L / */
    }
    else {
      *info = FINT_C(-12);
      return;
    }
  }
  else {
    if (A[1][1] != 0.0L) {
      const extended
        x = (A[1][0] / A[0][0]),
        y = (A[1][1] / A[0][0]);
      extended t2 = -((fabsl(x) <= y) ? (scalbnl(x, 1) * y) : (scalbnl(y, 1) * x));
      t2 /= fmal((x - y), (x + y), 1.0L);
      tu = (!(fabsl(t2) <= LDBL_MAX) ? copysignl(1.0L, t2) : (t2 / (1.0L + sqrtl(fmal(t2, t2, 1.0L)))));
      cu = sqrtl(fmal(tu, tu, 1.0L)); /* 1.0L / */
      tz = fmal(y, tu, -x);
      cz = sqrtl(fmal(tz, tz, 1.0L)); /* 1.0L / */
    }
    else {
      *info = FINT_C(1);
      tz = -(A[1][0] / A[0][0]);
      cz = sqrtl(fmal(tz, tz, 1.0L)); /* 1.0L / */
    }
  }

  if (*info)
    *info = FINT_C(0);
  else {
    const extended V[2][2] = { { cu, tu }, { -tu, cu } };
    CA2(V, U);
    CA2(V, A);
  }

  if (h) {
    const extended W[2][2] = { { cz, tz }, { tz, cz } };
    AC(A, W);
    AC(Z, W);
  }
  else {
    const extended W[2][2] = { { cz, -tz }, { tz, cz } };
    AC(A, W);
    AC(Z, W);
  }

  qhsvd2d(A, U, info);
  *info = FINT_C(1);
}

static inline void qhsvd2l(extended A[static 2][2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
{
  extended tu = 0.0L, cu = 1.0L, tz = 0.0L, cz = 1.0L;
  *info = FINT_C(0);

  if (A[0][0] != 0.0L) {
    const extended
      x = (A[0][0] / A[1][1]),
      y = (A[0][1] / A[1][1]);
    extended t2 = -((fabsl(y) <= x) ? (scalbnl(y, 1) * x) : (scalbnl(x, 1) * y));
    t2 /= fmal((x - y), (x + y), 1.0L);
    tu = (!(fabsl(t2) <= LDBL_MAX) ? copysignl(1.0L, t2) : (t2 / (1.0L + sqrtl(fmal(t2, t2, 1.0L)))));
    cu = sqrtl(fmal(tu, tu, 1.0L)); /* 1.0L / */
    tz = -fmal(x, tu, y);
    if (fabsl(tz) >= 1.0L) {
      *info = FINT_C(-13);
      return;
    }
    cz = sqrtl(fmal(-tz, tz, 1.0L)); /* 1.0L / */
  }
  else if (fabsl(A[0][1]) < fabsl(A[1][1])) {
    *info = FINT_C(1);
    tz = -(A[0][1] / A[1][1]);
    if (fabsl(tz) >= 1.0L) {
      *info = FINT_C(-15);
      return;
    }
    cz = sqrtl(fmal(-tz, tz, 1.0L)); /* 1.0L / */
  }
  else {
    *info = FINT_C(-17);
    return;
  }

  if (*info)
    *info = FINT_C(0);
  else {
    const extended V[2][2] = { { cu, tu }, { -tu, cu } };
    CA2(V, U);
    CA2(V, A);
  }

  const extended W[2][2] = { { cz, tz }, { tz, cz } };
  AC(A, W);
  AC(Z, W);

  qhsvd2d(A, U, info);
  *info = FINT_C(1);
}

static inline void qhsvd2g(const bool h, extended A[static 2][2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
{
  extended r = 0.0L,
    c = hypotl(A[0][0], A[0][1]),
    s = hypotl(A[1][0], A[1][1]);

  if (c < s) {
    r = s;
    qswp(&(A[0][0]), &(A[1][0]));
    qswp(&(A[0][1]), &(A[1][1]));
    if (!h) {
      qswp(&(Z[0][0]), &(Z[1][0]));
      qswp(&(Z[0][1]), &(Z[1][1]));
    }
    *info = FINT_C(1);
  }
  else {
    r = c;
    *info = FINT_C(0);
  }

  if (fabsl(A[0][0]) < fabsl(A[0][1])) {
    qswp(&(U[0][0]), &(U[0][1]));
    qswp(&(U[1][0]), &(U[1][1]));
    qswp(&(A[0][0]), &(A[0][1]));
    qswp(&(A[1][0]), &(A[1][1]));
  }

  r = copysignl(r, A[0][0]);
  s = A[0][1] / A[0][0];
  c = sqrtl(fmal(s, s, 1.0L)); /* 1.0L / */

  const extended Q[2][2] = { { c, -s }, { s, c } };
  CA1(Q, A[1]);
  CA2(Q, U);

  if (copysignl(1.0L, r) == -1.0L) {
    U[0][0] = -U[0][0];
    U[1][0] = -U[1][0];
    A[1][0] = -A[1][0];
    r = -r;
  }
  A[0][0] = r;
  if (copysignl(1.0L, A[1][1]) == -1.0L) {
    U[0][1] = -U[0][1];
    U[1][1] = -U[1][1];
    A[1][1] = -A[1][1];
  }
  A[0][1] = 0.0L;

  if (A[1][0] == 0.0L) {
    if (h) {
      if (*info == FINT_C(1)) {
        qswp(&(U[0][0]), &(U[0][1]));
        qswp(&(U[1][0]), &(U[1][1]));
        qswp(&(A[0][0]), &(A[1][1]));
      }
      else
        *info = FINT_C(1);
    }
    else
      *info = FINT_C(1);
  }
  else if (h && (*info == FINT_C(1))) {
    qswp(&(A[0][0]), &(A[1][0]));
    qswp(&(A[0][1]), &(A[1][1]));
    qswp(&(U[0][0]), &(U[0][1]));
    qswp(&(U[1][0]), &(U[1][1]));
    qswp(&(A[0][0]), &(A[0][1]));
    qswp(&(A[1][0]), &(A[1][1]));
    qhsvd2l(A, U, Z, info);
  }
  else
    qhsvd2u(h, A, U, Z, info);
}

static inline void qhsvd2s(const fint h, extended A[static 2][2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
{
  if (((h == FINT_C(2)) && (A[0][0] < A[1][1])) || ((h == FINT_C(-2)) && (A[0][0] > A[1][1]))) {
    qswp(&(U[0][0]), &(U[0][1]));
    qswp(&(U[1][0]), &(U[1][1]));
    qswp(&(A[0][0]), &(A[1][1]));
    qswp(&(Z[0][0]), &(Z[1][0]));
    qswp(&(Z[0][1]), &(Z[1][1]));
  }

  if ((U[0][0] != 1.0L) || (U[0][1] != 0.0L) || (U[1][0] != 0.0L) || (U[1][1] != 1.0L))
    *info += FINT_C(2);
  if ((Z[0][0] != 1.0L) || (Z[0][1] != 0.0L) || (Z[1][0] != 0.0L) || (Z[1][1] != 1.0L))
    *info += FINT_C(4);
}

static inline void qhsvd2_(extended A[static 2][2], const fint J[static 2], extended U[static 2][2], extended Z[static 2][2], fint info[static 1])
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
    dassgn2(A, (const extended (*)[2])A_);
    dassgn2(U, (const extended (*)[2])U_);
    dassgn2(Z, (const extended (*)[2])Z_);
  }
}
