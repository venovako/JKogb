#include "xtransf.h"
#ifdef USE_EXTENDED
#ifndef x0
#define x0 CMPLXL(0.0L, 0.0L)
#else /* x0 */
#error x0 already defined
#endif /* ?x0 */

#ifndef x1
#define x1 CMPLXL(1.0L, 0.0L)
#else /* x1 */
#error x1 already defined
#endif /* ?x1 */

static inline bool zisfinite(const dcomplex z)
{
  return (isfinite(creal(z)) && isfinite(cimag(z)));
}

static inline bool xisreal(const xcomplex x)
{
  return (cimagl(x) == 0.0L);
}

static inline bool xisimag(const xcomplex x)
{
  return (creall(x) == 0.0L);
}

static inline xcomplex xqfma(const extended a, const xcomplex b, const xcomplex c)
{
  return CMPLXL(fmal(a, creall(b), creall(c)), fmal(a, cimagl(b), cimagl(c)));
}

static inline xcomplex xjfma(const extended a, const xcomplex b, const xcomplex c)
{
  return CMPLXL(fmal(-a, cimagl(b), creall(c)), fmal(a, creall(b), cimagl(c)));
}

static inline xcomplex xxfma(const xcomplex a, const xcomplex b, const xcomplex c)
{
  return CMPLXL(fmal(creall(a), creall(b), fmal(-cimagl(a), cimagl(b), creall(c))), fmal(creall(a), cimagl(b), fmal(cimagl(a), creall(b), cimagl(c))));
}

static inline xcomplex xfma(const xcomplex a, const xcomplex b, const xcomplex c)
{
  if (xisreal(a))
    return xqfma(creall(a), b, c);
  if (xisimag(a))
    return xjfma(cimagl(a), b, c);
  return xxfma(a, b, c);
}

static inline xcomplex xxmul(const xcomplex a, const xcomplex b)
{
  return CMPLXL(fmal(creall(a), creall(b), -(cimagl(a) * cimagl(b))), fmal(creall(a), cimagl(b), (cimagl(a) * creall(b))));
}

static inline xcomplex xxdiv(const xcomplex a, const xcomplex b)
{
#ifdef NDEBUG
  return (a / b);
#else /* DEBUG */
  const extended f = cabsl(b);
  return (xxmul(conjl(b / f), a) / f);
#endif /* ?NDEBUG */
}

static inline void xswp(xcomplex a[static 1], xcomplex b[static 1])
{
  const xcomplex t = *a;
  *a = *b;
  *b = t;
}

static inline void xassgn2(xcomplex A[static 2][2], const xcomplex B[static 2][2])
{
  A[0][0] = B[0][0];
  A[0][1] = B[0][1];
  A[1][0] = B[1][0];
  A[1][1] = B[1][1];
}

static inline void zassgn2(dcomplex A[static 2][2], const xcomplex B[static 2][2])
{
  A[0][0] = (dcomplex)(B[0][0]);
  A[0][1] = (dcomplex)(B[0][1]);
  A[1][0] = (dcomplex)(B[1][0]);
  A[1][1] = (dcomplex)(B[1][1]);
}

/* unused
static inline void xassgn1(xcomplex A[static 2], const xcomplex B[static 2])
{
  A[0] = B[0];
  A[1] = B[1];
}

static inline void CA1(const xcomplex B[static 2][2], xcomplex A[static 2])
{
  const xcomplex C[2] =
    {
     (xfma(B[1][0], A[1], A[0]) / creall(B[0][0])),
     (xfma(B[0][1], A[0], A[1]) / creall(B[1][1]))
    };
  xassgn1(A, C);
}
*/

static inline void CA2(const xcomplex B[static 2][2], xcomplex A[static 2][2])
{
  const xcomplex C[2][2] =
    {
     {
      (xfma(B[1][0], A[0][1], A[0][0]) / creall(B[0][0])),
      (xfma(B[0][1], A[0][0], A[0][1]) / creall(B[1][1]))
     },
     {
      (xfma(B[1][0], A[1][1], A[1][0]) / creall(B[0][0])),
      (xfma(B[0][1], A[1][0], A[1][1]) / creall(B[1][1]))
     }
    };
  xassgn2(A, C);
}

static inline void AC(xcomplex A[static 2][2], const xcomplex B[static 2][2])
{
  const xcomplex C[2][2] =
    {
     {
      (xfma(A[1][0], B[0][1], A[0][0]) / creall(B[0][0])),
      (xfma(A[1][1], B[0][1], A[0][1]) / creall(B[0][0]))
     },
     {
      (xfma(A[0][0], B[1][0], A[1][0]) / creall(B[1][1])),
      (xfma(A[0][1], B[1][0], A[1][1]) / creall(B[1][1]))
     }
    };
  xassgn2(A, C);
}

static inline void xhsvd2d(xcomplex A[static 2][2], xcomplex U[static 2][2], fint info[static 1])
{
  *info = FINT_C(0);
  A[1][0] = A[0][1] = x0;

  if (xisreal(A[0][0])) {
    if (copysignl(1.0L, creall(A[0][0])) == -1.0L) {
      U[0][0] = -U[0][0];
      U[1][0] = -U[1][0];
      A[0][0] = -A[0][0];
    }
  }
  else if (xisimag(A[0][0])) {
    if (cimagl(A[0][0]) < 0.0L) {
      U[0][0] = CMPLXL(-cimagl(U[0][0]), creall(U[0][0]));
      U[1][0] = CMPLXL(-cimagl(U[1][0]), creall(U[1][0]));
      A[0][0] = CMPLXL(-cimagl(A[0][0]), creall(A[0][0]));
    }
    else {
      U[0][0] = CMPLXL(cimagl(U[0][0]), -creall(U[0][0]));
      U[1][0] = CMPLXL(cimagl(U[1][0]), -creall(U[1][0]));
      A[0][0] = CMPLXL(cimagl(A[0][0]), -creall(A[0][0]));
    }
  }
  else {
    const extended w = cabsl(A[0][0]);
    const xcomplex v = conjl(A[0][0] / w);
    U[0][0] = xxmul(v, U[0][0]);
    U[1][0] = xxmul(v, U[1][0]);
    A[0][0] = CMPLXL(w, 0.0L);
  }

  if (xisreal(A[1][1])) {
    if (copysignl(1.0L, creall(A[1][1])) == -1.0L) {
      U[0][1] = -U[0][1];
      U[1][1] = -U[1][1];
      A[1][1] = -A[1][1];
    }
  }
  else if (xisimag(A[1][1])) {
    if (cimagl(A[1][1]) < 0.0L) {
      U[0][1] = CMPLXL(-cimagl(U[0][1]), creall(U[0][1]));
      U[1][1] = CMPLXL(-cimagl(U[1][1]), creall(U[1][1]));
      A[1][1] = CMPLXL(-cimagl(A[1][1]), creall(A[1][1]));
    }
    else {
      U[0][1] = CMPLXL(cimagl(U[0][1]), -creall(U[0][1]));
      U[1][1] = CMPLXL(cimagl(U[1][1]), -creall(U[1][1]));
      A[1][1] = CMPLXL(cimagl(A[1][1]), -creall(A[1][1]));
    }
  }
  else {
    const extended w = cabsl(A[1][1]);
    const xcomplex v = conjl(A[1][1] / w);
    U[0][1] = xxmul(v, U[0][1]);
    U[1][1] = xxmul(v, U[1][1]);
    A[1][1] = CMPLXL(w, 0.0L);
  }
}

static inline void xhsvd2u(const bool h, xcomplex A[static 2][2], xcomplex U[static 2][2], xcomplex Z[static 2][2], fint info[static 1])
{
  extended tu = 0.0L, cu = 1.0L, tz = 0.0L, cz = 1.0L;
  xcomplex y_ = (A[1][0] * creall(A[1][1])), z_ = x0;
  const xcomplex x_ = conjl(A[1][0] / creall(A[0][0]));
  const extended x = cabsl(x_), y = (creall(A[1][1]) / creall(A[0][0]));

  if (h) {
    if (x == 1.0L) {
      *info = ((y == 0.0L) ? FINT_C(-8) : FINT_C(-9));
      return;
    }
    extended t2 = ((x < y) ? (scalbnl(x, 1) * y) : (scalbnl(y, 1) * x));
    t2 /= fmal((y - x), (y + x), 1.0L);
    tu = (!(fabsl(t2) < scalbnl(1.0, (LDBL_MANT_DIG + 1))) ? copysignl(1.0L, t2) : (t2 / (1.0L + sqrtl(fmal(t2, t2, 1.0L)))));
    cu = sqrtl(fmal(tu, tu, 1.0L));
    y_ = ((y_ == x0) ? CMPLXL(tu, 0.0L) : ((conjl(y_) / cabsl(y_)) * tu));
    z_ = xqfma(y, y_, -x_);
    tz = cabsl(z_);
    cz = sqrtl(fmal(-tz, tz, 1.0L));
  }
  else {
    extended t2 = -((x < y) ? (scalbnl(x, 1) * y) : (scalbnl(y, 1) * x));
    t2 /= fmal((x - y), (x + y), 1.0L);
    tu = (!(fabsl(t2) < scalbnl(1.0, (LDBL_MANT_DIG + 1))) ? copysignl(1.0L, t2) : (t2 / (1.0L + sqrtl(fmal(t2, t2, 1.0L)))));
    cu = sqrtl(fmal(tu, tu, 1.0L));
    y_ = ((y_ == x0) ? CMPLXL(tu, 0.0L) : ((conjl(y_) / cabsl(y_)) * tu));
    z_ = xqfma(y, y_, -x_);
    tz = cabsl(z_);
    cz = sqrtl(fmal(tz, tz, 1.0L));
  }

  const xcomplex V[2][2] = { { CMPLXL(cu, 0.0L), y_ }, { -conjl(y_), CMPLXL(cu, 0.0L) } };
  CA2(V, U);
  CA2(V, A);

  const xcomplex W[2][2] = { { CMPLXL(cz, 0.0L), (h ? z_ : -z_) }, { conjl(z_), CMPLXL(cz, 0.0L) } };
  AC(A, W);
  AC(Z, W);

  xhsvd2d(A, U, info);
  *info = FINT_C(1);
}

static inline void xhsvd2l(xcomplex A[static 2][2], xcomplex U[static 2][2], xcomplex Z[static 2][2], fint info[static 1])
{
  const xcomplex x_ = (A[0][1] / creall(A[1][1]));
  const extended x = cabsl(x_), y = (creall(A[0][0]) / creall(A[1][1]));

  if (x == 1.0L) {
    *info = ((y == 0.0L) ? FINT_C(-10) : FINT_C(-11));
    return;
  }

  xcomplex y_ = (conjl(A[0][1]) * creall(A[0][0]));
  extended t2 = -((x < y) ? (scalbnl(x, 1) * y) : (scalbnl(y, 1) * x));
  t2 /= fmal((y - x), (y + x), 1.0L);
  const extended tu = (!(fabsl(t2) < scalbnl(1.0, (LDBL_MANT_DIG + 1))) ? copysignl(1.0L, t2) : (t2 / (1.0L + sqrtl(fmal(t2, t2, 1.0L)))));
  const extended cu = sqrtl(fmal(tu, tu, 1.0L));
  y_ = ((y_ == x0) ? CMPLXL(tu, 0.0L) : ((conjl(y_) / cabsl(y_)) * tu));
  const xcomplex z_ = -xqfma(y, y_, x_);
  const extended tz = cabsl(z_);
  const extended cz = sqrtl(fmal(-tz, tz, 1.0L));

  const xcomplex V[2][2] = { { CMPLXL(cu, 0.0L), y_ }, { -conjl(y_), CMPLXL(cu, 0.0L) } };
  CA2(V, U);
  CA2(V, A);

  const xcomplex W[2][2] = { { CMPLXL(cz, 0.0L), z_ }, { conjl(z_), CMPLXL(cz, 0.0L) } };
  AC(A, W);
  AC(Z, W);

  xhsvd2d(A, U, info);
  *info = FINT_C(1);
}

static inline void xhsvd2g(const bool h, xcomplex A[static 2][2], xcomplex U[static 2][2], xcomplex Z[static 2][2], fint info[static 1])
{
  xcomplex
    r = CMPLXL(cabsl(A[0][0]), cabsl(A[0][1])),
    s = CMPLXL(cabsl(A[1][0]), cabsl(A[1][1]));
  extended
    c = cabsl(r),
    t = cabsl(s);

  if (c < t) {
    xswp(&(A[0][0]), &(A[1][0]));
    xswp(&(A[0][1]), &(A[1][1]));
    if (!h) {
      xswp(&(Z[0][0]), &(Z[1][0]));
      xswp(&(Z[0][1]), &(Z[1][1]));
    }
    *info = FINT_C(1);
  }
  else
    *info = FINT_C(0);

  if (cabsl(A[0][0]) < cabsl(A[0][1])) {
    xswp(&(U[0][0]), &(U[0][1]));
    xswp(&(U[1][0]), &(U[1][1]));
    xswp(&(A[0][0]), &(A[0][1]));
    xswp(&(A[1][0]), &(A[1][1]));
  }

  s = conjl((xisreal(A[0][0])) ? (A[0][1] / creall(A[0][0])) : xxdiv(A[0][1], A[0][0]));
  t = cabsl(s);
  c = sqrtl(fmal(t, t, 1.0L));

  const xcomplex Q[2][2] = { { c, -conjl(s) }, { s, c } };
  CA2(Q, U);
  CA2(Q, A);
  A[0][1] = x0;

  if (xisreal(A[0][0])) {
    if (copysignl(1.0L, creall(A[0][0])) == -1.0L) {
      U[0][0] = -U[0][0];
      U[1][0] = -U[1][0];
      A[0][0] = -A[0][0];
      A[1][0] = -A[1][0];
    }
  }
  else if (xisimag(A[0][0])) {
    if (cimagl(A[0][0]) < 0.0L) {
      U[0][0] = CMPLXL(-cimagl(U[0][0]), creall(U[0][0]));
      U[1][0] = CMPLXL(-cimagl(U[1][0]), creall(U[1][0]));
      A[0][0] = CMPLXL(-cimagl(A[0][0]), creall(A[0][0]));
      A[1][0] = CMPLXL(-cimagl(A[1][0]), creall(A[1][0]));
    }
    else {
      U[0][0] = CMPLXL(cimagl(U[0][0]), -creall(U[0][0]));
      U[1][0] = CMPLXL(cimagl(U[1][0]), -creall(U[1][0]));
      A[0][0] = CMPLXL(cimagl(A[0][0]), -creall(A[0][0]));
      A[1][0] = CMPLXL(cimagl(A[1][0]), -creall(A[1][0]));
    }
  }
  else {
    r = A[0][0];
    A[0][0] = CMPLXL(cabsl(r), 0.0L);
    r = conjl(r / creall(A[0][0]));
    U[0][0] = xxmul(r, U[0][0]);
    U[1][0] = xxmul(r, U[1][0]);
    A[1][0] = xxmul(r, A[1][0]);
  }

  if (xisreal(A[1][1])) {
    if (copysignl(1.0L, creall(A[1][1])) == -1.0L) {
      U[0][1] = -U[0][1];
      U[1][1] = -U[1][1];
      A[1][1] = -A[1][1];
    }
  }
  else if (xisimag(A[1][1])) {
    if (cimagl(A[1][1]) < 0.0L) {
      U[0][1] = CMPLXL(-cimagl(U[0][1]), creall(U[0][1]));
      U[1][1] = CMPLXL(-cimagl(U[1][1]), creall(U[1][1]));
      A[1][1] = CMPLXL(-cimagl(A[1][1]), creall(A[1][1]));
    }
    else {
      U[0][1] = CMPLXL(cimagl(U[0][1]), -creall(U[0][1]));
      U[1][1] = CMPLXL(cimagl(U[1][1]), -creall(U[1][1]));
      A[1][1] = CMPLXL(cimagl(A[1][1]), -creall(A[1][1]));
    }
  }
  else {
    r = A[1][1];
    A[1][1] = CMPLXL(cabsl(r), 0.0L);
    r = conjl(r / creall(A[1][1]));
    U[0][1] = xxmul(r, U[0][1]);
    U[1][1] = xxmul(r, U[1][1]);
  }

  if (A[1][0] == x0) {
    if (h) {
      if (*info == FINT_C(1)) {
        xswp(&(U[0][0]), &(U[0][1]));
        xswp(&(U[1][0]), &(U[1][1]));
        xswp(&(A[0][0]), &(A[1][1]));
      }
      else
        *info = FINT_C(1);
    }
    else
      *info = FINT_C(1);
  }
  else if (h && (*info == FINT_C(1))) {
    xswp(&(A[0][0]), &(A[1][0]));
    xswp(&(A[0][1]), &(A[1][1]));
    xswp(&(U[0][0]), &(U[0][1]));
    xswp(&(U[1][0]), &(U[1][1]));
    xswp(&(A[0][0]), &(A[0][1]));
    xswp(&(A[1][0]), &(A[1][1]));
    xhsvd2l(A, U, Z, info);
  }
  else
    xhsvd2u(h, A, U, Z, info);
}

static inline void xhsvd2s(const fint h, xcomplex A[static 2][2], xcomplex U[static 2][2], xcomplex Z[static 2][2], fint info[static 1])
{
  if (((h == FINT_C(2)) && (creall(A[0][0]) < creall(A[1][1]))) || ((h == FINT_C(-2)) && (creall(A[0][0]) > creall(A[1][1])))) {
    xswp(&(U[0][0]), &(U[0][1]));
    xswp(&(U[1][0]), &(U[1][1]));
    xswp(&(A[0][0]), &(A[1][1]));
    xswp(&(Z[0][0]), &(Z[1][0]));
    xswp(&(Z[0][1]), &(Z[1][1]));
  }

  if ((U[0][0] != x1) || (U[0][1] != x0) || (U[1][0] != x0) || (U[1][1] != x1))
    *info += FINT_C(2);
  if ((Z[0][0] != x1) || (Z[0][1] != x0) || (Z[1][0] != x0) || (Z[1][1] != x1))
    *info += FINT_C(4);
}

static inline void xhsvd2_(xcomplex A[static 2][2], const fint J[static 2], xcomplex U[static 2][2], xcomplex Z[static 2][2], fint info[static 1])
{
  if ((A[0][1] == x0) && (A[1][0] == x0))
    xhsvd2d(A, U, info);
  else /* A general */
    xhsvd2g((J[0] != J[1]), A, U, Z, info);

  if (*info >= FINT_C(0))
    xhsvd2s((J[0] + J[1]), A, U, Z, info);
}

void zhsvd2_(dcomplex A[static 2][2], const fint J[static 2], dcomplex U[static 2][2], dcomplex Z[static 2][2], fint info[static 1])
{
  if (!zisfinite(A[0][0]))
    *info = FINT_C(-1);
  else if (!zisfinite(A[0][1]))
    *info = FINT_C(-2);
  else if (!zisfinite(A[1][0]))
    *info = FINT_C(-3);
  else if (!zisfinite(A[1][1]))
    *info = FINT_C(-4);
  else /* A OK */
    *info = FINT_C(0);

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

  xcomplex A_[2][2] = { { A[0][0], A[0][1] }, { A[1][0], A[1][1] } };
  xcomplex U_[2][2] = { { x1, x0 }, { x0, x1 } };
  xcomplex Z_[2][2] = { { x1, x0 }, { x0, x1 } };

  xhsvd2_(A_, J, U_, Z_, info);

  if (*info >= FINT_C(0)) {
    zassgn2(A, (const xcomplex (*)[2])A_);
    zassgn2(U, (const xcomplex (*)[2])U_);
    zassgn2(Z, (const xcomplex (*)[2])Z_);
  }
}
#undef x1
#undef x0
#endif /* USE_EXTENDED */
