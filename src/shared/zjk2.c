#include "zjk2.h"

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
  if (cimagl(a) == 0.0L)
    return xqfma(creall(a), b, c);
  if (creall(a) == 0.0L)
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

#ifdef USE_EXTENDED
void zhsvd2_(dcomplex A[static 2][2], const fint J[static 2], dcomplex U[static 2][2], dcomplex Z[static 2][2], fint info[static 1])
{
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
  xcomplex U_[2][2] = { { CMPLXL(1.0L, 0.0L), CMPLXL(0.0L, 0.0L) }, { CMPLXL(0.0L, 0.0L), CMPLXL(1.0L, 0.0L) } };
  xcomplex Z_[2][2] = { { CMPLXL(1.0L, 0.0L), CMPLXL(0.0L, 0.0L) }, { CMPLXL(0.0L, 0.0L), CMPLXL(1.0L, 0.0L) } };
}
#endif /* USE_EXTENDED */
