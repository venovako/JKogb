#include "xoffpq.h"
#include "xoffsq.h"
#ifdef __cplusplus
#include <cassert>
#include <cstddef>
void
#ifdef _WIN32
XOFFPQ
#else /* !_WIN32 */
xoffpq_
#endif /* ?_WIN32 */
(const MKL_INT *const n, const double *const A, const MKL_INT *const ldA, const MKL_UINT *const p, const MKL_UINT *const q, const MKL_UINT *const b, long double *const x, MKL_INT *const info) throw()
#else /* !__cplusplus */
#include <assert.h>
#include <stddef.h>
void
#ifdef _WIN32
XOFFPQ
#else /* !_WIN32 */
xoffpq_
#endif /* ?_WIN32 */
(const MKL_INT n[static 1], const double A[static restrict 1], const MKL_INT ldA[static 1], const MKL_UINT p[static 1], const MKL_UINT q[static 1], const MKL_UINT b[static 1], long double x[static restrict 1], MKL_INT info[static restrict 1])
#endif /* ?__cplusplus */
{
  assert(n);
  assert(A);
  assert(ldA);
  assert(p);
  assert(q);
  assert(b);
  assert(x);
  assert(info);
  const size_t _ldA = (size_t)((*ldA < 0) ? -*ldA : *ldA);
  const MKL_UINT _n = (MKL_UINT)((*n < 0) ? -*n : *n);
  MKL_UINT _b = (*b << 1u);
  if (!_b || (_b > _n) || (_n % *b)) {
    *info = -6;
    return;
  }
  const MKL_UINT n_b = _n / *b;
  *info = 0;
  if (*q >= n_b)
    *info = -5;
  if (*p >= *q)
    *info = -4;
  if (_n > _ldA)
    *info = -1;
  if (*info)
    return;
  const MKL_UINT bp = *b * *p;
  const double *const App = A + bp * _ldA + bp;
  if (*q == (*p + 1u))
    xoffsq_((const MKL_INT*)&_b, App, ldA, x, info);
  else {
    const MKL_INT _b = -(MKL_INT)*b;
    const MKL_UINT bq = *b * *q;
    const double *const Aqp = A + bp * _ldA + bq;
    const double *const Apq = A + bq * _ldA + bp;
    const double *const Aqq = A + bq * _ldA + bq;
    long double spp, sqp, spq, sqq;
    xoffsq_((const MKL_INT*)b, App, ldA, &spp, info);
    if (*info)
      return;
    xoffsq_(&_b, Aqp, ldA, &sqp, info);
    if (*info)
      return;
    xoffsq_(&_b, Apq, ldA, &spq, info);
    if (*info)
      return;
    xoffsq_((const MKL_INT*)b, Aqq, ldA, &sqq, info);
    if (*info)
      return;
    if (sqp < spp) { *x = sqp; sqp = spp; spp = *x; }
    if (spq < sqp) { *x = spq; spq = sqp; sqp = *x; }
    if (sqq < spq) { *x = sqq; sqq = spq; spq = *x; }

    if (sqp < spp) { *x = sqp; sqp = spp; spp = *x; }
    if (spq < sqp) { *x = spq; spq = sqp; sqp = *x; }

    if (sqp < spp) { *x = sqp; sqp = spp; spp = *x; }
    *x = ((spp + sqp) + spq) + sqq;
  }
  char *const cx = (char*)x;
  char *const cp = (char*)p;
  cx[10u] = cp[0u];
  cx[11u] = cp[1u];
  cx[12u] = cp[2u];
  char *const cq = (char*)q;
  cx[13u] = cq[0u];
  cx[14u] = cq[1u];
  cx[15u] = cq[2u];
}
