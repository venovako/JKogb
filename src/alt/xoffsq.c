#ifdef NEG_LOOP
#error NEG_LOOP already defined
#else /* !NEG_LOOP */
#define NEG_LOOP                     \
for (size_t j = 0u; j < _n; ++j) {   \
  const size_t c = j * _ldA;         \
  for (size_t i = 0u; i < _n; ++i) { \
    const long double y = A[c + i];  \
    s += y * y;                      \
  }                                  \
}
#endif /* ?NEG_LOOP */

#ifdef POS_LOOP
#error POS_LOOP already defined
#else /* !POS_LOOP */
#define POS_LOOP                         \
for (size_t j = 0u; j < _n; ++j) {       \
  const size_t c = j * _ldA;             \
  for (size_t i = 0u; i < j; ++i) {      \
    const long double y = A[c + i];      \
    s += y * y;                          \
  }                                      \
  for (size_t i = j + 1u; i < _n; ++i) { \
    const long double y = A[c + i];      \
    s += y * y;                          \
  }                                      \
}
#endif /* ?POS_LOOP */

#include "xoffsq.h"
#ifdef __cplusplus
#include <cassert>
#include <cstddef>
void
#ifdef _WIN32
XOFFSQ
#else /* !_WIN32 */
#ifdef __powerpc__
xoffsq
#else /* !__powerpc__ */
xoffsq_
#endif /* ?__powerpc__ */
#endif /* ?_WIN32 */
(const F_INT *const n, const double *const A, const F_INT *const ldA, long double *const x, F_INT *const info) throw()
#else /* !__cplusplus */
#include <assert.h>
#include <stddef.h>
void
#ifdef _WIN32
XOFFSQ
#else /* !_WIN32 */
#ifdef __powerpc__
xoffsq
#else /* !__powerpc__ */
xoffsq_
#endif /* ?__powerpc__ */
#endif /* ?_WIN32 */
(const F_INT n[static 1], const double A[static restrict 1], const F_INT ldA[static 1], long double x[static restrict 1], F_INT info[static restrict 1])
#endif /* ?__cplusplus */
{
  assert(n);
  assert(A);
  assert(ldA);
  assert(x);
  assert(info);
  const size_t _n = (size_t)((*n < 0) ? -*n : *n);
  const size_t _ldA = (size_t)((*ldA < 0) ? -*ldA : *ldA);
  if ((*info = ((_n > _ldA) ? -1 : 0)))
    return;
  long double s = 0.0L;
  if (*n < 0) {
    if (*ldA <= 0) {
      NEG_LOOP
    }
    else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(_n,A,_ldA) reduction(+:s)
#endif /* _OPENMP */
      NEG_LOOP
    }
  }
  else if (*n > 0) {
    if (*ldA <= 0) {
      POS_LOOP
    }
    else {
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(_n,A,_ldA) reduction(+:s)
#endif /* _OPENMP */
      POS_LOOP
    }
  }
  else {
    *x = s;
    return;
  }
  if (s != s)
    *info = -2;
  else
    *x = s;
}
