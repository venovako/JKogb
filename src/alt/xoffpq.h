#ifndef XOFFPQ_H
#define XOFFPQ_H

#include <mkl.h>

#ifdef __cplusplus
extern "C" void
#ifdef _WIN32
XOFFPQ
#else /* !_WIN32 */
xoffpq_
#endif /* ?_WIN32 */
(const MKL_INT *const n, const double *const A, const MKL_INT *const ldA, const MKL_UINT *const p, const MKL_UINT *const q, const MKL_UINT *const b, long double *const x, MKL_INT *const info) throw();
#else /* !__cplusplus */
extern void
#ifdef _WIN32
XOFFPQ
#else /* !_WIN32 */
xoffpq_
#endif /* ?_WIN32 */
(const MKL_INT n[static 1], const double A[static restrict 1], const MKL_INT ldA[static 1], const MKL_UINT p[static 1], const MKL_UINT q[static 1], const MKL_UINT b[static 1], long double x[static restrict 1], MKL_INT info[static restrict 1]);
#endif /* ?__cplusplus */

#endif /* !XOFFPQ_H */
