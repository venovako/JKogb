#ifndef XOFFPQ_H
#define XOFFPQ_H

#ifndef F_INT
#ifdef MKL_ILP64
#ifdef _WIN32
#define F_INT long long
#else /* !_WIN32 */
#define F_INT long
#endif /* ?_WIN32 */
#else /* !MKL_ILP64 */
#define F_INT int
#endif /* ?MKL_ILP64 */
#endif /* !F_INT */

#ifndef F_UINT
#ifdef MKL_ILP64
#ifdef _WIN32
#define F_UINT unsigned long long
#else /* !_WIN32 */
#define F_UINT unsigned long
#endif /* ?_WIN32 */
#else /* !MKL_ILP64 */
#define F_UINT unsigned
#endif /* ?MKL_ILP64 */
#endif /* !F_UINT */

#ifdef __cplusplus
extern "C" void
#ifdef _WIN32
XOFFPQ
#else /* !_WIN32 */
xoffpq_
#endif /* ?_WIN32 */
(const F_INT *const n, const double *const A, const F_INT *const ldA, const F_UINT *const p, const F_UINT *const q, const F_UINT *const b, long double *const x, F_INT *const info) throw();
#else /* !__cplusplus */
extern void
#ifdef _WIN32
XOFFPQ
#else /* !_WIN32 */
xoffpq_
#endif /* ?_WIN32 */
(const F_INT n[static 1], const double A[static restrict 1], const F_INT ldA[static 1], const F_UINT p[static 1], const F_UINT q[static 1], const F_UINT b[static 1], long double x[static restrict 1], F_INT info[static restrict 1]);
#endif /* ?__cplusplus */

#endif /* !XOFFPQ_H */
