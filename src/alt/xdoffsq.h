#ifndef XDOFFSQ_H
#define XDOFFSQ_H

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

#ifdef __cplusplus
extern "C" void
#ifdef _WIN32
XDOFFSQ
#else /* !_WIN32 */
#ifdef __powerpc__
xdoffsq
#else /* !__powerpc__ */
xdoffsq_
#endif /* ?__powerpc__ */
#endif /* ?_WIN32 */
(const F_INT *const n, const double *const A, const F_INT *const ldA, long double *const x, F_INT *const info) throw();
#else /* !__cplusplus */
extern void
#ifdef _WIN32
XDOFFSQ
#else /* !_WIN32 */
#ifdef __powerpc__
xdoffsq
#else /* !__powerpc__ */
xdoffsq_
#endif /* ?__powerpc__ */
#endif /* ?_WIN32 */
(const F_INT n[static 1], const double A[static restrict 1], const F_INT ldA[static 1], long double x[static restrict 1], F_INT info[static restrict 1]);
#endif /* ?__cplusplus */

#endif /* !XDOFFSQ_H */
