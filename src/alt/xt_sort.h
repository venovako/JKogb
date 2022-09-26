#ifndef XT_SORT_H
#define XT_SORT_H

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
XT_SORT
#else /* !_WIN32 */
xt_sort_
#endif /* ?_WIN32 */
(const F_INT *const n, long double *const x) throw();
#else /* !__cplusplus */
extern void
#ifdef _WIN32
XT_SORT
#else /* !_WIN32 */
xt_sort_
#endif /* ?_WIN32 */
(const F_INT n[static restrict 1], long double x[static restrict 1]);
#endif /* ?__cplusplus */

#endif /* !XT_SORT_H */
