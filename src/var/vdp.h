#ifndef VDP_H
#define VDP_H
#include "common.h"

#ifdef USE_DOUBLE
#ifdef __SSE4_1__
#include <immintrin.h>
#endif /* __SSE4_1__ */
#endif /* USE_DOUBLE */

#ifdef USE_DOUBLE
static inline void wd(double d[static 1], const double apq, const double aqp)
{
#ifdef __SSE4_1__
  register const __m128d a = _mm_set_pd(aqp, apq);
  _mm_storel_pd(d, _mm_dp_pd(a, a, 0x31));
#else /* !__SSE4_1__ */
  *d = apq * apq + aqp * aqp;
#endif /* ?__SSE4_1__ */
}

static inline void wz(double d[static 1], const double apqr, const double apqi, const double aqpr, const double aqpi)
{
#ifdef __SSE4_1__
  register const __m128d pq = _mm_set_pd(apqi, apqr);
  register const __m128d qp = _mm_set_pd(aqpi, aqpr);
  _mm_storel_pd(d, _mm_add_pd(_mm_dp_pd(pq, pq, 0x31), _mm_dp_pd(qp, qp, 0x31)));
#else /* !__SSE4_1__ */
  *d = (apqr * apqr + apqi * apqi) + (aqpr * aqpr + aqpi * aqpi);
#endif /* ?__SSE4_1__ */
}
#else /* !USE_DOUBLE */
static inline void wd(long double d[static 1], const long double apq, const long double aqp)
{
  *d = apq * apq + aqp * aqp;
}

static inline void wz(long double d[static 1], const long double apqr, const long double apqi, const long double aqpr, const long double aqpi)
{
  *d = (apqr * apqr + apqi * apqi) + (aqpr * aqpr + aqpi * aqpi);
}
#endif /* ?USE_DOUBLE */
#endif /* !VDP_H */
