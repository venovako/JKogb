#ifndef VDP_H
#define VDP_H

#ifdef __SSE4_1__
#include <immintrin.h>
#else /* !__SEE4_1__ */
#error SSE4.1 instructions have to be available
#endif /* ?__SSE4_1__ */

static inline void dtw(double d[static 1], const double Apq, const double Aqp)
{
  register const __m128d a = _mm_set_pd(Aqp, Apq);
  _mm_storel_pd(d, _mm_dp_pd(a, a, 0x31));
}

static inline void ztw(double d[static 1], const double Apqr, const double Apqi, const double Aqpr, const double Aqpi)
{
  register const __m128d pq = _mm_set_pd(Apqi, Apqr);
  register const __m128d qp = _mm_set_pd(Aqpi, Aqpr);
  _mm_storel_pd(d, _mm_add_pd(_mm_dp_pd(pq, pq, 0x31), _mm_dp_pd(qp, qp, 0x31)));
}
#endif /* !VDP_H */
