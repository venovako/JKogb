#ifndef VDP_H
#define VDP_H
#include "common.h"

static inline void wd(double d[static 1], const double apq, const double aqp)
{
  register const __m128d a = _mm_set_pd(aqp, apq);
  _mm_storel_pd(d, _mm_dp_pd(a, a, 0x31));
}

static inline void wz(double d[static 1], const double apqr, const double apqi, const double aqpr, const double aqpi)
{
  register const __m128d pq = _mm_set_pd(apqi, apqr);
  register const __m128d qp = _mm_set_pd(aqpi, aqpr);
  _mm_storel_pd(d, _mm_add_pd(_mm_dp_pd(pq, pq, 0x31), _mm_dp_pd(qp, qp, 0x31)));
}

#endif /* !VDP_H */
