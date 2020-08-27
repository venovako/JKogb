#include "j2apq.h"

apq *apq_create(const uint16_t n, const fint *const j)
{
  apq *const o = (apq*)(n ? calloc(1u, sizeof(apq)) : NULL);
  if (o) {
    int e = 0;
    if (!(o->n_a = ((n * ((uint32_t)n - 1u)) >> 1u)))
      e = -1;
    else if (!(o->a = (wpqb*)malloc(o->n_a * sizeof(wpqb))))
      e = 1;
    else if (!(o->p = (uint16_t*)malloc(o->n_a * sizeof(uint16_t))))
      e = 2;
    else if (!(o->q = (uint16_t*)malloc(o->n_a * sizeof(uint16_t))))
      e = 3;

    if (e)
      return apq_free(o);

    const uint16_t n_1 = (uint16_t)(n - UINT16_C(1));
    if (j) {
      uint32_t ih = o->n_a;
      // row-cyclic pass
      for (uint16_t r = UINT16_C(0); r < n_1; ++r) {
        for (uint16_t c = (uint16_t)(r + UINT16_C(1)); c < n; ++c) {
          if (j[r] == j[c]) {
            (o->p)[o->n_t] = r;
            (o->q)[o->n_t] = c;
            ++(o->n_t);
          }
          else {
            --ih;
            (o->p)[ih] = r;
            (o->q)[ih] = c;
          }
        }
      }
    }
    else {
      // row-cyclic pass
      for (uint16_t r = UINT16_C(0); r < n_1; ++r) {
        for (uint16_t c = (uint16_t)(r + UINT16_C(1)); c < n; ++c) {
          (o->p)[o->n_t] = r;
          (o->q)[o->n_t] = c;
          ++(o->n_t);
        }
      }
    }
  }

  return o;
}

apq *apq_free(apq *const o)
{
  if (o) {
    if (o->q) {
      free(o->q);
      o->q = (uint16_t*)NULL;
    }
    if (o->p) {
      free(o->p);
      o->p = (uint16_t*)NULL;
    }
    if (o->a) {
      free(o->a);
      o->a = (wpqb*)NULL;
    }
    o->n_a = o->n_t = 0u;
    free(o);
  }
  return (apq*)NULL;
}
