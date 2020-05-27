#include "j2pq.h"

void pq_init(pq o[static 1], const uint16_t n, const int *const j)
{
  o->n_t = 0u;
  o->p = (uint16_t*)NULL;
  o->q = (uint16_t*)NULL;
  if (!(o->n_a = ((n * ((uint32_t)n - 1u)) >> 1u)))
    return;

  uint16_t *const x = (uint16_t*)malloc((o->n_a << 1u) * sizeof(uint16_t));
  if (!x)
    return;
  o->p = x;
  o->q = x + o->n_a;

  // row-cyclic pass
  const uint16_t n_1 = n - UINT16_C(1);
  if (j) {
    uint32_t ih = o->n_a;
    for (uint16_t r = UINT16_C(0); r < n_1; ++r) {
      for (uint16_t c = (r + UINT16_C(1)); c < n; ++c) {
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
    for (uint16_t r = UINT16_C(0); r < n_1; ++r) {
      for (uint16_t c = (r + UINT16_C(1)); c < n; ++c) {
        (o->p)[o->n_t] = r;
        (o->q)[o->n_t] = c;
        ++(o->n_t);
      }
    }
  }
}

void pq_free(pq o[static 1])
{
  free(o->p);
  o->n_t = o->n_a = 0u;
  o->q = o->p = (uint16_t*)NULL;
}
