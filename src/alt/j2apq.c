#include "j2apq.h"

int apq_init(apq o[static 1], const uint16_t n, const fint *const j)
{
  o->a = (wpqb*)NULL;
  o->q = o->p = (uint16_t*)NULL;

  int r = 0;
  if (!(o->n_a = ((n * ((uint32_t)n - 1u)) >> 1u)))
    r = -2;
  else if (!(o->a = (wpqb*)malloc(o->n_a * sizeof(wpqb))))
    r = 1;
  else if (!(o->p = (uint16_t*)malloc(o->n_a * sizeof(uint16_t))))
    r = 2;
  else if (!(o->q = (uint16_t*)malloc(o->n_a * sizeof(uint16_t))))
    r = 3;
  else
    o->n_t = 0u;

  if (r) {
    apq_free(o);
    return r;
  }

  const uint16_t n_1 = n - UINT16_C(1);
  if (j) {
    uint32_t ih = o->n_a;
    // row-cyclic pass
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
    // row-cyclic pass
    for (uint16_t r = UINT16_C(0); r < n_1; ++r) {
      for (uint16_t c = (r + UINT16_C(1)); c < n; ++c) {
        (o->p)[o->n_t] = r;
        (o->q)[o->n_t] = c;
        ++(o->n_t);
      }
    }
  }

  return r;
}

void apq_free(apq o[static 1])
{
  free(o->q);
  free(o->p);
  free(o->a);
  o->n_t = o->n_a = 0u;
  o->a = (wpqb*)NULL;
  o->q = o->p = (uint16_t*)NULL;
}
