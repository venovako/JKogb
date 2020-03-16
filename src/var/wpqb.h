#ifndef WPQB_H
#define WPQB_H
#ifdef USE_EXTENDED
#include "common.h"

typedef union {
  alignas(EXTENDED_ALIGN_B) extended w;
  alignas(EXTENDED_ALIGN_B) struct {
    uint8_t a[10];
    /* p, q, b are stored in the highest, unused 2+2+2=6 bytes of w */
    uint16_t p, q, b;
  } i;
} wpqb;

static inline void wpqb_init(wpqb a[static 1], const extended w, const uint16_t p, const uint16_t q)
{
  assert(p < q);
  a->w = w;
  a->i.p = p;
  a->i.q = q;
  a->i.b = (q - p);
}

extern wpqb *wpqb_alloc(const uint32_t n_a);
extern int wpqb_print(FILE f[static 1], wpqb a[static 1]);
extern int wpqb_dump(FILE f[static 1], wpqb *const a, const uint32_t n_a);
extern int wpqb_cmp(const wpqb a[static 1], const wpqb b[static 1]);
extern int wpqb_sort(wpqb *const a, const uint32_t n_a);
extern int wpqb_ncp0(const uint16_t n, wpqb *const a, const uint32_t n_a, uint16_t *const s, const uint16_t n_s, extended *const w);
extern int wpqb_ncp1(const uint16_t n, wpqb *const a, const uint32_t n_a, uint16_t *const s, const uint16_t n_s, extended *const w);
#endif /* USE_EXTENDED */
#endif /* !WPQB_H */
