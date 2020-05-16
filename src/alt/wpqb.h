#ifndef WPQB_H
#define WPQB_H
#include "common.h"

typedef union {
  long double w;
  struct {
    uint8_t a[10];
    // p, q, b are stored in the highest, unused 2+2+2=6 bytes of w
    uint16_t p, q, b;
  } i;
} wpqb;

typedef union {
  long double w;
  struct {
    uint8_t a[10];
    // s - step length
    uint16_t s;
    // f - starting index in a
    uint32_t f;
  } i;
} wpqb_info;

static inline void wpqb_init(wpqb a[static 1], const long double w, const uint16_t p, const uint16_t q)
{
  assert(p < q);
  a->w = w;
  a->i.p = p;
  a->i.q = q;
  a->i.b = (q - p);
}

static inline int wpqb_invalid(const wpqb a[static 1])
{
  if (!(a->w == a->w))
    return 1;
  if (a->w < -LDBL_MAX)
    return -1;
  return 0;
}

extern int wpqb_cmp(const wpqb a[static 1], const wpqb b[static 1]);
extern uint32_t wpqb_clean(wpqb *const a, const uint32_t n_a);

extern void wpqb_ncp0(const uint16_t n, const uint32_t n_a, wpqb *const restrict a, const uint32_t f, const uint16_t n_s, uint32_t *const restrict s, wpqb_info *const restrict w);
extern void wpqb_ncp1(const uint16_t n, const uint32_t n_a, wpqb *const restrict a, const uint32_t f, const uint16_t n_s, uint32_t *const restrict s, wpqb_info *const restrict w);
#endif /* !WPQB_H */
