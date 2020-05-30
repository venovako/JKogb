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

extern int wpqb_cmp(const wpqb a[static 1], const wpqb b[static 1]);
extern uint32_t wpqb_sort0(const uint32_t n_a, wpqb a[static 1]);
extern uint32_t wpqb_sort1(const uint32_t n_a, wpqb a[static 1]);
extern uint32_t wpqb_sort(const uint32_t n_a, wpqb a[static 1]);
extern void wpqb_ncp0(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1]);
extern void wpqb_ncp1(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1]);
extern void wpqb_ncp(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1]);
extern void wpqb_run0(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1]);
extern void wpqb_run1(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1]);
extern void wpqb_run(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1]);
#endif /* !WPQB_H */