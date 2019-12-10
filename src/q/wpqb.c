#include "wpqb.h"

wpqb *wpqb_alloc(const uint32_t n_a)
{
  return (n_a ? aligned_alloc((size_t)16, n_a * sizeof(wpqb)) : NULL);
}

int wpqb_free(wpqb *const a)
{
  if (a) {
    free(a);
    return 0;
  }
  return -1;
}

int wpqb_init(wpqb *const a, const uint16_t p, const uint16_t q, const extended *const w)
{
  if (a) {
    if (p >= q)
      return -2;
    a->i.p = p;
    a->i.q = q;
    a->i.b = (q - p);
    if (w)
      a->w = *w;
    else {
      union { uint16_t us; struct { uint8_t lo, hi; } lh; } x;
      x.us = a->i.p;
      a->i.a[0] = x.lh.lo;
      a->i.a[1] = x.lh.hi;
      x.us = a->i.q;
      a->i.a[2] = x.lh.lo;
      a->i.a[3] = x.lh.hi;
      x.us = a->i.b;
      a->i.a[4] = x.lh.lo;
      a->i.a[5] = x.lh.hi;
      a->i.a[6] = 0x00u;
      a->i.a[7] = 0xC0u;
      a->i.a[8] = 0xFFu;
      a->i.a[9] = 0x7Fu;
    }
    return 0;
  }
  return -1;
}

int wpqb_print(FILE *const f, wpqb *const a)
{
  return (f ? (a ? fprintf(f, "%# -30.21LE,%5hu,%5hu,%5hu", a->w, a->i.p, a->i.q, a->i.b) : -2) : -1);
}

int wpqb_dump(FILE *const f, wpqb *const a, const uint32_t n_a)
{
  if (!n_a)
    return 0;
  if (!a)
    return -2;
  if (!f)
    return -1;

  if (16 != fprintf(f, "\"w\",\"p\",\"q\",\"b\"\n"))
    return 1;
  for (uint32_t i = UINT32_C(0); i < n_a; ++i) {
    if (48 != wpqb_print(f, (a + i)))
      return (int)((i + UINT32_C(1)) << 1);
    if (1 != fprintf(f, "\n"))
      return (int)(((i + UINT32_C(1)) << 1) + UINT32_C(1));
  }

  return (fflush(f) ? (int)(n_a << 1) : 0);
}

int wpqb_cmp(const wpqb *const a, const wpqb *const b)
{
  assert(a);
  assert(b);
  if (a == b)
    return 0;
  if (a->w < b->w)
    return 1;
  if (a->w > b->w)
    return -1;
  if (a->w == b->w) {
    if (a->i.b < b->i.b)
      return 2;
    if (a->i.b > b->i.b)
      return -2;
    if (a->i.p < b->i.p)
      return 3;
    if (a->i.p > b->i.p)
      return -3;
    if (a->i.q < b->i.q)
      return 4;
    if (a->i.q > b->i.q)
      return -4;
    return 0;
  }
  if (b->w == b->w)
    return 5;
  return -5;
}

int wpqb_sort(wpqb *const a, const uint32_t n_a)
{
  if (n_a) {
    if (a)
      qsort(a, (size_t)n_a, sizeof(*a), (int (*)(const void*, const void*))wpqb_cmp);
    else /* NULL a invalid here */
      return -1;
  }
  return 0;
}

int wpqb_ncp(wpqb *const a, const uint32_t n_a, uint16_t *const s, const uint16_t n_s, extended *const w)
{
  if (!n_s)
    return 0;
  if (!s)
    return -3;
  if (!n_a)
    return 0;
  if (!a)
    return -1;

  uint16_t n_s_ = UINT16_C(0);
  uint32_t j = UINT32_C(0);
  for (uint16_t i = UINT16_C(0); i < n_s; ) {
    s[i] = j;
    n_s_ = ++i;
    for (++j; j < n_a; ) {
      bool c = !(a[j].w == a[j].w);
      if (c) {
        ++j;
        continue;
      }
      const extended ap = a[j].i.p;
      const extended aq = a[j].i.q;
      for (uint16_t k = i; k; ) {
        --k;
        const extended bp = a[s[k]].i.p;
        const extended bq = a[s[k]].i.q;
        if ((ap == bp) || (ap == bq) || (aq == bp) || (aq == bq)) {
          c = true;
          break;
        }
      }
      if (c)
        ++j;
      else
        break;
    }
    if (j >= n_a)
      break;
  }

  if (w) {
    uint16_t p = UINT16_C(0), n = n_s_;
    for (uint16_t i = UINT16_C(0); i < n_s_; ++i) {
      if (a[s[i]].w > 0.0L)
        p = i + UINT16_C(1);
      else if (a[s[i]].w < 0.0L) {
        n = i;
        break;
      }
    }
    w[0] =  0.0L;
    w[1] = -0.0L;
    for (uint16_t i = p; i; )
      w[0] += a[s[--i]].w;
    for (uint16_t i = n; i < n_s_; ++i)
      w[1] += a[s[i]].w;
  }

  return (int)n_s_;
}
