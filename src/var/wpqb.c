#include "wpqb.h"
#ifdef USE_EXTENDED
wpqb *wpqb_alloc(const uint32_t n_a)
{
  wpqb *p = (wpqb*)NULL;
  if (n_a) {
    const size_t s = n_a * sizeof(wpqb);
    p = (wpqb*)aligned_alloc((size_t)EXTENDED_ALIGN_B, s);
#ifndef NDEBUG
    if (p)
      (void)memset(p, 0, s);
#endif /* !NDEBUG */
  }
  return p;
}

int wpqb_print(FILE f[static 1], wpqb a[static 1])
{
  char s[31] = { '\0' };
  int l = sprintf(s, "%# -30.21LE", a->w);
  assert(l > 0);
  char *d = s + 30;
  for (--d; isblank(*d); --d)
    *d = '\0';
  char *e = strrchr(s, 'E');
  if (e) {
    e += 2;
    l = (int)(strchr(e, '\0') - e);
    if (l < 4) {
      d = s + 30;
      e += l;
      for (int i = 0; i < l; ++i)
        *--d = *--e;
      for (--d; isdigit(*d); --d)
        *d = '0';
    }
  }
  return fprintf(f, "%s,%5hu,%5hu,%5hu", s, a->i.p, a->i.q, a->i.b);
}

int wpqb_dump(FILE f[static 1], wpqb *const a, const uint32_t n_a)
{
  if (!n_a)
    return 0;
  if (!a)
    return -(int)(((n_a + 1u) << 1) + 1u);

  if (16 != fprintf(f, "\"w\",\"p\",\"q\",\"b\"\n"))
    return -1;
  for (uint32_t i = 0u; i < n_a; ++i) {
    if (48 != wpqb_print(f, (a + i)))
      return -(int)((i + 1u) << 1);
    if (1 != fprintf(f, "\n"))
      return -(int)(((i + 1u) << 1) + 1u);
  }

  return (fflush(f) ? -(int)((n_a + 1u) << 1) : (49 * (int)n_a + 16));
}

int wpqb_cmp(const wpqb a[static 1], const wpqb b[static 1])
{
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

#ifdef USE_OLD_NCP
int wpqb_ncp(const uint16_t n, wpqb *const a, const uint32_t n_a, uint16_t *const s, const uint16_t n_s, extended *const w)
{
  if (!n)
    return 0;
  if (!n_s)
    return 0;
  if (!s)
    return -4;
  if (!n_a)
    return 0;
  if (!a)
    return -2;

  uint16_t n_s_ = UINT16_C(0);
  uint32_t j = 0u;
  for (uint16_t i = UINT16_C(0); i < n_s; ) {
    s[i] = j;
    n_s_ = ++i;
    for (++j; j < n_a; ) {
      bool c = !(a[j].w == a[j].w);
      if (c) {
        ++j;
        continue;
      }
      const uint16_t ap = a[j].i.p;
      const uint16_t aq = a[j].i.q;
      for (uint16_t k = i; k; ) {
        --k;
        const uint16_t bp = a[s[k]].i.p;
        const uint16_t bq = a[s[k]].i.q;
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
    uint16_t r = UINT16_C(0), m = n_s_;
    for (uint16_t i = UINT16_C(0); i < n_s_; ++i) {
      if (a[s[i]].w > 0.0L)
        r = i + UINT16_C(1);
      else if (a[s[i]].w < 0.0L) {
        m = i;
        break;
      }
    }
    w[0] =  0.0L;
    w[1] = -0.0L;
    for (uint16_t i = r; i; )
      w[0] += a[s[--i]].w;
    for (uint16_t i = m; i < n_s_; ++i)
      w[1] += a[s[i]].w;
  }

  return (int)n_s_;
}
#else /* !USE_OLD_NCP */
int wpqb_ncp(const uint16_t n, wpqb *const a, const uint32_t n_a, uint16_t *const s, const uint16_t n_s, extended *const w)
{
  if (!n)
    return 0;
  if (!n_s)
    return 0;
  if (!s)
    return -4;
  if (!n_a)
    return 0;
  if (!a)
    return -2;

  bool p[n], q[n];
  (void)memset(p, 0, n * sizeof(bool));
  (void)memset(q, 0, n * sizeof(bool));

  uint16_t n_s_ = UINT16_C(0);
  uint32_t j = 0u;
  for (uint16_t i = UINT16_C(0); i < n_s; ) {
    s[i] = j;
    n_s_ = ++i;
    q[a[j].i.q] = p[a[j].i.p] = true;
    for (++j; j < n_a; ++j)
      if ((a[j].w == a[j].w) && !(p[a[j].i.p]) && !(q[a[j].i.q]))
        break;
    if (j >= n_a)
      break;
  }

  if (w) {
    uint16_t r = UINT16_C(0), m = n_s_;
    for (uint16_t i = UINT16_C(0); i < n_s_; ++i) {
      if (a[s[i]].w > 0.0L)
        r = i + UINT16_C(1);
      else if (a[s[i]].w < 0.0L) {
        m = i;
        break;
      }
    }
    w[0] =  0.0L;
    w[1] = -0.0L;
    for (uint16_t i = r; i; )
      w[0] += a[s[--i]].w;
    for (uint16_t i = m; i < n_s_; ++i)
      w[1] += a[s[i]].w;
  }

  return (int)n_s_;
}
#endif /* ?USE_OLD_NCP */
#endif /* USE_EXTENDED */
