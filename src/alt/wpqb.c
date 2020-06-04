#include "wpqb.h"

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
  if (a->w == a->w)
    return -5;
  // both weights are NaN
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

uint32_t wpqb_sort0(const uint32_t n_a, wpqb a[static 1])
{
  qsort(a, n_a, sizeof(wpqb), (int (*)(const void*, const void*))wpqb_cmp);
  for (uint32_t i = n_a; i; )
#ifdef USE_DOUBLE
    if (a[--i].w >= -DBL_MAX)
#else /* !USE_DOUBLE */
    if (a[--i].w >= -LDBL_MAX)
#endif /* ?USE_DOUBLE */
      return (i + 1u);
  return 0u;
}

uint32_t wpqb_sort1(const uint32_t n_a, wpqb a[static 1])
{
  if (!n_a)
    return 0u;
  const uint32_t na = (n_a - 1u);

  for (uint32_t swps = ~0u, s = 0u, oe = 0u; swps; (oe ^= 1u), (swps = s), (s = 0u)) {
    // odd-even sort
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(oe,na,a) reduction(+:s)
#endif /* _OPENMP */
    for (uint32_t i = oe; i < na; i += 2u) {
      const uint32_t j = (i + 1u);
      if (wpqb_cmp((a + i), (a + j)) > 0) {
        const wpqb t = a[i];
        a[i] = a[j];
        a[j] = t;
        ++s;
      }
    }
  }

  for (uint32_t i = n_a; i; )
#ifdef USE_DOUBLE
    if (a[--i].w >= -DBL_MAX)
#else /* !USE_DOUBLE */
    if (a[--i].w >= -LDBL_MAX)
#endif /* ?USE_DOUBLE */
      return (i + 1u);
  return 0u;
}

uint32_t wpqb_sort(const uint32_t n_a, wpqb a[static 1])
{
  return
#ifdef _OPENMP
    wpqb_sort1(n_a, a)
#else /* !_OPENMP */
    wpqb_sort0(n_a, a)
#endif /* ?_OPENMP */
    ;
}

void wpqb_ncp0(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  for (unsigned i = 0u; i < 10u; ++i)
    ((uint8_t*)w)[i] = UINT8_C(0xFF);
  w->i.s = UINT16_C(0);
  w->i.f = 0u;

  if (!n)
    return;
  if (!n_s)
    return;
  if (!n_a)
    return;

  bool r[n];
  for (uint16_t i = UINT16_C(0); i < n; ++i)
    r[i] = false;

  uint16_t ns = (n >> 1u);
  ns = ((ns <= n_a) ? ns : (uint16_t)n_a);
  ns = ((ns <= n_s) ? ns : n_s);

  for (uint32_t j = 0u, k; w->i.s < ns; j = k) {
    s[w->i.s] = j;
    r[a[j].i.p] = true;
    r[a[j].i.q] = true;
    for (k = j + 1u; k < n_a; ++k)
      if (!(r[a[k].i.p] || r[a[k].i.q]))
        break;
    w->i.s++;
    if (k >= n_a)
      break;
  }

#ifdef USE_DOUBLE
  w->w = 0.0;
#else /* !USE_DOUBLE */
  w->w = 0.0L;
#endif /* ?USE_DOUBLE */
  for (uint16_t i = w->i.s; i; )
    w->w += a[s[--i]].w;
}

static void wpqb_update(uint32_t s[static 1], const uint32_t my_s[static 1], wpqb_info w[static 1], const wpqb_info my_w[static 1])
{
  if (!(my_w->w <= w->w)) {
    for (uint16_t i = UINT16_C(0); i < my_w->i.s; ++i)
      s[i] = my_s[i] + my_w->i.f;
    *w = *my_w;
  }
}

static void wpqb_ncpt(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint32_t f, const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  uint32_t my_s[n_s];
  wpqb_info my_w;
  wpqb_ncp0(n, n_a, a, n_s, my_s, &my_w);
  if (my_w.i.s) {
    my_w.i.f = f;
#ifdef _OPENMP
#pragma omp critical
#endif /* _OPENMP */
    wpqb_update(s, my_s, w, &my_w);
  }
}

void wpqb_ncp1(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  for (unsigned i = 0u; i < 10u; ++i)
    ((uint8_t*)w)[i] = UINT8_C(0xFF);
  w->i.s = UINT16_C(0);
  w->i.f = 0u;

  if (!n)
    return;
  if (!n_s)
    return;
  if (!n_a)
    return;

#ifdef _OPENMP
#ifdef OLD_OMP
#pragma omp parallel for default(none) shared(n,n_a,a,n_s,s,w) schedule(dynamic,1)
#else /* !OLD_OMP */
#pragma omp parallel for default(none) shared(n,n_a,a,n_s,s,w) schedule(nonmonotonic:dynamic,1)
#endif /* ?OLD_OMP */
#endif /* _OPENMP */
  for (uint32_t i = 0u; i < n_a; ++i)
    wpqb_ncpt(n, n_a - i, a + i, i, n_s, s, w);
}

void wpqb_ncp(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
#ifdef _OPENMP
  wpqb_ncp1(n, n_a, a, n_s, s, w);
#else /* !_OPENMP */
  wpqb_ncp0(n, n_a, a, n_s, s, w);
#endif /* ?_OPENMP */
}

void wpqb_run0(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  for (unsigned i = 0u; i < 10u; ++i)
    ((uint8_t*)w)[i] = UINT8_C(0xFF);
  w->i.s = UINT16_C(0);

  if (!(w->i.f = wpqb_sort0(n_a, a)))
    return;
  wpqb_ncp0(n, w->i.f, a, n_s, s, w);
}

void wpqb_run1(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
  for (unsigned i = 0u; i < 10u; ++i)
    ((uint8_t*)w)[i] = UINT8_C(0xFF);
  w->i.s = UINT16_C(0);

  if (!(w->i.f = wpqb_sort1(n_a, a)))
    return;
  wpqb_ncp1(n, w->i.f, a, n_s, s, w);
}

void wpqb_run(const uint16_t n, const uint32_t n_a, wpqb a[static 1], const uint16_t n_s, uint32_t s[static 1], wpqb_info w[static 1])
{
#ifdef _OPENMP
  wpqb_run1(n, n_a, a, n_s, s, w);
#else /* !_OPENMP */
  wpqb_run0(n, n_a, a, n_s, s, w);
#endif /* ?_OPENMP */
}

#ifndef NDEBUG
int main(int argc, char *argv[])
{
  if (4 != argc) {
    (void)fprintf(stderr, "%s n n_s r\n", *argv);
    return EXIT_FAILURE;
  }

  const uint16_t n = (uint16_t)atoi(argv[1]);
  if (!n)
    return EXIT_FAILURE;
  const uint16_t n_s = (uint16_t)atoi(argv[2]);
  if (!n_s)
    return EXIT_FAILURE;
  const unsigned r = (unsigned)atoi(argv[3]);
  srand(r & ~1u);

  const uint32_t n_a = (((uint32_t)n * ((uint32_t)n - 1u)) >> 1u);
  wpqb *const a = (wpqb*)calloc(n_a, sizeof(wpqb));
  if (!a)
    return EXIT_FAILURE;

  uint32_t i = 0u;
  const uint16_t n_1 = n - UINT16_C(1);
  for (uint16_t p = UINT16_C(0); p < n_1; ++p)
    for (uint16_t q = p + UINT16_C(1); q < n; ++q, ++i)
#ifdef USE_DOUBLE
      wpqb_init((a + i), (((double)rand()) / rand()), p, q, 1, 1);
#else /* !USE_DOUBLE */
      wpqb_init((a + i), (((long double)rand()) / rand()), p, q);
#endif /* ?USE_DOUBLE */

  uint32_t *const s = (uint32_t*)calloc(n_s, sizeof(uint32_t));
  if (!s)
    return EXIT_FAILURE;

  wpqb_info w;
  ((r & 1u) ? wpqb_run1 : wpqb_run0)(n, n_a, a, n_s, s, &w);

#if (defined(_WIN32) || defined(USE_DOUBLE))
  (void)fprintf(stdout, "w=%# .17e\n", (double)(w.w));
#else /* !(_WIN32 || USE_DOUBLE) */
  (void)fprintf(stdout, "w=%# .21Le\n", w.w);
#endif /* ?(_WIN32 || USE_DOUBLE) */
  (void)fprintf(stdout, "s=%hu\n", w.i.s);
  (void)fprintf(stdout, "f=%u\n", w.i.f);
  for (i = 0u; i < w.i.s; ++i)
#if (defined(_WIN32) || defined(USE_DOUBLE))
    (void)fprintf(stdout, "%u=(%# .17e,%hu,%hu,%hu)\n", i, (double)(a[s[i]].w), a[s[i]].i.p, a[s[i]].i.q, a[s[i]].i.b);
#else /* !(_WIN32 || USE_DOUBLE) */
    (void)fprintf(stdout, "%u=(%# .21Le,%hu,%hu,%hu)\n", i, a[s[i]].w, a[s[i]].i.p, a[s[i]].i.q, a[s[i]].i.b);
#endif /* ?(_WIN32 || USE_DOUBLE) */

  free(s);
  free(a);
  return EXIT_SUCCESS;
}
#endif /* !NDEBUG */
