#include "wpqb.h"

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
  if (r)
    srand(r);

  const uint32_t n_a = ((n * ((uint32_t)n - 1u)) >> 1u);
  wpqb *const a = (wpqb*)calloc(n_a, sizeof(wpqb));
  if (!a)
    return EXIT_FAILURE;

  uint32_t i = 0u;
  const uint16_t n_1 = (uint16_t)(n - UINT16_C(1));
  for (uint16_t p = UINT16_C(0); p < n_1; ++p)
    for (uint16_t q = (uint16_t)(p + UINT16_C(1)); q < n; ++q, ++i)
      wpqb_init((a + i), (((long double)rand()) / rand()), p, q);

  uint32_t *const s = (uint32_t*)calloc(n_s, sizeof(uint32_t));
  if (!s)
    return EXIT_FAILURE;

  wpqb_info w;
  wpqb_ncp(n, n_a, a, n_s, s, &w);

#ifdef _WIN32
  (void)fprintf(stdout, "w=%# .17e, ", (double)(w.w));
#else /* !_WIN32 */
  (void)fprintf(stdout, "w=%# .21Le, ", w.w);
#endif /* ?_WIN32 */
  (void)fprintf(stdout, "s=%hu, f=%u\n", w.i.s, w.i.f);
  for (i = 0u; i < w.i.s; ++i)
#ifdef _WIN32
    (void)fprintf(stdout, "%u=(%# .17e,%hu,%hu,%hu)\n", i, (double)(a[s[i]].w), a[s[i]].i.p, a[s[i]].i.q, a[s[i]].i.b);
#else /* !_WIN32 */
    (void)fprintf(stdout, "%u=(%# .21Le,%hu,%hu,%hu)\n", i, a[s[i]].w, a[s[i]].i.p, a[s[i]].i.q, a[s[i]].i.b);
#endif /* ?_WIN32 */

  free(s);
  free(a);
  return EXIT_SUCCESS;
}
