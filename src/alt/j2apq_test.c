#include "j2apq.h"

int main(int argc, char *argv[])
{
  if (argc != 3) {
    (void)fprintf(stderr, "%s n r\n", *argv);
    return EXIT_FAILURE;
  }
  const uint16_t n = (uint16_t)atoi(argv[1]);
  const unsigned r = (unsigned)atoi(argv[2]);

  fint *j = (fint*)NULL;

  if (r) {
    if (!(j = calloc(n, sizeof(fint))))
      return EXIT_FAILURE;
    srand(r);
    (void)fprintf(stdout, "j= ");
    for (uint16_t i = UINT16_C(0); i < n; ++i) {
      j[i] = ((rand() & 1) ? FINT_C(-1) : FINT_C(1));
      (void)fprintf(stdout, "%3d", (int)(j[i]));
    }
    (void)fprintf(stdout, "\n");
  }

  apq *const o = apq_create(n, j);
  if (!o)
    return EXIT_FAILURE;

  if (r)
    for (uint32_t i = 0u; i < o->n_a; ++i)
      (void)fprintf(stdout, "%5u %5u %2d %2d\n", (unsigned)((o->p)[i]), (unsigned)((o->q)[i]), (int)(j[(o->p)[i]]), (int)(j[(o->q)[i]]));
  else
    for (uint32_t i = 0u; i < o->n_a; ++i)
      (void)fprintf(stdout, "%5u %5u %2d %2d\n", (unsigned)((o->p)[i]), (unsigned)((o->q)[i]), 1, 1);

  (void)apq_free(o);
  free(j);
  return EXIT_SUCCESS;
}
