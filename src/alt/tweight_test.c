#include "tweight.h"

int main(int argc, char *argv[])
{
  long double w = -0.0L;
  if (argc == 3) {
    const double Aqp = atof(argv[1]);
    const double Apq = atof(argv[2]);
    w = dtw(Aqp, Apq);
  }
  else if (argc == 5) {
    const double Aqp_r = atof(argv[1]);
    const double Aqp_i = atof(argv[2]);
    const double Apq_r = atof(argv[3]);
    const double Apq_i = atof(argv[4]);
    w = ztw(CMPLX(Aqp_r, Aqp_i), CMPLX(Apq_r, Apq_i));
  }
  else {
    (void)fprintf(stderr, "%s Aqp_r [Aqp_i] Apq_r [Apq_i]\n", *argv);
    return EXIT_FAILURE;
  }
#ifdef _WIN32
  (void)fprintf(stdout, "tw=%# .17e\n", (double)w);
#else /* !_WIN32 */
  (void)fprintf(stdout, "tw=%# .21Le\n", w);
#endif /* ?_WIN32 */
  return EXIT_SUCCESS;
}
