#include "common.h"

int main(int argc, char *argv[])
{
  if (argc != 1) {
    (void)fprintf(stderr, "%s\n", *argv);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
