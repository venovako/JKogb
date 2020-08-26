#include "common.h"

static alignas(long double) const uint64_t all_ones[2] = { UINT64_C(0xFFFFFFFFFFFFFFFF), UINT64_C(0xFFFFFFFFFFFFFFFF) };

long double qNaN()
{
  return *(const long double*)all_ones;
}
