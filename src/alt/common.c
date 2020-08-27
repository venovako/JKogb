#include "common.h"

static _Alignas(long double) const unsigned char all_ones[16] = { // 10 + 6
  UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF),
  UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF), UINT8_C(0xFF) };

long double qNaN()
{
  return *(const long double*)all_ones;
}
