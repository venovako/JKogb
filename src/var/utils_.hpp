#ifndef UTILS_HPP
#define UTILS_HPP

#include <cstdint>
#include <cstdlib>

struct aw {
  double w;
  int64_t p, q, b;
};

#ifdef _GNU_SOURCE
typedef int (*qcmp)(const aw *const, const aw *const, void *const);
#else // !_GNU_SOURCE
typedef int (*qcmp)(void *const, const aw *const, const aw *const);
#endif // ?_GNU_SOURCE

/* Fortran interface
  ! INTERFACE
  !    RECURSIVE FUNCTION PAR_SORT(T, A, N, C, X) BIND(C,NAME='par_sort')
  !      USE, INTRINSIC :: ISO_C_BINDING
  !      IMPLICIT NONE
  !      INTEGER(KIND=c_int), INTENT(IN), VALUE :: T
  !      TYPE(c_ptr), INTENT(IN), VALUE :: A, X
  !      INTEGER(KIND=c_size_t), INTENT(IN), VALUE :: N
  !      TYPE(c_funptr), INTENT(IN), VALUE :: C
  !      INTEGER(KIND=c_int) :: PAR_SORT
  !    END FUNCTION PAR_SORT
  ! END INTERFACE
 */
extern "C" int par_sort(const int t, aw *const a, const size_t n, const qcmp c, void *const x) throw();

#endif // !UTILS_HPP
