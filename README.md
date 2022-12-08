# JKogb
The J-Kogbetliantz algorithm for the hyperbolic singular value decomposition (HSVD)

This software is a supplementary material for the paper
doi:[10.1007/s11075-021-01197-4](https://doi.org/10.1007/s11075-021-01197-4 "A Kogbetliantz-type algorithm for the hyperbolic SVD")
(arXiv:[2003.06701](https://arxiv.org/abs/2003.06701 "A Kogbetliantz-type algorithm for the hyperbolic SVD") \[math.NA\]).

## Building

### Prerequisites

A recent 64-bit Linux (e.g., CentOS 7.9 with devtoolset-8), macOS (e.g., Big Sur), or Windows (e.g., 10) is needed.

### Make options

On Linux or macOS, run ``make`` (GNU make assumed) or ``gmake`` as follows:
```bash
cd src
make [COMPILER=gnu|x64|x200] [NDEBUG=0|1|2|3|4|5] [all|clean|help]
```
where ``COMPILER`` should be set for the Intel Fortran compiler to ``x64`` for Xeons, or to ``x200`` for Xeon Phi KNLs, respectively.
Building with a recent Intel Fortran is possible on Windows as well:
```bash
cd src
nmake [NDEBUG=d|1|2|3|4|5] [all|clean|help]
```

If ``COMPILER`` is not set or is ``gnu``, the GNU Fortran compiler will be used.
The major version of your GCC in this case should be 8, since the later ones will not work unless a fix is applied as noted in ``src/gnu.mk``.
Please see [this explanation](https://gcc.gnu.org/gcc-9/changes.html) regarding the new MIN and MAX intrinsics.

Here, ``NDEBUG`` should be set to the desired optimization level (``3`` is a sensible choice).
If unset, the predefined debug-mode build options will be used.

For example, ``make COMPILER=x200 NDEBUG=3 clean all`` will trigger a full, release-mode rebuild for the KNLs.

## Execution

### Command line

In the examples below, ``FN`` is the input and output file name prefix (without an extension), ``N`` for the matrix order, and ``N_2`` for the maximal number of transformations per step (up to ``N/2``, with the fail-safe but the slowest option being ``1``).

Assuming the top directory of the cloned repository has been made current, execute either
```bash
src/D/djk.exe FN N N_2
```
for the real, or
```bash
src/Z/zjk.exe FN N N_2
```
for the complex variant.
The standard OpenMP environment variables can be used to set the desired number and/or placement of threads.

### Data format

``FN.Y`` should be a binary (``DOUBLE PRECISION`` or ``DOUBLE COMPLEX``) file containing the input matrix in the Fortran (column-major) array order.
``FN.J`` should be a binary file containing a vector of ``N`` 8-byte ``INTEGER``s (either `1` or `-1`), representing the diagonal of the sign matrix.
The outputs will be stored in the binary files ``FN.SY`` (the hyperbolic singular values), ``FN.YU`` (for the left), and ``FN.ZZ`` (for the right hyperbolic singular vectors).

Makefiles for the old and unmaintained code are included for Windows and macOS; for Linux, please edit `src/old/GNUmakefile`.

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
