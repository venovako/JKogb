# JKogb
The J-Kogbetliantz algorithm for the hyperbolic singular value decomposition (HSVD)

This software is a supplementary material for the paper
arXiv:[2003.06701](https://arxiv.org/abs/2003.06701 "A Kogbetliantz-type algorithm for the hyperbolic SVD") \[math.NA\].

## Building

### Prerequisites

A recent 64-bit Linux (e.g., CentOS 7.8 with devtoolset-9) or macOS (e.g., Catalina) is needed.

Have the Intel MKL (Math Kernel Library) installed (only required as a JACSD's dependency).
Other (sequential) BLAS and LAPACK libraries might work with some makefile tweaking, if they support 8-byte INTEGERs.

Then, clone and build [JACSD](https://github.com/venovako/JACSD) in a directory parallel to this one.

### Make options

Run ``make`` as follows:
```bash
cd src
make [COMPILER=gnu|x64|x200|nvidia] [NDEBUG=0|1|2|3|4|5] [all|clean|help]
```
where ``COMPILER`` should be set for the Intel Fortran compilers (version 19.1+/2020+ recommended) to ``x64`` for Xeons, or to ``x200`` for Xeon Phi KNLs, respectively.
If ``COMPILER`` is not set or is ``gnu``, GNU Fortran compilers will be used (versions 8, 9, and 10 should work, but not extensively tested).
There is also a preliminary support for NVIDIA HPC SDK on Linux (building works, but the executables have not been tested).

Here, ``NDEBUG`` should be set to the desired optimization level (``3`` is a sensible choice).
If unset, the predefined debug-mode build options will be used.

For example, ``make COMPILER=x200 NDEBUG=3 clean all`` will trigger a full, release-mode rebuild for the KNLs.

## Execution

### Command line

In the examples below, ``FN`` is the input and output file name prefix (without an extension), ``N`` for the matrix order, and ``N_2`` for the maximal number of transformations per step (up to ``N/2``).

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
