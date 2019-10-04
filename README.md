# JKogb
J-Kogbetliantz algorithm for the hyperbolic singular value decomposition (HSVD)

...work in progress...

A recent 64-bit Linux (e.g., CentOS 7.7) or macOS (e.g., Mojave) is needed.

Have the Intel MKL (Math Kernel Library) installed.
Other (sequential) BLAS and LAPACK libraries might work with some makefile tweaking, if they support 8-byte INTEGERs.

Then, clone and build [JACSD](https://github.com/venovako/JACSD) in a directory parallel to this one.

Makefiles for the old and unmaintained code are included for Windows and macOS; for Linux, please edit `src/old/GNUmakefile`.

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
