# JKogb
J-Kogbetliantz algorithm for the hyperbolic singular value decomposition (HSVD)

This software is a supplementary material for the paper
arXiv:[2003.06701](https://arxiv.org/abs/2003.06701 "A Kogbetliantz-type algorithm for the hyperbolic SVD") \[math.NA\].

A recent 64-bit Linux (e.g., CentOS 7.8 with devtoolset-8) or macOS (e.g., Catalina) is needed.

Have the Intel MKL (Math Kernel Library) installed.
Other (sequential) BLAS and LAPACK libraries might work with some makefile tweaking, if they support 8-byte INTEGERs.
Intel C/C++ and Fortran compilers (version 19.1+/2020+) are recommended.

Then, clone and build [JACSD](https://github.com/venovako/JACSD) in a directory parallel to this one.

Makefiles for the old and unmaintained code are included for Windows and macOS; for Linux, please edit `src/old/GNUmakefile`.

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
