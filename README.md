# JKogb
J-Kogbetliantz algorithm for the hyperbolic singular value decomposition (HSVD)

...work in progress...

A recent 64-bit Linux (e.g., CentOS 7.6) or macOS (e.g., Mojave) is needed.

Have the Intel MKL (Math Kernel Library) and TBB (Thread Building Blocks) installed.

Then, clone and build [JACSD](https://github.com/venovako/JACSD) in a directory parallel to this one.

GNU Fortran 9 is *not* supported for the code in development!
Currently, only GPU Fortran *8* is fully supported.

Makefiles for the old and unmaintained code are included for Windows and macOS; for Linux, please edit `src/old/GNUmakefile`.

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
