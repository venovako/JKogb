SHELL=/bin/bash
ARCH=$(shell uname)
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
RM=rm -rfv
AR=ar
ARFLAGS=rsv
# comment out USE_X64 if not on Intel 64
CPUFLAGS=-DUSE_GNU -DUSE_X64 -fopenmp
FORFLAGS=-cpp -DMIND=C_FMIN -DMAXD=C_FMAX $(CPUFLAGS) -ffree-line-length-none -fstack-arrays
ifdef ANIMATE
FORFLAGS += -fdefault-integer-8
endif # ANIMATE
ifeq ($(ARCH),Darwin)
CPUFLAGS += -DUSE_MACOS
endif # Darwin
FPUFLAGS=-ffp-contract=fast
FC=gfortran
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller -fopt-info-optimized-vec #-DUSE_FAST
ifeq ($(ARCH),Darwin)
OPTFLAGS += -Wa,-q
endif # Darwin
DBGFLAGS=-DNDEBUG -pedantic -Wall -Wextra -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
else # DEBUG
OPTFLAGS=-O$(DEBUG) -march=native -fopt-info-optimized-vec
ifeq ($(ARCH),Darwin)
OPTFLAGS += -Wa,-q
endif # ?Darwin
DBGFLAGS=-$(DEBUG) -pedantic -Wall -Wextra -fcheck=all,no-recursion -finit-local-zero -finit-real=nan -finit-derived -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared -static-libgcc -static-libgfortran
ifdef ANIMATE
LIBFLAGS += -DUSE_MKL -DMKL_ILP64 -I../../../JACSD/vn -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
endif # ANIMATE
LDFLAGS=-L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
ifeq ($(ARCH),Darwin)
ifdef ANIMATE
LDFLAGS += -L../../../JACSD -lvn$(DEBUG) -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
endif # ANIMATE
else # Linux
ifdef ANIMATE
LDFLAGS += -L../../../JACSD -lvn$(DEBUG) -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core
endif # ANIMATE
endif # ?Darwin
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
