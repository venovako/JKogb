SHELL=/bin/bash
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
C11FLAGS=-std=gnu11 $(CPUFLAGS)
FORFLAGS=-cpp $(CPUFLAGS) -ffree-line-length-none -fstack-arrays
# uncomment the definition below if the GCC version is above 8
# FORFLAGS += -DMIND=C_FMIN -DMAXD=C_FMAX
FPUFLAGS=-ffp-contract=fast
CC=gcc$(GNU)
FC=gfortran$(GNU)
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller -fopt-info-optimized-vec #-DUSE_FAST
DBGFLAGS=-DNDEBUG -pedantic -Wall -Wextra #-DUSE_TEST
DBGFFLAGS=$(DBGFLAGS) -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
else # DEBUG
OPTFLAGS=-O$(DEBUG) -march=native -fopt-info-optimized-vec
DBGFLAGS=-$(DEBUG) -pedantic -Wall -Wextra
DBGFFLAGS=$(DBGFLAGS) -fcheck=all,no-recursion -finit-local-zero -finit-real=nan -finit-derived -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
endif # ?NDEBUG
OPTFFLAGS=$(OPTFLAGS)
FPUFFLAGS=$(FPUFLAGS)
LIBFLAGS=-I. -I../shared -static-libgcc
LDFLAGS=-L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(C11FLAGS) $(FPUFLAGS)
