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
ifndef MARCH
ifeq ($(shell uname -m),ppc64le)
MARCH=mcpu=native -mpower8-fusion -mtraceback=full
else # !ppc64le
MARCH=march=native
endif # ?ppc64le
endif # !MARCH
CPUFLAGS=-DUSE_GNU -DUSE_X64 -fPIC -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -fvect-cost-model=unlimited -fopenmp -$(MARCH)
ifeq ($(ARCH),Darwin)
CPUFLAGS += -Wa,-q
endif # Darwin
FORFLAGS=-cpp -DMIND=C_FMIN -DMAXD=C_FMAX $(CPUFLAGS) -ffree-line-length-none -fstack-arrays
ifeq ($(ARCH),Darwin)
CPUFLAGS += -DUSE_MACOS
endif # Darwin
FPUFLAGS=-ffp-contract=fast -fno-math-errno
FC=gfortran$(GNU)
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) #-DUSE_FAST
DBGFLAGS=-DNDEBUG -Wall -Wextra -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
else # DEBUG
OPTFLAGS=-O$(DEBUG)
DBGFLAGS=-$(DEBUG)gdb3 -Wall -Wextra -fcheck=all,no-recursion -finit-local-zero -finit-real=nan -finit-derived -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared
LDFLAGS=-rdynamic -L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
ifdef ANIMATE
LIBPVN=$(realpath ../../../libpvn)
LDFLAGS += -L$(LIBPVN)/src
ifeq ($(ARCH),Darwin)
LDFLAGS += -Wl,-rpath,$(LIBPVN)/src
else # !Darwin
LDFLAGS += -Wl,-rpath=$(LIBPVN)/src
endif # ?Darwin
LDFLAGS += -lpvn -ldl -lm
endif # ANIMATE
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
