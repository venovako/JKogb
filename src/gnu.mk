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
CPUFLAGS=-DUSE_GNU -DUSE_X64 -fPIC -fexceptions -fno-omit-frame-pointer -fopenmp -rdynamic
ifdef PROFILE
CPUFLAGS += -DVN_PROFILE=$(PROFILE) -fno-inline -finstrument-functions
endif # PROFILE
FORFLAGS=-cpp $(CPUFLAGS) -ffree-line-length-none -fstack-arrays
ifdef ANIMATE
FORFLAGS += -fdefault-integer-8
endif # ANIMATE
ifeq ($(ARCH),Darwin)
ifdef GNU
ifneq ($(GNU),-8)
FORFLAGS += -DMIND=C_FMIN -DMAXD=C_FMAX
endif # GNU > 8
else # !GNU
GNU=-8
endif # ?GNU
endif # Darwin
FC=gfortran$(GNU)
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller
DBGFLAGS=-DNDEBUG -fopt-info-optimized-vec -pedantic -Wall -Wextra
ifeq ($(ARCH),Darwin)
OPTFLAGS += -Wa,-q
endif # Darwin
OPTFFLAGS=$(OPTFLAGS)
DBGFFLAGS=$(DBGFLAGS) -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
FPUFLAGS=-ffp-contract=fast
FPUFFLAGS=$(FPUFLAGS)
else # DEBUG
OPTFLAGS=-O$(DEBUG) -march=native
DBGFLAGS=-$(DEBUG) -fsanitize=address -pedantic -Wall -Wextra
ifeq ($(ARCH),Darwin)
OPTFLAGS += -Wa,-q
else # Linux
DBGFLAGS += -fsanitize=leak
endif # ?Darwin
OPTFFLAGS=$(OPTFLAGS)
DBGFFLAGS=$(DBGFLAGS) -fcheck=array-temps -finit-local-zero -finit-real=snan -finit-derived -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all #-fcheck=all
FPUFLAGS=-ffp-contract=fast
FPUFFLAGS=$(FPUFLAGS)
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared
ifdef ANIMATE
LIBFLAGS += -DUSE_MKL -DMKL_ILP64 -I../../../JACSD/vn -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
endif # ANIMATE
LDFLAGS=-L../shared -ljk$(PROFILE)$(DEBUG)
ifeq ($(ARCH),Darwin)
ifdef ANIMATE
LDFLAGS += -L../../../JACSD -lvn$(PROFILE)$(DEBUG) -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
endif # ANIMATE
else # Linux
LIBFLAGS += -D_GNU_SOURCE
ifdef ANIMATE
LDFLAGS += -L../../../JACSD -lvn$(PROFILE)$(DEBUG) -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core
endif # ANIMATE
endif # ?Darwin
ifndef NDEBUG
LDFLAGS += -lubsan
endif # DEBUG
LDFLAGS += -lpthread -lm -ldl $(shell if [ -L /usr/lib64/libmemkind.so ]; then echo '-lmemkind'; fi)
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
