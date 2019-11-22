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
CPUFLAGS=-DUSE_GNU -DUSE_X64 -fexceptions -fopenmp
ifdef PROFILE
CPUFLAGS += -DVN_PROFILE=$(PROFILE) -fPIC -fno-inline -fno-omit-frame-pointer -finstrument-functions -rdynamic
endif # PROFILE
FORFLAGS=-cpp $(CPUFLAGS) -fdefault-integer-8 -ffree-line-length-none -fstack-arrays
C11FLAGS=$(CPUFLAGS) -fopenmp
ifeq ($(ARCH),Darwin)
FC=gfortran-9
CC=gcc-9
CXX=g++-9
else # Linux
FC=gfortran
CC=gcc
CXX=g++
endif # ?Darwin
CC += -std=gnu17
CXX += -std=gnu++17
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -march=native -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller
DBGFLAGS=-DNDEBUG -fopt-info-optimized-vec
ifeq ($(ARCH),Darwin)
OPTFLAGS += -Wa,-q
endif # Darwin
OPTFFLAGS=$(OPTFLAGS) 
OPTCFLAGS=$(OPTFLAGS)
DBGFFLAGS=$(DBGFLAGS) -pedantic -Wall -Wextra -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
DBGCFLAGS=$(DBGFLAGS)
FPUFLAGS=-ffp-contract=fast
FPUFFLAGS=$(FPUFLAGS)
FPUCFLAGS=$(FPUFLAGS)
else # DEBUG
OPTFLAGS=-O$(DEBUG) -march=native
ifeq ($(ARCH),Darwin)
OPTFLAGS += -Wa,-q
endif # Darwin
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-$(DEBUG)
DBGFFLAGS=$(DBGFLAGS) -fcheck=all -finit-local-zero -finit-real=snan -finit-derived -pedantic -Wall -Wextra -Wno-compare-reals -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
DBGCFLAGS=$(DBGFLAGS) -ftrapv
FPUFLAGS=-ffp-contract=fast
FPUFFLAGS=$(FPUFLAGS) -ffpe-trap=invalid,zero,overflow
FPUCFLAGS=$(FPUFLAGS)
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared -I../../../JACSD/vn
ifdef ANIMATE
LIBFLAGS += -DUSE_MKL -DMKL_ILP64 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
endif # ANIMATE
LDFLAGS=-L../shared -ljk$(PROFILE)$(DEBUG) -L../../../JACSD -lvn$(PROFILE)$(DEBUG)
ifeq ($(ARCH),Darwin)
ifdef ANIMATE
LDFLAGS += -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core #-lmkl_intel_thread -liomp5
endif # ANIMATE
else # Linux
LIBFLAGS += -D_GNU_SOURCE
ifdef ANIMATE
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core #-lmkl_gnu_thread -lgomp
endif # ANIMATE
endif # ?Darwin
LDFLAGS += -lpthread -lm -ldl $(shell if [ -L /usr/lib64/libmemkind.so ]; then echo '-lmemkind'; fi)
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C11FLAGS) $(FPUCFLAGS)
