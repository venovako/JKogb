SHELL=/bin/bash
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
ARCH=$(shell uname)
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
CC=gcc-9
FC=gfortran-9
CXX=g++-9
else # Linux
CC=gcc
FC=gfortran
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
LIBFLAGS=-DUSE_MKL -DMKL_ILP64 -I. -I../shared -I../../../JACSD/vn -I${TBBROOT}/include -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
LDFLAGS=-L../shared -ljk$(PROFILE)$(DEBUG) -L../../../JACSD -lvn$(PROFILE)$(DEBUG)
ifndef NDEBUG
LIBFLAGS += -DTBB_USE_DEBUG=1
endif # !NDEBUG
ifeq ($(ARCH),Darwin)
LDFLAGS += -L${TBBROOT}/lib -Wl,-rpath,${TBBROOT}/lib
ifdef NDEBUG
LDFLAGS += -ltbb -ltbbmalloc
else # DEBUG
LDFLAGS += -ltbb_debug -ltbbmalloc_debug
endif # ?NDEBUG
LDFLAGS += -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${MKLROOT}/../compiler/lib -Wl,-rpath,${MKLROOT}/../compiler/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core #-liomp5 -lmkl_intel_thread
else # Linux
LIBFLAGS += -D_GNU_SOURCE
LDFLAGS += -L${TBBROOT}/lib -Wl,-rpath=${TBBROOT}/lib
ifdef NDEBUG
LDFLAGS += -ltbb -ltbbmalloc
else # DEBUG
LDFLAGS += -ltbb_debug -ltbbmalloc_debug
endif # ?NDEBUG
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_ilp64 -lmkl_sequential -lmkl_core #-lgomp -lmkl_gnu_thread
endif # ?Darwin
LDFLAGS += -lstdc++ -lpthread -lm -ldl
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C11FLAGS) $(FPUCFLAGS)
CXXFLAGS=$(CFLAGS)
