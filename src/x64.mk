SHELL=/bin/bash
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
ARCH=$(shell uname)
RM=rm -rfv
AR=xiar
ARFLAGS=-qnoipo -lib rsv
CC=icc -std=c11
FC=ifort
CXX=icpc -std=c++17
CPUFLAGS=-DUSE_INTEL -DUSE_X64 -fexceptions -qopenmp
ifdef PROFILE
CPUFLAGS += -DVN_PROFILE=$(PROFILE) -fPIC -fno-inline -fno-omit-frame-pointer -finstrument-functions -rdynamic
endif # PROFILE
FORFLAGS=$(CPUFLAGS) -i8 -standard-semantics -cxxlib -threads
C11FLAGS=$(CPUFLAGS)
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -xHost
OPTFFLAGS=$(OPTFLAGS) -DMKL_DIRECT_CALL_SEQ_JIT
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-DNDEBUG -qopt-report=5 -traceback -diag-disable=10397
DBGFFLAGS=$(DBGFLAGS)
DBGCFLAGS=$(DBGFLAGS) -w3 -diag-disable=1572,2547
FPUFLAGS=-fma -fp-model source -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fimf-precision=high -fimf-use-svml=true
FPUFFLAGS=$(FPUFLAGS)
FPUCFLAGS=$(FPUFLAGS)
else # DEBUG
OPTFLAGS=-O0 -xHost
OPTFFLAGS=$(OPTFLAGS)
OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -traceback -diag-disable=10397
ifneq ($(ARCH),Darwin) # Linux
DBGFLAGS += -debug parallel
endif # ?Linux
DBGFFLAGS=$(DBGFLAGS) -debug-parameters all -check all -warn all
DBGCFLAGS=$(DBGFLAGS) -check=stack,uninit -w3 -diag-disable=1572,2547
FPUFLAGS=-fp-model strict -fp-stack-check -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fimf-precision=high
FPUFFLAGS=$(FPUFLAGS) -assume ieee_fpe_flags
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
LDFLAGS += -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core #-lmkl_intel_thread
else # Linux
LIBFLAGS += -D_GNU_SOURCE
LDFLAGS += -L${TBBROOT}/lib -Wl,-rpath=${TBBROOT}/lib
ifdef NDEBUG
LDFLAGS += -ltbb -ltbbmalloc
else # DEBUG
LDFLAGS += -ltbb_debug -ltbbmalloc_debug
endif # ?NDEBUG
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core #-lmkl_intel_thread
endif # ?Darwin
LDFLAGS += -lpthread -lm -ldl
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C11FLAGS) $(FPUCFLAGS)
CXXFLAGS=$(CFLAGS)
