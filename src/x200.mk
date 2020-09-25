SHELL=/bin/bash
ARCH=$(shell uname)
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
RM=rm -rfv
AR=xiar
ARFLAGS=-qnoipo -lib rsv
# CC=icc -std=c18
FC=ifort
# CXX=icpc -std=c++17
CPUFLAGS=-DUSE_INTEL -DUSE_X200 -fPIC -fexceptions -fno-omit-frame-pointer -qopenmp -rdynamic
ifdef PROFILE
CPUFLAGS += -DVN_PROFILE=$(PROFILE) -fno-inline -finstrument-functions
endif # PROFILE
FORFLAGS=$(CPUFLAGS) -i8 -standard-semantics -threads #-cxxlib
# C18FLAGS=$(CPUFLAGS)
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -xHost -qopt-zmm-usage=high
OPTFFLAGS=$(OPTFLAGS)
# OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-DNDEBUG -qopt-report=5 -traceback -diag-disable=10397
DBGFFLAGS=$(DBGFLAGS)
# DBGCFLAGS=$(DBGFLAGS) -w3 -diag-disable=1572,2547
FPUFLAGS=-fma -fp-model source -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fimf-precision=high -fimf-use-svml=true
FPUFFLAGS=$(FPUFLAGS)
# FPUCFLAGS=$(FPUFLAGS)
else # DEBUG
OPTFLAGS=-O0 -xHost -qopt-zmm-usage=high
OPTFFLAGS=$(OPTFLAGS)
# OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug parallel -debug pubnames -traceback -diag-disable=10397
DBGFFLAGS=$(DBGFLAGS) -debug-parameters all -check all -warn all
# DBGCFLAGS=$(DBGFLAGS) -check=stack,uninit -w3 -diag-disable=1572,2547
FPUFLAGS=-fp-model strict -fp-stack-check -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fimf-precision=high
FPUFFLAGS=$(FPUFLAGS) -assume ieee_fpe_flags
# FPUCFLAGS=$(FPUFLAGS)
endif # ?NDEBUG
LIBFLAGS=-D_GNU_SOURCE -I. -I../shared -I../../../JACSD/vn
ifdef ANIMATE
LIBFLAGS += -DUSE_MKL -DMKL_ILP64 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
endif # ANIMATE
LDFLAGS=-L../shared -ljk$(PROFILE)$(DEBUG) -L../../../JACSD -lvn$(PROFILE)$(DEBUG)
ifdef ANIMATE
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core #-lmkl_intel_thread
endif # ANIMATE
LDFLAGS += -lpthread -lm -ldl -lmemkind
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
# CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUCFLAGS)
