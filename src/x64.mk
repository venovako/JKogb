SHELL=/bin/bash
ARCH=$(shell uname)
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
ifndef FP
ifdef NDEBUG
FP=source
else # DEBUG
FP=source #strict
endif # ?NDEBUG
endif # !FP
RM=rm -rfv
AR=xiar
ARFLAGS=-qnoipo -lib rsv
# CC=icc -std=c18
FC=ifort
# CXX=icpc -std=c++17
CPUFLAGS=-DUSE_INTEL -DUSE_X64 -fPIC -fexceptions -fno-omit-frame-pointer -qopenmp -rdynamic
ifneq ($(ARCH),Darwin) # Linux
CPUFLAGS += -qopenmp-threadprivate=compat
endif # Linux
ifdef PROFILE
CPUFLAGS += -DVN_PROFILE=$(PROFILE) -fno-inline -finstrument-functions
endif # PROFILE
FORFLAGS=$(CPUFLAGS) -i8 -standard-semantics -threads #-cxxlib
# C18FLAGS=$(CPUFLAGS)
FPUFLAGS=-fp-model $(FP) -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fimf-precision=high
ifeq ($(FP),strict)
FPUFLAGS += -fp-stack-check -fimf-arch-consistency=true
else # !strict
FPUFLAGS += -fimf-use-svml=true
endif # ?strict
FPUFFLAGS=$(FPUFLAGS)
# FPUCFLAGS=$(FPUFLAGS)
ifeq ($(FP),strict)
FPUFFLAGS += -assume ieee_fpe_flags
endif # strict
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -xHost -qopt-zmm-usage=high
OPTFFLAGS=$(OPTFLAGS)
# OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-DNDEBUG -qopt-report=5 -traceback -diag-disable=10397
DBGFFLAGS=$(DBGFLAGS)
# DBGCFLAGS=$(DBGFLAGS) -w3 -diag-disable=1572,2547
else # DEBUG
OPTFLAGS=-O0 -xHost -qopt-zmm-usage=high
OPTFFLAGS=$(OPTFLAGS)
# OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -traceback -diag-disable=10397
ifneq ($(ARCH),Darwin) # Linux
DBGFLAGS += -debug parallel
endif # Linux
DBGFFLAGS=$(DBGFLAGS) -debug-parameters all -check all -warn all
# DBGCFLAGS=$(DBGFLAGS) -check=stack,uninit -w3 -diag-disable=1572,2547
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared -I../../../JACSD/vn
ifdef ANIMATE
LIBFLAGS += -DUSE_MKL -DMKL_ILP64 -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
endif # ANIMATE
LDFLAGS=-L../shared -ljk$(PROFILE)$(DEBUG) -L../../../JACSD -lvn$(PROFILE)$(DEBUG)
ifeq ($(ARCH),Darwin)
ifdef ANIMATE
LDFLAGS += -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core #-lmkl_intel_thread
endif # ANIMATE
else # Linux
LIBFLAGS += -D_GNU_SOURCE
ifdef ANIMATE
LDFLAGS += -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core #-lmkl_intel_thread
endif ANIMATE
endif # ?Darwin
LDFLAGS += -lpthread -lm -ldl
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
# CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUCFLAGS)
