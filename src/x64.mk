SHELL=/bin/bash
ARCH=$(shell uname)
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
ifndef FP
FP=precise
endif # !FP
RM=rm -rfv
AR=xiar
ARFLAGS=-qnoipo -lib rsv
FC=ifort
CPUFLAGS=-DUSE_INTEL -DUSE_X64 -qopenmp
FORFLAGS=$(CPUFLAGS) -standard-semantics -threads
ifdef ANIMATE
FORFLAGS += -i8
endif # ANIMATE
FPUFLAGS=-fp-model $(FP) -fprotect-parens -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt
ifeq ($(FP),strict)
FPUFLAGS += -fp-stack-check -assume ieee_fpe_flags
#else # !strict
#FPUFLAGS += -fimf-use-svml=true
endif # ?strict
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -xHost -qopt-multi-version-aggressive -vec-threshold0 -qopt-report=5 #-DUSE_FAST
DBGFLAGS=-DNDEBUG -diag-disable=10397
else # DEBUG
OPTFLAGS=-O0 -xHost -qopt-multi-version-aggressive #-qopt-report=5
DBGFLAGS=-$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -debug-parameters all -check all -warn all -diag-disable=10397
ifneq ($(ARCH),Darwin)
DBGFLAGS += -debug parallel
endif # Linux
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared
ifneq ($(ARCH),Darwin)
LIBFLAGS += -static-libgcc
endif # Linux
ifdef ANIMATE
LIBFLAGS += -DUSE_MKL -DMKL_ILP64 -I../../../JACSD/vn -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include
endif # ANIMATE
LDFLAGS=-L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
ifeq ($(ARCH),Darwin)
ifdef ANIMATE
LDFLAGS += -L../../../JACSD -lvn$(DEBUG) -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
endif # ANIMATE
else # Linux
ifdef ANIMATE
LDFLAGS += -L../../../JACSD -lvn$(DEBUG) -L${MKLROOT}/lib/intel64 -Wl,-rpath=${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
endif # ANIMATE
endif # ?Darwin
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
