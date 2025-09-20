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
ifndef MARCH
MARCH=Host
# COMMON-AVX512 for KNLs
endif # !MARCH
CPUFLAGS=-DUSE_INTEL -DUSE_X64 -fPIC -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -qopenmp -x$(MARCH) -qopt-multi-version-aggressive -qopt-zmm-usage=high -traceback -vec-threshold0
FORFLAGS=$(CPUFLAGS) -standard-semantics -threads
FPUFLAGS=-fp-model $(FP) -fprotect-parens -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt
ifeq ($(FP),strict)
FPUFLAGS += -fp-stack-check -assume ieee_fpe_flags
endif # ?strict
DBGFLAGS=-diag-disable=10448,10397
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -fno-math-errno -inline-level=2 -qopt-report=5 #-DUSE_FAST
DBGFLAGS += -DNDEBUG
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS += -$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -debug-parameters all -check all -warn all
ifneq ($(ARCH),Darwin)
DBGFLAGS += -debug parallel
endif # Linux
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared
LDFLAGS=-rdynamic
ifeq ($(ARCH),Darwin)
GCC=gcc-15
else # !Darwin
GCC=gcc
endif # !Darwin
LDFLAGS += -L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
ifdef ANIMATE
LDFLAGS += -L../../../libpvn/src
ifeq ($(ARCH),Darwin)
LDFLAGS += -Wl,-rpath,../../../libpvn/src
else # !Darwin
LDFLAGS += -Wl,-rpath=../../../libpvn/src
endif # ?Darwin
LDFLAGS += -lpvn -ldl -lm
endif # ANIMATE
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
