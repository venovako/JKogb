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
AR=ar
ARFLAGS=rsv
FC=ifx
ifndef MARCH
MARCH=Host
# common-avx512 for KNLs
endif # !MARCH
CPUFLAGS=-DUSE_INTEL -DUSE_X64 -fPIC -fexceptions -fasynchronous-unwind-tables -fno-omit-frame-pointer -mprefer-vector-width=512 -vec-threshold0 -qopenmp -x$(MARCH) -traceback
FORFLAGS=$(CPUFLAGS) -standard-semantics -threads
FPUFLAGS=-fp-model=$(FP) -fp-speculation=safe -fma -fprotect-parens -no-ftz -fimf-precision=high
ifeq ($(FP),strict)
FPUFLAGS += -assume ieee_fpe_flags
endif # ?strict
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -fno-math-errno -inline-level=2 -qopt-report=3 #-DUSE_FAST
DBGFLAGS=-DNDEBUG
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -debug parallel -debug-parameters all -check all -warn all
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared
LDFLAGS=-rdynamic -L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
ifdef ANIMATE
LIBPVN=$(realpath ../../../libpvn)
LDFLAGS += -L$(LIBPVN)/src -Wl,-rpath=$(LIBPVN)/src -lpvn -lquadmath -lgcc_s -ldl -lm
endif # ANIMATE
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
