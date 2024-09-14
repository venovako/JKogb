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
CPUFLAGS=-DUSE_INTEL -DUSE_X64 -qopenmp -xHost -qopt-multi-version-aggressive -vec-threshold0
FORFLAGS=$(CPUFLAGS) -standard-semantics -threads
FPUFLAGS=-fp-model $(FP) -fprotect-parens -fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt
ifeq ($(FP),strict)
FPUFLAGS += -fp-stack-check -assume ieee_fpe_flags
endif # ?strict
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -qopt-report=5 #-DUSE_FAST
DBGFLAGS=-DNDEBUG -diag-disable=10397
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -debug-parameters all -check all -warn all -diag-disable=10397
ifneq ($(ARCH),Darwin)
DBGFLAGS += -debug parallel
endif # Linux
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared
LDFLAGS=-rdynamic -static-libgcc -L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
ifdef ANIMATE
LDFLAGS += -L../../../libpvn/src -lpvn -ldl -lm
endif # ANIMATE
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
