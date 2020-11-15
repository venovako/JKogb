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
CPUFLAGS=-DUSE_INTEL -DUSE_X200 -fPIC -fexceptions -fno-omit-frame-pointer -qopenmp -qopenmp-threadprivate=compat -rdynamic
ifdef PROFILE
CPUFLAGS += -DVN_PROFILE=$(PROFILE) -fno-inline -finstrument-functions
endif # PROFILE
FORFLAGS=$(CPUFLAGS) -i8 -standard-semantics -threads #-cxxlib
# C18FLAGS=$(CPUFLAGS)
CPUFLAGS=-DUSE_INTEL -DUSE_X200 -DQX_WP=$(WP) -fPIC -fexceptions -fno-omit-frame-pointer -qopenmp -qopenmp-threadprivate=compat -rdynamic
ifdef PROFILE
CPUFLAGS += -DVN_PROFILE=$(PROFILE) -fno-inline -finstrument-functions
endif # PROFILE
CPUFLAGS += -DTSC_FREQ_HZ=$(shell if [ `if [ -r /etc/redhat-release ]; then grep -c 'release 7' /etc/redhat-release; else echo 0; fi` = 1 ]; then echo `dmesg | grep 'TSC clocksource calibration' | cut -d':' -f3 | cut -d' ' -f2 | sed 's/\.//g'`000; else echo 0; fi)ull
FORFLAGS=$(CPUFLAGS) -i8 -standard-semantics -threads
C18FLAGS=$(CPUFLAGS) -std=c18
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
DBGFLAGS=-$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug parallel -debug pubnames -traceback -diag-disable=10397
DBGFFLAGS=$(DBGFLAGS) -debug-parameters all -check all -warn all
# DBGCFLAGS=$(DBGFLAGS) -check=stack,uninit -w3 -diag-disable=1572,2547
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
