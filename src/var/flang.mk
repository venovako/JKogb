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
# comment out USE_X64 if not on Intel 64
ifndef MARCH
MARCH=native
endif # !MARCH
CPUFLAGS=-DUSE_FLANG -DUSE_X64 -fPIC -fimplicit-none -fno-omit-frame-pointer -march=$(MARCH) #-fopenmp
ifeq ($(ARCH),Darwin)
CPUFLAGS += -DUSE_MACOS -fintegrated-as
endif # Darwin
FORFLAGS=-cpp -DMIND=C_FMIN -DMAXD=C_FMAX $(CPUFLAGS)
FPUFLAGS=-ffp-contract=fast -fhonor-infinities -fhonor-nans
FC=flang-new
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) #-DUSE_FAST
DBGFLAGS=-DNDEBUG
else # DEBUG
OPTFLAGS=-O$(DEBUG)
DBGFLAGS=-$(DEBUG)
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared
LDFLAGS=-rdynamic -L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
