SHELL=/bin/bash
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
RM=rm -rfv
AR=ar
ARFLAGS=rsv
# comment out USE_X64 if not on Intel 64
CPUFLAGS=-DUSE_SUN -DUSE_X64 -m64 -xarch=native
FORFLAGS=$(CPUFLAGS) -u
FPUFLAGS=-fma=fused -ftrap=%none
FC=f95
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) -xopenmp=parallel #-DUSE_FAST
DBGFLAGS=-DNDEBUG
else # DEBUG
OPTFLAGS=-xopenmp=noopt
DBGFLAGS=-$(DEBUG) -C -xcheck=%all -xcommonchk=yes
endif # ?NDEBUG
LIBFLAGS=-I. -I../shared
LDFLAGS=-L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
