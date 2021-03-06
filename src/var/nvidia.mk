# BEWARE: the compiled executables might not work correctly!
SHELL=/bin/bash
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
RM=rm -rfv
AR=ar
ARFLAGS=rsv
CPUFLAGS=-DUSE_NVIDIA -DUSE_X64 -m64 -mp -Minfo
FORFLAGS=$(CPUFLAGS) -Mdclchk -Mstack_arrays
C11FLAGS=$(CPUFLAGS) -c11 -D_GNU_SOURCE
CC=nvc
FC=nvfortran
FPUFLAGS=-Kieee -Mfma -Mnodaz -Mnoflushz -Mnofpapprox -Mnofprelaxed
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) #-DUSE_FAST
DBGFLAGS=-DNDEBUG #-DUSE_TEST
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-g -Mbounds -Mchkstk
endif # ?NDEBUG
OPTFFLAGS=$(OPTFLAGS)
DBGFFLAGS=$(DBGFLAGS)
FPUFFLAGS=$(FPUFLAGS)
LIBFLAGS=-D_GNU_SOURCE -I. -I../shared
LDFLAGS=-L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(C11FLAGS) $(FPUFLAGS)
