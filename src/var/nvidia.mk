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
FORFLAGS=-DMIND=C_FMIN -DMAXD=C_FMAX $(CPUFLAGS) -Mdclchk -Mstack_arrays
FC=nvfortran
FPUFLAGS=-Kieee -Mfma -Mnodaz -Mnoflushz -Mnofpapprox -Mnofprelaxed
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) #-DUSE_FAST
DBGFLAGS=-DNDEBUG
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-g -Mbounds -Mchkstk
endif # ?NDEBUG
LIBFLAGS=-D_GNU_SOURCE -I. -I../shared
LDFLAGS=-L. -l$(TYPE)jk$(DEBUG) -L../shared -ljk$(DEBUG)
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
