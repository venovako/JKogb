SHELL=/bin/bash
ARCH=$(shell uname)
ifneq ($(ARCH),Linux)
$(error NVIDIA HPC SDK is supported on Linux only)
endif # !Linux
ifndef WP
WP=8
endif # !WP
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
RM=rm -rfv
AR=ar
ARFLAGS=rsv
CPUFLAGS=-DUSE_NVIDIA -DUSE_X64 -DOLD_OMP -m64 -mp -KPIC -Mframe -Meh_frame -Minfo
ifdef PROFILE
CPUFLAGS += -DVN_PROFILE=$(PROFILE) -Mnoinline -Minstrument=functions
endif # PROFILE
FORFLAGS=$(CPUFLAGS) -i8 -Mdclchk -Mlarge_arrays -Mrecursive -Mstack_arrays
# C11FLAGS=$(CPUFLAGS) -c11
# CC=nvc
FC=nvfortran
# CXX=nvc++
FPUFLAGS=-Kieee -Mfma -Mnodaz -Mnoflushz -Mnofpapprox -Mnofprelaxed
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
OPTFFLAGS=$(OPTFLAGS)
# OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-DNDEBUG
DBGFFLAGS=$(DBGFLAGS)
# DBGCFLAGS=$(DBGFLAGS)
FPUFFLAGS=$(FPUFLAGS)
# FPUCFLAGS=$(FPUFLAGS)
else # DEBUG
OPTFLAGS=-O0
OPTFFLAGS=$(OPTFLAGS)
# OPTCFLAGS=$(OPTFLAGS)
DBGFLAGS=-g -Mbounds -Mchkstk -traceback
DBGFFLAGS=$(DBGFLAGS)
# DBGCFLAGS=$(DBGFLAGS)
FPUFFLAGS=$(FPUFLAGS)
# FPUCFLAGS=$(FPUFLAGS)
endif # ?NDEBUG
LIBFLAGS=-D_GNU_SOURCE -I. -I../shared -I../../../JACSD/vn
LDFLAGS=-L../shared -ljk$(PROFILE)$(DEBUG) -L../../../JACSD -lvn$(PROFILE)$(DEBUG) -pgf90libs -lpthread -lm -ldl $(shell if [ -L /usr/lib64/libmemkind.so ]; then echo '-lmemkind'; fi)
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
# CFLAGS=$(OPTCFLAGS) $(DBGCFLAGS) $(LIBFLAGS) $(C11FLAGS) $(FPUCFLAGS)
