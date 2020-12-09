ifneq ($(ARCH),Linux)
$(error NVIDIA HPC SDK is supported on Linux only)
endif # !Linux
AR=ar
ARFLAGS=rsv
CC=nvc
C11FLAGS=-c11 -m64 -Minfo -pthread
FPUFLAGS=-Kieee -Mfma -Mnodaz -Mnoflushz -Mnofpapprox -Mnofprelaxed
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
DBGFLAGS=-DNDEBUG
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-g -Mchkstk
endif # ?NDEBUG
LIBFLAGS=-D_GNU_SOURCE -I.
LDFLAGS=-L. -lajk$(DEBUG)
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C11FLAGS) $(FPUFLAGS)
