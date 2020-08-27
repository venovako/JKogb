ifneq ($(ARCH),Linux)
$(error NVIDIA HPC SDK is supported on Linux only)
endif # !Linux
AR=ar
ARFLAGS=rsv
CC=nvc
C11FLAGS=-m64 -KPIC -Mframe -Meh_frame -Minfo -pthread -c11
FPUFLAGS=-Kieee -Mfma -Mnodaz -Mnoflushz -Mnofpapprox -Mnofprelaxed
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG)
DBGFLAGS=-DNDEBUG
else # DEBUG
OPTFLAGS=-O0
DBGFLAGS=-g -Mbounds -Mchkstk -traceback
endif # ?NDEBUG
LIBFLAGS=-D_GNU_SOURCE -I.
LDFLAGS=-L. -ljk$(PROFILE)$(DEBUG) -lm
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C11FLAGS) $(FPUFLAGS)
