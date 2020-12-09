AR=ar
ARFLAGS=rsv
CC=clang
C18FLAGS=-std=gnu18 -pthread
OPT=-march=native
DBG=-pedantic -Wall -Wextra
FPUFLAGS=-ffp-contract=fast
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) $(OPT)
DBGFLAGS=-DNDEBUG $(DBG)
ifneq ($(ARCH),Darwin)
FPUFLAGS += -fno-math-errno
endif # Linux
else # DEBUG
OPTFLAGS=-O$(DEBUG) $(OPT)
DBGFLAGS=-$(DEBUG) -ftrapv $(DBG)
endif # ?NDEBUG
LIBFLAGS=-I.
ifeq ($(ARCH),Darwin)
OPTFLAGS += -integrated-as
else # Linux
LIBFLAGS += -D_GNU_SOURCE
endif # ? Darwin
LDFLAGS=-L. -lajk$(DEBUG)
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUFLAGS)
