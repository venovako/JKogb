AR=ar
ARFLAGS=rsv
CC=clang
C18FLAGS=-std=gnu18 -fPIC -fexceptions -fno-omit-frame-pointer -pthread
OPT=-march=native
DBG=-pedantic -Wall -Wextra
FPU=-ffp-contract=fast
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) $(OPT)
DBGFLAGS=-DNDEBUG $(DBG)
FPUFLAGS=$(FPU)
ifneq ($(ARCH),Darwin)
FPUFLAGS += -fno-math-errno
endif # Linux
else # DEBUG
OPTFLAGS=-O$(DEBUG) $(OPT)
DBGFLAGS=-$(DEBUG) -fsanitize=address -fsanitize=undefined
ifneq ($(ARCH),Darwin)
DBGFLAGS += -fsanitize=leak
endif # Linux
DBGFLAGS += $(DBG)
FPUFLAGS=$(FPU)
endif # ?NDEBUG
ifeq ($(ARCH),Darwin)
OPTFLAGS += -integrated-as
endif # Darwin
LIBFLAGS=-I.
ifneq ($(ARCH),Darwin)
LIBFLAGS += -D_GNU_SOURCE
endif # Linux
LDFLAGS=-rdynamic -L. -ljk$(DEBUG)
LDFLAGS += -lm
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUFLAGS)
