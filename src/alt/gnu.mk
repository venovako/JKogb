ifeq ($(ARCH),Darwin)
ifndef GNU
GNU=-10
endif # !GNU
endif # Darwin
AR=ar
ARFLAGS=rsv
CC=gcc$(GNU)
C18FLAGS=-fPIC -fexceptions -fno-omit-frame-pointer -pthread -std=gnu18
OPT=-march=native
DBG=-pedantic -Wall -Wextra
FPUFLAGS=-ffp-contract=fast
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) $(OPT) -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller
DBGFLAGS=-DNDEBUG -fopt-info-optimized-vec $(DBG)
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
endif # ?NDEBUG
ifeq ($(ARCH),Darwin)
OPTFLAGS += -Wa,-q
endif # Darwin
LIBFLAGS=-I.
ifneq ($(ARCH),Darwin)
LIBFLAGS += -D_GNU_SOURCE
endif # Linux
LDFLAGS=-rdynamic -L. -ljk$(DEBUG)$(GNU)
ifndef NDEBUG
LDFLAGS += -lubsan
endif # DEBUG
LDFLAGS += -lm
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUFLAGS)
