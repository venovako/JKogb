ifeq ($(ARCH),Darwin)
ifndef GNU
GNU=-10
endif # !GNU
endif # Darwin
AR=ar
ARFLAGS=rsv
CC=gcc$(GNU)
C18FLAGS=-std=gnu18 -pthread
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
DBGFLAGS=-$(DEBUG) -ftrapv $(DBG)
endif # ?NDEBUG
LIBFLAGS=-I.
ifeq ($(ARCH),Darwin)
OPTFLAGS += -Wa,-q
else # Linux
LIBFLAGS += -D_GNU_SOURCE
endif # ?Darwin
LDFLAGS=-L. -lajk$(DEBUG)$(GNU)
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUFLAGS)
