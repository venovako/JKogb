AR=xiar
ARFLAGS=-qnoipo -lib rsv
CC=icx 
C18FLAGS=-std=c18 -qopenmp
OPT=-xHost -qopt-zmm-usage=high
FPU=-fprotect-parens -no-ftz
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) $(OPT)
DBGFLAGS=-DNDEBUG -qopt-report=3
FPUFLAGS=-fp-model precise $(FPU) -fma
else # DEBUG
OPTFLAGS=-O0 $(OPT)
DBGFLAGS=-$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -check=stack,uninit
ifneq ($(ARCH),Darwin)
DBGFLAGS += -debug parallel
endif # Linux
FPUFLAGS=-fp-model strict $(FPU)
endif # ?NDEBUG
LIBFLAGS=-I.
ifneq ($(ARCH),Darwin)
LIBFLAGS += -D_GNU_SOURCE
endif # Linux
LDFLAGS=-L. -lajk$(DEBUG)
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUCFLAGS)
