AR=xiar
ARFLAGS=-qnoipo -lib rsv
CC=icc 
C18FLAGS=-std=c18 -qopenmp
OPT=-xHost -qopt-zmm-usage=high
DBG=-w3 -diag-disable=1572,2547,10397
FPU=-fma -fprotect-parens -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) $(OPT)
DBGFLAGS=-DNDEBUG -qopt-report=5 $(DBG)
FPUFLAGS=-fp-model precise $(FPU) -fimf-use-svml=true
else # DEBUG
OPTFLAGS=-O0 $(OPT)
DBGFLAGS=-$(DEBUG) -debug emit_column -debug extended -debug inline-debug-info -debug pubnames -check=stack,uninit
ifneq ($(ARCH),Darwin)
DBGFLAGS += -debug parallel
endif # Linux
DBGFLAGS += $(DBG)
FPUFLAGS=-fp-model strict -fp-stack-check $(FPU)
endif # ?NDEBUG
LIBFLAGS=-I.
ifneq ($(ARCH),Darwin)
LIBFLAGS += -D_GNU_SOURCE
endif # Linux
LDFLAGS=-L. -lajk$(DEBUG)
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUCFLAGS)
