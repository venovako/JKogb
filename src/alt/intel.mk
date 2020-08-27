AR=xiar
ARFLAGS=-qnoipo -lib rsv
CC=icc 
C18FLAGS=-fPIC -fexceptions -fno-omit-frame-pointer -qopenmp -std=c18
OPT=-xHost -qopt-zmm-usage=high
DBG=-traceback -w3 -diag-disable=1572,2547,10397
FPU=-fma -no-ftz -no-complex-limited-range -no-fast-transcendentals -prec-div -prec-sqrt -fimf-precision=high -fimf-use-svml=true
ifdef NDEBUG
OPTFLAGS=-O$(NDEBUG) $(OPT)
DBGFLAGS=-DNDEBUG -qopt-report=5 $(DBG)
FPUFLAGS=-fp-model source $(FPU) -fimf-use-svml=true
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
LDFLAGS=-rdynamic -L. -ljk$(PROFILE)$(DEBUG) -lm
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUCFLAGS)
