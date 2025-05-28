TYPE=d
ifndef COMPILER
COMPILER=x64x
endif # !COMPILER
include ../$(COMPILER).mk
MKFS=GNUmakefile ../$(COMPILER).mk

.PHONY: all help clean

all: libdjk$(DEBUG).a djk.exe dtest.exe

help:
	@echo "gmake [COMPILER=x64x|x200|x64|gnu] [NDEBUG=optimization_level] [all|clean|help]"

libdjk$(DEBUG).a: dstep.o dtransf.o dtypes.o $(MKFS)
	$(AR) $(ARFLAGS) $@ dstep.o dtransf.o dtypes.o

ifdef ANIMATE
dstep.o dstep.mod: dstep.F90 dtransf.mod dtypes.mod ../shared/jstep.mod $(MKFS)
	$(FC) -DANIMATE=$(ANIMATE) $(FFLAGS) -c dstep.F90
else # !ANIMATE
dstep.o dstep.mod: dstep.F90 dtransf.mod dtypes.mod ../shared/jstep.mod $(MKFS)
	$(FC) $(FFLAGS) -c dstep.F90
endif # ?ANIMATE

dtransf.o dtransf.mod: dtransf.F90 ../shared/ktransf.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtransf.F90

dtypes.o dtypes.mod: dtypes.F90 ../shared/atypes.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtypes.F90

djk.exe: djk.o libdjk$(DEBUG).a ../shared/libjk$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) djk.o -o $@ $(LDFLAGS)

djk.o: djk.F90 bio.F90 ../shared/readcl.F90 dstep.mod ../shared/binio.mod $(MKFS)
	$(FC) $(FFLAGS) -c djk.F90

dtest.exe: dtest.o libdjk$(DEBUG).a ../shared/libjk$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) dtest.o -o $@ $(LDFLAGS)

dtest.o: dtest.F90 dtransf.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtest.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.dSYM
