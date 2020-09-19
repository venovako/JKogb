ifndef COMPILER
COMPILER=gnu
endif # !COMPILER
include ../$(COMPILER).mk
MKFS=GNUmakefile ../$(COMPILER).mk

.PHONY: all help clean

all: djk.exe dtest.exe

help:
	@echo "gmake [COMPILER=gnu|x64|x200|nvidia] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

dtransf.o dtransf.mod: dtransf.F90 ../shared/ktransf.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtransf.F90

dtypes.o dtypes.mod: dtypes.F90 ../shared/atypes.mod ../shared/timer.mod ../../../JACSD/vn/vn_sort_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtypes.F90
ifdef ANIMATE
dstep.o dstep.mod: dstep.F90 dtransf.mod dtypes.mod ../shared/utils.mod ../../../JACSD/vn/vn_mtxvis_f.mod $(MKFS)
	$(FC) -DANIMATE=$(ANIMATE) $(FFLAGS) -c dstep.F90
else # !ANIMATE
dstep.o dstep.mod: dstep.F90 dtransf.mod dtypes.mod ../shared/utils.mod $(MKFS)
	$(FC) $(FFLAGS) -c dstep.F90
endif # ?ANIMATE
djk.exe: djk.o dstep.o dtransf.o qtransf.o dtypes.o ../shared/libjk$(PROFILE)$(DEBUG).a ../../../JACSD/libvn$(PROFILE)$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) djk.o dstep.o dtransf.o qtransf.o dtypes.o -o$@ $(LDFLAGS)

djk.o: djk.F90 bio.F90 ../shared/readcl.F90 dstep.mod ../shared/binio.mod $(MKFS)
	$(FC) $(FFLAGS) -c djk.F90

qtransf.o: qtransf.c qtransf.h ../shared/common.h $(MKFS)
	$(CC) $(CFLAGS) -c qtransf.c

dtest.exe: dtest.o dtransf.o qtransf.o ../shared/libjk$(PROFILE)$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) dtest.o dtransf.o qtransf.o -o$@ $(LDFLAGS)

dtest.o: dtest.F90 dtransf.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtest.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.dSYM
