ifndef COMPILER
COMPILER=gnu
endif # !COMPILER
include ../$(COMPILER).mk
MKFS=GNUmakefile ../$(COMPILER).mk

.PHONY: all help clean

all: zjk.exe ztest.exe

help:
	@echo "gmake [COMPILER=gnu|x64|x200|nvidia] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

ztransf.o ztransf.mod: ztransf.F90 ../shared/ktransf.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztransf.F90

ztypes.o ztypes.mod: ztypes.F90 ../shared/atypes.mod ../shared/timer.mod ../../../JACSD/vn/vn_sort_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztypes.F90
ifdef ANIMATE
zstep.o zstep.mod: zstep.F90 ztransf.mod ztypes.mod ../shared/utils.mod ../../../JACSD/vn/vn_cmplxvis_f.mod $(MKFS)
	$(FC) -DANIMATE=$(ANIMATE) $(FFLAGS) -c zstep.F90
else # !ANIMATE
zstep.o zstep.mod: zstep.F90 ztransf.mod ztypes.mod ../shared/utils.mod $(MKFS)
	$(FC) $(FFLAGS) -c zstep.F90
endif # ?ANIMATE
zjk.exe: zjk.o zstep.o ztransf.o xtransf.o ztypes.o ../shared/libjk$(PROFILE)$(DEBUG).a ../../../JACSD/libvn$(PROFILE)$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) zjk.o zstep.o ztransf.o xtransf.o ztypes.o -o$@ $(LDFLAGS)

zjk.o: zjk.F90 bio.F90 ../shared/readcl.F90 zstep.mod ../shared/binio.mod $(MKFS)
	$(FC) $(FFLAGS) -c zjk.F90

xtransf.o: xtransf.c xtransf.h ../shared/common.h $(MKFS)
	$(CC) $(CFLAGS) -c xtransf.c

ztest.exe: ztest.o ztransf.o xtransf.o ../shared/libjk$(PROFILE)$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) ztest.o ztransf.o xtransf.o -o$@ $(LDFLAGS)

ztest.o: ztest.F90 ztransf.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztest.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.dSYM
