TYPE=z
ifndef COMPILER
COMPILER=gnu
endif # !COMPILER
include ../$(COMPILER).mk
MKFS=GNUmakefile ../$(COMPILER).mk

.PHONY: all help clean

all: libzjk$(DEBUG).a zjk.exe ztest.exe

help:
	@echo "gmake [COMPILER=gnu|x64|x64x|x200] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

libzjk$(DEBUG).a: zstep.o ztransf.o ztypes.o $(MKFS)
	$(AR) $(ARFLAGS) $@ zstep.o ztransf.o ztypes.o

ifdef ANIMATE
zstep.o zstep.mod: zstep.F90 ztransf.mod ztypes.mod ../shared/jstep.mod $(MKFS)
	$(FC) -DANIMATE=$(ANIMATE) $(FFLAGS) -c zstep.F90
else # !ANIMATE
zstep.o zstep.mod: zstep.F90 ztransf.mod ztypes.mod ../shared/jstep.mod $(MKFS)
	$(FC) $(FFLAGS) -c zstep.F90
endif # ?ANIMATE

ztransf.o ztransf.mod: ztransf.F90 ../shared/ktransf.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztransf.F90

ztypes.o ztypes.mod: ztypes.F90 ../shared/atypes.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztypes.F90

zjk.exe: zjk.o libzjk$(DEBUG).a ../shared/libjk$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) zjk.o -o $@ $(LDFLAGS)

zjk.o: zjk.F90 bio.F90 ../shared/readcl.F90 zstep.mod ../shared/binio.mod $(MKFS)
	$(FC) $(FFLAGS) -c zjk.F90

ztest.exe: ztest.o libzjk$(DEBUG).a ../shared/libjk$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) ztest.o -o $@ $(LDFLAGS)

ztest.o: ztest.F90 ztransf.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztest.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.opt.yaml
	-$(RM) *.dSYM
