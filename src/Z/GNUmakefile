ifndef CPU
CPU=gnu
endif # !CPU
include ../$(CPU).mk
MKFS=GNUmakefile ../$(CPU).mk

.PHONY: all help clean

all: zjk.exe

help:
	@echo "gmake [CPU=x64|x200|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

ztransf.o ztransf.mod: ztransf.F90 ../shared/params.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztransf.F90

ztypes.o ztypes.mod: ztypes.F90 ../shared/atypes.mod ../shared/timer.mod ../shared/utils.mod ../../../JACSD/vn/vn_sort_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztypes.F90
ifdef ANIMATE
zstep.o zstep.mod: zstep.F90 ztransf.mod ztypes.mod ../../../JACSD/vn/vn_cmplxvis_f.mod $(MKFS)
	$(FC) -DANIMATE $(FFLAGS) -c zstep.F90
else # !ANIMATE
zstep.o zstep.mod: zstep.F90 ztransf.mod ztypes.mod $(MKFS)
	$(FC) $(FFLAGS) -c zstep.F90
endif # ?ANIMATE
zjk.exe: zjk.o zstep.o ztypes.o ../shared/libjk$(PROFILE)$(DEBUG).a ../../../JACSD/libvn$(PROFILE)$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) zjk.o zstep.o ztypes.o -o$@ $(LDFLAGS)

zjk.o: zjk.F90 bio.F90 ../shared/readcl.F90 zstep.mod ../shared/binio.mod ../shared/blas_utils.mod $(MKFS)
	$(FC) $(FFLAGS) -c zjk.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.dSYM
