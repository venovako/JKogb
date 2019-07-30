ifndef CPU
CPU=gnu
endif # !CPU
include ../$(CPU).mk
MKFS=GNUmakefile ../$(CPU).mk

.PHONY: all help clean

all: zjk.exe

help:
	@echo "gmake [CPU=x64|x200|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

ztypes.o ztypes.mod: ztypes.F90 ../shared/atypes.mod ../shared/timer.mod ../shared/utils.mod ../../../JACSD/vn/vn_sort_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztypes.F90

zstep.o zstep.mod: zstep.F90 ztypes.mod $(MKFS)
	$(FC) $(FFLAGS) -c zstep.F90

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
