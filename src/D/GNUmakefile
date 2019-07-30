ifndef CPU
CPU=gnu
endif # !CPU
include ../$(CPU).mk
MKFS=GNUmakefile ../$(CPU).mk

.PHONY: all help clean

all: djk.exe

help:
	@echo "gmake [CPU=x64|x200|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

dtypes.o dtypes.mod: dtypes.F90 ../shared/atypes.mod ../shared/timer.mod ../shared/utils.mod ../../../JACSD/vn/vn_sort_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtypes.F90

dstep.o dstep.mod: dstep.F90 dtypes.mod $(MKFS)
	$(FC) $(FFLAGS) -c dstep.F90

djk.exe: djk.o dstep.o dtypes.o ../shared/libjk$(PROFILE)$(DEBUG).a ../../../JACSD/libvn$(PROFILE)$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) djk.o dstep.o dtypes.o -o$@ $(LDFLAGS)

djk.o: djk.F90 bio.F90 ../shared/readcl.F90 dstep.mod ../shared/binio.mod ../shared/blas_utils.mod $(MKFS)
	$(FC) $(FFLAGS) -c djk.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.dSYM
