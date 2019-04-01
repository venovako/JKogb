ifeq ($(CPU),x64) # Xeon / Intel Fortran
include x64.mk
MKFS=GNUmakefile x64.mk
else ifeq ($(CPU),x200) # Knights Landing / Intel Fortran
include x200.mk
MKFS=GNUmakefile x200.mk
else # !CPU
CPU=x86_64 # GNU Fortran
include gnu.mk
MKFS=GNUmakefile gnu.mk
endif # ?CPU

.PHONY: all help clean

all: libjkz.a jkz.exe

help:
	@echo "make [CPU=x64|x100|x200] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

libjkz.a: binio.o blas_utils.o dtypes.o jstep.o params.o timer.o $(MKFS)
	$(AR) $(ARFLAGS) $@ binio.o blas_utils.o dtypes.o jstep.o params.o timer.o

binio.o binio.mod: binio.F90 ../../../JACSD/vn/vn_binio_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c binio.F90

blas_utils.o blas_utils.mod: blas_utils.F90 params.mod ../../../JACSD/vn/vn_blas_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c blas_utils.F90

dtypes.o dtypes.mod: dtypes.F90 params.mod ../../../JACSD/vn/vn_sort_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c dtypes.F90

jstep.o jstep.mod: jstep.F90 dtypes.mod params.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c jstep.F90

params.o params.mod: params.F90 $(MKFS)
	$(FC) $(FFLAGS) -c params.F90

timer.o timer.mod: timer.F90 params.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c timer.F90

jkz.exe: jkz.o libjkz.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) jkz.o -o$@ -L. -ljkz -L../../../JACSD -lvn $(LDFLAGS)

jkz.o: jkz.F90 bio.F90 readcl.F90 binio.mod blas_utils.mod dtypes.mod jstep.mod params.mod timer.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c jkz.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *__genmod.f90
	-$(RM) *__genmod.mod
	-$(RM) *.dSYM
