ifndef CPU
CPU=gnu
endif # !CPU
include $(CPU).mk
MKFS=GNUmakefile $(CPU).mk

.PHONY: all help clean

all: libjkz.a jkz.exe

help:
	@echo "gmake [CPU=x64|x200|pgi|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

libjkz.a: binio.o blas_utils.o dtypes.o jstep.o params.o timer.o utils.o utils_.o $(MKFS)
	$(AR) $(ARFLAGS) $@ binio.o blas_utils.o dtypes.o jstep.o params.o timer.o utils.o utils_.o

binio.o binio.mod: binio.F90 ../../../JACSD/vn/vn_binio_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c binio.F90

blas_utils.o blas_utils.mod: blas_utils.F90 params.mod ../../../JACSD/vn/vn_blas_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c blas_utils.F90

dtypes.o dtypes.mod: dtypes.F90 timer.mod utils.mod ../../../JACSD/vn/vn_sort_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c dtypes.F90

jstep.o jstep.mod: jstep.F90 dtypes.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c jstep.F90

params.o params.mod: params.F90 $(MKFS)
	$(FC) $(FFLAGS) -c params.F90

timer.o timer.mod: timer.F90 params.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -I../../../JACSD/vn -c timer.F90

utils.o utils.mod: utils.F90 utils_.hpp $(MKFS)
	$(FC) $(FFLAGS) -c utils.F90

utils_.o: utils_.cpp utils_.hpp $(MKFS)
	$(CXX) $(CXXFLAGS) -c utils_.cpp

jkz.exe: jkz.o libjkz.a ../../../JACSD/libvn.a $(MKFS)
	$(FC) $(FFLAGS) jkz.o -o$@ -L. -ljkz -L../../../JACSD -lvn $(LDFLAGS)

jkz.o: jkz.F90 bio.F90 readcl.F90 binio.mod blas_utils.mod jstep.mod $(MKFS)
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
