ifndef CPU
CPU=gnu
endif # !CPU
include ../$(CPU).mk
MKFS=GNUmakefile ../$(CPU).mk

.PHONY: all help clean

all: libjk$(PROFILE)$(DEBUG).a

help:
	@echo "gmake [CPU=x64|x200|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

libjk$(PROFILE)$(DEBUG).a: binio.o blas_utils.o atypes.o jstep.o params.o timer.o utils.o utils_.o $(MKFS)
	$(AR) $(ARFLAGS) $@ binio.o blas_utils.o atypes.o jstep.o params.o timer.o utils.o utils_.o

binio.o binio.mod: binio.F90 ../../../JACSD/vn/vn_binio_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c binio.F90

blas_utils.o blas_utils.mod: blas_utils.F90 params.mod ../../../JACSD/vn/vn_blas_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c blas_utils.F90

atypes.o atypes.mod: atypes.F90 timer.mod utils.mod ../../../JACSD/vn/vn_sort_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c atypes.F90

jstep.o jstep.mod: jstep.F90 atypes.mod $(MKFS)
	$(FC) $(FFLAGS) -c jstep.F90

params.o params.mod: params.F90 ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) $(FFLAGS) -c params.F90

timer.o timer.mod: timer.F90 params.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c timer.F90

utils.o utils.mod: utils.F90 utils_.hpp $(MKFS)
	$(FC) $(FFLAGS) -c utils.F90

utils_.o: utils_.cpp utils_.hpp $(MKFS)
	$(CXX) $(CXXFLAGS) -c utils_.cpp

clean:
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.dSYM