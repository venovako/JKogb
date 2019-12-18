ifndef CPU
CPU=gnu
endif # !CPU
include ../$(CPU).mk
MKFS=GNUmakefile ../$(CPU).mk

.PHONY: all help clean

all: libjk$(PROFILE)$(DEBUG).a

help:
	@echo "gmake [CPU=x64|x200|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

libjk$(PROFILE)$(DEBUG).a: atypes.o binio.o common.o djk2.o jstep.o params.o timer.o utils.o wpqb.o $(MKFS)
	$(AR) $(ARFLAGS) $@  atypes.o binio.o common.o djk2.o jstep.o params.o timer.o utils.o wpqb.o

atypes.o atypes.mod: atypes.F90 timer.mod ../../../JACSD/vn/vn_sort_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c atypes.F90

binio.o binio.mod: binio.F90 ../../../JACSD/vn/vn_binio_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c binio.F90

common.o: common.c common.h $(MKFS)
	$(CC) $(CFLAGS) -c common.c

djk2.o: djk2.c djk2.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c djk2.c

jstep.o jstep.mod: jstep.F90 $(MKFS)
	$(FC) $(FFLAGS) -c jstep.F90

params.o params.mod: params.F90 ../../../JACSD/vn/vn_types_f.mod ../../../JACSD/vn/vn_params.F90 $(MKFS)
	$(FC) $(FFLAGS) -c params.F90

timer.o timer.mod: timer.F90 params.mod ../../../JACSD/vn/vn_timer_f.mod $(MKFS)
	$(FC) $(FFLAGS) -c timer.F90

utils.o utils.mod: utils.F90 $(MKFS)
	$(FC) $(FFLAGS) -c utils.F90

wpqb.o: wpqb.c wpqb.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c wpqb.c

clean:
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.dSYM
