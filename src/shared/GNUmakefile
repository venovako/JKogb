ifndef COMPILER
COMPILER=gnu
endif # !COMPILER
include ../$(COMPILER).mk
MKFS=GNUmakefile ../$(COMPILER).mk

.PHONY: all help clean

all: libjk$(DEBUG).a

help:
	@echo "gmake [COMPILER=gnu|x64|x64x|x200] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

libjk$(DEBUG).a: atypes.o binio.o jstep.o ktransf.o timer.o utils.o params.o $(MKFS) #cutils.o
	$(AR) $(ARFLAGS) $@ atypes.o binio.o jstep.o ktransf.o timer.o utils.o params.o #cutils.o

atypes.o atypes.mod: atypes.F90 timer.mod utils.mod $(MKFS)
	$(FC) $(FFLAGS) -c atypes.F90

binio.o binio.mod: binio.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) -c binio.F90

jstep.o jstep.mod: jstep.F90 $(MKFS)
	$(FC) $(FFLAGS) -c jstep.F90

ktransf.o ktransf.mod: ktransf.F90 utils.mod $(MKFS)
	$(FC) $(FFLAGS) -c ktransf.F90

params.o params.mod: params.F90 $(MKFS)
	$(FC) $(FFLAGS) -c params.F90

timer.o timer.mod: timer.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) -c timer.F90

utils.o utils.mod: utils.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) -c utils.F90

# cutils.o: cutils.c cutils.h $(MKFS)
# 	$(CC) $(CFLAGS) -c cutils.c

clean:
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.opt.yaml
	-$(RM) *.dSYM
