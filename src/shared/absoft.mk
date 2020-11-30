include ../af95.mk
MKFS=absoft.mk ../af95.mk

.PHONY: all help clean

all: libjk$(DEBUG).a

help:
	@echo "gmake [NDEBUG=1|2|3|fast|5] [all|clean|help]"

libjk$(DEBUG).a: atypes.o binio.o jstep.o ktransf.o timer.o utils.o params.o $(MKFS) #cutils.o
	$(AR) $(ARFLAGS) $@ atypes.o binio.o jstep.o ktransf.o timer.o utils.o params.o #cutils.o

atypes.o ATYPES.mod: atypes.F90 TIMER.mod UTILS.mod $(MKFS)
	$(FC) $(FFLAGS) -c atypes.F90

binio.o BINIO.mod: binio.F90 PARAMS.mod $(MKFS)
	$(FC) $(FFLAGS) -c binio.F90

jstep.o JSTEP.mod: jstep.F90 $(MKFS)
	$(FC) $(FFLAGS) -c jstep.F90

ktransf.o KTRANSF.mod: ktransf.F90 UTILS.mod $(MKFS)
	$(FC) $(FFLAGS) -c ktransf.F90

params.o PARAMS.mod: params.F90 $(MKFS)
	$(FC) $(FFLAGS) -c params.F90

timer.o TIMER.mod: timer.F90 PARAMS.mod $(MKFS)
	$(FC) $(FFLAGS) -c timer.F90

utils.o UTILS.mod: utils.F90 PARAMS.mod $(MKFS)
	$(FC) $(FFLAGS) -c utils.F90

# cutils.o: cutils.c cutils.h $(MKFS)
# 	$(FC) $(FFLAGS) -c cutils.c

clean:
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.dSYM
