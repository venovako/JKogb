TYPE=z
include ../af95.mk
MKFS=absoft.mk ../af95.mk

.PHONY: all help clean

all: libzjk$(DEBUG).a zjk.exe ztest.exe

help:
	@echo "gmake [NDEBUG=1|2|3|fast|5] [all|clean|help]"

libzjk$(DEBUG).a: zstep.o ztransf.o ztypes.o $(MKFS)
	$(AR) $(ARFLAGS) $@ zstep.o ztransf.o ztypes.o

zstep.o ZSTEP.mod: zstep.F90 ZTRANSF.mod ZTYPES.mod ../shared/JSTEP.mod $(MKFS)
	$(FC) $(FFLAGS) -c zstep.F90

ztransf.o ZTRANSF.mod: ztransf.F90 ../shared/KTRANSF.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztransf.F90

ztypes.o ZTYPES.mod: ztypes.F90 ../shared/ATYPES.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztypes.F90

zjk.exe: zjk.o libzjk$(DEBUG).a ../shared/libjk$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) zjk.o -o$@ $(LDFLAGS)

zjk.o: zjk.F90 bio.F90 ../shared/readcl.F90 ZSTEP.mod ../shared/BINIO.mod $(MKFS)
	$(FC) $(FFLAGS) -c zjk.F90

ztest.exe: ztest.o libzjk$(DEBUG).a ../shared/libjk$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) ztest.o -o$@ $(LDFLAGS)

ztest.o: ztest.F90 ZTRANSF.mod $(MKFS)
	$(FC) $(FFLAGS) -c ztest.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.dSYM
