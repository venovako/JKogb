TYPE=d
include ../af95.mk
MKFS=absoft.mk ../af95.mk

.PHONY: all help clean

all: libdjk$(DEBUG).a djk.exe dtest.exe

help:
	@echo "gmake [NDEBUG=1|2|3|fast|5] [all|clean|help]"

libdjk$(DEBUG).a: dstep.o dtransf.o dtypes.o $(MKFS)
	$(AR) $(ARFLAGS) $@ dstep.o dtransf.o dtypes.o

dstep.o DSTEP.mod: dstep.F90 DTRANSF.mod DTYPES.mod ../shared/JSTEP.mod $(MKFS)
	$(FC) $(FFLAGS) -c dstep.F90

dtransf.o DTRANSF.mod: dtransf.F90 ../shared/KTRANSF.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtransf.F90

dtypes.o DTYPES.mod: dtypes.F90 ../shared/ATYPES.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtypes.F90

djk.exe: djk.o libdjk$(DEBUG).a ../shared/libjk$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) djk.o -o$@ $(LDFLAGS)

djk.o: djk.F90 bio.F90 ../shared/readcl.F90 DSTEP.mod ../shared/BINIO.mod $(MKFS)
	$(FC) $(FFLAGS) -c djk.F90

dtest.exe: dtest.o libdjk$(DEBUG).a ../shared/libjk$(DEBUG).a $(MKFS)
	$(FC) $(FFLAGS) dtest.o -o$@ $(LDFLAGS)

dtest.o: dtest.F90 DTRANSF.mod $(MKFS)
	$(FC) $(FFLAGS) -c dtest.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.a
	-$(RM) *.optrpt
	-$(RM) *.dSYM
