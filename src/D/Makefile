TYPE=d
!INCLUDE ..\win.mk
MKFS=Makefile ..\win.mk

all: djk$(DEBUG).lib djk.exe dtest.exe

help:
	@echo "nmake.exe [NDEBUG=optimization_level] [all|clean|help]"

djk$(DEBUG).lib: dstep.obj dtransf.obj dtypes.obj $(MKFS)
	$(AR) $(ARFLAGS) /OUT:$@ dstep.obj dtransf.obj dtypes.obj

dstep.obj dstep.mod: dstep.F90 dtransf.mod dtypes.mod ..\shared\jstep.mod $(MKFS)
	$(FC) $(FFLAGS) /c dstep.F90

dtransf.obj dtransf.mod: dtransf.F90 ..\shared\ktransf.mod $(MKFS)
	$(FC) $(FFLAGS) /c dtransf.F90

dtypes.obj dtypes.mod: dtypes.F90 ..\shared\atypes.mod $(MKFS)
	$(FC) $(FFLAGS) /c dtypes.F90

djk.exe: djk.obj djk$(DEBUG).lib ..\shared\jk$(DEBUG).lib $(MKFS)
	$(FC) $(FFLAGS) djk.obj /Fe$@ $(LDFLAGS)

djk.obj: djk.F90 bio.F90 ..\shared\readcl.F90 dstep.mod ..\shared\binio.mod $(MKFS)
	$(FC) $(FFLAGS) /c djk.F90

dtest.exe: dtest.obj djk$(DEBUG).lib ..\shared\jk$(DEBUG).lib $(MKFS)
	$(FC) $(FFLAGS) dtest.obj /Fe$@ $(LDFLAGS)

dtest.obj: dtest.F90 dtransf.mod $(MKFS)
	$(FC) $(FFLAGS) /c dtest.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.obj
	-$(RM) *.lib
	-$(RM) *.optrpt
	-$(RM) *.opt.yaml
	-$(RM) *.pdb
