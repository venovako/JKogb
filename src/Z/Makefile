TYPE=z
!INCLUDE ..\win.mk
MKFS=Makefile ..\win.mk

all: zjk$(DEBUG).lib zjk.exe ztest.exe

help:
	@echo "nmake.exe [NDEBUG=d|1|2|3|4|5] [all|clean|help]"

zjk$(DEBUG).lib: zstep.obj ztransf.obj ztypes.obj $(MKFS)
	$(AR) $(ARFLAGS) /OUT:$@ zstep.obj ztransf.obj ztypes.obj

zstep.obj zstep.mod: zstep.F90 ztransf.mod ztypes.mod ..\shared\jstep.mod $(MKFS)
	$(FC) $(FFLAGS) /c zstep.F90

ztransf.obj ztransf.mod: ztransf.F90 ..\shared\ktransf.mod $(MKFS)
	$(FC) $(FFLAGS) /c ztransf.F90

ztypes.obj ztypes.mod: ztypes.F90 ..\shared\atypes.mod $(MKFS)
	$(FC) $(FFLAGS) /c ztypes.F90

zjk.exe: zjk.obj zjk$(DEBUG).lib ..\shared\jk$(DEBUG).lib $(MKFS)
	$(FC) $(FFLAGS) zjk.obj /Fe$@ $(LDFLAGS)

zjk.obj: zjk.F90 bio.F90 ..\shared\readcl.F90 zstep.mod ..\shared\binio.mod $(MKFS)
	$(FC) $(FFLAGS) /c zjk.F90

ztest.exe: ztest.obj zjk$(DEBUG).lib ..\shared\jk$(DEBUG).lib $(MKFS)
	$(FC) $(FFLAGS) ztest.obj /Fe$@ $(LDFLAGS)

ztest.obj: ztest.F90 ztransf.mod $(MKFS)
	$(FC) $(FFLAGS) /c ztest.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.obj
	-$(RM) *.lib
	-$(RM) *.optrpt
	-$(RM) *.opt.yaml
	-$(RM) *.pdb
