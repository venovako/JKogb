!INCLUDE ..\win.mk

all: zjk.exe ztest.exe

help:
	@echo "nmake.exe [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

ztransf.obj ztransf.mod: ztransf.F90 ..\shared\ktransf.mod Makefile
	$(FC) $(FFLAGS) /c ztransf.F90

ztypes.obj ztypes.mod: ztypes.F90 ..\shared\atypes.mod Makefile
	$(FC) $(FFLAGS) /c ztypes.F90

zstep.obj zstep.mod: zstep.F90 ztransf.mod ztypes.mod ..\shared\jstep.mod Makefile
	$(FC) $(FFLAGS) /c zstep.F90

zjk.exe: zjk.obj zstep.obj ztransf.obj ztypes.obj Makefile
	$(FC) $(FFLAGS) zjk.obj zstep.obj ztransf.obj ztypes.obj /Fe$@ $(LDFLAGS)

zjk.obj: zjk.F90 bio.F90 ..\shared\readcl.F90 zstep.mod ..\shared\binio.mod Makefile
	$(FC) $(FFLAGS) /c zjk.F90

ztest.exe: ztest.obj ztransf.obj ..\shared\jk$(DEBUG).lib Makefile
	$(FC) $(FFLAGS) ztest.obj ztransf.obj /Fe$@ $(LDFLAGS)

ztest.obj: ztest.F90 ztransf.mod Makefile
	$(FC) $(FFLAGS) /c ztest.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.obj
	-$(RM) *.lib
	-$(RM) *.optrpt
	-$(RM) *.pdb
