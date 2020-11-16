!INCLUDE ..\win.mk

all: djk.exe dtest.exe

help:
	@echo "nmake.exe [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

dtransf.obj dtransf.mod: dtransf.F90 ..\shared\ktransf.mod Makefile
	$(FC) $(FFLAGS) /c dtransf.F90

dtypes.obj dtypes.mod: dtypes.F90 ..\shared\atypes.mod Makefile
	$(FC) $(FFLAGS) /c dtypes.F90

dstep.obj dstep.mod: dstep.F90 dtransf.mod dtypes.mod ..\shared\jstep.mod Makefile
	$(FC) $(FFLAGS) /c dstep.F90

djk.exe: djk.obj dstep.obj dtransf.obj dtypes.obj Makefile
	$(FC) $(FFLAGS) djk.obj dstep.obj dtransf.obj dtypes.obj /Fe$@ $(LDFLAGS)

djk.obj: djk.F90 bio.F90 ..\shared\readcl.F90 dstep.mod ..\shared\binio.mod Makefile
	$(FC) $(FFLAGS) /c djk.F90

dtest.exe: dtest.obj dtransf.obj ..\shared\jk$(DEBUG).lib Makefile
	$(FC) $(FFLAGS) dtest.obj dtransf.obj /Fe$@ $(LDFLAGS)

dtest.obj: dtest.F90 dtransf.mod Makefile
	$(FC) $(FFLAGS) /c dtest.F90

clean:
	-$(RM) *.exe
	-$(RM) *.mod
	-$(RM) *.obj
	-$(RM) *.lib
	-$(RM) *.optrpt
	-$(RM) *.pdb
