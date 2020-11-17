!INCLUDE ..\win.mk

all: jk$(DEBUG).lib

help:
	@echo "nmake.exe [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

jk$(DEBUG).lib: atypes.obj binio.obj jstep.obj ktransf.obj params.obj timer.obj utils.obj Makefile
	$(AR) $(ARFLAGS) /OUT:$@  atypes.obj binio.obj jstep.obj ktransf.obj params.obj timer.obj utils.obj

atypes.obj atypes.mod: atypes.F90 timer.mod utils.mod Makefile
	$(FC) $(FFLAGS) /c atypes.F90

binio.obj binio.mod: binio.F90 params.mod Makefile
	$(FC) $(FFLAGS) /c binio.F90

jstep.obj jstep.mod: jstep.F90 Makefile
	$(FC) $(FFLAGS) /c jstep.F90

ktransf.obj ktransf.mod: ktransf.F90 utils.mod Makefile
	$(FC) $(FFLAGS) /c ktransf.F90

params.obj params.mod: params.F90 Makefile
	$(FC) $(FFLAGS) /c params.F90

timer.obj timer.mod: timer.F90 Makefile
	$(FC) $(FFLAGS) /c timer.F90

utils.obj utils.mod: utils.F90 params.mod Makefile
	$(FC) $(FFLAGS) /c utils.F90

clean:
	-$(RM) *.mod
	-$(RM) *.obj
	-$(RM) *.lib
	-$(RM) *.optrpt
	-$(RM) *.pdb
