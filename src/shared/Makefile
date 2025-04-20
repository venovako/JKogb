!INCLUDE ..\win.mk
MKFS=Makefile ..\win.mk

all: jk.lib

help:
	@echo "nmake.exe [NDEBUG=d|1|2|3|...] [all|clean|help]"

jk.lib: atypes.obj binio.obj jstep.obj ktransf.obj timer.obj utils.obj params.obj $(MKFS)
	$(AR) $(ARFLAGS) /OUT:$@ atypes.obj binio.obj jstep.obj ktransf.obj timer.obj utils.obj params.obj

atypes.obj atypes.mod: atypes.F90 timer.mod utils.mod $(MKFS)
	$(FC) $(FFLAGS) /c atypes.F90

binio.obj binio.mod: binio.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) /c binio.F90

jstep.obj jstep.mod: jstep.F90 $(MKFS)
	$(FC) $(FFLAGS) /c jstep.F90

ktransf.obj ktransf.mod: ktransf.F90 utils.mod $(MKFS)
	$(FC) $(FFLAGS) /c ktransf.F90

params.obj params.mod: params.F90 $(MKFS)
	$(FC) $(FFLAGS) /c params.F90

timer.obj timer.mod: timer.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) /c timer.F90

utils.obj utils.mod: utils.F90 params.mod $(MKFS)
	$(FC) $(FFLAGS) /c utils.F90

clean:
	-$(RM) *.mod
	-$(RM) *.obj
	-$(RM) *.lib
	-$(RM) *.optrpt
	-$(RM) *.pdb
