RM=del /F

HDF5=C:\HDF5
HDF5INC=/I$(HDF5)\include
HDF5LIB=/LIBPATH:$(HDF5)\lib libhdf5_hl_fortran.lib libhdf5_hl_f90cstub.lib libhdf5_hl.lib libhdf5_fortran.lib libhdf5_f90cstub.lib libhdf5.lib

MYINC=$(HDF5INC)
MYLIB=$(HDF5LIB)

F90=ifort.exe
F90FLAGS=/standard-semantics $(MYINC) /O3 /fpp /DNDEBUG /QxHost /fp:source /Qfma /Qcomplex-limited-range- /Qfast-transcendentals- /Qftz- /Qprec-div /Qprec-sqrt /Qpc80 /nologo /traceback /Qvec-report /Qvc12 /Qmkl:sequential /libs:dll /threads
LDFLAGS=/link /RELEASE $(MYLIB)

all: jk.exe

h5data.obj h5data.mod: h5data.f90 Makefile
	$(F90) $(F90FLAGS) /c h5data.f90

jk1.obj jk1.mod: jk1.f90 Makefile
	$(F90) $(F90FLAGS) /c jk1.f90

jk.exe: jk.f90 jk1.mod h5data.mod Makefile
	$(F90) $(F90FLAGS) jk.f90 jk1.obj h5data.obj /Fe$@ $(LDFLAGS)

clean:
	-$(RM) h5data.mod
	-$(RM) h5data.obj
	-$(RM) jk1.mod
	-$(RM) jk1.obj
	-$(RM) jk.obj
	-$(RM) jk.exe
