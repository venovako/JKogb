RM=rm -fv

HDF5=/usr/local
HDF5INC=-I$(HDF5)/include
HDF5LIB=-Wl,-rpath $(HDF5)/lib -L$(HDF5)/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5

MKLINC=-I${MKLROOT}/include
MKLLIB=-Wl,-rpath ${MKLROOT}/lib -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl
#MKLLIB=-Wl,-rpath ${MKLROOT}/lib -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

MYINC=$(HDF5INC) $(MKLINC)
MYLIB=$(HDF5LIB) $(MKLLIB)

F90=gfortran
ifdef NDEBUG
F90FLAGS=$(MYINC) -O$(NDEBUG) -cpp -DNDEBUG -march=native -Wa,-q -fopenmp -fstack-arrays -fgcse-las -fgcse-sm -fipa-pta -ftree-loop-distribution -ftree-loop-im -ftree-loop-ivcanon -fivopts -fvect-cost-model=unlimited -fvariable-expansion-in-unroller -Wall -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
else
F90FLAGS=$(MYINC) -Og -cpp -march=native -Wa,-q -fopenmp -fstack-arrays -fcheck=all -finit-local-zero -finit-real=snan -finit-derived -ffpe-trap=invalid,zero,overflow -Wall -Warray-temporaries -Wcharacter-truncation -Wimplicit-procedure -Wfunction-elimination -Wrealloc-lhs-all
endif
LDFLAGS=$(MYLIB)

all: jk.exe

h5data.o h5data.mod: h5data.f90 GNUmakefile
	$(F90) $(F90FLAGS) -c h5data.f90

jk1.o jk1.mod: jk1.f90 GNUmakefile
	$(F90) $(F90FLAGS) -c jk1.f90

jk.exe: jk.f90 jk1.mod h5data.mod GNUmakefile
	$(F90) $(F90FLAGS) jk.f90 jk1.o h5data.o -o $@ $(LDFLAGS)

clean:
	-$(RM) h5data.mod
	-$(RM) h5data.o
	-$(RM) jk1.mod
	-$(RM) jk1.o
	-$(RM) jk.exe
