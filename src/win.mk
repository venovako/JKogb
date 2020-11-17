RM=del /F
AR=xilib.exe
ARFLAGS=-qnoipo -lib /NOLOGO /VERBOSE
FC=ifort.exe
CPUFLAGS=/DUSE_INTEL /DUSE_X64 /Qopenmp
FORFLAGS=/nologo /fpp $(CPUFLAGS) /standard-semantics
FPUFLAGS=/fp:source /Qfma /Qftz- /Qcomplex-limited-range- /Qfast-transcendentals- /Qprec-div /Qprec-sqrt /Qimf-precision:high
!IFDEF NDEBUG
OPTFLAGS=/O$(NDEBUG) /QxHost /Qopt-multi-version-aggressive
OPTFFLAGS=$(OPTFLAGS)
DBGFLAGS=/DNDEBUG /Qopt-report:5 /traceback
DBGFFLAGS=$(DBGFLAGS)
FPUFFLAGS=$(FPUFLAGS) /Qimf-use-svml:true
LIBFLAGS=-I. -I..\shared /libs:dll /threads
LDFLAGS=/link /RELEASE /LIBPATH:..\shared jk.lib
!ELSE # DEBUG
OPTFLAGS=/O$(DEBUG) /QxHost /Qopt-multi-version-aggressive
OPTFFLAGS=$(OPTFLAGS)
DBGFLAGS=/debug:full /debug:inline-debug-info /traceback
DBGFFLAGS=$(DBGFLAGS) /debug-parameters:all /check:all /warn:all
FPUFFLAGS=$(FPUFLAGS) /Qfp-stack-check /Qimf-arch-consistency:true
LIBFLAGS=-I. -I..\shared /libs:dll /threads /dbglibs
LDFLAGS=/link /DEBUG /LIBPATH:..\shared jk$(DEBUG).lib
!ENDIF # ?NDEBUG
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
