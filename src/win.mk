RM=del /F
AR=xilib.exe
ARFLAGS=-qnoipo -lib /NOLOGO /VERBOSE
CC=icl.exe
FC=ifort.exe
CPUFLAGS=/DUSE_INTEL /DUSE_X64 /DUSE_WINDOWS /Qopenmp
C18FLAGS=/nologo /Qstd=c18 /Qlong-double
FORFLAGS=/nologo /fpp $(CPUFLAGS) /standard-semantics
OPTFLAGS=/QxHost /Qopt-multi-version-aggressive /Qopt-report:5
FPUFLAGS=/fp:precise /Qprotect-parens /Qfma /Qftz- /Qcomplex-limited-range- /Qfast-transcendentals- /Qprec-div /Qprec-sqrt
LIBFLAGS=-I. -I..\shared /libs:dll /threads
LDFLAGS=/link
!IFDEF NDEBUG
OPTFLAGS=/O$(NDEBUG) $(OPTFLAGS) #/DUSE_FAST
DBGFLAGS=/DNDEBUG #/DUSE_TEST
DBGFFLAGS=$(DBGFLAGS)
FPUFLAGS=$(FPUFLAGS) /Qimf-use-svml:true
C18FLAGS=$(C18FLAGS) /MD
LDFLAGS=$(LDFLAGS) /RELEASE /LIBPATH:. $(TYPE)jk.lib /LIBPATH:..\shared jk.lib
!ELSE # DEBUG
OPTFLAGS=/O$(DEBUG) $(OPTFLAGS)
DBGFLAGS=/debug:full /debug:inline-debug-info
DBGFFLAGS=$(DBGFLAGS) /debug-parameters:all /check:all /warn:all
FPUFLAGS=$(FPUFLAGS) /Qfp-stack-check
LIBFLAGS=$(LIBFLAGS) /dbglibs
C18FLAGS=$(C18FLAGS) /MDd
LDFLAGS=$(LDFLAGS) /DEBUG /LIBPATH:. $(TYPE)jk$(DEBUG).lib /LIBPATH:..\shared jk$(DEBUG).lib
!ENDIF # ?NDEBUG
OPTFFLAGS=$(OPTFLAGS)
FPUFFLAGS=$(FPUFLAGS)
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(C18FLAGS) $(FPUFLAGS)
