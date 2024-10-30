RM=del /F
AR=lib.exe
ARFLAGS=/NOLOGO /VERBOSE
FORT=ifx
FC=$(FORT).exe
CPUFLAGS=/DUSE_INTEL /DUSE_X64 /DUSE_WINDOWS /Qopenmp
FORFLAGS=/nologo /fpp $(CPUFLAGS) /standard-semantics
!IFNDEF MARCH
MARCH=Host
!ENDIF # !MARCH
OPTFLAGS=/Qx$(MARCH)
FPUFLAGS=/fp:precise /Qprotect-parens /Qfma /Qftz-
LIBFLAGS=-I. -I..\shared /libs:dll /threads
LDFLAGS=/link
!IFDEF NDEBUG
OPTFLAGS=/O$(NDEBUG) $(OPTFLAGS) /Qvec-threshold:0 #/DUSE_FAST
DBGFLAGS=/DNDEBUG
OPTFLAGS=$(OPTFLAGS) /Qopt-report:3
LDFLAGS=$(LDFLAGS) /RELEASE /LIBPATH:. $(TYPE)jk.lib /LIBPATH:..\shared jk.lib
!ELSE # DEBUG
OPTFLAGS=/O$(DEBUG) $(OPTFLAGS)
DBGFLAGS=/debug:full /debug:inline-debug-info /debug-parameters:all /check:all /warn:all
LIBFLAGS=$(LIBFLAGS) /dbglibs
LDFLAGS=$(LDFLAGS) /DEBUG /LIBPATH:. $(TYPE)jk$(DEBUG).lib /LIBPATH:..\shared jk$(DEBUG).lib
!ENDIF # ?NDEBUG
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
