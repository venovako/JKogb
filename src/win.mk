RM=del /F
AR=xilib.exe
ARFLAGS=-qnoipo -lib /NOLOGO /VERBOSE
!IFNDEF FORT
FORT=ifort
!ENDIF # !FORT
FC=$(FORT).exe
CPUFLAGS=/DUSE_INTEL /DUSE_X64 /DUSE_WINDOWS /Qopenmp
FORFLAGS=/nologo /fpp $(CPUFLAGS) /standard-semantics
OPTFLAGS=/QxHost
FPUFLAGS=/fp:precise /Qprotect-parens /Qfma /Qftz-
!IF "$(FORT)"=="ifort"
OPTFLAGS=$(OPTFLAGS) /Qopt-multi-version-aggressive
FPUFLAGS=$(FPUFLAGS) /Qcomplex-limited-range- /Qfast-transcendentals- /Qprec-div /Qprec-sqrt
!ENDIF # ifort
LIBFLAGS=-I. -I..\shared /libs:dll /threads
LDFLAGS=/link
!IFDEF NDEBUG
OPTFLAGS=/O$(NDEBUG) $(OPTFLAGS) /Qvec-threshold:0 #/DUSE_FAST
DBGFLAGS=/DNDEBUG
!IF "$(FORT)"=="ifort"
OPTFLAGS=$(OPTFLAGS) /Qopt-report:5
!ELSE # ifx
OPTFLAGS=$(OPTFLAGS) /Qopt-report:3
!ENDIF # ?FORT
LDFLAGS=$(LDFLAGS) /RELEASE /LIBPATH:. $(TYPE)jk.lib /LIBPATH:..\shared jk.lib
!ELSE # DEBUG
OPTFLAGS=/O$(DEBUG) $(OPTFLAGS)
DBGFLAGS=/debug:full /debug:inline-debug-info /debug-parameters:all /check:all /warn:all
!IF "$(FORT)"=="ifort"
FPUFLAGS=$(FPUFLAGS) /Qfp-stack-check
!ENDIF # ifort
LIBFLAGS=$(LIBFLAGS) /dbglibs
LDFLAGS=$(LDFLAGS) /DEBUG /LIBPATH:. $(TYPE)jk$(DEBUG).lib /LIBPATH:..\shared jk$(DEBUG).lib
!ENDIF # ?NDEBUG
FFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFLAGS)
