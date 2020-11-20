RM=del /F
AR=xilib.exe
ARFLAGS=-qnoipo -lib /NOLOGO /VERBOSE
CC=icl.exe
FC=ifort.exe
CPUFLAGS=/DUSE_INTEL /DUSE_X64 /DUSE_WINDOWS /Qopenmp
C18FLAGS=/nologo /Qstd=c18 /Qlong-double
FORFLAGS=/nologo /fpp $(CPUFLAGS) /standard-semantics
OPTFLAGS=/QxHost /Qopt-multi-version-aggressive
DBGFLAGS=/Qopt-report:5
FPUFLAGS=/fp:precise /Qfma /Qprotect-parens /Qftz- /Qcomplex-limited-range- /Qfast-transcendentals- /Qprec-div /Qprec-sqrt #/Qimf-use-svml:true
FPUFFLAGS=$(FPUFLAGS) /DHYPOT=HYPOTwX87 /DABSZ=ABSwX87
LIBFLAGS=-I. -I..\shared /libs:dll /threads
LDFLAGS=/link
!IFDEF NDEBUG
OPTFLAGS=/O$(NDEBUG) $(OPTFLAGS)
OPTFFLAGS=$(OPTFLAGS)
DBGFLAGS=$(DBGFLAGS) /DNDEBUG
DBGFFLAGS=$(DBGFLAGS)
FPUFFLAGS=$(FPUFFLAGS)
C18FLAGS=$(C18FLAGS) /MD
LDFLAGS=$(LDFLAGS) /RELEASE /LIBPATH:..\shared jk.lib
!ELSE # DEBUG
OPTFLAGS=/O$(DEBUG) $(OPTFLAGS)
OPTFFLAGS=$(OPTFLAGS)
DBGFLAGS=$(DBGFLAGS) /debug:full /debug:inline-debug-info
DBGFFLAGS=$(DBGFLAGS) /debug-parameters:all /check:all /warn:all
FPUFFLAGS=$(FPUFFLAGS) /Qfp-stack-check
LIBFLAGS=$(LIBFLAGS) /dbglibs
C18FLAGS=$(C18FLAGS) /MDd
LDFLAGS=$(LDFLAGS) /DEBUG /LIBPATH:..\shared jk$(DEBUG).lib
!ENDIF # ?NDEBUG
FFLAGS=$(OPTFFLAGS) $(DBGFFLAGS) $(LIBFLAGS) $(FORFLAGS) $(FPUFFLAGS)
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(C18FLAGS) $(FPUFLAGS)
