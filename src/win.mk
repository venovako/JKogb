RM=del /F
AR=lib.exe
ARFLAGS=/NOLOGO /VERBOSE
FC=ifx.exe
!IFNDEF MARCH
MARCH=Host
!ENDIF # !MARCH
!IFNDEF NDEBUG
NDEBUG=d
!ENDIF # !NDEBUG
FFLAGS=/nologo /fpp /Qopenmp /standard-semantics /traceback /DNDEBUG=$(NDEBUG) /DUSE_INTEL /DUSE_X64 /DUSE_WINDOWS /I. -I..\shared /MD /O$(NDEBUG) /Qx$(MARCH) /fp:precise /Qfma /Qftz- /Qprec-div /Qprotect-parens /Qopt-report:3 /Qvec-threshold:0
LDFLAGS=/link /RELEASE /LIBPATH:. $(TYPE)jk.lib /LIBPATH:..\shared jk.lib
