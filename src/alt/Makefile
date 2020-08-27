!IFDEF NDEBUG
DEBUG=
!ELSE # DEBUG
DEBUG=d
!ENDIF # ?NDEBUG
RM=del /F
AR=xilib.exe
ARFLAGS=-qnoipo -lib /NOLOGO /VERBOSE
CC=icl.exe
C11FLAGS=/nologo /Qstd=c11 /Qlong-double /Qopenmp
!IFDEF NDEBUG
OPTFLAGS=/O$(NDEBUG) /QxHost /Qopt-multi-version-aggressive
DBGFLAGS=/DNDEBUG /Qopt-report:5 /traceback
FPUFLAGS=/fp:source /Qimf-precision:high /Qfma /Qftz- /Qcomplex-limited-range- /Qfast-transcendentals- /Qprec-div /Qprec-sqrt
LIBFLAGS=/I. /MD
LDFLAGS=/link /RELEASE /LIBPATH:. jk$(DEBUG).lib
!ELSE # DEBUG
OPTFLAGS=/O$(DEBUG) /QxHost /Qopt-multi-version-aggressive
DBGFLAGS=/debug:full /debug:inline-debug-info /check:stack,uninit /traceback
FPUFLAGS=/fp:strict /Qfp-stack-check /Qimf-precision:high /Qfma /Qftz- /Qcomplex-limited-range- /Qfast-transcendentals- /Qprec-div /Qprec-sqrt
LIBFLAGS=/I. /MDd
LDFLAGS=/link /DEBUG /LIBPATH:. jk$(DEBUG).lib
!ENDIF # ?NDEBUG
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C11FLAGS) $(FPUFLAGS)

all: jk$(DEBUG).lib common_test.exe j2apq_test.exe tweight_test.exe wpqb_test.exe

help:
	@echo "nmake.exe [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

jk$(DEBUG).lib: common.obj j2apq.obj tweight.obj wpqb.obj Makefile
	$(AR) $(ARFLAGS) /OUT:$@ common.obj j2apq.obj tweight.obj wpqb.obj

common.obj: common.c common.h Makefile
	$(CC) $(CFLAGS) /c common.c /Fo$@

j2apq.obj: j2apq.c j2apq.h wpqb.h common.h Makefile
	$(CC) $(CFLAGS) /c j2apq.c /Fo$@

tweight.obj: tweight.c tweight.h common.h Makefile
	$(CC) $(CFLAGS) /c tweight.c /Fo$@

wpqb.obj: wpqb.c wpqb.h common.h Makefile
	$(CC) $(CFLAGS) /c wpqb.c /Fo$@

common_test.exe: common_test.c common.h jk$(DEBUG).lib Makefile
	$(CC) $(CFLAGS) common_test.c /Fe$@ $(LDFLAGS)

j2apq_test.exe: j2apq_test.c j2apq.h wpqb.h common.h jk$(DEBUG).lib Makefile
	$(CC) $(CFLAGS) j2apq_test.c /Fe$@ $(LDFLAGS)

tweight_test.exe: tweight_test.c tweight.h common.h jk$(DEBUG).lib Makefile
	$(CC) $(CFLAGS) tweight_test.c /Fe$@ $(LDFLAGS)

wpqb_test.exe: wpqb_test.c wpqb.h common.h jk$(DEBUG).lib Makefile
	$(CC) $(CFLAGS) wpqb_test.c /Fe$@ $(LDFLAGS)

clean:
!IFDEF NDEBUG
	-$(RM) wpqb_test.optrpt
	-$(RM) tweight_test.optrpt
	-$(RM) j2apq_test.optrpt
	-$(RM) common_test.optrpt
	-$(RM) wpqb.optrpt
	-$(RM) tweight.optrpt
	-$(RM) j2apq.optrpt
	-$(RM) common.optrpt
!ELSE # DEBUG
	-$(RM) *.pdb
!ENDIF # ?NDEBUG
	-$(RM) wpqb_test.exe
	-$(RM) wpqb_test.obj
	-$(RM) tweight_test.exe
	-$(RM) tweight_test.obj
	-$(RM) j2apq_test.exe
	-$(RM) j2apq_test.obj
	-$(RM) common_test.exe
	-$(RM) common_test.obj
	-$(RM) jk$(DEBUG).lib
	-$(RM) wpqb.obj
	-$(RM) tweight.obj
	-$(RM) j2apq.obj
	-$(RM) common.obj
