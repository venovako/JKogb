!IFDEF NDEBUG
DEBUG=
!ELSE # DEBUG
DEBUG=d
!ENDIF # ?NDEBUG
RM=del /F
AR=xilib.exe
ARFLAGS=-qnoipo -lib /NOLOGO /VERBOSE
!IFNDEF ICC
ICC=icl
!ENDIF # !ICC
CC=$(ICC).exe
C18FLAGS=/nologo /Qstd=c18 /Qlong-double /Qopenmp /Qvla
OPT=/QxHost /Qprotect-parens
DBG=/Wall
FPU=/Qfma /Qftz-
!IF "$(ICC)"=="icl"
C18FLAGS=$(C18FLAGS) /Qrestrict /Quse-intel-optimized-headers
OPT=$(OPT) /Qopt-multi-version-aggressive
DBG=$(DBG) /Wcheck /Qdiag-disable:10441
FPU=$(FPU) /Qcomplex-limited-range- /Qfast-transcendentals- /Qprec-div /Qprec-sqrt
!ENDIF # icl
LNK=/LIBPATH:. ajk$(DEBUG).lib
!IFDEF NDEBUG
!IF "$(ICC)"=="icl"
DBG=$(DBG) /Qopt-report:5
!ELSE # icx
DBG=$(DBG) /Qopt-report:3
!ENDIF # ?ICC
OPTFLAGS=/O$(NDEBUG) $(OPT)
DBGFLAGS=/DNDEBUG $(DBG)
FPUFLAGS=/fp:precise $(FPU) /Qimf-use-svml:true
LDFLAGS=/link /RELEASE $(LNK)
!ELSE # DEBUG
!IF "$(ICC)"=="icl"
DBG=/check:stack,uninit $(DBG)
FPUFLAGS=/fp:strict /Qfp-stack-check $(FPU)
!ELSE # icx
FPUFLAGS=/fp:precise $(FPU)
!ENDIF # ?ICC
OPTFLAGS=/O$(DEBUG) $(OPT)
DBGFLAGS=/debug:full /debug:inline-debug-info /debug:expr-source-pos $(DBG)
LDFLAGS=/link /DEBUG $(LNK)
!ENDIF # ?NDEBUG
LIBFLAGS=/I. /MD$(DEBUG)
CFLAGS=$(OPTFLAGS) $(DBGFLAGS) $(LIBFLAGS) $(C18FLAGS) $(FPUFLAGS)

all: ajk$(DEBUG).lib common_test.exe j2apq_test.exe tweight_test.exe wpqb_test.exe

help:
	@echo "nmake.exe [ICC=icl|icx] [NDEBUG=d|1|2|3|4|5] [all|clean|help]"

ajk$(DEBUG).lib: common.obj j2apq.obj tweight.obj wpqb.obj Makefile
	$(AR) $(ARFLAGS) /OUT:$@ common.obj j2apq.obj tweight.obj wpqb.obj

common.obj: common.c common.h Makefile
	$(CC) $(CFLAGS) /c common.c /Fo$@

j2apq.obj: j2apq.c j2apq.h wpqb.h common.h Makefile
	$(CC) $(CFLAGS) /c j2apq.c /Fo$@

tweight.obj: tweight.c tweight.h common.h Makefile
	$(CC) $(CFLAGS) /c tweight.c /Fo$@

wpqb.obj: wpqb.c wpqb.h common.h Makefile
	$(CC) $(CFLAGS) /c wpqb.c /Fo$@

common_test.exe: common_test.c common.h ajk$(DEBUG).lib Makefile
	$(CC) $(CFLAGS) common_test.c /Fe$@ $(LDFLAGS)

j2apq_test.exe: j2apq_test.c j2apq.h wpqb.h common.h ajk$(DEBUG).lib Makefile
	$(CC) $(CFLAGS) j2apq_test.c /Fe$@ $(LDFLAGS)

tweight_test.exe: tweight_test.c tweight.h common.h ajk$(DEBUG).lib Makefile
	$(CC) $(CFLAGS) tweight_test.c /Fe$@ $(LDFLAGS)

wpqb_test.exe: wpqb_test.c wpqb.h common.h ajk$(DEBUG).lib Makefile
	$(CC) $(CFLAGS) wpqb_test.c /Fe$@ $(LDFLAGS)

clean:
!IFDEF NDEBUG
	-$(RM) *.optrpt
!ELSE # DEBUG
	-$(RM) *.pdb
!ENDIF # ?NDEBUG
	-$(RM) ajk$(DEBUG).lib
	-$(RM) *.exe
	-$(RM) *.obj
