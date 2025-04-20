RM=del /F
AR=lib.exe
ARFLAGS=/NOLOGO /VERBOSE
CC=icx.exe
!IFNDEF MARCH
MARCH=Host
!ENDIF # !MARCH
!IFNDEF NDEBUG
NDEBUG=d
!ENDIF # NDEBUG
CFLAGS=/nologo /Qstd:c18 /Qlong-double /Qoption,cpp,--extended_float_types /DNDEBUG=$(NDEBUG) /DPVN_OS="\"Windows\"" /DPVN_ARCH="\"x64\"" /DPVN_COMPILER="\"icx\"" /DPVN_CR_MATH="\"bundled_noerrno\"" /DPVN_QUADMATH="\"Intel\"" /I. /EHsc /MT /O$(NDEBUG) /Qx$(MARCH) /Ob2 /fp:precise /Qfma /Qftz- /Qopt-report:3 /Qvec-threshold:0 /traceback /Wno-deprecated-declarations /Wno-overriding-option
LDFLAGS=/link /RELEASE /LIBPATH:. ajk.lib

all: ajk.lib common_test.exe j2apq_test.exe tweight_test.exe wpqb_test.exe

help:
	@echo "nmake.exe [MARCH=Host|...] [NDEBUG=d|1|2|3|...] [all|clean|help]"

ajk.lib: common.obj j2apq.obj tweight.obj wpqb.obj Makefile
	$(AR) $(ARFLAGS) /OUT:$@ common.obj j2apq.obj tweight.obj wpqb.obj

common.obj: common.c common.h Makefile
	$(CC) $(CFLAGS) /c common.c /Fo$@

j2apq.obj: j2apq.c j2apq.h wpqb.h common.h Makefile
	$(CC) $(CFLAGS) /c j2apq.c /Fo$@

tweight.obj: tweight.c tweight.h common.h Makefile
	$(CC) $(CFLAGS) /c tweight.c /Fo$@

wpqb.obj: wpqb.c wpqb.h common.h Makefile
	$(CC) $(CFLAGS) /c wpqb.c /Fo$@

common_test.exe: common_test.c common.h ajk.lib Makefile
	$(CC) $(CFLAGS) common_test.c /Fe$@ $(LDFLAGS)

j2apq_test.exe: j2apq_test.c j2apq.h wpqb.h common.h ajk.lib Makefile
	$(CC) $(CFLAGS) j2apq_test.c /Fe$@ $(LDFLAGS)

tweight_test.exe: tweight_test.c tweight.h common.h ajk.lib Makefile
	$(CC) $(CFLAGS) tweight_test.c /Fe$@ $(LDFLAGS)

wpqb_test.exe: wpqb_test.c wpqb.h common.h ajk.lib Makefile
	$(CC) $(CFLAGS) wpqb_test.c /Fe$@ $(LDFLAGS)

clean:
	-$(RM) *.exe
	-$(RM) *.obj
	-$(RM) *.lib
	-$(RM) *.optrpt
	-$(RM) *.pdb
