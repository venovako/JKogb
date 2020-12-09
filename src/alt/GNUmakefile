ifndef COMPILER
COMPILER=gnu
endif # !COMPILER
ARCH=$(shell uname)
ifdef NDEBUG
DEBUG=
else # DEBUG
DEBUG=g
endif # ?NDEBUG
RM=rm -rfv
include $(COMPILER).mk
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean

all: libajk$(DEBUG)$(GNU).a common_test.exe j2apq_test.exe tweight_test.exe wpqb_test.exe

help:
	@echo "[g]make [COMPILER=gnu|clang|nvidia|intel] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

libajk$(DEBUG)$(GNU).a: common.o j2apq.o tweight.o wpqb.o $(MKFS)
	$(AR) $(ARFLAGS) $@ common.o j2apq.o tweight.o wpqb.o

common.o: common.c common.h $(MKFS)
	$(CC) $(CFLAGS) -c common.c -o $@

j2apq.o: j2apq.c j2apq.h wpqb.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c j2apq.c -o $@

tweight.o: tweight.c tweight.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c tweight.c -o $@

wpqb.o: wpqb.c wpqb.h common.h $(MKFS)
	$(CC) $(CFLAGS) -c wpqb.c -o $@

common_test.exe: common_test.c common.h libajk$(DEBUG)$(GNU).a $(MKFS)
	$(CC) $(CFLAGS) common_test.c -o $@ $(LDFLAGS)

j2apq_test.exe: j2apq_test.c j2apq.h wpqb.h common.h libajk$(DEBUG)$(GNU).a $(MKFS)
	$(CC) $(CFLAGS) j2apq_test.c -o $@ $(LDFLAGS)

tweight_test.exe: tweight_test.c tweight.h common.h libajk$(DEBUG)$(GNU).a $(MKFS)
	$(CC) $(CFLAGS) tweight_test.c -o $@ $(LDFLAGS)

wpqb_test.exe: wpqb_test.c wpqb.h common.h libajk$(DEBUG)$(GNU).a $(MKFS)
	$(CC) $(CFLAGS) wpqb_test.c -o $@ $(LDFLAGS)

clean:
ifdef NDEBUG
ifeq ($(COMPILER),intel)
	-$(RM) *.optrpt
endif # intel
else # DEBUG
ifeq ($(ARCH),Darwin)
	-$(RM) *.exe.dSYM
endif # Darwin
endif # ?NDEBUG
	-$(RM) libajk$(DEBUG)$(GNU).a
	-$(RM) *.exe
	-$(RM) *.o
