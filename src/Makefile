MKFS=Makefile win.mk

all: SHARED REAL8 COMPLEX8

help:
	@echo "nmake.exe [NDEBUG=d|1|2|3|...] [all|clean|help]"

SHARED: $(MKFS)
	pushd shared && $(MAKE) /NOLOGO NDEBUG=$(NDEBUG) && popd

REAL8: SHARED $(MKFS)
	pushd D && $(MAKE) /NOLOGO NDEBUG=$(NDEBUG) && popd

COMPLEX8: SHARED $(MKFS)
	pushd Z && $(MAKE) /NOLOGO NDEBUG=$(NDEBUG) && popd

clean:
	pushd Z && $(MAKE) /NOLOGO NDEBUG=$(NDEBUG) clean && popd
	pushd D && $(MAKE) /NOLOGO NDEBUG=$(NDEBUG) clean && popd
	pushd shared && $(MAKE) /NOLOGO NDEBUG=$(NDEBUG) clean && popd
