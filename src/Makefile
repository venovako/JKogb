MKFS=Makefile win.mk

all: SHARED REAL8 COMPLEX8

help:
	@echo nmake.exe [all|clean|help]

SHARED: $(MKFS)
	pushd shared && $(MAKE) /NOLOGO && popd

REAL8: SHARED $(MKFS)
	pushd D && $(MAKE) /NOLOGO && popd

COMPLEX8: SHARED $(MKFS)
	pushd Z && $(MAKE) /NOLOGO && popd

clean:
	pushd Z && $(MAKE) /NOLOGO clean && popd
	pushd D && $(MAKE) /NOLOGO clean && popd
	pushd shared && $(MAKE) /NOLOGO clean && popd
