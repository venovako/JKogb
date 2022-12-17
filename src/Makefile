MKFS=Makefile win.mk

all: SHARED REAL8 COMPLEX8

help:
	@echo "nmake.exe [FORT=ifort|ifx] [NDEBUG=d|1|2|3|4|5] [all|clean|help]"

SHARED: $(MKFS)
!IFDEF NDEBUG
	pushd shared && $(MAKE) NDEBUG=$(NDEBUG) && popd
!ELSE # DEBUG
	pushd shared && $(MAKE) DEBUG=d && popd
!ENDIF # ?NDEBUG

REAL8: SHARED $(MKFS)
!IFDEF NDEBUG
	pushd D && $(MAKE) NDEBUG=$(NDEBUG) && popd
!ELSE # DEBUG
	pushd D && $(MAKE) DEBUG=d && popd
!ENDIF # ?NDEBUG

COMPLEX8: SHARED $(MKFS)
!IFDEF NDEBUG
	pushd Z && $(MAKE) NDEBUG=$(NDEBUG) && popd
!ELSE # DEBUG
	pushd Z && $(MAKE) DEBUG=d && popd
!ENDIF # ?NDEBUG

clean:
!IFDEF NDEBUG
	pushd Z && $(MAKE) NDEBUG=$(NDEBUG) clean && popd
	pushd D && $(MAKE) NDEBUG=$(NDEBUG) clean && popd
	pushd shared && $(MAKE) NDEBUG=$(NDEBUG) clean && popd
!ELSE # DEBUG
	pushd Z && $(MAKE) DEBUG=d clean && popd
	pushd D && $(MAKE) DEBUG=d clean && popd
	pushd shared && $(MAKE) DEBUG=d clean && popd
!ENDIF # ?NDEBUG
