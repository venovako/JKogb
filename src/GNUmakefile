SHELL=/bin/bash
ifndef COMPILER
COMPILER=x64
endif # !COMPILER
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean SHARED REAL8 COMPLEX8

all: SHARED REAL8 COMPLEX8

help:
	@echo "gmake [COMPILER=x64|gnu] [ANIMATE=ppe] [all|clean|help]"

SHARED: $(MKFS)
	pushd shared && $(MAKE) COMPILER=$(COMPILER) && popd

REAL8: SHARED $(MKFS)
ifdef ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?ANIMATE

COMPLEX8: SHARED $(MKFS)
ifdef ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?ANIMATE

clean:
	pushd Z && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd D && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd shared && $(MAKE) COMPILER=$(COMPILER) clean && popd
