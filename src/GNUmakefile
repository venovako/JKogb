SHELL=/bin/bash
ifndef COMPILER
COMPILER=x64x
endif # !COMPILER
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean SHARED REAL8 COMPLEX8

all: SHARED REAL8 COMPLEX8

help:
	@echo "gmake [COMPILER=x64x|x200|x64|gnu] [MARCH=...] [NDEBUG=optimization_level] [ANIMATE=ppe] [all|clean|help]"

SHARED: $(MKFS)
ifdef NDEBUG
	pushd shared && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd shared && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?NDEBUG

REAL8: SHARED $(MKFS)
ifdef NDEBUG
ifdef ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
endif # ?ANIMATE
else # DEBUG
ifdef ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?ANIMATE
endif # ?NDEBUG

COMPLEX8: SHARED $(MKFS)
ifdef NDEBUG
ifdef ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
endif # ?ANIMATE
else # DEBUG
ifdef ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?ANIMATE
endif # ?NDEBUG

clean:
ifdef NDEBUG
	pushd Z && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
	pushd D && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
	pushd shared && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
else # DEBUG
	pushd Z && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd D && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd shared && $(MAKE) COMPILER=$(COMPILER) clean && popd
endif # ?NDEBUG
