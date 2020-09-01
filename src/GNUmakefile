ifndef COMPILER
COMPILER=gnu
endif # !COMPILER
include $(COMPILER).mk
MKFS=GNUmakefile $(COMPILER).mk

.PHONY: all help clean SHARED REAL8 COMPLEX8

all: SHARED REAL8 COMPLEX8

help:
	@echo "gmake [COMPILER=gnu|x64|x200|nvidia] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

SHARED: $(MKFS)
ifdef NDEBUG
ifdef PROFILE
	pushd shared && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd shared && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd shared && $(MAKE) COMPILER=$(COMPILER) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd shared && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?PROFILE
endif # ?NDEBUG

REAL8: SHARED $(MKFS)
ifdef NDEBUG
ifdef PROFILE
ifdef ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
endif # ?ANIMATE
else # !PROFILE
ifdef ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
endif # ?ANIMATE
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
ifdef ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) PROFILE=$(PROFILE) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) PROFILE=$(PROFILE) && popd
endif # ?ANIMATE
else # !PROFILE
ifdef ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?ANIMATE
endif # ?PROFILE
endif # ?NDEBUG

COMPLEX8: SHARED $(MKFS)
ifdef NDEBUG
ifdef PROFILE
ifdef ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
endif # ?ANIMATE
else # !PROFILE
ifdef ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) && popd
endif # ?ANIMATE
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
ifdef ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) PROFILE=$(PROFILE) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) PROFILE=$(PROFILE) && popd
endif # ?ANIMATE
else # !PROFILE
ifdef ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) && popd
endif # ?ANIMATE
endif # ?PROFILE
endif # ?NDEBUG

clean:
ifdef NDEBUG
ifdef PROFILE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
	pushd D && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
	pushd shared && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
else # !PROFILE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
	pushd D && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
	pushd shared && $(MAKE) COMPILER=$(COMPILER) NDEBUG=$(NDEBUG) clean && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) PROFILE=$(PROFILE) clean && popd
	pushd D && $(MAKE) COMPILER=$(COMPILER) PROFILE=$(PROFILE) clean && popd
	pushd shared && $(MAKE) COMPILER=$(COMPILER) PROFILE=$(PROFILE) clean && popd
else # !PROFILE
	pushd Z && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd D && $(MAKE) COMPILER=$(COMPILER) clean && popd
	pushd shared && $(MAKE) COMPILER=$(COMPILER) clean && popd
endif # ?PROFILE
endif # ?NDEBUG
