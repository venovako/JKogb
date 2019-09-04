ifndef CPU
CPU=gnu
endif # !CPU
include $(CPU).mk
MKFS=GNUmakefile $(CPU).mk

.PHONY: all help clean SHARED REAL8 COMPLEX8

all: SHARED REAL8 COMPLEX8

help:
	@echo "gmake [CPU=x64|x200|gnu] [NDEBUG=0|1|2|3|4|5] [all|clean|help]"

SHARED: $(MKFS)
ifdef NDEBUG
ifdef PROFILE
	pushd shared && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd shared && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd shared && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) && popd
else # !PROFILE
	pushd shared && $(MAKE) CPU=$(CPU) && popd
endif # ?PROFILE
endif # ?NDEBUG

REAL8: SHARED $(MKFS)
ifdef NDEBUG
ifdef PROFILE
ifdef ANIMATE
	pushd D && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
endif # ?ANIMATE
else # !PROFILE
ifdef ANIMATE
	pushd D && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) && popd
endif # ?ANIMATE
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
ifdef ANIMATE
	pushd D && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) && popd
endif # ?ANIMATE
else # !PROFILE
ifdef ANIMATE
	pushd D && $(MAKE) CPU=$(CPU) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd D && $(MAKE) CPU=$(CPU) && popd
endif # ?ANIMATE
endif # ?PROFILE
endif # ?NDEBUG

COMPLEX8: SHARED $(MKFS)
ifdef NDEBUG
ifdef PROFILE
ifdef ANIMATE
	pushd Z && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) && popd
endif # ?ANIMATE
else # !PROFILE
ifdef ANIMATE
	pushd Z && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) && popd
endif # ?ANIMATE
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
ifdef ANIMATE
	pushd Z && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) && popd
endif # ?ANIMATE
else # !PROFILE
ifdef ANIMATE
	pushd Z && $(MAKE) CPU=$(CPU) ANIMATE=$(ANIMATE) && popd
else # !ANIMATE
	pushd Z && $(MAKE) CPU=$(CPU) && popd
endif # ?ANIMATE
endif # ?PROFILE
endif # ?NDEBUG

clean:
ifdef NDEBUG
ifdef PROFILE
	pushd Z && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
	pushd D && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
	pushd shared && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) PROFILE=$(PROFILE) clean && popd
else # !PROFILE
	pushd Z && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) clean && popd
	pushd D && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) clean && popd
	pushd shared && $(MAKE) CPU=$(CPU) NDEBUG=$(NDEBUG) clean && popd
endif # ?PROFILE
else # DEBUG
ifdef PROFILE
	pushd Z && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) clean && popd
	pushd D && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) clean && popd
	pushd shared && $(MAKE) CPU=$(CPU) PROFILE=$(PROFILE) clean && popd
else # !PROFILE
	pushd Z && $(MAKE) CPU=$(CPU) clean && popd
	pushd D && $(MAKE) CPU=$(CPU) clean && popd
	pushd shared && $(MAKE) CPU=$(CPU) clean && popd
endif # ?PROFILE
endif # ?NDEBUG
