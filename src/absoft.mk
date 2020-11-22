MKFS=absoft.mk af95.mk

.PHONY: all help clean SHARED REAL8 COMPLEX8

all: SHARED REAL8 COMPLEX8

help:
	@echo "gmake [NDEBUG=1|2|3|fast|5] [all|clean|help]"

SHARED: $(MKFS)
ifdef NDEBUG
	pushd shared && $(MAKE) -f absoft.mk NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd shared && $(MAKE) -f absoft.mk && popd
endif # ?NDEBUG

REAL8: SHARED $(MKFS)
ifdef NDEBUG
	pushd D && $(MAKE) -f absoft.mk NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd D && $(MAKE) -f absoft.mk && popd
endif # ?NDEBUG

COMPLEX8: SHARED $(MKFS)
ifdef NDEBUG
	pushd Z && $(MAKE) -f absoft.mk NDEBUG=$(NDEBUG) && popd
else # DEBUG
	pushd Z && $(MAKE) -f absoft.mk && popd
endif # ?NDEBUG

clean:
ifdef NDEBUG
	pushd Z && $(MAKE) -f absoft.mk NDEBUG=$(NDEBUG) clean && popd
	pushd D && $(MAKE) -f absoft.mk NDEBUG=$(NDEBUG) clean && popd
	pushd shared && $(MAKE) -f absoft.mk NDEBUG=$(NDEBUG) clean && popd
else # DEBUG
	pushd Z && $(MAKE) -f absoft.mk clean && popd
	pushd D && $(MAKE) -f absoft.mk clean && popd
	pushd shared && $(MAKE) -f absoft.mk clean && popd
endif # ?NDEBUG
