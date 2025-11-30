SHELL=/bin/bash
ARCH=$(shell uname)
RM=rm -rfv
AR=ar
ARFLAGS=rsv
include ../../../libpvn/src/pvn.mk
FC=$(PVN_FC)
FFLAGS=$(PVN_FCFLAGS) $(PVN_CPPFLAGS) -DUSE_INTEL -DUSE_X64 -I. -I../shared #-DUSE_FAST
LDFLAGS=$(PVN_LDFLAGS) -L. -l$(TYPE)jk -L../shared -ljk $(PVN_LIBS)
