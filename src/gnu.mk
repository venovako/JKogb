SHELL=/bin/bash
ARCH=$(shell uname)
RM=rm -rfv
AR=ar
ARFLAGS=rsv
include ../../../libpvn/src/pvn.mk
FC=$(PVN_FC)
# comment out USE_X64 if not on Intel 64
FFLAGS=$(PVN_FCFLAGS) $(PVN_CPPFLAGS) -DUSE_GNU -DUSE_X64 -DMIND=C_FMIN -DMAXD=C_FMAX -I. -I../shared #-DUSE_FAST
LDFLAGS=$(PVN_LDFLAGS) -L. -l$(TYPE)jk -L../shared -ljk $(PVN_LIBS)
