#
# Makefile for Miscellaneous routines
#

LIBS 	= g
PROGS	= $(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

# 7/1/99 johnt - G_DLL defined under Visual C++ to export globals from DLL
CFLAGS += $(SHLIB_CFLAGS) $(G_DLL)

# Optimise IO
#COPTDEBUG=$(COPT)

GBIN=$(O)

#
# Objects
#

GSIO = \
	$(GBIN)/g-files.o \
	$(GBIN)/g-db.o \
	$(GBIN)/g-struct.o

GSLIB = \
	$(GBIN)/freetree.o

GSMISC = \
	$(GBIN)/g-error.o \
	$(GBIN)/g-io.o \
	$(GBIN)/g-debug.o

GSREQ = \
	$(GBIN)/g-connect.o \
	$(GBIN)/g-request.o

OBJS = $(GSIO) $(GSLIB) $(GSMISC) $(GSREQ)


DEPS = $(MISC_DEP)

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

# 1/7/99 johnt - added abstractions to support Visual C++ linker
$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(DEPS)

# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for WINNT
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)


DEPEND_OBJ = $(OBJS)

backup: FORCE
	tar cvf - *.[chf] README Makefile \
	    | gzip > backup/`date +%d_%m_%y`.tar.gz

FORCE::

include dependencies

install:
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)
