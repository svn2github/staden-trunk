#
# Makefile for Miscellaneous routines
#

LIBS 	= g
PROGS	= $(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

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


G_LIBS = $(MISC_LIB)

#
# Main dependency
#
$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(G_LIBS)

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
