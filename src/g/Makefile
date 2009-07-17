#
# Makefile for Miscellaneous routines
#

LIBS 	= g
PROGS	= $(LIBS)
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

INCLUDES_S += $(MISC_INC)

# CACHE_FREETREE makes Gap4 append the freetree data structure to the end of
# the .aux file when a database is shutdown. Similarly it's loaded on database
# open.
# For now this is disabled for several reasons:
# 1. The format has changed (freelists, AVL balance factors) and the code has
#    not yet been updated.
# 2. The existing trees are unbalanced, so we want to force recomputing.
# 3. The speed of recomputing from scratch when using a balanced tree seems
#    to be only barely slower than loading a cached copy.

#CFLAGS += $(SHLIB_CFLAGS) $(G_DLL) -DCACHE_FREETREE
CFLAGS += $(SHLIB_CFLAGS) $(G_DLL) 

# Optimise IO
#COPTDEBUG=$(COPT)

#
# Objects
#

GSIO = \
	g-files.o \
	g-db.o \
	g-struct.o

GSLIB = \
	freetree.o

GSMISC = \
	g-error.o \
	g-io.o

GSREQ = \
	g-connect.o \
	g-request.o

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
	cp $(PROGLIBS) $(INSTALLLIB)
