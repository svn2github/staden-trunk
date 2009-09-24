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

distsrc: distsrc_dirs
	cp $(S)/*.[ch] $(S)/Makefile $(DIRNAME)

install:
	cp $(PROGLIBS) $(INSTALLLIB)

# DO NOT DELETE THIS LINE -- make depend depends on it.

freetree.o: $(SRCROOT)/Misc/misc.h
freetree.o: $(SRCROOT)/Misc/os.h
freetree.o: $(SRCROOT)/Misc/xalloc.h
freetree.o: $(SRCROOT)/Misc/xerror.h
freetree.o: $(PWD)/staden_config.h
freetree.o: freetree.h g-error.h
freetree.o: g-os.h
g-connect.o: $(SRCROOT)/Misc/array.h
g-connect.o: $(SRCROOT)/Misc/os.h
g-connect.o: $(SRCROOT)/Misc/xerror.h
g-connect.o: $(PWD)/staden_config.h
g-connect.o: freetree.h g-error.h
g-connect.o: g-connect.h g-struct.h
g-connect.o: g-defs.h g-db.h
g-connect.o: g-filedefs.h
g-connect.o: g-os.h
g-db.o: $(SRCROOT)/Misc/array.h
g-db.o: $(SRCROOT)/Misc/os.h
g-db.o: $(SRCROOT)/Misc/xerror.h
g-db.o: $(PWD)/staden_config.h
g-db.o: g-db.h g-struct.h freetree.h
g-db.o: g-error.h g-os.h
g-db.o: g-filedefs.h g-defs.h
g-db.o: g-files.h g-io.h g-request.h
g-error.o: $(SRCROOT)/Misc/xerror.h
g-error.o: g-error.h
g-error.o: g-misc.h
g-files.o: $(SRCROOT)/Misc/array.h
g-files.o: $(SRCROOT)/Misc/os.h
g-files.o: $(SRCROOT)/Misc/xalloc.h
g-files.o: $(SRCROOT)/Misc/xerror.h
g-files.o: $(PWD)/staden_config.h
g-files.o: freetree.h g-error.h
g-files.o: g-db.h g-defs.h
g-files.o: g-filedefs.h g-io.h
g-files.o: g-files.h g-struct.h
g-files.o: g-os.h
g-io.o: $(SRCROOT)/Misc/array.h
g-io.o: $(SRCROOT)/Misc/os.h
g-io.o: $(SRCROOT)/Misc/xerror.h
g-io.o: $(PWD)/staden_config.h
g-io.o: freetree.h g-error.h
g-io.o: g-filedefs.h
g-io.o: g-files.h g-struct.h
g-io.o: g-io.h
g-io.o: g-os.h
g-request.o: $(SRCROOT)/Misc/array.h
g-request.o: $(SRCROOT)/Misc/os.h
g-request.o: $(SRCROOT)/Misc/xalloc.h
g-request.o: $(SRCROOT)/Misc/xerror.h
g-request.o: $(PWD)/staden_config.h
g-request.o: freetree.h g-error.h
g-request.o: g-defs.h g-db.h
g-request.o: g-filedefs.h g-files.h
g-request.o: g-misc.h
g-request.o: g-os.h
g-request.o: g-request.h g-struct.h
g-struct.o: $(SRCROOT)/Misc/array.h
g-struct.o: $(SRCROOT)/Misc/os.h
g-struct.o: $(SRCROOT)/Misc/xalloc.h
g-struct.o: $(SRCROOT)/Misc/xerror.h
g-struct.o: $(PWD)/staden_config.h
g-struct.o: freetree.h g-error.h
g-struct.o: g-filedefs.h g-io.h
g-struct.o: g-files.h g-struct.h
g-struct.o: g-os.h
