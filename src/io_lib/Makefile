#-----------------------------------------------------------------------------
# Makefile for Read library.
#

LIBS	= read
PROGS	= $(LIBS) progs
PROGLIBS= $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk
CFLAGS += $(SHLIB_CFLAGS) -Iinclude
#COPTDEBUG=-O -g3

%.o:	../%.c
	$(CC) $(CFLAGS) $(COBJFLAG)$@ -c $<
	$(MKDEFC) $(MKFLAGS) $@

include options.mk

#-----------------------------------------------------------------------------
# Mandatory object files
READOBJS=\
	read/$(O)/Read.o \
	read/$(O)/translate.o \
	read/$(O)/scf_extras.o

UTILSOBJS=\
	utils/$(O)/find.o \
	utils/$(O)/mach-io.o \
	utils/$(O)/traceType.o \
	utils/$(O)/read_alloc.o \
	utils/$(O)/compress.o \
	utils/$(O)/open_trace_file.o

# Deleted from Staden Package Makefile only due to duplicates in Misc
#	utils/$(O)/error.o 
#	utils/$(O)/files.o 
#	utils/$(O)/array.o 
#	utils/$(O)/xalloc.o

OBJS	= $(READOBJS) $(UTILSOBJS)

#-----------------------------------------------------------------------------
# Optional objects, depending on above IOLIB_* definitions
ifdef IOLIB_SCF
OBJS	+=\
	scf/$(O)/read_scf.o \
	scf/$(O)/write_scf.o \
	scf/$(O)/misc_scf.o
endif

ifdef IOLIB_EXP
OBJS	+=\
	exp_file/$(O)/expFileIO.o
endif

ifdef IOLIB_PLN
OBJS	+=\
	plain/$(O)/seqIOPlain.o
endif

ifdef IOLIB_ABI
OBJS	+=\
	abi/$(O)/fpoint.o \
	abi/$(O)/seqIOABI.o
endif

ifdef IOLIB_ALF
OBJS	+=\
	alf/$(O)/seqIOALF.o
endif

ifdef IOLIB_CTF
OBJS	+=\
	ctf/$(O)/ctfCompress.o \
	ctf/$(O)/seqIOCTF.o
endif

ifdef IOLIB_ZTR
OBJS	+=\
	ztr/$(O)/compression.o\
	ztr/$(O)/ztr_translate.o\
	ztr/$(O)/ztr.o
endif

#-----------------------------------------------------------------------------
# Build rules

RLIBS += $(MISC_LIB)

$(LIBS): $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)
	@

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX):	$(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(RLIBS)

# Rule used when $(DEF_FILE) defined - currently only for windows
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)


# Platform specific clean rules
ifeq ($(MACHINE),windows)
CLEANCOMMAND = -rm -f */$(O)/*.o */$(O)/*.def */$(O)/*.ilk */$(O)/*.pdb
else
CLEANCOMMAND = -rm -f $(OBJS); cd progs && $(MAKE) clean
endif


# To rebuild the dependencies, type "gmake depend"
DEPEND_OBJ = $(OBJS)
include dependencies

.PHONY:	progs
include/.links:
	cd include; \
	for i in abi alf ctf exp_file plain read scf utils ztr;	\
	do \
	    ln -s ../$$i/*.h .;\
	done
	touch include/.links

progs: include/.links
	cd $@ && $(MAKE)

#clean:
#	$(CLEANRULE)


spotless:	clean
	-rm -f $(PROGLIBS) $(L)/so_locations
	cd progs && $(MAKE) spotless

install:
	-mkdir -p $(INSTALLBIN)
	-mkdir -p $(INSTALLLIB)
	-mkdir -p $(INSTALLLIB)/$(O)
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)
	cd progs && $(MAKE) install

FORCE:
distsrc: FORCE
	echo DIRNAME=$(DIRNAME)
	-rmdir $(DIRNAME)/alpha-binaries
	-rmdir $(DIRNAME)/solaris-binaries
	-rmdir $(DIRNAME)/sgi-binaries
	-rmdir $(DIRNAME)/linux-binaries
	cpd() { \
	    mkdir $(DIRNAME)/$$1; \
	    mkdir $(DIRNAME)/$$1/alpha-binaries; \
	    mkdir $(DIRNAME)/$$1/solaris-binaries; \
	    mkdir $(DIRNAME)/$$1/sgi-binaries; \
	    mkdir $(DIRNAME)/$$1/linux-binaries; \
	    cp -R $$1/*.[ch] $(DIRNAME)/$$1; \
	}; \
	cpd abi; \
	cpd alf; \
	cpd ctf; \
	cpd exp_file; \
	cpd plain; \
	cpd progs; \
	cpd read; \
	cpd scf; \
	cpd utils; \
	cpd ztr
	-cp -R CHANGES COPYRIGHT Makefile README options.mk include man mk \
	    $(DIRNAME)
	-cp -R ztr/FORMAT $(DIRNAME)/ztr
	-cp -R progs/Makefile $(DIRNAME)/progs

