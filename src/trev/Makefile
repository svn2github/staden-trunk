#PROGS = $(O)/trevsh

# Uncomment this for compiling of remote sources
#REMOTESRC=$(TREVSRC)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(IOLIB_INC) $(TK_INC)

#TREVBIN=$(O)

#$(O)/hex2c:	$(O)/hex2c.o
#	$(CLD) -o $@ $< $(LIBSC)

#$(O)/tcl2c:	$(O)/tcl2c.o
#	$(CLD) -o $@ $< $(LIBSC)

#$(O)/c2tcl:	$(O)/c2tcl.o
#	$(CLD) -o $@ $< $(LIBSC)

#%.c:	%.tcl $(O)/hex2c $(O)/tcl2c
#	tcl_cruncher-1.11/tcl_cruncher $< | $(O)/tcl2c | $(O)/hex2c > $@

#TREV_OBJS = \
#	$(TREVBIN)/tkAppInit.o

#DEPS=\
#	$(IOLIBBIN)/libread.a \
#	$(IOLIBBIN)/libio-utils.a \
#	$(MISCBIN)/libmisc.a

#TREV_LIBS=\
#	$(TKUTILS_LIB) \
#	$(IOLIB_LIB) \
#	$(MISC_LIB) \
#	$(TK_LIB)

#$(O)/trevsh:	$(TREV_OBJS)
#	$(CLD) -o $@ $(TREV_OBJS) $(TREV_LIBS) $(LIBSC) $(DEPS)

distsrc: distsrc_dirs
	-cp -R *.tcl Makefile trev trev.bat $(DIRNAME)

install:
#	cp $(PROGS) $(INSTALLBIN)
	cp trev $(INSTALLSCRIPT)
	-mkdir $(INSTALLLIB)/trev
	cp trev.tcl $(INSTALLLIB)/trev

DEPEND_OBJ = $(TREV_OBJS)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
