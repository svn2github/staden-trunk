#
# Makefile for hetins
#

PROGS = $(O)/hetins

#SRCROOT=$(STADENROOT)/src
SRCROOT=..

include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

HETINSBIN=$(O)


INCLUDES_E += $(IOLIB_INC) $(HETINS_INC)


HLIBS=\
	$(MISC_LIB) \
	$(IOLIB_LIB)

OBJSD	= \
	$(HETINSBIN)/hetins.o


$(O)/hetins: $(OBJSD)
	$(CXXLD) $(CXXLDFLAGS) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJSD) $(HLIBS) $(LIBSC)


distsrc: distsrc_dirs
	-cp -R *.[ch] *.[ch]pp Makefile dependencies $(DIRNAME)

install:
	cp $(PROGS) $(INSTALLBIN)

DEPEND_OBJ = $(OBJS) $(OBJS) $(OBJSA)

include dependencies

