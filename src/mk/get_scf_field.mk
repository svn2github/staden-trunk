#
# Makefile for get_scf_field
#

PROGS= 	$(O)/get_scf_field

SRCROOT=..

include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(SCF_INC)

SCFFIELDSRC=$(O)

all: $(O)/get_scf_field

OBJS	= $(SCFFIELDSRC)/get_scf_field.o

# 7/1/99 johnt - added MISC_DEP (for getopt with Visual C++)
DEPS=\
	$(IOLIB_DEP) \
	$(IOUTILS_DEP) \
	$(MISC_DEP)

# 7/1/99 johnt - added abstractions to support Visual C++
$(O)/get_scf_field:	$(OBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUBSYSTEMCONSOLE) $(OBJS) $(IOLIB_LIB) $(LIBSC) $(DEPS)

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
