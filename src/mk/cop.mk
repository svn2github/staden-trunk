#
# Makefile for COP (Check Out Project)
#

PROGS = $(O)/cop $(O)/cop-bap

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

COPBIN=$(O)
INCLUDES_E += -I$(BAPSRC) -I$(STADENSRC) -I$(TEDSRC) -I$(CONVERTSRC) $(XAW_INC)

#
# Main and miscellaneous routines
#
COPOBJ = \
	$(COPBIN)/cop.o

COPBAPOBJ = \
	$(COPBIN)/cop-bap.o

#
# For reading in xdap database
#
XDAPIO = \
	$(COPBIN)/dapIO.o

XBAPIO = \
	$(COPBIN)/bapIO.o

#
# For reading in trace files
#
TRACES = \
	$(COPBIN)/getSeq.o \
	$(COPBIN)/seqIOABI.o \
	$(COPBIN)/seqIOALF.o \
	$(COPBIN)/seqIOSCF.o \
	$(COPBIN)/seqIOPlain.o \
	$(COPBIN)/opp.o \
	$(COPBIN)/seq.o \
	$(COPBIN)/fpoint.o\
	$(COPBIN)/mach-io.o

#
# For aligning sequences
#
ALIGN = \
	$(COPBIN)/llin.o

#
# The whole lot
#
OBJ = \
	$(COPOBJ) \
	$(XDAPIO) \
	$(TRACES) \
	$(ALIGN)

OBJBAP = \
	$(COPBAPOBJ) \
	$(XBAPIO) \
	$(TRACES) \
	$(ALIGN)


#
# Cop
#
$(O)/cop : $(OBJ)
	$(CLD) -o $@ $(OBJ) $(MISC_LIB) $(LIBSC)

$(O)/cop-bap : $(OBJBAP)
	$(CLD) -o $@ $(OBJBAP) $(MISC_LIB) $(LIBSC)

$(O)/cop-bap.o: $(COPSRC)/cop.c
	$(CC) -DBAP_VERSION $(CFLAGS) -o $@ -c $(COPSRC)/cop.c

#
# Special rule needed for seq*.o from the ted directory.
#
$(O)/%.o: $(TEDSRC)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

DEPEND_OBJ = $(OBJ) $(OBJBAP)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
