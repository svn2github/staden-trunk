#
# Makefile for getMCH
#

PROGS = $(O)/trace2seq

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

GETMCHBIN=$(O)

#
#Turning on the AUTO_CLIP define allows ted to automatically
#select a left (using the -enzyme option or the .enzyme in the
#Xted file) and right cutoff on your sequence (using 2 out of
#5 N's)
#
# Turning on the SAVE_EDITS define allows the user to maintain
# copies of their edits, i.e. when you open up ted on a file
# that you have already edited, the old edits show up in the
# ted edit window.  The user may also call up any of their
# past editing sessions.  
#
#Turning on the DEF_OUT  define makes
# trace_name.seq the default output file name
#
#DEFINES += -DAUTO_CLIP  -DDEF_OUT  -DSAVE_EDITS
INCLUDES_E += -I$(TEDSRC) $(XAW_INC)

#
# Object files
#
OBJS=\
	$(GETMCHBIN)/getMCH.o\
	$(GETMCHBIN)/seq.o\
	$(GETMCHBIN)/seqIOPlain.o\
	$(GETMCHBIN)/seqIOABI.o\
	$(GETMCHBIN)/seqIOEdit.o\
	$(GETMCHBIN)/seqIOALF.o\
	$(GETMCHBIN)/seqIOSCF.o\
	$(GETMCHBIN)/seqOutput.o\
	$(GETMCHBIN)/opp.o\
	$(GETMCHBIN)/fpoint.o\
	$(GETMCHBIN)/match.o\
	$(GETMCHBIN)/mach-io.o


#
# Linking
#
$(O)/trace2seq: $(OBJS)
	$(CLD) -o $@ $(OBJS) $(EXP_LIB) $(LIBSC)

#
# Special rule needed for Graph.o, seq*.o from the ted directory.
#
$(O)/%.o: $(TEDSRC)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

DEPEND_OBJ = $(OBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
