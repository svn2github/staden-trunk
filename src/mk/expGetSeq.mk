#
# Makefile for expGetSeq
#

PROGS = $(O)/expGetSeq

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += -I$(TEDSRC) $(XAW_INC)
EXPGETSEQBIN=$(O)

#
# Turning on the AUTO_CLIP switch allows ted to automatically
# select a left (using the -enzyme option or the .enzyme in the
# Xted file) and right cutoff on your sequence (using 2 out of
# 5 N's)
#
# Turning on the SAVE_EDITS switch allows the user to maintain
# copies of their edits, i.e. when you open up ted on a file
# that you have already edited, the old edits show up in the
# ted edit window.  The user may also call up any of their
# past editing sessions.  
#
# Turning on the DEF_OUT  switch makes
# trace_name.seq the default output file name
#
# Turning on the QUAL_CODE switch enables LaDeana's trace quality
# clipping code.
#
# Turning on the QUAL_CHECK switch (at the same time as QUAL_CODE)
# enables the overall trace quality check.
#
#DEFINES += -DAUTO_CLIP  -DDEF_OUT  -DSAVE_EDITS  -DQUAL_CODE -DQUAL_CHECK
DEFINES += -DAUTO_CLIP -DQUAL_CODE

#
# Object files
#
OBJS=\
	$(EXPGETSEQBIN)/getMCH.o\
	$(EXPGETSEQBIN)/seq.o\
	$(EXPGETSEQBIN)/seqIOPlain.o\
	$(EXPGETSEQBIN)/seqIOABI.o\
	$(EXPGETSEQBIN)/seqIOEdit.o\
	$(EXPGETSEQBIN)/seqIOALF.o\
	$(EXPGETSEQBIN)/seqIOSCF.o\
	$(EXPGETSEQBIN)/opp.o\
	$(EXPGETSEQBIN)/fpoint.o\
	$(EXPGETSEQBIN)/match.o\
	$(EXPGETSEQBIN)/mach-io.o\
	$(EXPGETSEQBIN)/seqRead.o\
	$(EXPGETSEQBIN)/traceType.o


#
# Linking
#
$(O)/expGetSeq: $(OBJS)
	$(CLD) -o $@ $(OBJS) $(XAW_LIB) $(LIBSC)

#
# Special rule needed for match.o, seq*.o from the newted directory.
#
$(O)/%.o: $(TEDSRC)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

DEPEND_OBJ = $(OBJS)

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
