#
# Makefile for ted (trace editor)
#

PROGS = $(O)/ted $(O)/autoted $(O)/autotede

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

NEWTEDBIN=$(O)
INCLUDES_E += $(XAW_INC)

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
# Parameters for QUAL_CODE and QUAL_CHECK are defined at the top of
# seqIOEdit.c. You may wish to tune these yourself.
#
#DEFINES += -DAUTO_CLIP  -DDEF_OUT  -DSAVE_EDITS  -DQUAL_CODE -DQUAL_CHECK
DEFINES += -DAUTO_CLIP -DQUAL_CODE -DQUAL_CHECK

#
# Object files
#
TEDOBJS=\
	$(NEWTEDBIN)/ted.o\
	$(NEWTEDBIN)/dialogues.o\
	$(NEWTEDBIN)/seq.o\
	$(NEWTEDBIN)/seqIOPlain.o\
	$(NEWTEDBIN)/seqIOABI.o\
	$(NEWTEDBIN)/help.o\
	$(NEWTEDBIN)/display.o\
	$(NEWTEDBIN)/traceDisplay.o\
	$(NEWTEDBIN)/basesDisplay.o\
	$(NEWTEDBIN)/Graph.o\
	$(NEWTEDBIN)/seqIOEdit.o\
	$(NEWTEDBIN)/seqIOALF.o\
	$(NEWTEDBIN)/seqIOSCF.o\
	$(NEWTEDBIN)/seqOutput.o\
	$(NEWTEDBIN)/opp.o\
	$(NEWTEDBIN)/info.o\
	$(NEWTEDBIN)/fpoint.o\
	$(NEWTEDBIN)/mach-io.o\
	$(NEWTEDBIN)/seqRead.o\
	$(NEWTEDBIN)/traceType.o\
	$(NEWTEDBIN)/match.o\
	$(NEWTEDBIN)/expFileIO.o

AUTOTEDCOBJS=\
	$(NEWTEDBIN)/seq.o\
	$(NEWTEDBIN)/seqIOPlain.o\
	$(NEWTEDBIN)/seqIOABI.o\
	$(NEWTEDBIN)/seqIOEdit.o\
	$(NEWTEDBIN)/seqIOALF.o\
	$(NEWTEDBIN)/seqIOSCF.o\
	$(NEWTEDBIN)/seqOutput.o\
	$(NEWTEDBIN)/opp.o\
	$(NEWTEDBIN)/traceType.o\
	$(NEWTEDBIN)/fpoint.o\
	$(NEWTEDBIN)/mach-io.o\
	$(NEWTEDBIN)/match.o\
	$(NEWTEDBIN)/seqRead.o\
	$(NEWTEDBIN)/expFileIO.o

AUTOTEDOBJS=\
        $(NEWTEDBIN)/autoted.o\
        $(AUTOTEDCOBJS)

AUTOTEDEOBJS=\
        $(NEWTEDBIN)/autotede.o\
        $(AUTOTEDCOBJS)

#
# Linking
#
$(O)/ted: $(TEDOBJS)
	$(CLD) -o $@ $(TEDOBJS) $(XAW_LIB) $(MISC_LIB) $(LIBSC)

$(O)/autoted: $(AUTOTEDOBJS)
	$(CLD) -o $@ $(AUTOTEDOBJS) $(XAW_LIB) $(MISC_LIB) $(LIBSC)

$(O)/autotede: $(AUTOTEDEOBJS)
	$(CLD) -o $@ $(AUTOTEDEOBJS) $(XAW_LIB) $(MISC_LIB) $(LIBSC)

DEPEND_OBJ = $(TEDOBJS) $(AUTOTEDOBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
