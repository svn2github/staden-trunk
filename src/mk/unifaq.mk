# FAK Makefile

PROGS=	$(O)/create_graph $(O)/add_to_graph $(O)/show_graph $(O)/show_multi \
	$(O)/show_layout $(O)/assemble $(O)/write_exp_file \
	$(O)/create_exp_constraints $(O)/caf_multi

FLIB=KERNEL

SRCROOT=../..
include $(SRCROOT)/mk/global.mk
include $(FLIB)/$(MACHINE).mk

FAKBIN=$(O)
FAKLIBS=-L$(FLIB)/$(O) -lfa -lm

INCLUDES_E += $(SEQUTILS_INC) -I$(GAP4SRC) $(G_INC) $(IOLIB_INC) -I$(FLIB)

CGOBJS=\
	$(FAKBIN)/library.o \
	$(FAKBIN)/create_graph.o \
	$(FAKBIN)/fak_exp_utils.o \
	$(FAKBIN)/exp_utils.o \
	$(GAP4BIN)/seqInfo.o \
	$(GAP4BIN)/scf_extras.o \
	$(SEQUTILSBIN)/dna_utils.o

COOBJS=\
	$(FAKBIN)/library.o \
	$(FAKBIN)/create_overlap.o

AGOBJS=\
	$(FAKBIN)/add_to_graph.o \
	$(FAKBIN)/library.o \
	$(FAKBIN)/fak_exp_utils.o \
	$(FAKBIN)/exp_utils.o \
	$(GAP4BIN)/seqInfo.o \
	$(GAP4BIN)/scf_extras.o \
	$(SEQUTILSBIN)/dna_utils.o

SGOBJS=\
	$(FAKBIN)/show_graph.o \
	$(FAKBIN)/regexp.o \
	$(FAKBIN)/library.o

ASOBJS=\
	$(FAKBIN)/assemble.o \
	$(FAKBIN)/library.o

SMOBJS=\
	$(FAKBIN)/show_multi.o \
	$(FAKBIN)/library.o

SLOBJS=\
	$(FAKBIN)/show_layout.o \
	$(FAKBIN)/library.o

ECOBJS=\
	$(FAKBIN)/create_exp_constraints.o\
	$(FAKBIN)/fak_exp_utils.o \
	$(FAKBIN)/exp_utils.o \
	$(GAP4BIN)/seqInfo.o \
	$(GAP4BIN)/scf_extras.o	\
	$(SEQUTILSBIN)/dna_utils.o \
	$(FAKBIN)/library.o

EFOBJS= \
	$(FAKBIN)/library.o \
	$(FAKBIN)/write_exp_file.o \
	$(FAKBIN)/fak_exp_utils.o \
	$(FAKBIN)/exp_utils.o \
	$(GAP4BIN)/seqInfo.o \
	$(GAP4BIN)/scf_extras.o	\
	$(SEQUTILSBIN)/dna_utils.o

CMOBJS=\
	$(FAKBIN)/library.o \
	$(FAKBIN)/caf_multi.o

FAK_LIBS = \
	$(IOLIB_LIB) \
	$(MISC_LIB) \
	$(TEXTUTILS_LIB)

$(O)/create_graph:	$(CGOBJS)
	$(CLD) -o $@ $(CGOBJS) $(FAKLIBS) $(FAK_LIBS) $(LIBSC)

$(O)/add_to_graph:	$(AGOBJS)
	$(CLD) -o $@ $(AGOBJS) $(FAKLIBS) $(FAK_LIBS) $(LIBSC)

$(O)/show_graph:	$(SGOBJS)
	$(CLD) -o $@ $(SGOBJS) $(FAKLIBS) $(LIBSC)

$(O)/assemble:		$(ASOBJS)
	$(CLD) -o $@ $(ASOBJS) $(FAKLIBS) $(LIBSC)

$(O)/show_multi:	$(SMOBJS)
	$(CLD) -o $@ $(SMOBJS) $(FAKLIBS) $(LIBSC)

$(O)/show_layout:	$(SLOBJS)
	$(CLD) -o $@ $(SLOBJS) $(FAKLIBS) $(LIBSC)

$(O)/write_exp_file:	$(EFOBJS)
	$(CLD) -o $@ $(EFOBJS) $(FAKLIBS) $(FAK_LIBS) $(LIBSC)

$(O)/create_exp_constraints:	$(ECOBJS)
	$(CLD) -o $@ $(ECOBJS) $(FAKLIBS) $(FAK_LIBS) $(LIBSC)

$(O)/caf_multi:		$(CMOBJS)
	$(CLD) -o $@ $(CMOBJS) $(FAKLIBS) $(LIBSC)



