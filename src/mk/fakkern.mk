# FAKLIB Makefile

LIBS = $(O)/libfa.a
PROGS= $(LIBS) $(O)/fa_scoretab.i

SRCROOT=../..
include $(SRCROOT)/mk/global.mk
include $(MACHINE).mk

FAKLIBBIN=$(O)

OBJS = \
	$(FAKLIBBIN)/fa_assemble.o \
	$(FAKLIBBIN)/fa_lists.o \
	$(FAKLIBBIN)/fa_cnst.o \
	$(FAKLIBBIN)/fa_unionfind.o \
	$(FAKLIBBIN)/fa_multi.o \
	$(FAKLIBBIN)/fa_edit.o \
	$(FAKLIBBIN)/fa_system.o \
	$(FAKLIBBIN)/fa_graph.o \
	$(FAKLIBBIN)/fa_screen.o \
	$(FAKLIBBIN)/fa_overlap.o \
	$(FAKLIBBIN)/fa_score.o \
	$(FAKLIBBIN)/fa_reduce.o \
	$(FAKLIBBIN)/fa_asmrw.o \
	$(FAKLIBBIN)/fa_ultra.o \
	$(FAKLIBBIN)/fa_cnstrplc.o

#
# Main dependency
#
$(LIBS) : $(OBJS)
	/bin/rm -f $(LIBS) ;\
	$(AR) $(ARFLAGS) $(LIBS) $(OBJS) ;\
	$(RANLIB) $(LIBS)

SOBJ = \
	$(FAKLIBBIN)/fa_system.o \
	$(FAKLIBBIN)/fa_scoregen.o

$(O)/fa_scoretab.i : $(SOBJ)
	$(CLD) $(SOBJ) -o $(O)/fa_scoretab_gen $(MATH_LIB) $(LIBSC)
	cd $(O); ./fa_scoretab_gen

DEPEND_OBJ = $(OBJS)

install:

include dependencies
