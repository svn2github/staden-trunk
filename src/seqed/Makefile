#
# Makefile for seqed routines
#

LIBS = seqed
PROGS = x$(LIBS) 
PROGLIBS=$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX)	

SRCROOT=..
include	$(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(TK_INC) $(TKUTILS_INC) $(SEQUTILS_INC) $(IOLIB_INC) $(SPIN_INC) $(ITK_INC)

CFLAGS += $(SHLIB_CFLAGS)

SEQEDBIN=$(O)


OBJS = \
	$(SEQEDBIN)/read_sequence.o \
	$(SEQEDBIN)/parse_feature.o \
	$(SEQEDBIN)/feature_table.o \
	$(SEQEDBIN)/feature_colour.o \
	$(SEQEDBIN)/parse_db.o \
	$(SEQEDBIN)/renzyme_box.o \
	$(SEQEDBIN)/renzyme_search.o \
	$(SEQEDBIN)/graphic_editor.o \
	$(SEQEDBIN)/renzyme_map_canvas.o \
	$(SEQEDBIN)/text_editor.o \
	$(SEQEDBIN)/editor.o \
	$(SEQEDBIN)/end_editor.o \
	$(SEQEDBIN)/editor_reg.o \
	$(SEQEDBIN)/feature_editor.o 


#
# Main dependency
#
x$(LIBS) : $(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX) 
	@

DEPS = $(TCL_LIB)\
       $(TK_LIB) \
       $(IOUTILS_LIB) \
       $(TKUTILS_LIB) \
       $(SEQUTILS_LIB) \
       $(MISC_LIB) \
       $(SPIN_LIB) \
       $(ITCL_LIB)

$(L)/$(SHLIB_PREFIX)$(LIBS)$(SHLIB_SUFFIX): $(OBJS) $(DEF_FILE)
	$(SHLIB_LD) $(SHLIB_LDFLAGS) $(SHLIB_OUTFLAG)$@ $(SHLIB_SONAME) $(OBJS) $(DEPS) 

# 7/1/99 johnt - Rule used when $(DEF_FILE) defined - currently only for windows
$(L)/$(SHLIB_PREFIX)$(LIBS).def: $(OBJS)
	$(MKDEFL) $@ $(OBJS)


DEPEND_OBJ = $(OBJS)

distsrc: distsrc_dirs
	-cp -R *.[ch] *.tcl tclIndex Makefile seqed seqed.bat dependencies \
		$(DIRNAME)

install:
#	cp seqed /home/yc/seqed/linux-bin/
#	cp .seqedrc /home/yc/
#	-mkdir /home/yc/seqed/lib
#	cp *.itcl *.gif *.tcl *.itk tclIndex /home/yc/seqed/lib	
#	-mkdir /home/yc/seqed/lib/linux-binaries

	cp seqed $(INSTALLSCRIPT)
	-mkdir $(INSTALLLIB)/seqed
	cp *.itcl *.gif *.tcl *.itk tclIndex $(INSTALLLIB)/seqed
	cp $(PROGLIBS) $(INSTALLLIB)/$(O)

# DO NOT DELETE THIS LINE -- make depend depends on it.
