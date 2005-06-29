#
# Makefile for example io_lib programs
#

include ../options.mk

PROGS = $(O)/convert_trace $(O)/index_tar $(O)/get_comment $(O)/trace_dump \
	$(O)/hash_tar

ifdef IOLIB_PLN
PROGS += $(O)/extract_seq
endif

ifdef IOLIB_SCF
PROGS += $(O)/scf_update $(O)/scf_info $(O)/scf_dump $(O)/makeSCF
endif

ifdef IOLIB_ZTR
PROGS += $(O)/ztr_dump
endif

SRCROOT=../..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk
CFLAGS += -I../include

IOLIBBIN=$(L)
INCLUDES_X += -I$(IOLIBSRC)/include
CLDFLAGS   += $(LINK_PATHFLAG)$(IOLIBBIN)

ESOBJS=$(O)/extract_seq.o
SUOBJS=$(O)/scf_update.o
SIOBJS=$(O)/scf_info.o
STOBJS=$(O)/scf_dump.o
TDOBJS=$(O)/trace_dump.o
ZDOBJS=$(O)/ztr_dump.o
MSOBJS=$(O)/makeSCF.o
CTOBJS=$(O)/convert_trace.o
ITOBJS=$(O)/index_tar.o
HTOBJS=$(O)/hash_tar.o
GCOBJS=$(O)/get_comment.o

PLIBS=$(LINK_LIBFLAG)read$(LIB_EXT) $(LIBSC) $(LINK_LIBFLAG)misc$(LIB_EXT)

$(O)/extract_seq:	$(ESOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(ESOBJS) $(PLIBS)

$(O)/get_comment:	$(GCOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(GCOBJS) $(PLIBS)

$(O)/scf_update:	$(SUOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SUOBJS) $(PLIBS)

$(O)/scf_info: 		$(SIOBJS) 
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(SIOBJS) $(PLIBS)

$(O)/scf_dump: 		$(STOBJS) 
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(STOBJS) $(PLIBS)

$(O)/trace_dump: 	$(TDOBJS) 
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(TDOBJS) $(PLIBS)

$(O)/ztr_dump: 	$(ZDOBJS) 
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(ZDOBJS) $(PLIBS) $(MATH_LIB)

$(O)/makeSCF: 		$(MSOBJS) 
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(MSOBJS) $(PLIBS)

$(O)/convert_trace: 	$(CTOBJS) 
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(CTOBJS) $(PLIBS)

$(O)/index_tar:		$(ITOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(ITOBJS) $(PLIBS)

$(O)/hash_tar:		$(HTOBJS)
	$(CLD) $(LDEXEFLAG)$@$(EXE_SUFFIX) $(HTOBJS) $(PLIBS)

DEPEND_OBJ = $(ESOBJS) $(SUOBJS) $(SIOBJS) $(STOBJS) $(TDOBJS) $(MSOBJS) $(CTOBJS) $(ITOBJS) $(GCOBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

#clean:
#	-rm -f $(O)/*.o

spotless:	clean
	-rm -f $(PROGS)

include dependencies
# DO NOT DELETE THIS LINE -- make depend depends on it.
