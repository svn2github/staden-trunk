PROGS= 	$(O)/gap	$(O)/xgap

# Local comment: Uncomment next line for remote compilation
# REMOTESRC=$(GAPSRC)

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

# Local comment: Comment out next line for remote compilation
GAPBIN=$(O)

# What a mess gap is! The order is important due to include files with
# identical names. We must have staden before ted, and ted before io_lib.
INCLUDES_E += $(G_INC) $(XAW_INC) -I$(OSPSRC) \
	      -I$(STADENSRC) -I$(TEDSRC) $(EXP_INC) $(TT_INC)

#
# The C objects, needed by every X program
#
CCORE=\
	$(STADENBIN)/postscript.o\
	$(GAPBIN)/Graph.o\
	$(STADENBIN)/plotLog.o\
	$(STADENBIN)/help.o\
	$(STADENBIN)/dialogues.o\
	$(STADENBIN)/userfacecom.o\
	$(STADENBIN)/xmenu.o\
	$(STADENBIN)/mcspec.o\
	$(STADENBIN)/helpnmenu.o\
	$(STADENBIN)/FtoC.o

#COBJS=\
#	$(STADENBIN)/main.o\
#	$(CCORE)

COBJS2=\
	$(GAPBIN)/xgapmain.o\
	$(CCORE)


#
# Common objects, needed by every program
#
COMMONOBJSB=\
	$(STADENBIN)/seeme-$(MACHINE).o\
	$(STADENBIN)/Cme-$(MACHINE).o\
	$(STADENBIN)/nxspec.o\
	$(STADENBIN)/userface.o\
	$(STADENBIN)/userfacecom.o\
	$(STADENBIN)/nxhelpmenu.o\
	$(STADENBIN)/helpnmenu.o\
	$(STADENBIN)/FtoC.o\
	$(STADENBIN)/subs89.a

COMMONOBJS=\
	$(COMMONOBJSB)\
	$(STADENBIN)/postscript.o

#SCOMMONOBJS=\
#	$(COMMONOBJSB)

#XCOMMONOBJS=\
#	$(STADENBIN)/seeme-$(MACHINE).o\
#	$(STADENBIN)/Cme-$(MACHINE).o\
#	$(STADENBIN)/xspec.o\
#	$(STADENBIN)/subs89.a\
#	$(COBJS)

XCOMMONOBJS2=\
	$(STADENBIN)/seeme-$(MACHINE).o\
	$(STADENBIN)/Cme-$(MACHINE).o\
	$(STADENBIN)/xspec.o\
	$(STADENBIN)/subs89.a\
	$(COBJS2)


#
# Building the programs
# This should be just a linking phase because all of the object
# files and library files are generated using implicit rules.
# We use the fortran compiler to do linking.
#

#
# Trace manager objects
#
TEDFILES=\
	$(STADENBIN)/seq.o\
	$(GAPBIN)/seqIOSCF.o\
	$(STADENBIN)/seqIOABI.o\
	$(STADENBIN)/seqIOALF.o\
	$(STADENBIN)/seqIOPlain.o\
	$(STADENBIN)/mach-io.o\
	$(STADENBIN)/fpoint.o

XTEDFILES=\
	$(TEDFILES)\
	$(STADENBIN)/opp.o\
	$(STADENBIN)/info.o

TMANOBJS=\
	$(GAPBIN)/tman_main.o\
	$(GAPBIN)/tman_display.o\
	$(GAPBIN)/tman_traceDisplay.o\
	$(GAPBIN)/tman_basesDisplay.o\
	$(GAPBIN)/tman_context.o\
	$(GAPBIN)/tman_gadgets.o\
	$(GAPBIN)/tman_interface.o\
	$(XTEDFILES)

OSPOBJS=\
	$(OSPBIN)/analysis.o\
	$(OSPBIN)/our_allo.o\
	$(OSPBIN)/paramIOX.o\
	$(OSPBIN)/paramIO.o\
	$(OSPBIN)/get_scores.o\
	$(OSPBIN)/utils.o\
	$(GAPBIN)/mess.o\
	$(GNULIB)

XHACK=\
	$(STADENBIN)/SmeLine.o

#	$(GAPBIN)/oligo-needs-work.o
EDITOR=\
	$(XHACK)\
	$(GAPBIN)/xsapConEdit.o\
	$(GAPBIN)/contigEditor.o\
	$(GAPBIN)/join.o\
	$(GAPBIN)/edUtils.o\
	$(GAPBIN)/undo.o\
	$(GAPBIN)/Sheet.o\
	$(GAPBIN)/select.o\
	$(GAPBIN)/extend.o\
	$(GAPBIN)/searchUtils.o\
	$(GAPBIN)/edMenu.o\
	$(GAPBIN)/oligocom.o\
	$(GAPBIN)/oligo.o\
	$(GAPBIN)/subclone.o\
	$(GAPBIN)/myparams.o\
	$(STADENBIN)/trans.o

TAGEDITOR=\
	$(GAPBIN)/tagEditor.o\
	$(GAPBIN)/tagU2.o\
	$(GAPBIN)/tagU1.o


GAP=\
	$(GAPBIN)/gap-local.o\
	$(GAPBIN)/gap-remote.o\
	$(GAPBIN)/gap-if.o\
	$(GAPBIN)/gap-init.o\
	$(GAPBIN)/gap-dbstruct.o\
	$(GAPBIN)/gap-create.o\
	$(GAPBIN)/gap-error.o\
	$(GAPBIN)/gap-io.o


NEWDB=\
	$(GAP)\
	$(GAPBIN)/IO.o\
	$(GAPBIN)/IO2.o\
	$(GAPBIN)/IO3.o\
	$(GAPBIN)/seqInfo.o

#	$(GAPBIN)/cop.o
NEWOPTS=\
	$(GAPBIN)/dbfix.o\
	$(GAPBIN)/dbcheck.o\
	$(GAPBIN)/clones.o\
	$(GAPBIN)/trace.o\
	$(GAPBIN)/readpair.o\
	$(GAPBIN)/extract.o




BAP=\
	$(GAPBIN)/cli.o\
        $(GAPBIN)/align.o\
        $(GAPBIN)/align_lib.o\
        $(GAPBIN)/align_ss.o\
        $(GAPBIN)/align_sv.o\
	$(GAPBIN)/gap.o\
	$(GAPBIN)/dbsysnew.o\
	$(GAPBIN)/dbsyscommon.o\
	$(GAPBIN)/consen.o\
	$(GAPBIN)/active_tags.o\
	$(GAPBIN)/dbsysplus.o\
	$(NEWDB)\
	$(NEWOPTS)\
	$(GAPBIN)/actf.o\
	$(GAPBIN)/dstrand.o\
	$(GAPBIN)/oligo_sel.o\
	$(GAPBIN)/qual.o\
	$(GAPBIN)/reactions.o\
	$(GAPBIN)/list.o\
	$(GAPBIN)/utils.o\
	$(GAPBIN)/preass.o\
	$(STADENBIN)/asubs89.a\
	$(STADENBIN)/plot92.a


# A bit of a lie in places
GAPOBJS=\
	$(BAP)\
	$(COMMONOBJS)\
	$(STADENBIN)/pl4010.o\
	$(STADENBIN)/conEdit.o\
	$(GAPBIN)/tagU2.o\
	$(OSPBIN)/analysis.o\
	$(OSPBIN)/our_allo.o\
	$(OSPBIN)/get_scores.o\
	$(OSPBIN)/utils.o\
	$(GAPBIN)/mess.o\
	$(GAPBIN)/gaponly.o\
	$(GAPBIN)/oligocom.o\
	$(TEDFILES)\
	$(GAPLIBS)


XGAPOBJS=\
	$(BAP)\
	$(XCOMMONOBJS2)\
	$(GAPBIN)/plX.o\
	$(STADENBIN)/textOutput.o\
	$(EDITOR)\
	$(TMANOBJS)\
	$(TT_OBJS) \
	$(TAGEDITOR)\
	$(OSPOBJS)

GAP_LIBS=\
	$(G_LIB) \
	$(EXP_LIB) \
	$(MISC_LIB) \
	$(TK_LIB)

GAP_DEP=\
	$(G_DEP) \
	$(EXP_DEP) \
	$(MISC_DEP)

$(O)/gap: $(GAPOBJS)
	$(FLD) -o $@ $(GAPOBJS) $(GAP_LIBS) $(LIBSF) $(GAP_DEP)

$(O)/xgap: $(XGAPOBJS)
	$(FLD) -o $@ $(XGAPOBJS) $(TT_LIB) $(GAP_LIBS) $(XAW_LIB) $(LIBSF) $(GAP_DEP)

DEPEND_OBJ =\
	$(GAPOBJS) $(XGAPOBJS)

backup: FORCE
	tar cvf - *.[chf] README TODO help/GAP.RNO tables/Xgap tables/TAGDB \
	    Makefile ../ChangeLog dependencies\
	    | gzip > backup/`date +%d_%m_%y`.tar.gz

FORCE::

include dependencies

install:
	cp $(PROGS) $(INSTALLBIN)
