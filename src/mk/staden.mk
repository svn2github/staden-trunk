NPROGS=	$(O)/mep	$(O)/nip	$(O)/pip	\
	$(O)/sap	$(O)/sip	\
	$(O)/splitp3	$(O)/gip	$(O)/sethelp	$(O)/rep	\
	$(O)/lip	$(O)/nipf	$(O)/vep	\
	$(O)/pipf	$(O)/splitseq	$(O)/splitseqf

#	$(O)/convert_project $(O)/splitp1 $(O)/splitp2 $(O)/sipl
#	$(O)/xsap $(O)/xdap $(O)/dap $(O)/sapf

LPROGS=	$(O)/nipl	$(O)/pipl

XPROGS=	$(O)/xmep	$(O)/xnip	$(O)/xpip	$(O)/xsip	

PROGS= $(NPROGS) $(LPROGS) $(XPROGS)

# Redefine PROGS - for the time being we only want to compile and
# install splitseq
PROGS = $(O)/splitseq

#
# For 'remote' compilation, uncomment the following two lines and comment out
# the current SRCROOT definition. Also, don't forget to copy the dependencies
# file and to create the $MACHINE-binaries directory.
#
# REMOTESRC=$(STADENSRC)
# SRCROOT=$(STADENROOT)/src
#

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

INCLUDES_E += $(XAW_INC) -I$(TEDSRC) $(TT_INC)

STADENBIN=$(O)

#
# Sequence library handling routines
#
SEQLIB=\
	$(STADENBIN)/seqlibsubs.o\
	$(STADENBIN)/bit.o

#
# The C objects, needed by every X program
#
CCORE=\
	$(STADENBIN)/postscript.o\
	$(STADENBIN)/Graph.o\
	$(STADENBIN)/plotLog.o\
	$(STADENBIN)/help.o\
	$(STADENBIN)/dialogues.o\
	$(STADENBIN)/userfacecom.o\
	$(STADENBIN)/helpnmenu.o\
	$(STADENBIN)/xmenu.o\
	$(STADENBIN)/FtoC.o\
	$(STADENBIN)/mcspec.o

COBJS=\
	$(STADENBIN)/main.o\
	$(CCORE)

XDAPCOBJS=\
	$(STADENBIN)/xdapmain.o\
	$(CCORE)


#
# Common objects, needed by every program
#
COMMONOBJSB=\
	$(STADENBIN)/fmain.o\
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

SCOMMONOBJS=\
	$(COMMONOBJSB)

XCOMMONOBJS=\
	$(STADENBIN)/seeme-$(MACHINE).o\
	$(STADENBIN)/Cme-$(MACHINE).o\
	$(STADENBIN)/xspec.o\
	$(STADENBIN)/subs89.a\
	$(COBJS)

XDAPCOMMONOBJS=\
	$(STADENBIN)/seeme-$(MACHINE).o\
	$(STADENBIN)/Cme-$(MACHINE).o\
	$(STADENBIN)/xspec.o\
	$(STADENBIN)/subs89.a\
	$(XDAPCOBJS)


#
# Building the programs
# This should be just a linking phase because all of the object
# files and library files are generated using implicit rules.
# We use the fortran compiler to do linking.
#
GIP=\
	$(STADENBIN)/gip.o

GIPOBJS=\
	$(GIP)\
	$(SCOMMONOBJS)

$(O)/gip: $(GIPOBJS)
	$(FLD) -o $@ $(GIPOBJS) $(LIBSF)


LIP=\
	$(STADENBIN)/lip.o
 
LIPOBJS=\
	$(LIP)\
	$(SEQLIB)\
	$(COMMONOBJS)
 
$(O)/lip: $(LIPOBJS)
	$(FLD) -o $@ $(LIPOBJS) $(LIBSF)


MEP=\
	$(STADENBIN)/mep.o\
	$(STADENBIN)/mepsub.o\
	$(STADENBIN)/asubs89.a\
	$(STADENBIN)/plot92.a

MEPOBJS=\
	$(MEP)\
	$(COMMONOBJS)\
	$(STADENBIN)/pl4010.o

XMEPOBJS=\
	$(MEP)\
	$(XCOMMONOBJS)\
	$(STADENBIN)/plX.o\
	$(STADENBIN)/textOutput.o 

$(O)/mep: $(MEPOBJS)
	$(FLD) -o $@ $(MEPOBJS) $(LIBSF)

$(O)/xmep: $(XMEPOBJS)
	$(FLD) -o $@ $(XMEPOBJS) $(XAW_LIB) $(MISC_LIB) $(LIBSF)




NIP=\
	$(STADENBIN)/nip.o\
	$(STADENBIN)/patternn.a\
	$(STADENBIN)/patternnc.a\
	$(STADENBIN)/anals89.a\
	$(STADENBIN)/asubs89.a\
	$(STADENBIN)/plot92.a \
	$(SEQLIB)

NIPOBJS=\
	$(NIP)\
	$(COMMONOBJS)\
	$(STADENBIN)/pl4010.o

XNIPOBJS=\
	$(NIP)\
	$(XCOMMONOBJS)\
	$(STADENBIN)/plX.o\
	$(STADENBIN)/textOutput.o 

$(O)/nip: $(NIPOBJS)
	$(FLD) -o $@ $(NIPOBJS) $(LIBSF)

$(O)/xnip: $(XNIPOBJS)
	$(FLD) -o $@ $(XNIPOBJS) $(XAW_LIB) $(MISC_LIB) $(LIBSF)




NIPL=\
	$(STADENBIN)/nipl.o\
	$(STADENBIN)/patternnc.a\
	$(STADENBIN)/anals89.a\
	$(STADENBIN)/asubs89.a\
	$(SEQLIB)

NIPLOBJS=\
	$(NIPL)\
	$(SCOMMONOBJS)

$(O)/nipl: $(NIPLOBJS)
	$(FLD) -o $@ $(NIPLOBJS) $(LIBSF)



NIPF=\
	$(STADENBIN)/nipf.o\
	$(STADENBIN)/asubs89.a\
	$(STADENBIN)/plot92.a

NIPFOBJS=\
	$(NIPF)\
	$(COMMONOBJS)\
	$(STADENBIN)/pl4010.o


$(O)/nipf: $(NIPFOBJS)
	$(FLD) -o $@ $(NIPFOBJS) $(LIBSF)



PIPF=\
	$(STADENBIN)/pipf.o\
	$(STADENBIN)/asubs89.a\
	$(STADENBIN)/plot92.a

PIPFOBJS=\
	$(PIPF)\
	$(COMMONOBJS)\
	$(STADENBIN)/pl4010.o


$(O)/pipf: $(PIPFOBJS)
	$(FLD) -o $@ $(PIPFOBJS) $(LIBSF)





PIP=\
	$(STADENBIN)/pip.o\
	$(STADENBIN)/analps89.a\
	$(STADENBIN)/patternp.a\
	$(STADENBIN)/patternpc.a\
	$(STADENBIN)/asubs89.a\
	$(STADENBIN)/plot92.a\
	$(SEQLIB)

PIPOBJS=\
	$(PIP)\
	$(COMMONOBJS)\
	$(STADENBIN)/pl4010.o

XPIPOBJS=\
	$(PIP)\
	$(XCOMMONOBJS)\
	$(STADENBIN)/plX.o\
	$(STADENBIN)/textOutput.o 

$(O)/pip: $(PIPOBJS)
	$(FLD) -o $@ $(PIPOBJS) $(LIBSF)

$(O)/xpip:$(XPIPOBJS)
	$(FLD) -o $@ $(XPIPOBJS) $(XAW_LIB) $(MISC_LIB) $(LIBSF)




PIPL=\
	$(STADENBIN)/pipl.o\
	$(STADENBIN)/patternpc.a\
	$(STADENBIN)/analps89.a\
	$(STADENBIN)/asubs89.a\
	$(SEQLIB)

PIPLOBJS=\
	$(PIPL)\
	$(SCOMMONOBJS)

$(O)/pipl: $(PIPLOBJS)
	$(FLD) -o $@ $(PIPLOBJS) $(LIBSF)



#
# Trace manager objects
#
STDTEDFILES=\
	$(STADENBIN)/seq.o\
	$(STADENBIN)/seqIOABI.o\
	$(STADENBIN)/seqIOALF.o\
	$(STADENBIN)/seqIOSCF.o\
	$(STADENBIN)/seqIOPlain.o\
	$(STADENBIN)/opp.o\
	$(STADENBIN)/fpoint.o\
	$(STADENBIN)/mach-io.o\
	$(STADENBIN)/info.o

TMANOBJS=\
	$(STADENBIN)/tman_main.o\
	$(STADENBIN)/tman_display.o\
	$(STADENBIN)/tman_traceDisplay.o\
	$(STADENBIN)/tman_basesDisplay.o\
	$(STADENBIN)/tman_context.o\
	$(STADENBIN)/tman_gadgets.o\
	$(STDTEDFILES)

# Some versions of X11R4 may have a bug in SmeLine.c
XHACK=\
	$(STADENBIN)/SmeLine.o

EDITOR=\
	$(XHACK)\
	$(STADENBIN)/xsapConEdit.o\
	$(STADENBIN)/contigEditor.o\
	$(STADENBIN)/edUtils.o\
	$(STADENBIN)/undo.o\
	$(STADENBIN)/Sheet.o\
	$(STADENBIN)/select.o\
	$(STADENBIN)/extend.o\
	$(STADENBIN)/searchUtils.o\
	$(STADENBIN)/edMenu.o\
	$(STADENBIN)/trans.o

TAGEDITOR=\
	$(STADENBIN)/tagEditor.o\
	$(STADENBIN)/tagdbparse.o\
	$(STADENBIN)/tagU2.o\
	$(STADENBIN)/tagU1.o

DAP=\
	$(STADENBIN)/dap.o\
	$(STADENBIN)/dbsysnew.o\
	$(STADENBIN)/dbsyscommon.o\
	$(STADENBIN)/asubs89.a\
	$(STADENBIN)/plot92.a

DAPOBJS=\
	$(DAP)\
	$(COMMONOBJS)\
	$(STADENBIN)/pl4010.o\
	$(STADENBIN)/conEdit.o\
	$(STADENBIN)/tagU2.o

XDAPOBJS=\
	$(DAP)\
	$(XDAPCOMMONOBJS)\
	$(STADENBIN)/plX.o\
	$(STADENBIN)/textOutput.o\
	$(EDITOR)\
	$(TMANOBJS)\
	$(TT_OBJS) \
	$(TAGEDITOR)

$(O)/dap: $(DAPOBJS)
	$(FLD) -o $@ $(DAPOBJS) $(LIBSF)

$(O)/xdap: $(XDAPOBJS)
	$(FLD) -o $@ $(XDAPOBJS) $(TT_LIB) $(XAW_LIB) $(MISC_LIB) $(LIBSF)

$(O)/convert_project: cvt.o
	$(CLD) -o $@ cvt.o $(LIBSC)

SAP=\
	$(STADENBIN)/sap.o\
	$(STADENBIN)/dbsyscommon.o\
	$(STADENBIN)/dbsysold.o\
	$(STADENBIN)/asubs89.a\
	$(STADENBIN)/plot92.a

SAPOBJS=\
	$(SAP)\
	$(COMMONOBJS)\
	$(STADENBIN)/pl4010.o

XSAPOBJS=\
	$(SAP)\
	$(XCOMMONOBJS)\
	$(STADENBIN)/plX.o\
	$(STADENBIN)/textOutput.o 

$(O)/sap: $(SAPOBJS)
	$(FLD) -o $@ $(SAPOBJS) $(MISC_LIB) $(LIBSF)

$(O)/xsap: $(XSAPOBJS)
	$(FLD) -o $@ $(XSAPOBJS) $(XAW_LIB) $(MISC_LIB) $(LIBSF)

SAPF=\
	$(STADENBIN)/sapf.o\
	$(STADENBIN)/dbsyscommon.o\
	$(STADENBIN)/dbsysold.o\
	$(STADENBIN)/asubs89.a\
	$(STADENBIN)/plot92.a

SAPFOBJS=\
	$(SAPF)\
	$(STADENBIN)/sapmen.o\
	$(COMMONOBJS)\
	$(STADENBIN)/pl4010.o

$(O)/sapf: $(SAPFOBJS)
	$(FLD) -o $@ $(SAPFOBJS) $(LIBSF)

SIP=\
	$(STADENBIN)/sip.o\
	$(STADENBIN)/dias89.a\
	$(STADENBIN)/plot92.a\
	$(SEQLIB)

SIPOBJS=\
	$(SIP)\
	$(COMMONOBJS)\
	$(STADENBIN)/pl4010.o

XSIPOBJS=\
	$(SIP)\
	$(XCOMMONOBJS)\
	$(STADENBIN)/plX.o\
	$(STADENBIN)/textOutput.o 

$(O)/sip: $(SIPOBJS)
	$(FLD) -o $@ $(SIPOBJS) $(LIBSF)

$(O)/xsip: $(XSIPOBJS)
	$(FLD) -o $@ $(XSIPOBJS) $(XAW_LIB) $(MISC_LIB) $(LIBSF)


SIPL=\
	$(STADENBIN)/sipl.o\
	$(STADENBIN)/dias89.a\
	$(SEQLIB)

SIPLOBJS=\
	$(SIPL)\
	$(SCOMMONOBJS)

$(O)/sipl: $(SIPLOBJS)
	$(FLD) -o $@ $(SIPLOBJS) $(LIBSF)


SETHELP=\
	$(STADENBIN)/sethelp.o

SETHELPOBJS=\
	$(SETHELP)

$(O)/sethelp: $(SETHELPOBJS)
	$(CLD) -o $@ $(SETHELPOBJS) $(LIBSC)


SPLITP1=\
	$(STADENBIN)/splitp1.o
SPLITP1OBJS=\
	$(SPLITP1)\
	$(SCOMMONOBJS)

$(O)/splitp1: $(SPLITP1OBJS)
	$(FLD) -o $@ $(SPLITP1OBJS) $(LIBSF)


SPLITP2=\
	$(STADENBIN)/splitp2.o
SPLITP2OBJS=\
	$(SPLITP2)\
	$(SCOMMONOBJS)

$(O)/splitp2: $(SPLITP2OBJS)
	$(FLD) -o $@ $(SPLITP2OBJS) $(LIBSF)


SPLITP3=\
	$(STADENBIN)/splitp3.o
SPLITP3OBJS=\
	$(SPLITP3)\
	$(SCOMMONOBJS)

$(O)/splitp3: $(SPLITP3OBJS)
	$(FLD) -o $@ $(SPLITP3OBJS) $(LIBSF)


REP=\
	$(STADENBIN)/rep.o\
	$(STADENBIN)/dias89.a\
	$(STADENBIN)/subs89.a

REPOBJS=\
	$(REP)\
	$(COMMONOBJS)
 
$(O)/rep:	$(REPOBJS)
	$(FLD) -o $@ $(REPOBJS) $(LIBSF)


VEP=\
	$(STADENBIN)/vep.o\
	$(STADENBIN)/dias89.a\
	$(STADENBIN)/subs89.a


VEPOBJS=\
	$(VEP)\
	$(SCOMMONOBJS)


$(O)/vep: $(VEPOBJS)
	$(FLD) -o $@ $(VEPOBJS) $(LIBSF)

SPLITSEQ=\
	$(STADENBIN)/splitseq.o
 
SPLITSEQOBJS=\
	$(SPLITSEQ)\
	$(COMMONOBJS)
 
$(O)/splitseq: $(SPLITSEQOBJS)
	$(FLD) -o $@ $(SPLITSEQOBJS) $(LIBSF)

SPLITSEQF=\
	$(STADENBIN)/splitseqf.o
 

SPLITSEQFOBJS=\
	$(SPLITSEQF)\
	$(COMMONOBJS)
 
$(O)/splitseqf: $(SPLITSEQFOBJS)
	$(FLD) -o $@ $(SPLITSEQFOBJS) $(LIBSF)

#
# Special rule needed for Graph.o, seq*.o from the ted directory.
#
#$(O)/%.o: $(TEDSRC)/%.c
#	$(CC) $(CFLAGS) -c $< -o $@

DEPEND_OBJ =\
	$(GIPOBJS) $(LIPOBJS) $(MEPOBJS) $(XMEPOBJS) $(NIPOBJS) $(XNIPOBJS)\
	$(NIPLOBJS) $(NIPFOBJS) $(PIPOBJS) $(XPIPOBJS) $(PIPLOBJS) $(DAPOBJS)\
	$(XDAPOBJS) $(SAPOBJS) $(XSAPOBJS) $(SAPFOBJS) $(SIPOBJS) $(XSIPOBJS)\
	$(SIPLOBJS) $(SETHELPOBJS) $(SPLITP1OBJS) $(SPLITP2OBJS)\
	$(SPLITP3OBJS) $(REPOBJS) $(VEPOBJS)

#
#
#
nprogs: $(NPROGS)

xprogs: $(XPROGS)

lprogs: $(LPROGS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies
