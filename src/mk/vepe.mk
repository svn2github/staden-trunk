#
# Makefile for vepe and repe
#

PROGS= $(O)/vepe $(O)/repe

SRCROOT=..
include $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

VEPEBIN=$(O)
REPEBIN=$(O)

#
# Common objects, needed by every program
#
COMMONOBJS=\
	$(STADENBIN)/seeme-$(MACHINE).o\
	$(STADENBIN)/Cme-$(MACHINE).o\
	$(STADENBIN)/nxspec.o\
	$(STADENBIN)/userface.o\
	$(STADENBIN)/userfacecom.o\
	$(STADENBIN)/helpnmenu.o\
	$(STADENBIN)/nxhelpmenu.o\
	$(STADENBIN)/FtoC.o\
	$(STADENBIN)/subs89.a

#DEPS=\
#	$(IOLIBBIN)/libexp.a \
#	$(MISCBIN)/libmisc.a
DEPS=\
	$(EXP_DEP) \
	$(MISC_DEP)

VEPE=\
	$(VEPEBIN)/vepe.o\
	$(STADENBIN)/dias89.a\
	$(VEPEBIN)/expio.o


VEPEOBJS=\
	$(VEPE)\
	$(COMMONOBJS)

VRLIBS=\
	$(EXP_LIB) \
	$(MISC_LIB)

$(O)/vepe: $(VEPEOBJS)
	$(FLD) -o $@ $(VEPEOBJS) $(VRLIBS) $(LIBSF) $(DEPS)


REPE=\
	$(VEPEBIN)/repe.o\
	$(STADENBIN)/dias89.a\
	$(VEPEBIN)/expio.o


REPEOBJS=\
	$(REPE)\
	$(COMMONOBJS)


$(O)/repe: $(REPEOBJS) $(DEPS)
	$(FLD) -o $@ $(REPEOBJS) $(VRLIBS) $(LIBSF)



DEPEND_OBJ = $(VEPEOBJS) $(REPEOBJS)

install:
	cp $(PROGS) $(INSTALLBIN)

include dependencies

