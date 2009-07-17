SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

all:
	@echo Nothing to build for spin_emboss.

clean:
	@echo Nothing to clean for spin_emboss.

cleanprogs:
cleanlibs:
spotless: clean cleanprogs cleanlibs

.PHONY: distsrc install acdtcl

ACDDIR=parser/acd
acdtcl:
	stash create_emboss_files.tcl

distsrc:
	-mkdir $(DIRNAME)
	-cp Makefile *.tcl tclIndex create_emboss_files $(DIRNAME)
	-cp -R acdtcl $(DIRNAME)

S=$(SRCROOT)/spin_emboss
install:
	-mkdir $(INSTALLTCL)/spin_emboss
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/spin_emboss
	cp $(S)/create_emboss_files $(INSTALLBIN)
	-mkdir $(INSTALLTCL)/spin_emboss/acdtcl
	-cp $(S)/acdtcl/* $(INSTALLTCL)/spin_emboss/acdtcl
	cp $(S)/emboss_menu $(S)/spin_embossrc $(INSTALLETC)
