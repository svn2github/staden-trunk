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

distsrc: distsrc_dirs
	cp $(S)/Makefile $(S)/*.tcl $(S)/tclIndex $(S)/create_emboss_files \
	   $(S)/emboss_menu $(S)/spin_embossrc $(DIRNAME)

install:
	-mkdir $(INSTALLTCL)/spin_emboss
	cp $(S)/*.tcl $(S)/tclIndex $(INSTALLTCL)/spin_emboss
	cp $(S)/create_emboss_files $(INSTALLBIN)
	-mkdir $(INSTALLTCL)/spin_emboss/acdtcl
	-cp $(S)/acdtcl/* $(INSTALLTCL)/spin_emboss/acdtcl
	cp $(S)/emboss_menu $(S)/spin_embossrc $(INSTALLETC)
