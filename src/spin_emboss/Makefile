SRCROOT=..
include	 $(SRCROOT)/mk/global.mk
include $(SRCROOT)/mk/$(MACHINE).mk

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

install:
	-mkdir $(INSTALLLIB)/spin_emboss
	cp *.tcl tclIndex $(INSTALLLIB)/spin_emboss
	cp create_emboss_files $(INSTALLBIN)
	-mkdir $(INSTALLLIB)/spin_emboss/acdtcl
	cp acdtcl/* $(INSTALLLIB)/spin_emboss/acdtcl
