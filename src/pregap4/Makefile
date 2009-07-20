SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

all:
	@echo Nothing to build for Pregap4.

clean:
	@echo Nothing to clean for Pregap4.


cleanprogs:


cleanlibs:


spotless: clean cleanprogs cleanlibs



.PHONY: distsrc install
distsrc: distsrc_dirs
	cp $(S)/Makefile $(S)/*.tcl $(S)/*.xbm $(S)/tclIndex $(S)/pregap4 \
	   $(S)/pregap4.bat $(S)/widget_licence $(S)/pregap4rc $(DIRNAME)
	-mkdir $(DIRNAME)/modules
	cp $(S)/modules/*.p4m $(DIRNAME)/modules
	-mkdir $(DIRNAME)/templates
	-cp $(S)/templates/*.p4t $(DIRNAME)/templates
	-mkdir $(DIRNAME)/naming_schemes
	cp $(S)/naming_schemes/*.p4t $(DIRNAME)/naming_schemes

S=$(SRCROOT)/pregap4
install:
	cp $(S)/pregap4 $(INSTALLBIN)/pregap4
	-mkdir $(INSTALLTCL)/pregap4
	cp $(S)/*.tcl $(S)/widget_licence $(S)/tclIndex $(INSTALLTCL)/pregap4
	-mkdir $(INSTALLTCL)/pregap4/modules
	-mkdir $(INSTALLTCL)/pregap4/templates
	-mkdir $(INSTALLTCL)/pregap4/naming_schemes
	cp $(S)/modules/*.p4m $(INSTALLTCL)/pregap4/modules
	-cp $(S)/templates/*.p4t $(INSTALLTCL)/pregap4/templates
	cp $(S)/naming_schemes/*.p4t $(INSTALLTCL)/pregap4/naming_schemes
	cp $(S)/pregap4rc $(INSTALLETC)