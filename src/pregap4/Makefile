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
distsrc:
	-mkdir $(DIRNAME)
	-cp Makefile *.tcl *.xbm tclIndex pregap4 pregap4.bat \
	    widget_licence $(DIRNAME)
	-mkdir $(DIRNAME)/modules
	-mkdir $(DIRNAME)/templates
	-mkdir $(DIRNAME)/naming_schemes
	cp modules/*.p4m $(DIRNAME)/modules
	-cp templates/*.p4t $(DIRNAME)/templates
	cp naming_schemes/*.p4t $(DIRNAME)/naming_schemes

install:
	-cp pregap4 $(INSTALLBIN)/pregap4
	-mkdir $(INSTALLLIB)/pregap4
	cp *.tcl widget_licence tclIndex pregap4 $(INSTALLLIB)/pregap4
	-mkdir $(INSTALLLIB)/pregap4/modules
	-mkdir $(INSTALLLIB)/pregap4/templates
	-mkdir $(INSTALLLIB)/pregap4/naming_schemes
	cp modules/*.p4m $(INSTALLLIB)/pregap4/modules
	-cp templates/*.p4t $(INSTALLLIB)/pregap4/templates
	-cp naming_schemes/*.p4t $(INSTALLLIB)/pregap4/naming_schemes
