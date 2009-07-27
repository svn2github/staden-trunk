SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

all:
	@echo Nothing to build for data.

clean:
	@echo Nothing to clean for data.


cleanprogs:


cleanlibs:


spotless: clean cleanprogs cleanlibs


.PHONY: distsrc cleanprogs cleanlibs

distsrc: distsrc_dirs
	cp $(S)/RENZYM.* $(S)/renzyme_bairoch $(S)/m13mp18_primers \
	   $(S)/vector_primer $(S)/pam250 $(S)/Makefile $(DIRNAME)
	-mkdir $(DIRNAME)/alus
	cp $(S)/alus/* $(DIRNAME)/alus
	-mkdir $(DIRNAME)/blast_data
	cp $(S)/blast_data/* $(DIRNAME)/blast_data
	-mkdir $(DIRNAME)/gcodes
	cp $(S)/gcodes/* $(DIRNAME)/gcodes
	-mkdir $(DIRNAME)/vectors
	cp $(S)/vectors/* $(DIRNAME)/vectors

install:
	-mkdir $(INSTALLDATA)/alus
	cp $(S)/alus/* $(INSTALLDATA)/alus
	-mkdir $(INSTALLDATA)/blast_data
	cp $(S)/blast_data/* $(INSTALLDATA)/blast_data
	-mkdir $(INSTALLDATA)/gcodes
	cp $(S)/gcodes/* $(INSTALLDATA)/gcodes
	-mkdir $(INSTALLDATA)/vectors
	cp $(S)/vectors/* $(INSTALLDATA)/vectors
	cp $(S)/RENZYM.* $(INSTALLDATA)
	cp $(S)/renzyme_bairoch $(INSTALLDATA)
	cp $(S)/m13mp18_primers $(INSTALLDATA)
	cp $(S)/vector_primer $(INSTALLDATA)
	cp $(S)/pam250 $(INSTALLDATA)

