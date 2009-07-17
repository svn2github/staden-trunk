SRCROOT=$(SRC)/..
include $(SRCROOT)/global.mk
include ../system.mk

all:
	@echo Nothing to build for scripts.

clean:
	@echo Nothing to clean for scripts.


cleanprogs:


cleanlibs:


spotless: clean cleanprogs cleanlibs


.PHONY: distsrc cleanprogs cleanlibs
distsrc:
	-mkdir $(DIRNAME)
	cp gzip_trace bzip_trace bzip2_trace assemble4 fasta-split \
	   fasta2exp make_fasta staden2exp splitseq_da staden_help \
	   Makefile $(DIRNAME)
	-rmdir *-binaries

# Not all scripts are particularly useful, so don't bother to install these.

S=$(SRCROOT)/scripts
install:
#	cp $(S)/gzip_trace $(INSTALLSCRIPTS)
#	cp $(S)/bzip_trace $(INSTALLSCRIPTS)
#	cp $(S)/bzip2_trace $(INSTALLSCRIPTS)
	cp $(S)/assemble4 $(INSTALLSCRIPTS)
	cp $(S)/fasta-split $(INSTALLSCRIPTS)
	cp $(S)/fasta2exp $(INSTALLSCRIPTS)
#	cp $(S)/make_fasta $(INSTALLSCRIPTS)
#	cp $(S)/staden2exp $(INSTALLSCRIPTS)
	cp $(S)/splitseq_da $(INSTALLSCRIPTS)
	cp $(S)/staden_help $(INSTALLSCRIPTS)
