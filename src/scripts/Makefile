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
distsrc: distsrc_dirs
	cp $(S)/gzip_trace $(S)/bzip_trace $(S)/bzip2_trace $(S)/assemble4 \
	   $(S)/fasta-split $(S)/fasta2exp $(S)/make_fasta $(S)/staden2exp \
	   $(S)/splitseq_da $(S)/staden_help $(S)/Makefile $(DIRNAME)

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
