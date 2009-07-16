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

install:
#	cp gzip_trace $(INSTALLBIN)
#	cp bzip_trace $(INSTALLBIN)
#	cp bzip2_trace $(INSTALLBIN)
	cp assemble4 $(INSTALLBIN)
	cp fasta-split $(INSTALLBIN)
	cp fasta2exp $(INSTALLBIN)
#	cp make_fasta $(INSTALLBIN)
#	cp staden2exp $(INSTALLBIN)
	cp splitseq_da $(INSTALLBIN)
	cp staden_help $(INSTALLBIN)
