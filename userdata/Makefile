.PHONY: dist install

dist:
	-mkdir $(DIST)/userdata
	-cp *.seq *.aln *.ab1 $(DIST)/userdata
	-cp ecoli.* $(DIST)/userdata

distsrc:
	$(MAKE) $(SUBFLAGS) DIST=$(DISTSRC) dist
	cp Makefile $(DISTSRC)/userdata

install:

