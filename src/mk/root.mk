#
# Root level Makefile for Staden package.
# MACHINE define must be set before using this makefile. Set either using
# "gmake MACHINE=alpha" or as an environment variable.
#
# Valid values are "alpha", "al", "dec", "sgi", "solaris", and "sun".
#
# Invoking this makefile with JOB set will pass that as the item to make onto
# all sublevel makefiles. eg "gmake JOB=clean".
# 

# gap/osp-bits
# gap
DIRS =	Misc\
	io_lib\
	seq_utils\
	tk_utils\
	spin\
	text_utils\
	seqlib\
	abi\
	alf\
	ted\
	trev\
	staden\
	g\
	gap4/osp-bits\
	gap4\
	convert\
	frog\
	getMCH\
	indexseqlibs\
	vepe\
	clip\
	eba\
	init_exp\
	get_scf_field\
	lookup\
	trace_clip\
	sip4\
	vector_clip\
	trace_diff\
	screen_seq\
	nip4\
	qclip\
	pregap4\
	scripts

JOB = all

all:	$(DIRS)

#
# Sun's make can't handle redefining targets, so we need to rename the target
# in io_lib from (eg) clean to cleanx.
#
clean:
	cd io_lib; $(MAKE) $(FLAGS) cleanx
	$(MAKE) $(FLAGS) JOB=clean

spotless:
	cd io_lib; $(MAKE) $(FLAGS) spotlessx
	$(MAKE) $(FLAGS) JOB=spotless
	-rm lib/$(O)/so_locations

relink:
	-rm lib/$(O)/so_locations
	$(MAKE) $(FLAGS) JOB=relink

depend:
	$(MAKE) $(FLAGS) JOB=depend

install:
	$(MAKE) $(FLAGS) JOB=install

$(DIRS): FORCE
	if [ -d $@ ]; then cd $@; $(MAKE) $(FLAGS) $(MFLAGS) $(JOB); fi

FORCE:
# DO NOT DELETE THIS LINE -- make depend depends on it.
