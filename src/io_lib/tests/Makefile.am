## Makefile.am -- Process this file with automake to produce Makefile.in

EXTRA_DIST              = $(TESTS) data
MAINTAINERCLEANFILES    = Makefile.in

test_outdir              = test.out

TESTS_ENVIRONMENT       = \
			top_builddir=$(top_builddir) \
			top_srcdir=$(top_srcdir) \
			outdir=$(test_outdir)

TESTS                   = \
                        illumina2srf.test \
			srf2illumina.test \
                        srf2fasta.test \
                        srf2fastq.test \
			srf_info.test \
			srf_filter.test \
			srf_index.test

distclean-local:
	-rm -rf $(test_outdir)

