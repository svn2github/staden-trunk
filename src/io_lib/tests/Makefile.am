## Makefile.am -- Process this file with automake to produce Makefile.in

EXTRA_DIST              = (TESTS)
MAINTAINERCLEANFILES    = Makefile.in

TESTS_ENVIRONMENT       = \
			top_builddir=$(top_builddir) \
			top_srcdir=$(top_srcdir)
TESTS                   = \
                        illumina2srf.test \
			srf2illumina.test \
                        srf2fasta.test \
                        srf2fastq.test \
			srf_info.test \
			srf_filter.test \
			srf_index.test
