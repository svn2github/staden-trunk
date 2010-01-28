/*  Last edited: Jun 15 11:43 2009 (badger) */
/*
 * Author:         James Bonfield, Feb 2007
 *                 Wellcome Trust Sanger Institute
 *
 * g_index: A program to index Tony Cox's .dat file format using a
 * recursive relative binning strategy. The data is stored in a g-library
 * format DB, for use with g_view for visualisation.
 *
 * The .dat alignment format consists of one line per sequence in order:
 *
 * seq_name ref_name start_pos end_pos direction(+1 or -1) clip_left \
 *    clip_right 0(?) 0(?) ? ? sequence confidence_values
 * 
 * Coordinates all start with 1 being the first base.
 */

#include <staden_config.h>
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <assert.h>

#include <staden_config.h>

#include "maqmap.h" /* lh3-rev */

#include "array.h"
#include "tg_gio.h"

#include "maq.h"
#include "ace.h"
#include "baf.h"
#include "tg_index_common.h"
#include "zfio.h"

#ifdef HAVE_SAMTOOLS
#include "sam.h"
#endif


void usage(void) {
    fprintf(stderr, "Usage: g_index [options] data_file ...\n");
    fprintf(stderr, "      -o output            Specify ouput filename (g_db)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -m                   Input is MAQ format\n");
    fprintf(stderr, "      -M                   Input is MAQ-long format\n");
    fprintf(stderr, "      -A                   Input is ACE format\n");
    fprintf(stderr, "      -B                   Input is BAF format\n");
#ifdef HAVE_SAMTOOLS
    fprintf(stderr, "      -b                   Input is BAM format\n");
    fprintf(stderr, "      -s                   Input is SAM format (with @SQ headers)\n");
#endif
    fprintf(stderr, "\n");
    fprintf(stderr, "      -u                   Also store unmapped reads (from SAM/BAM only)\n");
    fprintf(stderr, "      -x                   Also store auxillary records (from SAM/BAM only)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -p                   Link read-pairs together (default on)\n");
    fprintf(stderr, "      -P                   Do not link read-pairs together\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -a                   Append to existing db\n");
    fprintf(stderr, "      -n                   New contigs always (relevant if appending)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -g                   When appending to an existing db, assume the\n");
    fprintf(stderr, "                           alignment was performed against an ungapped copy\n");
    fprintf(stderr, "                           of the existing consensus. Add gaps back in to\n");
    fprintf(stderr, "                           reads and/or consensus as needed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -t                   Index sequence names\n");
    fprintf(stderr, "      -T                   Do not index sequence names (default on)\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "      -z value             Specify minimum bin size (default is '4k')\n"); 
    fprintf(stderr, "\n");
    fprintf(stderr, "      -f                   Fast mode: read-pair links are unidirectional\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -r nseq              Reserve space. Only necessary for exceptionally\n");
    fprintf(stderr, "                           large databases, eg n.seq > 100 million.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -d data_types        Only copy over certain data types. This is a comma\n"
	            "                           separated list containing one or more words\n"
	            "                           from: seq, qual, anno, name, all or none\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -c method            Specifies the compression method. This shold be\n");
    fprintf(stderr, "                           one of 'none', 'zlib' or 'lzma'.\n");
    fprintf(stderr, "                           Zlib is the default.\n");
}

//#include <malloc.h>

int main(int argc, char **argv) {
    tg_args a;
    GapIO *io;
    int opt, err = 0;
    char *cp;

    a.fmt            = 'a'; /* auto */
    a.out_fn         = "g_db";
    a.no_tree        = 1;
    a.pair_reads     = 1;
    a.append         = 0;
    a.merge_contigs  = -1;
    a.min_bin_size   = MIN_BIN_SIZE;
    a.fast_mode      = 0;
    a.reserved_seqs  = 0;
    a.data_type      = DATA_ALL;
    a.comp_mode      = COMP_MODE_ZLIB;
    a.repad          = 0;
    a.store_unmapped = 0;
    a.sam_aux        = 0;

    printf("\n\tg_index:\tShort Read Alignment Indexer, version 1.2.6\n");
    printf("\n\tAuthor: \tJames Bonfield (jkb@sanger.ac.uk)\n");
    printf("\t        \t2007-2010, Wellcome Trust Sanger Institute\n\n");

    //mallopt(M_TRIM_THRESHOLD, 100000);
    //mallopt(M_MMAP_MAX, 0);

    /* Arg parsing */
#ifdef HAVE_SAMTOOLS
    while ((opt = getopt(argc, argv, "aBsbtThAmMo:pPnz:fr:d:c:gux")) != -1) {
#else
    while ((opt = getopt(argc, argv, "aBstThAmMo:pPnz:fr:d:c:gux")) != -1) {
#endif
	switch(opt) {
	case 'g':
	    a.repad = 1;
	    break;

	case 'a':
	    a.append = 1;
	    if (a.merge_contigs == -1)
		a.merge_contigs = 1;
	    break;

	case 't':
	    a.no_tree = 0;
	    break;

	case 'T':
	    a.no_tree = 1;
	    break;

	case 'm':
	case 'M':
	case 'A':
	case 'B':
	case 's':
	case 'b':
	    a.fmt = opt;
	    break;

	case 'o':
	    a.out_fn = optarg;
	    break;

	case 'p':
	    a.pair_reads = 1;
	    break;

	case 'P':
	    a.pair_reads = 0;
	    break;

	case 'h':
	    usage();
	    return 0;

	case 'n':
	    a.merge_contigs = 0;
	    break;

	case 'z':
	    a.min_bin_size = strtol(optarg, &cp, 10);
	    if (*cp == 'k' || *cp == 'K') a.min_bin_size *= 1024;
	    if (*cp == 'm' || *cp == 'M') a.min_bin_size *= 1024*1024;
	    if (*cp == 'g' || *cp == 'G') a.min_bin_size *= 1024*1024*1024;
	    break;

	case 'f':
	    a.fast_mode = 1;
	    break;

	case 'r':
	    a.reserved_seqs = strtol(optarg, &cp, 10);
	    if (*cp == 'k' || *cp == 'K') a.reserved_seqs *= 1024;
	    if (*cp == 'm' || *cp == 'M') a.reserved_seqs *= 1024*1024;
	    if (*cp == 'g' || *cp == 'G') a.reserved_seqs *= 1024*1024*1024;
	    break;

	case 'd':
	    a.data_type = parse_data_type(optarg);
	    break;

	case 'c':
	    if (0 == strcmp(optarg, "none")) {
		a.comp_mode = COMP_MODE_NONE;
	    } else if (0 == strcmp(optarg, "zlib")) {
		a.comp_mode = COMP_MODE_ZLIB;
	    } else if (0 == strcmp(optarg, "lzma")) {
		a.comp_mode = COMP_MODE_LZMA;
	    } else {
		fprintf(stderr, "Unknown compression mode '%s'\n", optarg);
		usage();
		return 1;
	    }
	    break;

	case 'u':
	    a.store_unmapped = 1;
	    break;
	    
	case 'x':
	    a.sam_aux = 1;
	    break;
	    
	default:
	    if (opt == ':')
		fprintf(stderr, "Missing parameter\n");
	    else
		fprintf(stderr, "Unknown option '%c'\n", opt);
	    usage();
	    return 1;
	}
    }
    if (a.merge_contigs == -1)
	a.merge_contigs = 0;

    if (a.reserved_seqs)
	set_reserved_seqs(a.reserved_seqs);

    if (optind == argc) {
	usage();
	return 1;
    }

    /* Open the DB */
    io = gio_open(a.out_fn, 0, a.append ? 0 : 1);
    if (NULL == io) {
	perror("gio_open");
	return 1;
    }
    io->iface->setopt(io->dbh, OPT_COMP_MODE, a.comp_mode);

    if (a.no_tree) {
	io->db = cache_rw(io, io->db);
	io->db->seq_name_index = 0;
    }
    io->min_bin_size = a.min_bin_size;

    
    /* Open a temporary file for B+Tree indexing if needed */
    a.tmp = a.no_tree ? NULL : bttmp_file_open();



    /* File processing loop */
    while (optind < argc) {
	int fmt = a.fmt;

	/* Auto detect file type if appropriate */
	if (fmt == 'a')
	    fmt = tg_index_file_type(argv[optind]);

	switch (fmt) {
	case 'm':
	case 'M':
	    printf("Processing MAQ file %s\n", argv[optind]);
	    parse_maqmap(io, argv[optind++], &a);
	    break;

	case 'A':
	    printf("Processing ACE file %s\n", argv[optind]);
	    parse_ace(io, argv[optind++], &a);
	    break;

	case 'B':
	    printf("Processing BAF file %s\n", argv[optind]);
	    parse_baf(io, argv[optind++], &a);
	    break;

#ifdef HAVE_SAMTOOLS
	case 'b':
	    printf("Processing BAM file %s\n", argv[optind]);
	    parse_bam(io, argv[optind++], &a);
	    break;
	case 's':	
	    printf("Processing SAM file %s\n", argv[optind]);
	    parse_sam(io, argv[optind++], &a);
	    break;
#endif

	default:
	    fprintf(stderr, "Unknown file type for '%s' - skipping\n",
		    argv[optind++]);
	    err = 1;
	    break;
	}
    }

    /* Add to our sequence name B+Tree */
    if (a.tmp) {
	char *name;
	int rec;

	puts("Sorting sequence name index");
	bttmp_file_sort(a.tmp);

	puts("Building index");
	while (name = bttmp_file_get(a.tmp, &rec)) {
	    sequence_index_update(io, name, strlen(name), rec);
	}
	
	bttmp_file_close(a.tmp);
    }


    /* system("ps lx"); */

    gio_close(io);

    return err;
}
