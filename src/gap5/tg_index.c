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
#include "caf.h"
#include "tg_index_common.h"
#include "zfio.h"
#include "fasta.h"

#include "sam_index.h"
#include "bam.h"
#include "afg.h"


void usage(void) {
    fprintf(stderr, "Usage: g_index [options] data_file ...\n");
    fprintf(stderr, "      -o output            Specify ouput filename (g_db)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -m                   Input is MAQ format\n");
    fprintf(stderr, "      -M                   Input is MAQ-long format\n");
    fprintf(stderr, "      -A                   Input is ACE format\n");
    fprintf(stderr, "      -B                   Input is BAF format\n");
    fprintf(stderr, "      -C                   Input is CAF format\n");
    fprintf(stderr, "      -f                   Input is FASTA format\n");
    fprintf(stderr, "      -F                   Input is FASTQ format\n");
    fprintf(stderr, "      -b                   Input is BAM format\n");
    fprintf(stderr, "      -s                   Input is SAM format (with @SQ headers)\n");
/*    fprintf(stderr, "      -V                   Input is AFG (Velvet)format\n"); */
    fprintf(stderr, "\n");
    fprintf(stderr, "      -u                   Also store unmapped reads           (SAM/BAM only)\n");
    fprintf(stderr, "      -x                   Also store auxillary records        (SAM/BAM only)\n");
    fprintf(stderr, "      -r                   Store reference-position data (on)  (SAM/BAM only)\n");
    fprintf(stderr, "      -R                   Don't store reference-position data (SAM/BAM only)\n");
    fprintf(stderr, "      -D                   Do not remove duplicates (SAM/BAM only)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -p                   Link read-pairs together (default on)\n");
    fprintf(stderr, "      -P                   Do not link read-pairs together\n");
    fprintf(stderr, "      -L                   Do not link pairs spanning bam files or with '-a'.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -q value             Number of reads to queue in memory while waiting\n");                   
    fprintf(stderr, "                           for pairing.  Use to reduce memory requirements\n");
    fprintf(stderr, "                           for assemblies with lots of single reads at the\n");
    fprintf(stderr, "                           expense of running time.  0 for all in memory,\n");
    fprintf(stderr, "                           suggest 1000000 if used (default 0)\n");           
    fprintf(stderr, "\n");
    fprintf(stderr, "      -a                   Append to existing db\n");
    fprintf(stderr, "      -n                   New contigs always (relevant if appending)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -g                   When appending to an existing db, assume the\n");
    fprintf(stderr, "                           alignment was performed against an ungapped copy\n");
    fprintf(stderr, "                           of the existing consensus. Add gaps back in to\n");
    fprintf(stderr, "                           reads and/or consensus as needed.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -t                   Index sequence names (default)\n");
    fprintf(stderr, "      -T                   Do not index sequence names\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "      -z value             Specify minimum bin size (default is '4k')\n"); 
    fprintf(stderr, "\n");
    fprintf(stderr, "      -f                   Fast mode: read-pair links are unidirectional\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "                           large databases, eg n.seq > 100 million.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -d data_types        Only copy over certain data types. This is a comma\n"
	            "                           separated list containing one or more words\n"
	            "                           from: seq, qual, anno, name, all or none\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -c method            Specifies the compression method. This shold be\n");
    fprintf(stderr, "                           one of 'none', 'zlib' or 'lzma'.\n");
    fprintf(stderr, "                           Zlib is the default.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "      -[1-9]               Use a fixed compression level from 1 to 9\n");
    fprintf(stderr, "      -v version_num       Request a specific database formation version\n");
}

//#include <malloc.h>

int main(int argc, char **argv) {
    tg_args a;
    GapIO *io;
    int opt, err = 0;
    char *cp;

    a.fmt            = 'a'; /* auto */
    a.out_fn         = "";
    a.no_tree        = 0;
    a.pair_reads     = 1;
    a.append         = 0;
    a.merge_contigs  = -1;
    a.min_bin_size   = MIN_BIN_SIZE;
    a.fast_mode      = 0;
    a.data_type      = DATA_ALL;
    a.comp_mode      = COMP_MODE_ZLIB;
    a.repad          = 0;
    a.store_unmapped = 0;
    a.sam_aux        = 0;
    a.pair_queue     = 0;
    a.store_refpos   = 1;
    a.remove_dups    = 1;
    a.version        = DB_VERSION;
    a.link_pairs     = 1;

    printf("\n\tg_index:\tShort Read Alignment Indexer, version 1.2.13"SVN_VERS"\n");
    printf("\n\tAuthor: \tJames Bonfield (jkb@sanger.ac.uk)\n");
    printf("\t        \t2007-2011, Wellcome Trust Sanger Institute\n\n");

    //mallopt(M_TRIM_THRESHOLD, 100000);
    //mallopt(M_MMAP_MAX, 0);

    /* Arg parsing */
    while ((opt = getopt(argc, argv, "aBCsVbtThAmMo:pPq:nz:fd:c:"
			 "gux123456789rRDv:L")) != -1) {
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
	case 'C':
	case 'V':
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
	    
	case 'q':
	    a.pair_queue = strtol(optarg, &cp, 10);
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

	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
	    set_tg_compression_level(opt-'0');
	    break;

	case 'r':
	    a.store_refpos = 1;
	    break;

	case 'R':
	    a.store_refpos = 0;
	    break;

	case 'D':
	    a.remove_dups = 0;
	    break;

	case 'v':
	    a.version = atoi(optarg);
	    break;
	    
	case 'L':
	    a.link_pairs = 0;
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
    if (!a.pair_reads)
	a.link_pairs = 0;
    if (a.merge_contigs == -1)
	a.merge_contigs = 0;

    if (optind == argc) {
	usage();
	return 1;
    }

    /* Pick a default name for the DB */
    if (!*a.out_fn) {
	char *cp;

	a.out_fn = malloc(strlen(argv[optind])+3);
	strcpy(a.out_fn, argv[optind]);
	if ((cp = strrchr(a.out_fn, '.')))
	    *cp = 0;
	strcat(a.out_fn, ".0");

	printf("Selecting output database filename %s\n", a.out_fn);
    }

    /* Open the DB */
    if (a.version != DB_VERSION)
	gio_set_db_version(a.version);

    io = gio_open(a.out_fn, 0, a.append ? 0 : 1);
    if (NULL == io) {
	perror("gio_open");
	return 1;
    }
    io->iface->setopt(io->dbh, OPT_COMP_MODE, a.comp_mode);

    if (a.no_tree || (a.data_type & DATA_NAME) == 0) {
	io->db = cache_rw(io, io->db);
	io->db->seq_name_index = 0;
    }
    io->min_bin_size = a.min_bin_size;

    
    /* File processing loop */
    while (optind < argc) {
	int fmt = a.fmt;

	/* Open a temporary file for B+Tree indexing if needed */
	a.tmp = a.no_tree ? NULL : bttmp_file_open();

	/* Auto detect file type if appropriate */
	if (fmt == 'a')
	    fmt = tg_index_file_type(argv[optind]);

	switch (fmt) {
	case 'm':
	case 'M':
	    printf("Processing MAQ file %s\n", argv[optind]);
	    err = parse_maqmap(io, argv[optind++], &a) ? 1 : 0;
	    break;

	case 'A':
	    printf("Processing ACE file %s\n", argv[optind]);
	    err = parse_ace(io, argv[optind++], &a) ? 1 : 0;
	    break;

	case 'B':
	    printf("Processing BAF file %s\n", argv[optind]);
	    err = parse_baf(io, argv[optind++], &a) ? 1 : 0;
	    break;

	case 'C':
	    printf("Processing CAF file %s\n", argv[optind]);
	    parse_caf(io, argv[optind++], &a);
	    break;

	case 'b':
	    printf("Processing BAM file %s\n", argv[optind]);
	    err = parse_bam(io, argv[optind++], &a) ? 1 : 0;
	    break;

	case 's':	
	    printf("Processing SAM file %s\n", argv[optind]);
	    err = parse_sam(io, argv[optind++], &a) ? 1 : 0;
	    break;

	case 'F':
	    printf("Processing FASTA file %s\n", argv[optind]);
	    err = parse_fasta_or_fastq(io, argv[optind++], &a, 'a') ? 1 : 0;
	    break;

	case 'Q':
	    printf("Processing FASTQ file %s\n", argv[optind]);
	    err = parse_fasta_or_fastq(io, argv[optind++], &a, 'q') ? 1 : 0;
	    break;

	case 'V':
	    printf("Processing AFG file %s\n", argv[optind]);
	    err = parse_afg(io, argv[optind++], &a) ? 1 : 0;
	    break;

	default:
	    fprintf(stderr, "Unknown file type for '%s' - skipping\n",
		    argv[optind++]);
	    err = 1;
	    break;
	}

	/* Force final update of cached bin nseq */
	bin_add_range(io, NULL, NULL, NULL, NULL, -1);
	cache_flush(io);

	if (err) {
	    fprintf(stderr, "\nERROR: Failed to parse input data file - exiting\n");
	    exit(1);
	}

	/* Add to our sequence name B+Tree */
	if (a.tmp) {
	    char *name;
	    tg_rec rec;
	    int cnt = 0;

	    puts("Sorting sequence name index");
	    bttmp_file_sort(a.tmp);

	    puts("Building index: one dot per 10k reads");
	    while (name = bttmp_file_get(a.tmp, &rec)) {
		sequence_index_update(io, name, strlen(name), rec);
		if (++cnt == 10000) {
		    putchar('.'); fflush(stdout);
		    cnt = 0;
		    cache_flush(io);
		}
	    }
	    putchar('\n');
	
	    bttmp_file_close(a.tmp);
	}
    }

    /* system("ps lx"); */

    gio_close(io);

    return err;
}
