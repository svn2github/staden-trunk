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

#ifdef HAVE_SAMTOOLS
#include "sam.h"
#endif

/* .dat version format to handle */
#define DAT_VERSION 3

#define MAX_LINE_LEN 10240


/* ------------------------------------------------------------------------ */
/*
 * File format reading code. We should add in here solexa eland support and
 * similar.
 */

/*
 * Parses a line of text and stores the data in seq.
 *
 * Returns 0 for success
 *        -1 for failure
  */
int parse_line(seq_t *s, char *line, char *contig) {
    char name[MAX_LINE_LEN], seq[MAX_LINE_LEN], conf[MAX_LINE_LEN];
    int start, end, dir, cleft, cright, ind, i, qb;
    size_t len;
    char *cp;

#if DAT_VERSION==1
#    define TEMPLATE "%s %s %d %d %d %d %d %*d %*d %s %n"
#else
#    define TEMPLATE "%s %s %d %d %d %d %d %*d %*d %*d %*d %s %n"
#endif

    /* Prevent valgrind from complaining about portions of seq being unset */
    memset(s, 0, sizeof(*s));

    /* Decode the line */
    if (sscanf(line, TEMPLATE,
               name, contig,
	       &start, &end, &dir, &cleft, &cright, seq, &ind) == 8) {

        len = strlen(seq);
        line = &line[ind];

        for (i = 0; i < len; i++) {
            conf[i] = strtol(line, &line, 10);
        }

	/* possibly 4 confidence values per base instead */
	strtol(line, &cp, 10);
	if (line != cp) {
	    for (; i < len*4; i++) {
		conf[i] = strtol(line, &line, 10);
	    }

	    qb = 4;
	    s->format = SEQ_FORMAT_CNF4;
	} else {
	    qb = 1;
	    s->format = SEQ_FORMAT_MAQ; /* pack bytes */

	    /* In MAQ style padding conf 0 implies N */
	    for (i = 0; i < len; i++) {
		if (conf[i] == 0)
		    conf[i] = 1;
		if (seq[i] != 'A' && seq[i] != 'C' &&
		    seq[i] != 'G' && seq[i] != 'T')
		    conf[i] = 0;
	    }
	}

        s->pos = start;
        s->len = dir * len;
	s->flags = s->len < 0 ? SEQ_COMPLEMENTED : 0;
        s->left = cleft;
        s->right = cright;
	s->mapping_qual = 50; /* arbitrary default value */
	if (dir == 1)
	    s->pos -= s->left-1;
	else
	    s->pos -= -s->len - s->right;
#if DAT_VERSION >= 2
        s->left = cleft < cright ? cleft : cright;
        s->right = cright > cleft ? cright : cleft;
#endif

	/*
	printf("%s: s->pos=%d+%d st%d en%d\n", name, s->pos, s->len,
	       s->left, s->right);
	*/

	/* FIXME: use commands in tg_sequence.c instead */
	s->name_len = strlen(name);
	s->name = s->data = (char *)malloc(s->name_len+3+len+qb*len);
        strcpy(s->name, name);

	s->trace_name = s->name + s->name_len+1;
	*s->trace_name = 0;
	s->trace_name_len = 0;

	s->alignment = s->trace_name + s->trace_name_len+1;
	*s->alignment = 0;
	s->alignment_len = 0;

	s->seq = s->alignment + s->alignment_len+1;
        memcpy(s->seq, seq, len);

	s->conf = s->seq + len;
        memcpy(s->conf, conf,
	       (s->format == SEQ_FORMAT_CNF4 ? 4 : 1) * len);
        return 0;
    }

    return -1;
}

/*
 * Called from pair Hache when querying a name that is not present.
 */
static HacheData *pair_load(void *clientdata, char *key, int key_len,
			    HacheItem *hi) {
    GapIO *io = (GapIO *)clientdata;
    int rec;
    static HacheData hd;

    rec = sequence_index_query(io, key);
    printf("Query %s => %d\n", key, rec);

    if (!rec)
	return NULL;

    hd.i = rec;

    return &hd;
}

typedef struct {
    int rec;
    int bin;
    int idx;
    int crec;
} pair_loc_t;

/*
 * Parses the .dat file passed in, writing sequences to 'io' and adding to
 * bins/ranges pointed to by 'index'.
 *
 * Returns 0 on success
 *	  -1 on error
 */
int parse_file(GapIO *io, char *dat_fn, tg_args *a) {
    int nseqs = 0;
    struct stat sb;
    FILE *dat_fp;
    unsigned char line[MAX_LINE_LEN];
    off_t pos;
    contig_t *c;
    char last_contig[MAX_LINE_LEN];
    HacheTable *pair = NULL;
    int ncontigs = 0;
	
    *last_contig = 0;
    
    printf("Loading %s...\n", dat_fn);
    if (-1 == stat(dat_fn, &sb) ||
	NULL == (dat_fp = fopen(dat_fn, "r"))) {
	perror(dat_fn);
	return -1;
    }

    if (a->pair_reads) {
	pair = HacheTableCreate(1024, HASH_DYNAMIC_SIZE);
	pair->name = "pair";
	pair->load = pair_load;
	pair->del  = NULL;
    }

    /* Loop:
     * Read 1 sequence
     * Parse it & create a range
     * Save sequence
     * Insert to index
     */
    pos = 0;
    while (fgets((char *)line, MAX_LINE_LEN, dat_fp)) {
	seq_t seq;
	range_t r, *r_out;
	int recno;
	bin_index_t *bin;
	char contig[MAX_LINE_LEN];
	HacheItem *hi;
	pair_loc_t *pl = NULL;

	pos += strlen((char *)line);

	if (-1 == parse_line(&seq, (char *)line, contig)) {
	    fprintf(stderr, "Parser error on line: %s", line);
	    seq.len = 0;
	}

	if (strcmp(last_contig, contig)) {
	    /* printf("Processing contig %s\n", contig); */
	    ncontigs++;
	    strcpy(last_contig, contig);
	    if (!a->merge_contigs ||
		NULL == (c = find_contig_by_name(io, contig))) {
		c = contig_new(io, contig);
		contig_index_update(io, contig, strlen(contig), c->rec);
	    }
	    cache_incr(io, c);
	}

	/* Create range */
	r.start = seq.pos;
	r.end   = seq.pos + (seq.len > 0 ? seq.len : -seq.len) - 1;
	r.rec   = 0;
	r.pair_rec = 0;
	r.mqual = seq.mapping_qual;
	r.flags = GRANGE_FLAG_TYPE_SINGLE;
	/* Guess work here. For now all <--- are rev, all ---> are fwd */
	r.flags|= seq.len > 0
	    ? GRANGE_FLAG_END_FWD
	    : GRANGE_FLAG_END_REV;
	if (seq.len < 0)
	    r.flags |= GRANGE_FLAG_COMP1;

	/* Add the range to a bin, and see which bin it was */
	bin = bin_add_range(io, &c, &r, &r_out, NULL);

	/* Save sequence */
	seq.bin = bin->rec;
	seq.bin_index = r_out - ArrayBase(range_t, bin->rng);
	recno = sequence_new_from(io, &seq);

	/* Find pair if appropriate */
	if (pair) {
	    int new = 0;
	    HacheData hd;
	    
	    /* Add data for this end */
	    pl = (pair_loc_t *)malloc(sizeof(*pl));
	    pl->rec  = recno;
	    pl->bin  = bin->rec;
	    pl->crec = c->rec;
	    pl->idx  = seq.bin_index;
	    hd.p = pl;

	    hi = HacheTableAdd(pair, seq.name, seq.name_len, hd, &new);

	    /* Pair existed already */
	    if (!new) {
		pair_loc_t *po = (pair_loc_t *)hi->data.p;
		bin_index_t *bo;
		range_t *ro;

		/* We found one so update r_out now, before flush */
		r_out->flags &= ~GRANGE_FLAG_TYPE_MASK;
		r_out->flags |=  GRANGE_FLAG_TYPE_PAIRED;
		r_out->pair_rec = po->rec;

		if (!a->fast_mode) {
		    /* Link other end to 'us' too */
		    bo = (bin_index_t *)cache_search(io, GT_Bin, po->bin);
		    bo = cache_rw(io, bo);
		    bo->flags |= BIN_RANGE_UPDATED;
		    ro = arrp(range_t, bo->rng, po->idx);
		    ro->flags &= ~GRANGE_FLAG_TYPE_MASK;
		    ro->flags |=  GRANGE_FLAG_TYPE_PAIRED;
		    ro->pair_rec = pl->rec;
		}

		/* And, making an assumption, remove from hache */
		HacheTableDel(pair, hi, 1);
		free(pl);
	    }
	}

	if (!a->no_tree)
	    sequence_index_update(io, seq.name, seq.name_len, recno);
	free(seq.data);

	/* Link bin back to sequence too before it gets flushed */
	r_out->rec = recno;

	nseqs++;
	//cache_flush(io);

	/*
	if (nseqs == 7)
	    exit(0);
	*/

	if ((nseqs & 0xffff) == 0) {
	    static int perc = 0;
	    if (perc < 100.0 * pos / sb.st_size) {
		perc = 100.0 * pos / sb.st_size;
		printf("\r%d%%", perc);
		//HacheTableStats(io->cache, stdout);
		//HacheTableStats(((GDB **)io->dbh)[0]->gfile->idx_hash, stdout);
		{
		    static struct timeval last, curr;
		    static int first = 1;
		    static int last_contig = 0;
		    long delta;

		    gettimeofday(&curr, NULL);
		    if (first) {
			last = curr;
			first = 0;
		    }

		    delta = (curr.tv_sec - last.tv_sec) * 1000000
			+ (curr.tv_usec - last.tv_usec);
		    printf(" - %g sec %d contigs\n", delta/1000000.0,
			   ncontigs - last_contig);
		    last = curr;
		    last_contig = ncontigs;
		}
		fflush(stdout);
	    }
	}

	/*
	 * If we rely on the input being sorted then we know that no
	 * subsequent sequences can be leftwards of where we are now,
	 * but they may possibly be in a lower bin (in our left child) if
	 * we just inserted a very long sequence and the next sequence is
	 * very short.
	 */
	if ((nseqs & 0xffff) == 0) {
	    cache_flush(io);
	}
    }

    cache_flush(io);
    fclose(dat_fp);

    if (pair)
	HacheTableDestroy(pair, 0);

    return 0;
}

/* ------------------------------------------------------------------------ */
void usage(void) {
    fprintf(stderr, "Usage: g_index [options] [-m] [-T] dat_file ...\n");
    fprintf(stderr, "      -o output            Specify ouput filename (g_db)\n");
    fprintf(stderr, "      -m                   Input is MAQ format\n");
    fprintf(stderr, "      -M                   Input is MAQ-long format\n");
    fprintf(stderr, "      -A                   Input is ACE format\n");
    fprintf(stderr, "      -B                   Input is BAF format\n");
#ifdef HAVE_SAMTOOLS
    fprintf(stderr, "      -b                   Input is BAM format\n");
#endif
    fprintf(stderr, "      -p                   Link read-pairs together (default on)\n");
    fprintf(stderr, "      -P                   Do not link read-pairs together\n");
    fprintf(stderr, "      -a                   Append to existing db\n");
    fprintf(stderr, "      -n                   New contigs always (relevant if appending)\n");
    fprintf(stderr, "      -t                   Index sequence names\n");
    fprintf(stderr, "      -T                   Do not index sequence names (default on)\n");
    fprintf(stderr, "      -z value             Specify minimum bin size (default is '4k')\n"); 
    fprintf(stderr, "      -f                   Fast mode: read-pair links are unidirectional\n");
}

#include <malloc.h>

int main(int argc, char **argv) {
    tg_args a;
    GapIO *io;
    int opt;
    char *cp;

    a.fmt           = 'a'; /* aln */
    a.out_fn        = "g_db";
    a.no_tree       = 1;
    a.pair_reads    = 1;
    a.append        = 0;
    a.merge_contigs = -1;
    a.min_bin_size  = MIN_BIN_SIZE;
    a.fast_mode     = 0;

    printf("\n\tg_index:\tShort Read Alignment Indexer, version 1.2.3\n");
    printf("\n\tAuthor: \tJames Bonfield (jkb@sanger.ac.uk)\n");
    printf("\t        \t2007-2009, Wellcome Trust Sanger Institute\n\n");

    //mallopt(M_MMAP_MAX, 0);

    /* Arg parsing */
#ifdef HAVE_SAMTOOLS
    while ((opt = getopt(argc, argv, "aBsbtThAmMo:pPnz:f")) != -1) {
#else
    while ((opt = getopt(argc, argv, "aBstThAmMo:pPnz:f")) != -1) {
#endif
	switch(opt) {
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

    if (optind == argc) {
	usage();
	return 1;
    }

    /* Open the DB */
    io = gio_open(a.out_fn, 0, a.append ? 0 : 1);
    if (a.no_tree) {
	io->db = cache_rw(io, io->db);
	io->db->seq_name_index = 0;
    }
    io->min_bin_size = a.min_bin_size;

    
    /* Open a temporary file for B+Tree indexing if needed */
    a.tmp = a.no_tree ? NULL : bttmp_file_open();


    /* File processing loop */
    while (optind < argc) {
	switch (a.fmt) {
	case 'm':
	case 'M':
	    parse_maqmap(io, argv[optind++], &a);
	    break;

	case 'A':
	    parse_ace(io, argv[optind++], &a);
	    break;

	case 'B':
	    parse_baf(io, argv[optind++], &a);
	    break;

#ifdef HAVE_SAMTOOLS
	case 'b':
	    parse_bam(io, argv[optind++], &a);
	    break;
#endif

	default:
	    parse_file(io, argv[optind++], &a);
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

    return 0;
}
