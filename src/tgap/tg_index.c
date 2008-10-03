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

#include "maqmap.h" /* lh3-rev */

#include "array.h"
#include "tg_gio.h"

#include "maq.h"
#include "ace.h"

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
    int start, end, dir, cleft, cright, ind, i;
    size_t len;

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
	    if (conf[i] == 0)
		conf[i] = 1; /* conf 0 implies N */
	    if (seq[i] != 'A' && seq[i] != 'C' &&
		seq[i] != 'G' && seq[i] != 'T')
		conf[i] = 0;
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

	s->name_len = strlen(name);
	s->name = s->data = (char *)malloc(s->name_len+2*len);
	s->seq = s->data + s->name_len;
	s->conf = s->seq + len;
        strcpy(s->name, name);
        memcpy(s->seq, seq, len);
        memcpy(s->conf, conf, len);
        return 0;
    }

    return -1;
}

/*
 * Parses the .dat file passed in, writing sequences to 'io' and adding to
 * bins/ranges pointed to by 'index'.
 *
 * Returns 0 on success
 *	  -1 on error
 */
int parse_file(GapIO *io, int max_size, char *dat_fn, int no_tree,
	       int pair_reads, int merge_contigs) {
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

    if (pair_reads) {
	pair = HacheTableCreate(32768, HASH_DYNAMIC_SIZE);
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
	HacheItem *hi, *other_hi;

	pos += strlen((char *)line);

	if (-1 == parse_line(&seq, (char *)line, contig)) {
	    fprintf(stderr, "Parser error on line: %s", line);
	    seq.len = 0;
	}

	if (strcmp(last_contig, contig)) {
	    /* printf("Processing contig %s\n", contig); */
	    ncontigs++;
	    strcpy(last_contig, contig);
	    if (!merge_contigs ||
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

	/* Add the range to a bin, and see which bin it was */
	bin = bin_add_range(io, &c, &r, &r_out);

	/* Find pair if appropriate */
	if (pair) {
	    int new = 0;
	    HacheData hd;
	    
	    hd.i = 0;
	    hi = HacheTableAdd(pair, seq.name, seq.name_len, hd, &new);

	    /* Pair existed already */
	    if (!new) {
		seq.other_end = hi->data.i;

		/* Link other end to 'us' too */
		other_hi = hi;
		hi = NULL;
	    } else {
		other_hi = NULL;
	    }
	} else {
	    hi = NULL;
	    other_hi = NULL;
	}

	/* Save sequence */
	seq.bin = bin->rec;
	recno = sequence_new_from(io, &seq);
	if (hi)
	    hi->data.i = recno;

	if (!no_tree)
	    sequence_index_update(io, seq.name, seq.name_len, recno);
	free(seq.data);

	/* Link bin back to sequence too before it gets flushed */
	r_out->rec = recno;

	/* Link other end of pair back to this recno if appropriate */
	if (other_hi) {
	    seq_t *other = (seq_t *)cache_search(io, GT_Seq, other_hi->data.i);
	    sequence_set_other_end(io, &other, recno);
	    HacheTableDel(pair, other_hi, 0);
	}


	nseqs++;
	//cache_flush(io);

	/*
	if (nseqs == 7)
	    exit(0);
	*/

	if ((nseqs & 0x3fff) == 0) {
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
	if ((nseqs & 0x3fff) == 0) {
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
void usage() {
    fprintf(stderr, "Usage: g_index [options] [-m] [-T] dat_file ...\n");
    fprintf(stderr, "      -o output		Specify ouput filename (g_db)\n");
    fprintf(stderr, "      -m			Input is MAQ format (off)\n");
    fprintf(stderr, "      -A			Input is ACE format (off)\n");
    fprintf(stderr, "      -p			Link read-pairs together\n");
    fprintf(stderr, "      -a			Append to existing db\n");
    fprintf(stderr, "      -n			New contigs always (relevant if appending)\n");
}

int main(int argc, char **argv) {
    unsigned int max_size = 63;
    GapIO *io;
    int opt, fmt = 'a' /*aln */;
    char *out_fn = "g_db";
    int no_tree=0, pair_reads=0, append=0, merge_contigs=1;

    printf("\n\tg_index:\tShort Read Alignment Indexer, version 0.11\n");
    printf("\n\tAuthor: \tJames Bonfield (jkb@sanger.ac.uk)\n");
    printf("\t        \t2007-2008, Wellcome Trust Sanger Institute\n\n");

    /* Arg parsing */
    while ((opt = getopt(argc, argv, "aThAmo:pn")) != -1) {
	switch(opt) {
	case 'a':
	    append = 1;
	    break;

	case 'T':
	    no_tree = 1;
	    break;

	case 'm':
	case 'A':
	    fmt = opt;
	    break;

	case 'o':
	    out_fn = optarg;
	    break;

	case 'p':
	    pair_reads = 1;
	    break;

	case 'h':
	    usage();
	    return 0;

	case 'n':
	    merge_contigs = 0;
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

    if (optind == argc) {
	usage();
	return 1;
    }

    /* Open the DB */
    io = gio_open(out_fn, 0, append ? 0 : 1);
    if (no_tree) {
	io->db = cache_rw(io, io->db);
	io->db->seq_name_index = 0;
    }

    /* File processing loop */
    while (optind < argc) {
	switch (fmt) {
	case 'm':
	    parse_maqmap(io, max_size, argv[optind++], no_tree, pair_reads,
			 merge_contigs);
	    break;

	case 'A':
	    parse_ace(io, max_size, argv[optind++], no_tree, pair_reads,
		      merge_contigs);
	    break;

	default:
	    parse_file(io, max_size, argv[optind++], no_tree, pair_reads,
		       merge_contigs);
	}
    }

    /* system("ps lx"); */

    gio_close(io);

    return 0;
}
