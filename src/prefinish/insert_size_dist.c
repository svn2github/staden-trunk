#include <stdio.h>
#include <stdlib.h>

#include "IO.h"
#include "template.h"
#include "xalloc.h"
#include "qualIO.h"
#include "qual.h"

#define MAX_LENGTH 5000

/*
 * Stub functions to avoid link errors from IO.c.
 * This is the 'note execution' code, which we wish to avoid anyway.
 */
void execute_database_notes(GapIO *io, char *type) {}
void process_rawdata_note(GapIO *io) {}
void fix_notes(GapIO *io) {}

/*
 * Count the number of insert sizes of each length for templates claiming to
 * have size ranges between 'from' and 'to'.
 * lengths should be allocated as size 'lengths_size'.
 *
 * Returns the total count of templates processed.
 */
int insert_size_dist(GapIO *io, int from, int to,
		     int *lengths, int lengths_size) {
    template_c **tarr;
    int i, j, cnum, count = 0;

    for (cnum = 1; cnum <= NumContigs(io); cnum++) {
	char *cons;
	int *trans;
	/* Compute padded position to unpadded translation buffer */
	cons = malloc(io_clength(io, cnum));
	trans = calloc(io_clength(io, cnum), sizeof(int));
	calc_consensus(cnum, 1, io_clength(io, cnum), CON_SUM,
		       cons, NULL, NULL, NULL, 0.01, 0,
		       database_info, (void *)io);
	for (j = i = 0; i < io_clength(io, cnum); i++) {
	    trans[i] = j;
	    if (cons[i] != '*')
		j++;
	}
	xfree(cons);

	/* Compute template_c structures */
	tarr = init_template_checks(io, 1, &cnum, 1);
	check_all_templates(io, tarr);

	/*
	 * Scan through templates finding ones with both ends known and within
	 * the same contig.
	 */
	for (i = 1; i < Ntemplates(io); i++) {
	    item_t *clist;
	    int err, contig;
	    GTemplates t;
	    int len;

	    if (!tarr[i])
		continue;

	    if (tarr[i]->consistency)
		continue;

	    if (tarr[i]->flags &
		(TEMP_FLAG_GUESSED_START |
		 TEMP_FLAG_GUESSED_END |
		 TEMP_FLAG_SPANNING))
		continue;
	
	    err = 0;
	    contig = 0;

	    for (clist = head(tarr[i]->gel_cont); clist; clist = clist->next) {
		gel_cont_t *gc = (gel_cont_t *)clist->data;
		if (contig && gc->contig != contig) {
		    err = 1;
		    break;
		} else {
		    contig = gc->contig;
		}
	    }
	    if (err)
		continue;

	    template_read(io, i, t);
	    if (t.insert_length_min != from ||
		t.insert_length_max != to)
		continue;

	    if (tarr[i]->end >= 1 &&
		tarr[i]->end <= io_clength(io, cnum) &&
		tarr[i]->start >= 1 &&
		tarr[i]->start <= io_clength(io, cnum));
	    len = abs(trans[tarr[i]->end-1] - trans[tarr[i]->start-1] + 1);
	    if (len >= lengths_size) {
		fprintf(stderr,
			"Abnormally long template %d spotted, len %d\n",
			i, len);
		continue;
	    }
	    lengths[len]++;
	    
	    count++;
	}

	xfree(trans);
    }

    return count;
}

typedef struct {
    int min_size;
    int max_size;
    int freq;
} tsize_stats;

void dump_template_sizes(GapIO *io) {
    tsize_stats *stats = NULL;
    int nstats = 0;
    int i, j;
    GTemplates t;

    for (i = 1; i <= Ntemplates(io); i++) {
	template_read(io, i, t);
	for (j = 0; j < nstats; j++) {
	    if (stats[j].min_size == t.insert_length_min &&
		stats[j].max_size == t.insert_length_max)
		break;
	}

	if (j == nstats) {
	    nstats++;
	    stats = realloc(stats, nstats * sizeof(*stats));
	    stats[j].min_size = t.insert_length_min;
	    stats[j].max_size = t.insert_length_max;
	    stats[j].freq = 0;
	}

	stats[j].freq++;
    }

    for (i = 0; i < nstats; i++) {
	fprintf(stderr, "Size %5d - %5d count %d\n",
		stats[j].min_size, stats[j].max_size, stats[j].freq);
    }
}

static void usage(void) {
    fprintf(stderr, "insert_size_dist min_size max_size GAP4_D.B ...\n");
    exit(1);
}

int main(int argc, char **argv) {
    GapIO *io;
    char *project, *version;
    int status;
    int argn, i;
    int lengths[MAX_LENGTH];
    int count = 0;
    int min_size, max_size;

    if (argc < 4) {
	usage();
    }

    min_size = atoi(argv[1]);
    max_size = atoi(argv[2]);

    memset(lengths, 0, MAX_LENGTH * sizeof(int));

    for (argn = 3; argn < argc; argn++) {
	project = argv[argn];
	version = strchr(project, '.');
	if (!version) {
	    fprintf(stderr, "Malformed project name\n");
	    return 1;
	}
	*version++ = 0;

	fprintf(stderr, "Opening project '%s' version '%s'\n",
		project, version);
	io = open_db(project, version, &status, 0, 1);
	if (!io) {
	    fprintf(stderr, "Couldn't open database\n");
	    return 1;
	}

	dump_template_sizes(io);
	count += insert_size_dist(io, min_size, max_size, lengths, MAX_LENGTH);
	
	close_db(io);
    }

    fprintf(stderr, "Averaged over %d templates\n", count);
    for (i = 0; i < MAX_LENGTH; i++) {
	printf("%4d %3d\n", i, lengths[i]);
    }

    return 0;
}
