#include <string.h>

#include "probe.h"
#include "gap_globals.h"
#include "assemble_direct.h"
#include "dna_utils.h"
#include "text_output.h"
#include "sequence_formats.h"
#include "misc.h"
#include "fort.h"
#include "primlib.h"

#define MAX_OLIGO_LEN 100
#define MAX_SCAN_REGION 1024
#define MAX_MATCHES 10
#define MAX_VECTOR_SEQ 100000

#define FORWARD 0
#define REVERSE 1

typedef struct {
    int start;
    int length;
    int distance;
    float primer_score;
    float best_match;
    int matches;
    char sequence[MAX_OLIGO_LEN+1];
    float tm;
} oligo_t;

/*
 * Scan through an array of oligo_t structures checking for uniqueness of
 * each oligo by checking the vector sequences. The matches
 * field of the oligo_t structure is updated using this routine.
 *
 * Returns  0 for success,
 *         -1 for failure
 */
static int find_uniques_vector(GapIO *io, int contig, float match_perc,
			       char *vectors, oligo_t *ol, int num_ol)
{
    FILE *fp;
    char *file;
    /* char seq[MAX_VECTOR_SEQ+1]; */
    char *seq = NULL;
    int i, nmatches;
    char oligo[MAX_OLIGO_LEN+1];
    int max_olen, vector_len, min_match;
    int match_arr[MAX_MATCHES], score_arr[MAX_MATCHES];
    float bi, best;

    if (0 == num_ol)
	return 0;

    if (NULL == (fp = open_fofn(vectors)))
	return -1;

    /*
     * Find the longest oligo - we need this for our later fudge on making
     * the vector a circle.
     */
    for (max_olen = i = 0; i < num_ol; i++) {
	if (ol[i].length > max_olen)
	    max_olen = ol[i].length;
    }


    /* Read the vector sequences */
    while (file = read_fofn(fp)) {
	if (0 != get_seq(&seq, MAX_VECTOR_SEQ-max_olen, &vector_len,
			 file, NULL)) {
	    verror(ERR_WARN, "find_probes", "Couldn't load file '%s'", file);
	    continue;
	}

	if (1 != get_seq_type(seq, vector_len)) {
	    verror(ERR_WARN, "find_probes", "File '%s' read is not DNA", file);
	}

	/* The vector is a circle, so duplicate the ends */
	seq = xrealloc(seq, vector_len + max_olen + 1);
	strncpy(&seq[vector_len], seq, max_olen);
	vector_len += max_olen;


	for (i = 0; i < num_ol; i++) {
	    UpdateTextOutput();
	    min_match = match_perc * ol[i].length;
	    best = 0;

	    /* Scan in the forward sense */
	    strcpy(oligo, ol[i].sequence);
	    nmatches = inexact_match(seq, vector_len, oligo, ol[i].length,
				     min_match, match_arr, score_arr,
				     MAX_MATCHES);
	    if (nmatches < 0) {
		ol[i].matches += MAX_MATCHES+1;
	    } else {
		ol[i].matches += nmatches;
	    }
	    bi = best_inexact_match(seq, vector_len, oligo, ol[i].length, NULL)
		/ (float)ol[i].length;
	    if (bi > best)
		best = bi;
	    ol[i].best_match = best;

	    /* Find matches in the reverse sense */
	    complement_seq(oligo, ol[i].length);
	    nmatches = inexact_match(seq, vector_len, oligo, ol[i].length,
				     min_match, match_arr, score_arr,
				     MAX_MATCHES);
	    if (nmatches < 0) {
		ol[i].matches += MAX_MATCHES+1;
	    } else {
		ol[i].matches += nmatches;
	    }
	    bi = best_inexact_match(seq, vector_len, oligo, ol[i].length, NULL)
		/ (float)ol[i].length;
	    if (bi > best)
		best = bi;
	    ol[i].best_match = best;
	}

	xfree(seq);
    }

    close_fofn(fp);

    return 0;
}

/*
 * Scan through an array of oligo_t structures checking for uniqueness of
 * each oligo by checking the existing consensus sequence. The matches
 * field of the oligo_t structure is updated using this routine.
 *
 * Returns  0 for success,
 *         -1 for failure
 */
static int find_uniques_con(GapIO *io, int contig, float match_perc,
			    consen_info *ci, oligo_t *ol, int num_ol) {
    int nmatches, min_match, i, j, oligo_len;
    char oligo[MAX_OLIGO_LEN+1], con_tmp[MAX_OLIGO_LEN];
    int match_arr[MAX_MATCHES], score_arr[MAX_MATCHES];
    char *cp;
    float bi, best;

    /*
     * Scan the oligos for matches in the consensus sequence. We expect 1
     * match (the oligo itself). Any more than this is a failure.
     */
    for (i = 0; i < num_ol; i++) {
	UpdateTextOutput();
	cp = &ci->con_item[contig-1][ol[i].start];
	best = 0;

	oligo_len = ol[i].length;

	strncpy(oligo, ol[i].sequence, oligo_len);
	oligo[oligo_len] = '\0';

	min_match = match_perc * oligo_len;


	/* Mask the consensus for this oligo with dashes */
	strncpy(con_tmp, cp, oligo_len);
	for (j=0; j < oligo_len; j++)
	    cp[j] = '-';

	/* Find matches in the forward sense */
	nmatches = inexact_match(ci->con_all, ci->con_len, oligo, oligo_len,
				 min_match, match_arr, score_arr, MAX_MATCHES);
	if (nmatches < 0) {
	    ol[i].matches += MAX_MATCHES+1;
	} else {
	    ol[i].matches += nmatches;
	}
	bi = (float)best_inexact_match(ci->con_all, ci->con_len, oligo,
					 oligo_len, NULL) / (float)oligo_len;
	if (bi > best)
	    best = bi;

	/* Find matches in the reverse sense */
	complement_seq(oligo, oligo_len);
	nmatches = inexact_match(ci->con_all, ci->con_len, oligo, oligo_len,
				 min_match, match_arr, score_arr, MAX_MATCHES);
	if (nmatches < 0) {
	    ol[i].matches += MAX_MATCHES+1;
	} else {
	    ol[i].matches += nmatches;
	}
	bi = (float)best_inexact_match(ci->con_all, ci->con_len, oligo,
					 oligo_len, NULL) / (float)oligo_len;
	if (bi > best)
	    best = bi;
	ol[i].best_match = best;


	/* Unmask the consensus for this oligo with dashes */
	strncpy(cp, con_tmp, oligo_len);
    }

    return 0;
}


/*
 * Get_probes uses primlib to find oligos within a given range of the consensus
 * sequence.
 *
 * Returns oligo_t* for success (maybe NULL *num_oligos is 0) and
 *         stores the size of this result in *num_oligos.
 * Returns NULL for failure with *num_oligos == -1.
 */
static oligo_t *get_probes(GapIO *io, int contig, int from, int to,
			   int sense, int min_size, int max_size,
			   float match_perc, consen_info *ci,
			   int *num_oligos, char *primer_defs) {
    int i, ok, len;
    char scan[MAX_SCAN_REGION+1];
    oligo_t *ol;
    primlib_state *state;
    primlib_args *args;

    *num_oligos = 0;
    UpdateTextOutput();

    /* Determine range to search for oligos within */
    if (to < 0) {
	int tmp = to;
	to = io_clength(io, contig) + from - 1;
	from = io_clength(io, contig) + tmp - 1;
    }

    if (to < 0)
	to = 0;
    if (from < 0)
	from = 0;
    if (to >= io_clength(io, contig))
	to = io_clength(io, contig)-1;
    if (from >= io_clength(io, contig))
	from = io_clength(io, contig)-1;

    if (to - from > MAX_SCAN_REGION)
	to = from + MAX_SCAN_REGION;

    if (to - from < min_size)
	return NULL;

    len = to-from;


    /* Initialise primlib */
    state = primlib_create();
    args = primlib_str2args(primer_defs);
    if (!args) {
	verror(ERR_WARN, "get_probes",
	       "Failed to parse primer arguments");
	return NULL;
    }
    args->min_len = min_size;
    args->max_len = max_size;
    args->opt_len = (min_size + max_size) / 2;
    primlib_set_args(state, args);
    free(args);

    memcpy(scan, &ci->con_item[contig-1][from], len);
    scan[len] = '\0';

    if (sense == REVERSE) {
	complement_seq(scan, len);
    }


    /*
     * Find the oligos. Fix the position to be relative to the original
     * consensus again rather than our 'scan' region.
     */
    ok = primlib_choose(state, scan);
    if (-1 == ok || 0 == state->nprimers) {
	/* *num_oligos = -1; */
	return NULL;
    }

    if (NULL == (ol = (oligo_t *)xmalloc(state->nprimers * sizeof(oligo_t)))) {
	*num_oligos = -1;
	primlib_destroy(state);
	return NULL;
    }

    for (i = 0; i < state->nprimers; i++) {
	if (sense == REVERSE) {
	    ol[i].start = to - 1 -
		(state->primers[i].start + state->primers[i].length - 1);
	} else {
	    ol[i].start = state->primers[i].start + from;
	}
	ol[i].length = state->primers[i].length;

	if (ol[i].length > MAX_OLIGO_LEN)
	    ol[i].length = MAX_OLIGO_LEN;

	ol[i].tm = state->primers[i].temp;
	ol[i].primer_score = state->primers[i].quality;
	strncpy(ol[i].sequence, &ci->con_item[contig-1][ol[i].start],
		ol[i].length);
	ol[i].sequence[ol[i].length] = '\0';
	if (sense == REVERSE) {
	    complement_seq(ol[i].sequence, ol[i].length);
	}

	ol[i].matches = 0;

	if (sense == FORWARD) {
	    ol[i].distance = ol[i].start;
	} else {
	    ol[i].distance = io_clength(io, contig) -
		(ol[i].start + ol[i].length - 1);
	}
    }

    *num_oligos = state->nprimers;

    primlib_destroy(state);
    return ol;
}


/*
 * Output information about the matches found. In time this may be replace
 * by a proper GUI, but only when we work out what information is required.
 */
static void list_probes(GapIO *io, int contig, int end,
			oligo_t *ol, int num_ol, Tcl_DString *dstr) {
    int i, rejected;
    char buf[1024], *name;

    Tcl_DStringStartSublist(dstr);

    vmessage("Contig %s(%d): %s\n",
	     name = get_read_name(io, io_clnbr(io, contig)),
	     io_clnbr(io, contig),
	     end==FORWARD ? "Start" : "End");

    sprintf(buf, "\"%s\" %d %s",
	    name,
	    io_clnbr(io, contig),
	    end==FORWARD ? "Start" : "End");
    Tcl_DStringAppend(dstr, buf, -1);

    if (0 == num_ol) {
	vmessage("    No oligos found\n");
	Tcl_DStringAppend(dstr, " {} 0}", -1);
	return;
    }

    rejected = 0;
    Tcl_DStringStartSublist(dstr);
    for (i = 0; i < num_ol; i++) {
	int pos;

	if (ol[i].matches > 0) {
	    rejected++;
	    continue;
	}

	pos = end == FORWARD ? ol[i].distance + 1 :
	    io_clength(io, contig) - ol[i].distance + 1;
	vmessage("    Pos %6d, Dist %3d, primer=%2.0f, Tm=%2.0f, match=%2.0f%%, %s\n",
		 pos,
		 ol[i].distance,
		 ol[i].primer_score,
		 ol[i].tm,
		 ol[i].best_match * 100,
		 ol[i].sequence);

	Tcl_DStringStartSublist(dstr);

	sprintf(buf, "%6d %3d %2.0f %2.0f %2.0f \"%s\"",
		pos,
		ol[i].distance,
		ol[i].primer_score,
		ol[i].tm,
		ol[i].best_match * 100,
		ol[i].sequence);
	Tcl_DStringAppend(dstr, buf, -1);

	Tcl_DStringEndSublist(dstr);
    }
    Tcl_DStringEndSublist(dstr);

    if (rejected) {
	vmessage("    Rejected %d oligo%s due to non uniqueness\n",
		 rejected,
		 rejected == 1 ? "" : "s");
    }

    sprintf(buf, " %d", rejected);
    Tcl_DStringAppend(dstr, buf, -1);

    Tcl_DStringEndSublist(dstr);
}


/*
 * As find_probes() routine, but for one end of one contig only.
 */
static void find_probes_end(GapIO *io, int contig, int sense, consen_info *ci,
			    int min_size, int max_size, float match_perc,
			    int from, int to, char *vectors,
			    Tcl_DString *dstr, char *primer_defs) {
    oligo_t *ol;
    int num_oligos;

    ol = get_probes(io, contig, from-1, to-1, sense,
		    min_size, max_size, match_perc, ci, &num_oligos,
		    primer_defs);
    if (-1 == num_oligos) {
	verror(ERR_WARN, "find_probes", "couldn't find oligos");
	return;
    }

    find_uniques_con(io, contig, match_perc, ci, ol, num_oligos);

    if (vectors)
	find_uniques_vector(io, contig, match_perc, vectors, ol, num_oligos);

    list_probes(io, contig, sense, ol, num_oligos, dstr);

    if (ol)
	xfree(ol);
}


/*
 * Find oligos suitable for the 'probe' sequencing strategy. The strategy
 * is to find oligos near the end of the contigs pointing towards the centre
 * of the contig. These oligos are then used in a screening process to find
 * templates containing these oligos. The forward readings for each template
 * is then sequenced and assembled. Note that the oligos chosen therefore need
 * to be unique in the consensus sequence and not found in the vector
 * sequences.
 *
 * This routine lists suitable oligos. Arguments are:
 *
 * min_size, max_size
 *	The range of lengths of suitable oligos
 *
 * match_perc
 *	Oligos matching other sequence with >= match_perc are considered
 *	non-unique.
 *
 * from, to
 *	Which area of contig consensus to search for oligos in. Specified as
 *	offsets from the contig ends.
 *
 * vectors
 *	A file of filenames of vector sequences to screen against.
 *
 * Returns:  -1 for success
 *           0 for failure
 */
int find_probes(GapIO *io, int num_contigs, int *contig_arr,
		int min_size, int max_size, float match_perc,
		int from, int to, char *vectors, char *primer_defs,
		Tcl_DString *dstr) {
    consen_info *ci;
    int i;

    /* Calculate consensus */
    if (NULL == (ci = all_consensus(io, consensus_cutoff)))
	return -1;

    for (i = 0; i < num_contigs; i++) {
	find_probes_end(io, contig_arr[i], FORWARD, ci, min_size, max_size,
			match_perc, from-1, to-1, vectors, dstr,
			primer_defs);

	find_probes_end(io, contig_arr[i], REVERSE, ci, min_size, max_size,
			match_perc, -(from-1), -(to-1), vectors, dstr,
			primer_defs);
    }

    free_all_consensus(ci);

    return 0;
}
