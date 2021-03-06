#ifndef _CONSENSUS_H_
#define _CONSENSUS_H_

#include <tg_gio.h>

#define CONS_NO_END_N   1
#define CONS_SCORES     2
#define CONS_DISCREP    4
#define CONS_COUNTS     8

#define CONS_ALL        15

typedef struct {
    /* the most likely base call - we never call N here */
    /* A=0, C=1, G=2, T=3, *=4 */
    int call;

    /* The most likely heterozygous base call */
    /* Use "ACGT*"[het / 5] vs "ACGT*"[het % 5] for the combination */
    int het_call;

    /* Log-odds values for A, C, G and T and gap. 5th is N, with 0 prob */
    /* scores[6] is score for het_call above */
    float scores[7];

    /* Single phred style call */
    unsigned char phred;

    /* Sequence depth */
    int depth;

    /* Individual base type counts */
    int counts[6];

    /* Discrepancy search score */
    float discrep;
} consensus_t;


/*
 * The consensus calculation function - rewritten for tgap style
 * databases.
 *
 * Here con/qual holds the consensus sequence and a single phred style
 * quality value. Either can be passed in as NULL if that array is not needed,
 * but if supplied it is up to the caller to ensure they are large enough.
 *
 * start and end specify the inclusive range, counting from base 0.
 *
 * Returns 0 on success,
 *        -1 on failure
 */
int calculate_consensus_simple(GapIO *io, tg_rec contig, int start, int end,
			       char *con, float *qual);

/*
 * The consensus calculation function - rewritten for tgap style
 * databases.
 *
 * Similar to calculate_consensus_simple, but we present the full
 * probabilities instead of just the called base.
 * "cons" is filled out with consensus_t structs. It should be passed in
 * by the caller and is expected to be the correct length.
 *
 * start and end specify the inclusive range, counting from base 0.
 *
 * Returns 0 on success,
 *        -1 on failure
 *
 */
int calculate_consensus(GapIO *io, tg_rec contig, int start, int end,
			consensus_t *cons);
int calculate_consensus_fast(GapIO *io, tg_rec contig, int start, int end,
			     consensus_t *cons);

/*
 * Internal function, called by calculate_consensus(). Exposed here for when
 * we already have a rangec_t array loaded so we can avoid recomputing the
 * same data. (Used in quality_plot.c of the template display).
 */
int calculate_consensus_bit_het(GapIO *io, tg_rec contig,
				int start, int end,
				int end_N,
				rangec_t *r,
				int nr,
				consensus_t *cons);
/*
 * Finds the portion of a contig that has non-clipped data.
 * This is a somewhat crude method by just computing consensus at the ends
 * and trimming back the zero-depth regions.
 *
 * Specify start and end as pointers for the results. Passing over NULL
 * indicates that you are not interested in that end.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
int consensus_valid_range(GapIO *io, tg_rec contig, int *start, int *end);

/*
 * As per consensus_valid_range(), but we return the unclipped size of the
 * contig instead. This is somewhat easier and faster to compute.
 * Passing NULL as start or end implies you are uninterested in that end.
 *
 * Returns 0 for success
 *        -1 for failure.
 */
int consensus_unclipped_range(GapIO *io, tg_rec contig, int *start, int *end);

/*
 * Converts a padded position into an unpadded position.
 * Returns 0 for success and writes to upos
 *        -1 for error
 */
int consensus_unpadded_pos(GapIO *io, tg_rec contig, int pos, int *upos);

/*
 * Converts an unpadded position into a padded position.
 * Returns 0 for success and writes to pos
 *        -1 for error
 */
int consensus_padded_pos(GapIO *io, tg_rec contig, int upos, int *pos);

#endif /* _CONSENSUS_H_ */
