#ifndef _CONSENSUS_H_
#define _CONSENSUS_H_

#include <tg_gio.h>

typedef struct {
    /* the most likely base call - we never call N here */
    /* A=0, C=1, G=2, T=3, *=4 */
    int call;

    /* Log-odds values for A, C, G and T and gap. 5th is N, with 0 prob */
    float scores[6];

    /* Single phred style call */
    unsigned char phred;
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
int calculate_consensus_simple(GapIO *io, int contig, int start, int end,
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
int calculate_consensus(GapIO *io, int contig, int start, int end,
			consensus_t *cons);

#endif /* _CONSENSUS_H_ */
