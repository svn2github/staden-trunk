#include <io_lib/scf.h>

/*
 * MODULE    SeqQual12
 *
 * area of the called divided by the area of the max non-called
 * for range from 1/2 between this base and the lastbase to 1/2
 * between this base and the next base
 *
 * So find the area for each base....divide area of the called
 * by other max area
 *
 * phred_scale is a boolean flag to indicate whether the values should be
 * rescaled to be the same as used by phred.
 *
 * average_qual indicates whether the quality values should be averaged.
 *
 * cosarea indicates whether to use a raised cosine (width = +/- 5 samples)
 * instead of a square for computing the trace area.
 */
void calc_conf_values(Read *r, int phred_scale, int average_qual, int offset);

/*
 * Averages the confidence arrays in r over a window of length
 * WINDOW_SIZE (#define at top of conf.c).
 */
int average_conf(Read *r);

/*
 * Rescales eba scores to phred-style scores.
 */
void rescale_scores(Read *r, int table);

/*
 * Applies a simple filter to the trace data, producing new traces.
 * Equivalent to 1st-order differentiation followed by fitting a +-1 square
 * wave of width 11.
 */
void filter_trace(Read *r);

