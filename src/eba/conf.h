#include "scf.h"

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
 */
void calc_conf_values(Read *r, int phred_scale);
