#ifndef _TMAN_DIFF_H
#define _TMAN_DIFF_H

/*
 * This function produces an array of trace sample positions for a given
 * sequence and range using its trace data and original position information.
 * Bases without original positions are assigned new ones by looking at the
 * positions of neighbouring bases.
 */
int *get_trace_pos(Read *r, EdStruct *xx, int seq, int pos, int from, int to,
		   char *s, int do_comp);

#endif
