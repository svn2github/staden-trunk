/*
 * File: oligo.h
 * Version:
 *
 * Author:
 *         MRC Laboratory of Molecular Biology
 *	   Hills Road
 *	   Cambridge CB2 2QH
 *	   United Kingdom
 *
 * Description:
 *
 * Created:
 * Updated:
 *
 */

#ifndef _OLIGO_H_
#define _OLIGO_H_

#include "edUtils.h"

/*
 * Find suitable oligos using Primlib with current parameter settings.
 * Return number of oligos found (or -1 for error).
 */
int edSelectOligoGenerate(EdStruct *xx, int sense, int bkwd_width,
			  int fwd_width, int avg_read_len, char *pdefs);

/*
 * We cycle through the oligo list
 * curr-oligo gives the current oligo entry under consideration.
 *
 * We return a list of strings in the form:
 * default_tname tname ...
 *
 * The string returned is malloced and so is expected to be freed by the
 * calling routine.
 */
char *edSelectOligoNext(EdStruct *xx);
char *edSelectOligoPrev(EdStruct *xx);

/*
 * Select the current oligo.
 * Creats a tag for it and everything else...
 *
 * Returns a status line containing the template and sequence
 */
char *edSelectOligoAccept(EdStruct *xx, char *template_name);

/*
 * Frees up the temporary structures used by the oligo selection code.
 */
void edSelectOligoQuit(EdStruct *xx);

/*
 * Creates a temporary annotation view, much like the suggest primers and the
 * pick pcr primer functions.
 * It bypasses the undo mechanism and will not be saved.
 *
 * NB: Currently this ignores many of the arguments.
 */
void createTmpAnnotation(EdStruct *xx, int seq, int pos, int len,
			 char *type, char *text, int strand);

#endif /* _OLIGO_H_ */
