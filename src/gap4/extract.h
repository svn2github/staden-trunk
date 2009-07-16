#include <io_lib/expFileIO.h>

/*
 * C interface
 *
 * Returns:
 * -1 failure
 *  0 success
 */
int extract_readings(GapIO *io,
		     int num_readings,
		     char **reading_array,
		     char *dir,
		     int format);
 
/*
 * Output TG lines.
 * Common code used by both the gel and consensus experiment file routines
 *
 * Offset is a value to add to all annotations. Used by the consensus output
 * to shift gel annotations onto the consensus coordinates.
 *
 * If orig is set, then we output the tags relative to their original positions
 * in the sequences. Otherwise we swap around so that they are in their
 * positions as seen in the contig editor. g_sense and g_len need only be sent
 * over when orig is 0.
 * 
 * clip_left and clip_right, when non zero values, are used to clip annotations
 * outside of a given range. This currently only works (as it's only needed)
 * when not shuffling the tags to their original senses.
 *
 * If comment != NULL, a CC line is outputed before the first annotation.
 *
 * We don't output SVEC and CVEC tags as these are included in the SL, SR
 * and CS lines.
 *
 * Returns 0 for success, -1 for failure.
 */
int output_annotations(GapIO *io, Exp_info *e, int anno, int offset,
		       int orig, int g_sense, int g_len, int consensus,
		       int clip_left, int clip_right, char *cc_line,
		       int *pads, int npads);

/*
 * Output NT lines (notes).
 * Common code used by both the gel and consensus experiment file routines
 *
 * Returns 0 for success, -1 for failure.
 */
int output_notes(GapIO *io, Exp_info *e, int note,
		 int source_type, int source_num);
