#include <string.h>

#include "xalloc.h"
#include "seq_results.h"
#include "sip_results.h"
#include "readpam.h"
#include "dna_utils.h"

void SipRescanMatches(Tcl_Interp *interp,
		      seq_result *result,
		      int id,
		      int min_score) 
{
    out_raster *output = result->output;
    d_plot *data = result->data;
    int num_pts = data->n_pts;
    int i, j;
    double coords[2];
    Tk_Raster *raster;
    Tcl_CmdInfo info;
    char *opts[3];
    int index;
    char *seq1, *seq2;
    int seq1_len, seq2_len;
    int score;
    int start;
    int x, y;
    int index1, index2;
    double wx0, wy0, wx1, wy1;

    if (output->hidden) {
	return;
    }

    index1 = GetSeqNum(result->seq_id[0]);
    index2 = GetSeqNum(result->seq_id[1]);

    /* if either index is -1; the corresponding seq must have been deleted */
    if ((index1 == -1) || (index2 == -1)) {
	return;
    }
    seq1 = GetSeqSequence(index1);
    seq2 = GetSeqSequence(index2);
    seq1_len = GetSeqLength(index1);
    seq2_len = GetSeqLength(index2);

    Tcl_GetCommandInfo(interp, output->raster_win, &info);
    raster = (Tk_Raster*)info.clientData;
    opts[0] = "-fg";
    opts[1] = "purple";
    opts[2] = NULL;
    index = CreateDrawEnviron(interp, raster, 2, opts);

    SetDrawEnviron(output->interp, raster, index);

    RasterGetWorldScroll(raster, &wx0, &wy0, &wx1, &wy1);

    start = (int)(data->win_len / 2);

    /* printing */
    for (i = 0; i < num_pts; i++) {

	x = data->p_array[i].x - start;
	y = data->p_array[i].y - start;

	for (j = 0; j < data->win_len; j++, x++, y++) {
	    
	    if ((x < 1) || (y < 1) || (x > seq1_len) || (y > seq2_len)) {
		continue;
	    }
	    score = score_matrix[char_lookup[seq1[x - 1]]]
	      [char_lookup[seq2[y - 1]]];
    
	    if (score >= min_score) {
		/* printf("j %d x %d y %d score %d \n", j, x, y, score); */
		coords[0] =  x;
		coords[1] =  (int)wy1 - y + wy0;
		RasterDrawPoints(raster, coords, 1);
	    }
	}
    }
    tk_RasterRefresh(raster);
}
