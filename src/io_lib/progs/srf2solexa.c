/*
 * ======================================================================
 * This software has been created by Genome Research Limited (GRL).
 *
 * GRL hereby grants permission to use, copy, modify and distribute
 * this software and its documentation for non-commercial purposes
 * without fee at the user's own risk on the basis set out below.
 *
 * GRL neither undertakes nor accepts any duty whether contractual or
 * otherwise in connection with the software, its use or the use of
 * any derivative, and makes no representations or warranties, express
 * or implied, concerning the software, its suitability, fitness for
 * a particular purpose or non-infringement.
 *
 * In no event shall the authors of the software or GRL be responsible
 * or liable for any loss or damage whatsoever arising in any way
 * directly or indirectly out of the use of this software or its
 * derivatives, even if advised of the possibility of such damage.
 *
 * Our software can be freely distributed under the conditions set out
 * above, and must contain this copyright notice.
 * ======================================================================
 */

/*
 * This converts an SRF file back to the solexa text format. It's has
 * minimal checking and is designed purely for debugging to test the validity
 * of SRF files.
 */

#include <stdio.h>
#include "Read.h"
#include "misc.h"
#include "ztr.h"
#include "srf.h"

/* ------------------------------------------------------------------------ */

/*
 * Ripped out of io_lib's trace_dump program.
 * It reformats a trace to as printable ASCII.
 */
void dump(mFILE *mf, char *name, FILE *fp_seq, FILE *fp_sig2, FILE *fp_prb) {
    Read *read = mfread_reading(mf, name, TT_ZTR);
    int i;
    int lane=1, tile=1, x=1, y=1;

    if (read == NULL) {
	fprintf(stderr, "Tracedump was unable to open file %s\n", name );
	return;
    }

    sscanf(name, "s_%d_%d_%d_%d", &lane, &tile, &x, &y);

    /* Sequence */
    fprintf(fp_seq, "%d\t%d\t%d\t%d\t%.*s\n",
	    lane, tile, x, y, read->NBases, read->base);

    /* Traces */
    fprintf(fp_sig2, "%d\t%d\t%d\t%d", lane, tile, x, y);
    for (i = 0; i < read->NBases; i++) {
	fprintf(fp_sig2, "\t%d %d %d %d",
		(int)read->traceA[i] - 32768,
		(int)read->traceC[i] - 32768,
		(int)read->traceG[i] - 32768,
		(int)read->traceT[i] - 32768);
    }
    fprintf(fp_sig2, "\n");

    /* Confidence */
    for (i = 0; i < read->NBases; i++) {
	if (i)
	    fputc('\t', fp_prb);
	fprintf(fp_prb, "%4d %4d %4d %4d",
		(int)read->prob_A[i],
		(int)read->prob_C[i],
		(int)read->prob_G[i],
		(int)read->prob_T[i]);
    }
    fprintf(fp_prb, "\n");

    read_deallocate(read);
}

/* ------------------------------------------------------------------------ */
int main(int argc, char **argv) {
    char *ar_name;
    mFILE *mf;
    srf_t *srf;
    char name[512];
    char *prefix = "s_lane_tile_";
    FILE *fp_seq, *fp_sig2, *fp_prb;

    if (argc != 2 && argc != 3) {
	fprintf(stderr, "Usage: srf_dump archive_name [output_prefix]\n");
	return 1;
    }
    ar_name = argv[1];
    if (argc == 3)
	prefix = argv[2];

    if (NULL == (srf = srf_open(ar_name, "r"))) {
	perror(ar_name);
	return 4;
    }

    sprintf(name, "%sseq.txt", prefix);
    if (NULL == (fp_seq = fopen(name, "w"))) {
	perror(name);
	return 5;
    }
    sprintf(name, "%sprb.txt", prefix);
    if (NULL == (fp_prb = fopen(name, "w"))) {
	perror(name);
	return 5;
    }
    sprintf(name, "%ssig2.txt", prefix);
    if (NULL == (fp_sig2 = fopen(name, "w"))) {
	perror(name);
	return 5;
    }

    while (NULL != (mf = srf_next_trace(srf, name))) {
	dump(mf, name, fp_seq, fp_sig2, fp_prb);
	mfclose(mf);
    }

    fclose(fp_seq);
    fclose(fp_sig2);
    fclose(fp_prb);

    srf_destroy(srf, 1);
    
    return 0;
}
