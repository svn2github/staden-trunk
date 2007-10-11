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
 * This performs a linear (non-indexed) search for a trace in an SRF archive.
 *
 * It's not intended as a suitable production program or as a library of code
 * to use, but as a test and benchmark statistic.
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
void dump(mFILE *mf, char *name) {
    Read *read = mfread_reading(mf, name, TT_ZTR);
    int i;

    if (read == NULL) {
	fprintf(stderr, "Tracedump was unable to open file %s\n", name );
	return;
    }

    printf("[Trace]\n");
    printf("%s\n", read->trace_name);

    printf("\n[Header]\n");
    printf("%d\t\t# format\n",          read->format);
    printf("%d\t\t# NPoints\n",         read->NPoints);
    printf("%d\t\t# NBases\n",          read->NBases);
    printf("%d\t\t# NFlows\n",          read->nflows);
    printf("%d\t\t# maxTraceVal\n",     (int)read->maxTraceVal);
    printf("%d\t\t# baseline\n",        read->baseline);
    printf("%d\t\t# leftCutoff\n",      read->leftCutoff);
    printf("%d\t\t# rightCutoff\n",     read->rightCutoff);

    puts("\n[Bases]");
    for (i = 0; i < read->NBases; i++) {
    printf("%c %05d %+03d %+03d %+03d %+03d #%3d\n",
           read->base[i],
           read->basePos ? read->basePos[i] : 0,
           (int)read->prob_A[i],
           (int)read->prob_C[i],
           (int)read->prob_G[i],
           (int)read->prob_T[i],
           i);
    }

    if (read->NPoints) {
	puts("\n[A_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceA[i]-32768, i);

	puts("\n[C_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceC[i]-32768, i);

	puts("\n[G_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceG[i]-32768, i);

	puts("\n[T_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceT[i]-32768, i);
    }

    if (read->flow_order) {
	puts("\n[Flows]");
	for (i = 0; i < read->nflows; i++) {
	    printf("%c %5.2f  %u\t#%5d\n",
		   read->flow_order[i],
		   read->flow ? read->flow[i] : 0,
		   read->flow_raw ? read->flow_raw[i] : 0,
		   i);
	}
    }

    if (read->info) {
	puts("\n[Info]");
	printf("%s\n", read->info);
    }

    read_deallocate(read);
}

/* ------------------------------------------------------------------------ */
int main(int argc, char **argv) {
    char *ar_name;
    mFILE *mf;
    srf_t *srf;
    char name[512];

    if (argc != 2) {
	fprintf(stderr, "Usage: srf_dump archive_name\n");
	return 1;
    }
    ar_name = argv[1];

    if (NULL == (srf = srf_open(ar_name, "r"))) {
	perror(ar_name);
	return 4;
    }

    while (NULL != (mf = srf_next_trace(srf, name))) {
	dump(mf, name);
	mfclose(mf);
    }

    srf_destroy(srf, 1);

    return 0;
}
