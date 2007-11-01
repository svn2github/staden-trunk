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
 * If this program or something akin to it is every seriously used in earnest
 * rather than just an example/debugging tool then note that printf is indeed
 * something like 78% of the entire CPU time.
 *
 * Not worth the hassle probably, but here's an example of a faster printf
 * style function to replace printf("%5d", i) with fast_printf(buf, i, 5, ' ').
 */
#if 0
static char *fast_printf(char *buf, int val, int prec, char padding) {
    int digits[100];
    int i;

    if (val < 0) {
	val = -val;
	*buf++ = '-';
    }

    if (val) {
	for (i = 0; val; val /= 10)
	    digits[i++] = val % 10;
    } else {
	digits[0] = 0;
	i = 1;
    }

    for (; prec > i; prec--)
	*buf++ = padding;

    for (; i-- > 0;)
	*buf++ = digits[i] + '0';
    *buf++ = 0;

    return buf-1;
}
#endif

/*
 * Ripped out of io_lib's trace_dump program.
 * It reformats a trace to as printable ASCII.
 */
void dump(ztr_t *z, char *name) {
    Read *read;
    int i;

    uncompress_ztr(z);
    read = ztr2read(z); /* Inefficient; can do direct */

    if (read == NULL) {
	fprintf(stderr, "Tracedump was unable to open file %s\n", name );
	return;
    }

    printf("[Trace]\n");
    printf("%s\n", name);

    printf("\n[Header]\n");
    printf("%d\t\t# format\n",          read->format);
    printf("%d\t\t# NPoints\n",         read->NPoints);
    printf("%d\t\t# NBases\n",          read->NBases);
    printf("%d\t\t# NFlows\n",          read->nflows);
    printf("%d\t\t# maxTraceVal\n",     (int)read->maxTraceVal-read->baseline);
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
	    printf("%d\t#%5d\n", (int)read->traceA[i] - read->baseline, i);

	puts("\n[C_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceC[i] - read->baseline, i);

	puts("\n[G_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceG[i] - read->baseline, i);

	puts("\n[T_Trace]");
	for(i = 0; i < read->NPoints; i++)
	    printf("%d\t#%5d\n", (int)read->traceT[i] - read->baseline, i);
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
    ztr_t *ztr;

    if (argc != 2) {
	fprintf(stderr, "Usage: srf_dump archive_name\n");
	return 1;
    }
    ar_name = argv[1];

    if (NULL == (srf = srf_open(ar_name, "r"))) {
	perror(ar_name);
	return 4;
    }

    while (NULL != (ztr = srf_next_ztr(srf, name))) {
	dump(ztr, name);
	delete_ztr(ztr);
    }

    srf_destroy(srf, 1);

    return 0;
}
