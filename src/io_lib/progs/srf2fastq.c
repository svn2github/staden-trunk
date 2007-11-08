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
#include <math.h>
#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/srf.h>

static int qlookup[256];
void init_qlookup(void) {
    int i;
    for (i = -128; i < 128; i++) {
        qlookup[i+128] = '!' + (int)((10*log(1+pow(10, i/10.0))/log(10)+.499));
    }
}

/* ------------------------------------------------------------------------ */

#define MAX_READ_LEN 1024
void ztr2fastq(ztr_t *z, char *name) {
    Read *read;
    int i;
    char buf[MAX_READ_LEN*2 + 512 + 6];
    char *seq, *qual;

    read = ztr2read(z); /* Inefficient; can do direct */

    if (read == NULL) {
	fprintf(stderr, "Tracedump was unable to open file %s\n", name );
	return;
    }
    if (read->NBases > MAX_READ_LEN) {
	fprintf(stderr, "Recompiler with larger MAX_READ_LEN (%d)\n",
		read->NBases);
	return;
    }

    /* Copy name and leave seq ptr at place to store sequence */
    seq = buf;
    *seq++ = '@';
    while (*name)
	*seq++ = *name++;
    *seq++ = '\n';

    /* Set qual at place to store quality */
    qual = seq + read->NBases;
    *qual++ = '\n';
    *qual++ = '+';
    *qual++ = '\n';
    qual[read->NBases] = '\n';

    /* Decode and copy */
    for (i = 0; i < read->NBases; i++) {
	switch (read->base[i]) {
	case 'A': case 'a':
	    qual[i] = qlookup[read->prob_A[i]+128];
	    seq[i] = 'A';
	    break;

	case 'C': case 'c':
	    qual[i] = qlookup[read->prob_C[i]+128];
	    seq[i] = 'C';
	    break;

	case 'G': case 'g':
	    qual[i] = qlookup[read->prob_G[i]+128];
	    seq[i] = 'G';
	    break;

	case 'T': case 't':
	    qual[i] = qlookup[read->prob_T[i]+128];
	    seq[i] = 'T';
	    break;

	default:
	    qual[i] = qlookup[read->prob_A[i]+128];
	    seq[i] = 'N';
	    break;
	}
    }

    fwrite(buf, 1, qual+read->NBases+1 - buf, stdout);

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
	fprintf(stderr, "Usage: srf2fastq archive_name\n");
	return 1;
    }
    ar_name = argv[1];

    if (NULL == (srf = srf_open(ar_name, "r"))) {
	perror(ar_name);
	return 4;
    }

    read_sections(READ_BASES);
    init_qlookup();

    while (NULL != (ztr = srf_next_ztr(srf, name))) {
	ztr2fastq(ztr, name);
	delete_ztr(ztr);
    }

    srf_destroy(srf, 1);

    return 0;
}
