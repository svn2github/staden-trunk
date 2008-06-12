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

static char qlookup[256];
void init_qlookup(void) {
    int i;
    for (i = -128; i < 128; i++) {
        qlookup[i+128] = '!' + (int)((10*log(1+pow(10, i/10.0))/log(10)+.499));
    }
}

/* ------------------------------------------------------------------------ */

#define MAX_READ_LEN 1024
void ztr2fastq(ztr_t *z, char *name, int calibrated) {
    int i, nc, seq_len;
    char buf[MAX_READ_LEN*2 + 512 + 6];
    char *seq, *qual, *sdata, *qdata;
    ztr_chunk_t **chunks;

    /* Extract the sequence only */
    chunks = ztr_find_chunks(z, ZTR_TYPE_BASE, &nc);
    if (nc != 1) {
	fprintf(stderr, "Zero or greater than one BASE chunks found.\n");
	if (chunks)
	    free(chunks);
	return;
    }
    uncompress_chunk(z, chunks[0]);
    sdata = chunks[0]->data+1;
    seq_len = chunks[0]->dlength-1;

    /* Extract the quality */
    free(chunks);
    if (calibrated)
	chunks = ztr_find_chunks(z, ZTR_TYPE_CNF1, &nc);
    else
	chunks = ztr_find_chunks(z, ZTR_TYPE_CNF4, &nc);

    if (nc != 1) {
	fprintf(stderr, "Zero or greater than one CNF chunks found.\n");
	if (chunks)
	    free(chunks);
	return;
    }
    uncompress_chunk(z, chunks[0]);
    qdata = chunks[0]->data+1;

    /* Construct fastq entry */
    seq = buf;
    *seq++ = '@';
    while (*name)
	*seq++ = *name++;
    *seq++ = '\n';

    qual = seq + seq_len;
    *qual++ = '\n';
    *qual++ = '+';
    *qual++ = '\n';

    for (i = 0; i < seq_len; i++) {
	if (*sdata != '.') {
	    *seq++ = *sdata++;
	} else {
	    *seq++ = 'N';
	    sdata++;
	}
	*qual++ = calibrated
	    ? *qdata++ + '!'
	    : qlookup[*qdata++ + 128];
    }
    *qual++ = '\n';

    fwrite(buf, 1, qual - buf, stdout);
    free(chunks);
}

/* ------------------------------------------------------------------------ */
void usage(void) {
    fprintf(stderr, "Usage: srf2fastq [-c] [-C] archive_name ...\n");
    exit(0);
}

int main(int argc, char **argv) {
    int calibrated = 0;
    int mask = 0, i;

    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;
	} else if (!strcmp(argv[i], "-C")) {
	    mask = SRF_READ_FLAG_BAD_MASK;
	} else if (!strcmp(argv[i], "-c")) {
	    calibrated = 1;
	} else {
	    usage();
	}
    }    

    if (i == argc) {
	usage();
    }

    read_sections(READ_BASES);
    init_qlookup();

    for (; i < argc; i++) {
	char *ar_name;
	srf_t *srf;
	char name[512];
	ztr_t *ztr;

	ar_name = argv[i];

	if (NULL == (srf = srf_open(ar_name, "r"))) {
	    perror(ar_name);
	    return 4;
	}


	while (NULL != (ztr = srf_next_ztr(srf, name, mask))) {
	    ztr2fastq(ztr, name, calibrated);
	    delete_ztr(ztr);
	}

	srf_destroy(srf, 1);
    }

    return 0;
}
