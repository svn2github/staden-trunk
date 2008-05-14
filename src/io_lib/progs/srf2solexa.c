/* 
 * IMPORTANT NOTE: This file differs from the original software created by
 * Genome Research Limited (GRL). It has been modified by Illumina for
 * integration with Illumina's Genome Analyzer data analysis pipeline.
 */

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

#ifdef HAVE_CONFIG_H
#include "io_lib_config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/srf.h>

#define RAW       (1<<0)
#define PROCESSED (1<<1)

/* ------------------------------------------------------------------------ */

/*
 * Parses a name assuming it consists of:
 * prefix<separator><lane><separator><tile><separator><x><separator><y>
 *
 * We fill out the supplied lane, tile, x and y parameters, or set them to
 * 0, 0, 0 and 0 if unknown.
 *
 * Returns 0 for success
 *        -1 for failure (unknown)
 */
int parse_name(char *name, int *lane, int *tile, int *x, int *y) {
    size_t len = strlen(name), i;

    /* FIXME: lack of error checking */

    i = len-1;
    while (i >= 0 && isdigit(name[i]))
	i--;
    *y = atoi(&name[i+1]);
    i--;
    
    while (i >= 0 && isdigit(name[i]))
	i--;
    *x = atoi(&name[i+1]);
    i--;

    while (i >= 0 && isdigit(name[i]))
	i--;
    *tile = atoi(&name[i+1]);
    i--;

    while (i >= 0 && isdigit(name[i]))
	i--;
    *lane = atoi(&name[i+1]);
    i--;

    return 0;
}

#define MAX_CYCLES 500
void dump_samples4(FILE *fp, int lane, int tile, int x, int y,
		   int baseline, unsigned char *bytes, int nbytes) {
    int i, j, nb;
    int A[MAX_CYCLES], C[MAX_CYCLES], G[MAX_CYCLES], T[MAX_CYCLES];

    nb = nbytes/8;

    for (i = j = 0; i < nb; i++, j+= 2)
	A[i] = (bytes[j] << 8) | bytes[j+1];

    for (i = 0; i < nb; i++, j+= 2)
	C[i] = (bytes[j] << 8) | bytes[j+1];

    for (i = 0; i < nb; i++, j+= 2)
	G[i] = (bytes[j] << 8) | bytes[j+1];

    for (i = 0; i < nb; i++, j+= 2)
	T[i] = (bytes[j] << 8) | bytes[j+1];

    fprintf(fp, "%d\t%d\t%d\t%d", lane, tile, x, y);
    for (i = 0; i < nb; i++)
	fprintf(fp, "\t%d %d %d %d",
		A[i] - baseline,
		C[i] - baseline,
		G[i] - baseline,
		T[i] - baseline);
    fprintf(fp, "\n");
}

void dump_conf4(FILE *fp, char *seq, signed char *bytes, int nbytes) {
    int i, j, nb;
    int A[MAX_CYCLES], C[MAX_CYCLES], G[MAX_CYCLES], T[MAX_CYCLES];

    nb = nbytes/4;

    for (i = 0, j = nb; i < nb; i++) {
	switch(seq[i]) {
	case 'A': case 'a':
	    A[i] = bytes[i];
	    C[i] = bytes[j++];
	    G[i] = bytes[j++];
	    T[i] = bytes[j++];
	    break;
	case 'C': case 'c':
	    A[i] = bytes[j++];
	    C[i] = bytes[i];
	    G[i] = bytes[j++];
	    T[i] = bytes[j++];
	    break;
	case 'G': case 'g':
	    A[i] = bytes[j++];
	    C[i] = bytes[j++];
	    G[i] = bytes[i];
	    T[i] = bytes[j++];
	    break;
	default:
	    A[i] = bytes[j++];
	    C[i] = bytes[j++];
	    G[i] = bytes[j++];
	    T[i] = bytes[i];
	    break;
	}
    }

    for (i = 0; i < nb; i++) {
	if (i)
	    fputc('\t', fp);
	fprintf(fp, "%4d %4d %4d %4d", A[i], C[i], G[i], T[i]);
    }

    fprintf(fp, "\n");
}

void dump_qcal(FILE *fp[], signed char *bytes, int nbytes) {
    int read_1 = strlen(bytes);
    int i = 0;
    while (read_1 > i) {
        fprintf(fp[0], "%c", (char)(*(bytes+i)));
        ++i;
    }
    fprintf(fp[0], "\n");
    if (NULL != fp[1]) {
        ++i; /* skip the 0 */
        while (nbytes > i) {
            fprintf(fp[1], "%c", (char)(*(bytes+i)));
            ++i;
        }
        fprintf(fp[1], "\n");
    }
}

/*
 * Ripped out of io_lib's trace_dump program.
 * It reformats a trace to as printable ASCII.
 */
void dump(ztr_t *z,
	  int lane, int tile, int x, int y,
	  FILE *fp_int, FILE *fp_nse,
	  FILE *fp_seq, FILE *fp_sig2, FILE *fp_prb, FILE *fp_qcal[]) {
    int i, nc;
    ztr_chunk_t **chunks;
    char *seq;

    uncompress_ztr(z);

    /* Sequence */
    if (fp_seq || fp_prb) {
	chunks = ztr_find_chunks(z, ZTR_TYPE_BASE, &nc);
	if (nc != 1) {
	    fprintf(stderr, "Zero or greater than one BASE chunks found.\n");
	    return;
	}
	seq = chunks[0]->data+1;

	if (fp_seq)
	    fprintf(fp_seq, "%d\t%d\t%d\t%d\t%.*s\n",
		    lane, tile, x, y,
		    chunks[0]->dlength-1,
		    chunks[0]->data+1);
    }

    /* Confidence */
    if (fp_prb) {
	chunks = ztr_find_chunks(z, ZTR_TYPE_CNF4, &nc);
	if (nc != 1) {
	    fprintf(stderr, "Zero or greater than one CNF chunks found.\n");
	    return;
	}

	dump_conf4(fp_prb, seq, chunks[0]->data+1, chunks[0]->dlength-1);
    }

    /* QCal */
    if (fp_qcal[0]) {
	chunks = ztr_find_chunks(z, ZTR_TYPE_CNF1, &nc);
	if (nc != 1) {
	    fprintf(stderr, "Zero or greater than one CNF chunks found.\n");
	    return;
	}

	dump_qcal(fp_qcal, chunks[0]->data+1, chunks[0]->dlength-1);
    }

    /* Traces */
    if (fp_sig2) {
	chunks = ztr_find_chunks(z, ZTR_TYPE_SMP4, &nc);
	for (i = 0; i < nc; i++) {
	    char *key = ztr_lookup_mdata_value(z, chunks[i], "TYPE");
	    if (!key || 0 == strcmp(key, "PROC")) {
		key = ztr_lookup_mdata_value(z, chunks[i], "OFFS");
		dump_samples4(fp_sig2, lane, tile, x, y, key ? atoi(key) : 0,
			      chunks[i]->data+2, chunks[i]->dlength-2);
		break;
	    }
	}
    }

    if (fp_int) {
	chunks = ztr_find_chunks(z, ZTR_TYPE_SMP4, &nc);
	for (i = 0; i < nc; i++) {
	    char *key = ztr_lookup_mdata_value(z, chunks[i], "TYPE");
	    if (key && 0 == strcmp(key, "SLXI")) {
		key = ztr_lookup_mdata_value(z, chunks[i], "OFFS");
		dump_samples4(fp_int, lane, tile, x, y, key ? atoi(key) : 0,
			      chunks[i]->data+2, chunks[i]->dlength-2);
		break;
	    }
	}
    }

    if (fp_nse) {
	chunks = ztr_find_chunks(z, ZTR_TYPE_SMP4, &nc);
	for (i = 0; i < nc; i++) {
	    char *key = ztr_lookup_mdata_value(z, chunks[i], "TYPE");
	    if (key && 0 == strcmp(key, "SLXN")) {
		key = ztr_lookup_mdata_value(z, chunks[i], "OFFS");
		dump_samples4(fp_nse, lane, tile, x, y, key ? atoi(key) : 0,
			      chunks[i]->data+2, chunks[i]->dlength-2);
		break;
	    }
	}
    }
}


/* ------------------------------------------------------------------------ */
FILE *fopen_slx(char *dir, char *type, int lane, int tile) {
    char fn[1024];
    FILE *fp;

    sprintf(fn, "%s/s_%d_%04d_%s.txt", dir, lane, tile, type);
    if (NULL == (fp = fopen(fn, "w+"))) {
	perror(fn);
	return NULL;
    }

    return fp;
}

int dump_text(char *dir, int lane, ztr_chunk_t **chunks)
{
    char fn[1024];
    FILE *fp;
    int j=0, k=0, key=0;
    for (j = 0; j < chunks[0]->dlength; j++) {
        if(chunks[0]->data[j] == 0) {
            ++k;
            if(!(k % 2)) {
                sprintf(fn, "%s/%d_%s.txt", dir, lane, chunks[0]->data + key + 1);
                if (NULL == (fp = fopen(fn, "w+"))) {
                    perror(fn);
                    return -1;
                } else {
                    fprintf(fp, "%s", chunks[0]->data + j + 1);
                    fclose(fp);
                }
            } else {
                key = j;
            }
        }
    }

    return 0;
}

int process_srf(char *file, char *dir, int mode, int qcal,
		int filter_mask) {
    srf_t *srf;
    char name[1024], dir2[1024];
    ztr_t *ztr;
    int last_lane = 0, last_tile = 0;
    FILE *fp_int = NULL;
    FILE *fp_nse = NULL;
    FILE *fp_seq = NULL;
    FILE *fp_prb = NULL;
    FILE *fp_sig = NULL;
    FILE *fp_qcal[2] = {NULL, NULL};

    if (NULL == (srf = srf_open(file, "rb"))) {
	perror(file);
	return 1;
    }

    if (dir) {
	strcpy(dir2, dir);
    } else {
	char *cp = strrchr(file, '.');
	if (cp) *cp = 0;
	sprintf(dir2, "%s.run", file);
    }
    mkdir(dir2, 0777);

    while (NULL != (ztr = srf_next_ztr(srf, name, filter_mask))) {
	int lane, tile, x, y;
	parse_name(name, &lane, &tile, &x, &y);

	if (last_lane != lane) {
	    int nc;
	    ztr_chunk_t **chunks;

	    /* Text chunks */
	    uncompress_ztr(ztr);
	    chunks = ztr_find_chunks(ztr, ZTR_TYPE_TEXT, &nc);
	    dump_text(dir2, lane, chunks);
	}

	if (last_lane != lane || last_tile != tile) {
	    fprintf(stderr, "New tile: %d/%d\n", lane, tile);
	    last_lane = lane;
	    last_tile = tile;

	    if (fp_int) { fclose(fp_int); fp_int = NULL; }
	    if (fp_nse) { fclose(fp_nse); fp_nse = NULL; }
	    if (fp_seq) { fclose(fp_seq); fp_seq = NULL; }
	    if (fp_prb) { fclose(fp_prb); fp_prb = NULL; }
	    if (fp_sig) { fclose(fp_sig); fp_sig = NULL; }
	    if (fp_qcal[0]) { fclose(fp_qcal[0]); fp_qcal[0] = NULL; }
	    if (fp_qcal[1]) { fclose(fp_qcal[1]); fp_qcal[1] = NULL; }

	    if (mode & RAW) {
		fp_int = fopen_slx(dir2, "int", lane, tile);
		fp_nse = fopen_slx(dir2, "nse", lane, tile);
	    }
	    if (mode & PROCESSED) {
		fp_seq = fopen_slx(dir2, "seq", lane, tile);
		fp_prb = fopen_slx(dir2, "prb", lane, tile);
		fp_sig = fopen_slx(dir2, "sig2", lane, tile);
	    }
            if (qcal == 1) {
		fp_qcal[0] = fopen_slx(dir2, "qcal", lane, tile);
            }
            if (qcal ==2) {
		fp_qcal[0] = fopen_slx(dir2, "1_qcal", lane, tile);
		fp_qcal[1] = fopen_slx(dir2, "2_qcal", lane, tile);
            }
	}
	dump(ztr, lane, tile, x, y, fp_int, fp_nse, fp_seq, fp_sig, fp_prb, fp_qcal);
	delete_ztr(ztr);
    }

    if (fp_int) fclose(fp_int);
    if (fp_nse) fclose(fp_nse);
    if (fp_seq) fclose(fp_seq);
    if (fp_prb) fclose(fp_prb);
    if (fp_sig) fclose(fp_sig);
    if (fp_qcal[0]) fclose(fp_qcal[0]);
    if (fp_qcal[1]) fclose(fp_qcal[1]);

    srf_destroy(srf, 1);
    return 0;
}

/* ------------------------------------------------------------------------ */
void usage(void) {
    printf("Usage: srf2illumina [options]\n\n");
    printf(" Options:\n");
    printf("    -r       Output raw (.int/.nse) data only\n");
    printf("    -p       Output processed (.sig2/.seq/.prb) data only\n");
    printf("    -q num   QCal values, with num (1 or 2) reads\n");
    printf("    -d dir   Set the output directory to 'dir'\n");
    printf("    -C       Filter out reads marked as bad.\n");
    exit(0);
}

int main(int argc, char **argv) {
    int mode = RAW | PROCESSED;
    int qcal = 0;
    int filter_mask = 0;
    char *dir = NULL;
    int i;

    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;
	} else if (!strcmp(argv[i], "-r")) {
	    mode = RAW;
	} else if (!strcmp(argv[i], "-p")) {
	    mode = PROCESSED;
	} else if (!strcmp(argv[i], "-q")) {
	    qcal = atoi(argv[++i]);
	} else if (!strcmp(argv[i], "-d")) {
	    dir = argv[++i];
	} else if (!strcmp(argv[i], "-C")) {
	    filter_mask = SRF_READ_FLAG_BAD_MASK;
	} else {
	    usage();
	}
    }

    if (i == argc)
	usage();

    for (; i < argc; i++) {
	if (process_srf(argv[i], dir, mode, qcal, filter_mask)) {
	    perror(argv[i]);
	    return 1;
	}
    }
    
    return 0;
}
