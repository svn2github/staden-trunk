/* TODO:
 * Define a skip-pattern (e.g) 0/1/1/1 to indicate how many huffman codes
 * we need and what order they get used in. In the above case for each set of
 * 4 bytes the first byte uses the 1st stored set of huffman codes and
 * the last 3 bytes use the 2nd set. For 16-bit trace data we'd use "0/1".
 *
 * This can be passed into calc_bit_lengths, generate_codes, etc. So the
 * huffman_codes_t struct needs multiple code arrays.
 *
 * Also they get serialised into a single stream for Deflate, so we can
 * inline multiple codes too. The only difference I think is a wrapper
 * around the outside to indicate a skip-pattern. It defaults "0" (each
 * byte using the same set).
 *
 * Skip EOF symbol in code_sets that do not need it? No difference I'd guess
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
 * Author: James Bonfield, March 2006
 * Updated April 2007 to add a format string for naming samples.
 * Updated June 2007 to wrap the output format in SRF (example).
 *
 * This code converts _sig, _prb and _seq files to ZTR files held within
 * an SRF container.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <fcntl.h>

#include "Read.h"
#include "misc.h"
#include "ztr.h"
#include "array.h"
#include <zlib.h>

#include "srf.h"

//#define SINGLE_HUFF

/* A model saves approx 2% data size, vs 1.3% saving for delta */
//#define USE_MODEL

#define N_TR_CODE 8
#define SHIFT_BY 32768

#define BASE_CODE (CODE_USER)
#define SMP4_CODE (CODE_USER+1)
#define CNF4_CODE (CODE_USER+2)

void delta_called(ztr_t *ztr, ztr_chunk_t *chunk) {
    int nc, nbases, i;
    ztr_chunk_t **base = ztr_find_chunks(ztr, ZTR_TYPE_BASE, &nc);
    char *bases;
    unsigned short *told, last;

    if (nc == 0)
	return;

    bases  = base[nc-1]->data+1;
    nbases = base[nc-1]->dlength-1;

    told = ((unsigned short *)chunk->data) + 4;
    last = be_int2(*told);
    *told += 4;
    for (i = 1; i < nbases; i++) {
	unsigned short curr = be_int2(*told);
	*told = be_int2(curr-last);
	last = curr;
	told += 4;
    }
}

#define SYM_EOF 256
static void output_code_set(FILE *fp, huffman_codes_t *cds) {
    int i, j;
    int nbits = 0;
    huffman_code_t *codes = cds->codes;
    int ncodes = cds->ncodes;

    fprintf(fp, "static huffman_code_t codes_?[] = {\n");
    for (i = j = 0; i < ncodes; i++) {
	nbits += codes[i].nbits * codes[i].freq;
	if (j == 0)
	    fprintf(fp, "    ");
	if (codes[i].symbol == SYM_EOF) {
	    fprintf(fp, "{SYM_EOF,%3d}, ", codes[i].nbits);
	    j = 10;
	} else {
	    if (isalnum(codes[i].symbol)) {
		fprintf(fp, "{'%c',%3d}, ", codes[i].symbol, codes[i].nbits);
	    } else {
		fprintf(fp, "{%3d,%3d}, ", codes[i].symbol, codes[i].nbits);
	    }
	}
	j++;
	
	if (j >= 6) {
	    fputc('\n', fp);
	    j = 0;
	}
    }
    if (j)
	fputc('\n', fp);
    fprintf(fp, "};\n");
    fprintf(fp, "/* Tested on %d bits of compressed data */\n", nbits);
}

/* Internal struct used during parsing */
typedef struct {
    int x, y;
    int tile;
    int lane;
} loc_t;

#define MAX_CYCLES 100

/*
 * ztr_mwrite_header
 *
 * Writes a ZTR file header.
 *
 * Arguments:
 * 	fp		A mFILE pointer
 *	h		A pointer to the header to write
 *
 * Returns:
 *	Success:  0
 *	Failure: -1
 */
static int ztr_mwrite_header(mFILE *fp, ztr_header_t *h) {
    if (1 != mfwrite(h, sizeof(*h), 1, fp))
	return -1;

    return 0;
}

/*
 * ztr_mwrite_chunk
 *
 * Writes a ZTR chunk including chunk header and data
 *
 * Arguments:
 * 	fp		A mFILE pointer
 *	chunk		A pointer to the chunk to write
 *
 * Returns:
 *	Success:  0
 *	Failure: -1
 */
static int ztr_mwrite_chunk(mFILE *fp, ztr_chunk_t *chunk) {
    int4 bei4;

    /*
    {
	char str[5];
	fprintf(stderr, "Write chunk %.4s %08x length %d\n",
		ZTR_BE2STR(chunk->type, str), chunk->type, chunk->dlength);
    }
    */

    /* type */
    bei4 = be_int4(chunk->type);
    if (1 != mfwrite(&bei4, 4, 1, fp))
	return -1;

    /* metadata length */
    bei4 = be_int4(chunk->mdlength);
    if (1 != mfwrite(&bei4, 4, 1, fp))
	return -1;

    /* metadata */
    if (chunk->mdlength)
	if (chunk->mdlength != mfwrite(chunk->mdata, 1, chunk->mdlength, fp))
	    return -1;

    /* data length */
    bei4 = be_int4(chunk->dlength);
    if (1 != mfwrite(&bei4, 4, 1, fp))
	return -1;

    /* data */
    if (chunk->dlength != mfwrite(chunk->data, 1, chunk->dlength, fp))
	return -1;

    return 0;
}

/*
 * Fetches a sequence from a *_*_*_seq.txt file. The format is lines of:
 *
 * 3       56      5       818     TACACACAGAGTGAAAGAAAATCT.T.ATCTACACA
 *
 * We don't care about any field here except the last; the sequence.
 * 'trim' designates how many bases to clip off from the left end. Eg in the
 * above the first "T" is present on all sequences and is not part of the
 * result returned.
 *
 * Returns: Seq on success (static, non-reentrant)
 *          NULL on fialure
 */
char *get_seq(FILE *fp, int trim, int *lane, int *tile, int *x, int *y) {
    static char line[1024];
    char *cp;

    if (NULL == fgets(line, 1023, fp))
	return NULL;

    if ((cp = strchr(line, '\n')))
	*cp = 0;
    cp = strrchr(line, '\t');

    if (4 != sscanf(line, "%d %d %d %d", lane, tile, x, y))
	return NULL;

    return cp ? cp + 1 + trim : NULL;
}

/*
 * Fetches the probability arrays as an Nx4 array from the *_prb.txt file.
 * The format of the lines is a series of 4-tuples (4 integer numbers,
 * themselves separated by spaced) which are separated by tabs.
 *
 * The data returned is statically allocated, so it should not be freed and
 * the function is not rentrant.
 *
 * Returns point to Nx4 array of ints holding the log-odds scores.
 *         NULL on failure
 */
int (*get_prb(FILE *fp, int trim))[4] {
    char line[MAX_CYCLES*20 +1];
    static int prb[MAX_CYCLES][4];
    char *cp;
    int c; /* cycle */

    if (NULL == fgets(line, MAX_CYCLES*20, fp))
	return NULL;

    for (c = -trim, cp = strtok(line, "\t"); cp;
	 cp = strtok(NULL, "\t"), c++) {
	if (c < 0)
	    continue;
	if (4 != sscanf(cp, "%d %d %d %d",
			&prb[c][0], &prb[c][1], &prb[c][2], &prb[c][3]))
	    return NULL;
    }
    
    return prb;
}

/*
 * Fetches the signal strength arrays as an Nx4 array from the *_sig.txt file.
 * The format of the lines is a series of 4-tuples (4 floating point numbers,
 * themselves separated by spaced) which are separated by tabs.
 *
 * The data returned is statically allocated, so it should not be freed and
 * the function is not rentrant.
 *
 * Returns point to Nx4 array of ints holding the signal strengths.
 *         NULL on failure
 */
float (*get_sig(FILE *fp, int trim))[4] {
    char line[MAX_CYCLES*30 +1];
    static float sig[MAX_CYCLES][4];
    char *cp;
    int c, i, j; /* cycle */

    if (NULL == fgets(line, MAX_CYCLES*30, fp))
	return NULL;

    for (c = -4-trim, cp = strtok(line, "\t"); cp;
	 cp = strtok(NULL, "\t"), c++) {
	if (c < 0)
	    continue;

	if (4 != sscanf(cp, "%f %f %f %f",
			&sig[c][0], &sig[c][1], &sig[c][2], &sig[c][3]))
	    return NULL;
    }

    return sig;
}

/*
 * Creates an io_lib Read object from sequence, probabilities and signal
 * strengths.
 *
 * Returns: allocated Read on success
 *	    NULL on failure
 */
Read *create_read(char *seq, int (*prb)[4], float (*sig)[4]) {
    size_t nbases = strlen(seq), i;
    Read *r;
    int max = 0;
    
    if (NULL == (r = read_allocate(nbases, nbases)))
	return NULL;

    for (i = 0; i < nbases; i++) {
	/* Confidence values */
#if 0
	/* FIXME: for now we don't have -ve confidences */
	r->prob_A[i] = MAX(prb[i][0], 0);
	r->prob_C[i] = MAX(prb[i][1], 0);
	r->prob_G[i] = MAX(prb[i][2], 0);
	r->prob_T[i] = MAX(prb[i][3], 0);
#else
	r->prob_A[i] = prb[i][0];
	r->prob_C[i] = prb[i][1];
	r->prob_G[i] = prb[i][2];
	r->prob_T[i] = prb[i][3];
#endif

	/* Traces */
	if ((r->traceA[i] = (int)sig[i][0] + SHIFT_BY) > max) max = sig[i][0];
	if ((r->traceC[i] = (int)sig[i][1] + SHIFT_BY) > max) max = sig[i][1];
	if ((r->traceG[i] = (int)sig[i][2] + SHIFT_BY) > max) max = sig[i][2];
	if ((r->traceT[i] = (int)sig[i][3] + SHIFT_BY) > max) max = sig[i][3];

	/* Sequence & position */
	r->base[i] = seq[i];
	r->basePos[i] = i;
    }

    r->maxTraceVal = max;
    r->leftCutoff = 0;
    r->rightCutoff = nbases+1;

    /* Info: none for present */

    return r;
}

void format_name(char *name, char *fmt, int lane, int tile, int x, int y,
		 int count) {
    int n;

    for (; *fmt; fmt++) {
	switch(*fmt) {

	/* A format specifier */
	case '%': {
	    char *endp;
	    int l1 = 0, val = 0;

	    /* Precision specifier */
	    l1 = strtol(++fmt, &endp, 10);
	    fmt = endp;

	    switch(tolower(*fmt)) {
	    case 'l':
		val = lane;
		break;
	    case 't':
		val = tile;
		break;
	    case 'x':
		val = x;
		break;
	    case 'y':
		val = y;
		break;
	    case 'c':
		val = count;
		break;
	    default:
		fprintf(stderr, "Invalid %% rule '%%%c'\n", *fmt);
	    }
	    n = sprintf(name, isupper(*fmt) ? "%0*x" : "%0*d", l1, val);
	    name += n;
	    break;
	}

	default:
	    *name++ = *fmt;
	}
    }

    *name++ = '\0';
}

/*
 * Given one or two ZTR chunk types this aggregates all together and
 * computes a set of huffman codes matching data from those chunk types.
 *
 * Returns huffman_code_t on success
 *         NULL on failure
 */
huffman_codeset_t *ztr2codes(Array za, int code_set, int nc, int type) {
    unsigned char *buf = NULL;
    size_t sz = 0, alloc = 0;
    int nreads = ArrayMax(za);
    int i, j, skip;
    huffman_codeset_t *cds;

    /* Accumulated concatenated blocks of data from all ZTR files */
    for (i = 0; i < nreads; i++) {
	ztr_t *z = arr(ztr_t *, za, i);
	for (j = 0; j < z->nchunks; j++) {
	    int pos, len;

	    if (z->chunk[j].type == type) {
#if 0
		/*
		 * Optimisation for base calls: only use huffman
		 * compression on data that is not pure A,C,G,T.
		 *
		 * For pure ACGT we could use an additional ZTR format that
		 * utilises 2-bit encodings. I predict this would reduce the
		 * space taken up by BASE data blocks by 22% (although this
		 * is overall only 1% of the total file size).
		 */
		if (z->chunk[j].type == ZTR_TYPE_BASE) {
		    int k, l = z->chunk[j].dlength;
		    for (k = 1; k < l; k++) {
			if (z->chunk[j].data[k] != 'A' &&
			    z->chunk[j].data[k] != 'C' &&
			    z->chunk[j].data[k] != 'G' &&
			    z->chunk[j].data[k] != 'T')
			    break;
		    }

		    if (k == l)
			continue;
		}
#endif
		len = (int)((z->chunk[j].dlength + nc-1)/nc) * nc;
		while (len + sz > alloc) {
		    alloc += 65536;
		    if (!(buf = (unsigned char *)realloc(buf, alloc)))
			return NULL;
		}

		memcpy(&buf[sz], z->chunk[j].data, z->chunk[j].dlength);
		sz += z->chunk[j].dlength;

		/* Pad out to ensure data copied is a multiple of nc */
		if (len > z->chunk[j].dlength) {
		    memset(&buf[sz], 0, len - z->chunk[j].dlength);
		    sz += len - z->chunk[j].dlength;
		}
	    }
	}
    }

    /* Analyse concatenated data to compute codes */
    cds = generate_code_set(code_set, nc, buf, sz, nreads, MAX_CODE_LEN, 0);

#if 0
    for (i = 0; i < nc; i++) {
	output_code_set(stderr, cds->codes[i]);
    }
#endif

    /* Tidy up */
    free(buf);
    return cds;
}

/*
 * Reorders a ZTR file so that the chunks we wish to be in the common
 * header are at the start. The purpose of this is to allow us to write
 * out a header and footer separately for data saving reasons.
 * Header:
 *     ZTR magic
 *     HUFF chunks
 *     BPOS chunk
 *     CLIP chunk
 *     SMP4 chunk header (type codes, but not data itself)
 * Footer:
 *     SMP4 data
 *     BASE chunk
 *     CNF4 chunk
 *     ...
 */
static void reorder_ztr(ztr_t *ztr) {
    int i, j, cnum;
    ztr_chunk_t tmp;
    int headers[] = {ZTR_TYPE_HUFF, ZTR_TYPE_BPOS, ZTR_TYPE_CLIP,
		     ZTR_TYPE_SMP4};

    /* Reorder */
    for (cnum = j = 0; j < sizeof(headers)/sizeof(*headers); j++) {
	for (i = 0; i < ztr->nchunks; i++) {
	    if (ztr->chunk[i].type == headers[j]) {
		/* Swap */
		tmp = ztr->chunk[cnum];
		ztr->chunk[cnum++] = ztr->chunk[i];
		ztr->chunk[i] = tmp;
	    }
	}
    }
}

#ifdef USE_MODEL

#define MAX_BASES 74
static double model[MAX_BASES][4];
int model_count = 0;

void init_model() {
    int i, j;
    for (i = 0; i < MAX_BASES; i++) {
	for (j = 0; j < 4; j++) {
	    model[i][j] = 0;
	}
    }
}

/*
 * Normalises a trace to a fixed average intensity and builds up a
 * model trace from it.
 */
int add_to_model(ztr_t *ztr, ztr_chunk_t *chunk) {
    int i, j, sum[4];
    unsigned short *tr;
    int nb;
    double r[4];

    nb = (chunk->dlength-2)/8;

#if 1
    /* Compute average intensity over entire channel */
    tr = ((unsigned short *)chunk->data)+4;
    sum[0] = sum[1] = sum[2] = sum[3] = 0;
    for (i = 0; i < nb; i++) {
	for (j = 0; j < 4; j++) {
	    sum[j] += be_int2(*tr)-SHIFT_BY;
	    tr++;
	}
    }

    /* Turn this into a scale factor so our model averages 10000 */
    for (j = 0; j < 4; j++) {
	sum[j]/=nb;
	r[j] = 10000.0 / sum[j];

	/* Skip outliers */
	if (r[j] > 2000)
	    r[j] = 2000;
	if (r[j] < -2000)
	    r[j] = -2000;
    }
#else
    r[0] = r[1] = r[2] = r[3] = 1.0;
#endif

    /* Now sum this trace to the model */
    tr = ((unsigned short *)chunk->data)+4;
    for (i = 0; i < MAX_BASES && i < nb; i++) {
	for (j = 0; j < 4; j++) {
	    double val = (be_int2(*tr)-SHIFT_BY);
	    val *= r[j];
	    model[i][j] += val;
	    tr++;
	}
    }
    model_count++;
}

/*
 * The model so far is just a summation, so scale it down to make it a
 * true average of our trace data.
 */
void average_model(void) {
    int i, j;
    for (i = 0; i < MAX_BASES; i++) {
	for (j = 0; j < 4; j++) {
	    model[i][j] /= model_count;
	}
	/*
	printf("%2d %7.1f %7.1f %7.1f %7.1f\n",
	       i, model[i][0], model[i][1], model[i][2], model[i][3]);
	*/
    }
}

/*
 * Subtracts an appropriate scaled model trace from this one.
 */
int subtract_model(ztr_t *ztr, ztr_chunk_t *chunk) {
    int i, j, sum[4];
    unsigned short *tr;
    int nb;
    float r[4];

    nb = (chunk->dlength-2)/8;

#if 1
    /* Compute average intensity over entire channel */
    tr = ((unsigned short *)chunk->data)+4;
    sum[0] = sum[1] = sum[2] = sum[3] = 0;
    for (i = 0; i < nb; i++) {
	for (j = 0; j < 4; j++) {
	    sum[j] += be_int2(*tr)-SHIFT_BY;
	    tr++;
	}
    }

    /* Compute the scale factor for the model */
    for (j = 0; j < 4; j++) {
	sum[j]/=nb;
	r[j] = 10000.0 / sum[j];

	/* Skip outliers */
	if (r[j] > 2000)
	    r[j] = 2000;
	if (r[j] < -2000)
	    r[j] = -2000;
    }
#else
    r[0] = r[1] = r[2] = r[3] = 1.0;
#endif

    /* Quantise r[0] down to 8 bit data so we can estimate overhead */
    r[0] = (int)(log(r[0]) * 33 + 0.5);
    if (r[0] < 0)   r[0] = 0;
    if (r[0] > 255) r[0] = 255;
    r[0] = exp((r[0] - 0.5)/33.0);

    /*
     * FIXME: to reverse this we also need to store the 4 r[j] values
     * in the trace. Skipping for now to see the impact on compression.
     */

    /* And finally subtract the model */
    tr = ((unsigned short *)chunk->data)+4;
    for (i = 0; i < MAX_BASES && i < nb; i++) {
	for (j = 0; j < 4; j++) {
	    double val = (be_int2(*tr)-SHIFT_BY);
	    val -= model[i][j] / r[j];
	    if (j == 0) /* experiment - better model only of called peak? */
		*tr = be_int2((int)(val+SHIFT_BY));
	    tr++;
	}
    }
}

void subtract_models(Array za) {
    int i, j;
    int nreads = ArrayMax(za);

    for (i = 0; i < nreads; i++) {
	ztr_t *z = arr(ztr_t *, za, i);
	for (j = 0; j < z->nchunks; j++) {
	    if (z->chunk[j].type == ZTR_TYPE_SMP4) {
		subtract_model(z, &z->chunk[j]);
	    }
	}
    }
}
#endif

#define EBASE 65536
static double entropy16(unsigned short *data, int len) {
    double E[EBASE];
    double P[EBASE];
    double e;
    int i;
    
    for (i = 0; i < EBASE; i++)
        P[i] = 0;

    for (i = 0; i < len; i++)
        P[data[i]]++;

    for (i = 0; i < EBASE; i++) {
        if (P[i]) {
            P[i] /= len;
            E[i] = -(log(P[i])/log(EBASE));
        } else {
            E[i] = 0;
        }
    }

    for (e = i = 0; i < len; i++)
        e += E[data[i]];

    return e * log(EBASE)/log(256);
}

#define EBASE8 256
static double entropy8(unsigned char *data, int len) {
    double E[EBASE8];
    double P[EBASE8];
    double e;
    int i;
    
    for (i = 0; i < EBASE8; i++)
        P[i] = 0;

    for (i = 0; i < len; i++)
        P[data[i]]++;

    for (i = 0; i < EBASE8; i++) {
        if (P[i]) {
            P[i] /= len;
            E[i] = -(log(P[i])/log(EBASE8));
        } else {
            E[i] = 0;
        }
    }

    for (e = i = 0; i < len; i++)
        e += E[data[i]];

    return e * log(EBASE8)/log(256);
}

#define FD_TRACE_RAW 0
#define FD_TRACE_ORD 1
#define FD_TRACE_CMP 2
#define FD_QUAL_RAW 3
#define FD_QUAL_ORD 4
#define FD_QUAL_RLE 5
#define FD_QUAL_CMP 6
#define FD_BASE_RAW 7
#define FD_BASE_CMP 8
int fds[10];

/*
 * Applies compression to ZTR chunks in a Solexa appropriate manner.
 * There are two levels of compression here. The first is before entropy
 * encoding (level 0) and the second is to apply the final entropy encoding
 * (level 1).
 */
static void srf_compress_ztr(ztr_t *ztr, int level) {
    int i;
    for (i = 0; i < ztr->nchunks; i++) {
	switch(ztr->chunk[i].type) {
	case ZTR_TYPE_SMP4:
	    if (level == 0) {
		write(fds[FD_TRACE_RAW],
		      ztr->chunk[i].data+2, ztr->chunk[i].dlength-2);
		compress_chunk(ztr, &ztr->chunk[i], ZTR_FORM_TSHIFT, 0, 0);
		/* delta_called(ztr, &ztr->chunk[i]); */
#ifdef USE_MODEL
		add_to_model(ztr, &ztr->chunk[i]);
#endif
	    }
	    if (level == 1) {
		write(fds[FD_TRACE_ORD],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
		compress_chunk(ztr, &ztr->chunk[i],
			       ZTR_FORM_STHUFF, SMP4_CODE, 0);
		write(fds[FD_TRACE_CMP],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
	    }
	    break;

	case ZTR_TYPE_CNF4:
	    if (level == 0) {
		write(fds[FD_QUAL_RAW],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
		compress_chunk(ztr, &ztr->chunk[i], ZTR_FORM_QSHIFT, 0, 0);
		write(fds[FD_QUAL_ORD],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
		compress_chunk(ztr, &ztr->chunk[i], ZTR_FORM_XRLE2, 4, 0);
		write(fds[FD_QUAL_RLE],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
	    } else {
		compress_chunk(ztr, &ztr->chunk[i],
			       ZTR_FORM_STHUFF, CNF4_CODE, 0);
		write(fds[FD_QUAL_CMP],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
	    }
	    break;

	case ZTR_TYPE_BPOS:
	    if (level == 0) {
		compress_chunk(ztr, &ztr->chunk[i], ZTR_FORM_DELTA4, 1, 0);
		compress_chunk(ztr, &ztr->chunk[i], ZTR_FORM_32TO8,  0, 0);
		compress_chunk(ztr, &ztr->chunk[i], ZTR_FORM_RLE,   -1, 0);
	    } else {
		compress_chunk(ztr, &ztr->chunk[i], ZTR_FORM_ZLIB, Z_HUFFMAN_ONLY, 0);
	    }
	    break;

	case ZTR_TYPE_BASE:
	    if (level == 1) {
		write(fds[FD_BASE_RAW],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
		compress_chunk(ztr, &ztr->chunk[i],
		               ZTR_FORM_STHUFF, BASE_CODE, 0);
		write(fds[FD_BASE_CMP],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
	    }
	    break;
	}
    }
}

/*
 * Encodes a ztr file into an mFILE structure.
 * The 'footer' parameter, if supplied as non-NULL, is filled out with
 * the location 14 bytes into the SMP4 chunk or 0 if not found.
 *
 * returns NULL on failure
 */
mFILE *encode_ztr(ztr_t *ztr, int *footer) {
    mFILE *mf = mfcreate(NULL, 0);
    int pos, i;

    ztr_store_hcodes(ztr);
    reorder_ztr(ztr);
    ztr_mwrite_header(mf, &ztr->header);

    if (footer)
	*footer = 0;

    /* Write out chunks */
    for (i = 0; i < ztr->nchunks; i++) {
	pos = mftell(mf);
	ztr_mwrite_chunk(mf, &ztr->chunk[i]);
	if (ztr->chunk[i].type == ZTR_TYPE_SMP4 && footer) {
	    *footer = pos + 10; /* allows traces up to 64k */
	}
    }

    return mf;
}

/*
 * Takes a s_*_*_seq.txt file as input and creates a ZTR file as output.
 * It uses the associated prb and sig files too to do this.
 *
 * By default it uses the _sig.txt file, but if raw_mode is defined it'll
 * use _int.txt from the parent directory and if phased if defined it'll
 * use _sig2.txt for the phase-correct data. Specifying both raw and
 * phased mode has undefined behaviour.
 *
 * Phased indicates whether to use the _sig.txt or _sig2.txt file (with
 * _sig2.txt being the phase-corrected data).
 *
 * Assumption: within a single input file, all reads are the same length and
 * we're using unclipped data. We then utilise this to construct a common
 * header to reduce overhead of data.
 *
 * Skip indicates how many bases to remove from the front of the trace.
 * FIXME: Maybe this should not be remove, but mark as clipped instead.
 *
 * Returns: 0 for success
 *	   -1 for failure
 */
int append(srf_t *srf, char *seq_file, int raw_mode, int skip,
	   int phased, char *name_fmt, char *prefix_fmt) {
    char *cp;
    char prb_file[1024], sig_file[1024];
    FILE *fp_seq, *fp_prb, *fp_sig;
    char *seq;
    int seq_num = 0;
    char last_prefix[1024] = {'\0'};
    Array za, la;
    loc_t l;
    int nreads = 0;
    int t;
    huffman_codeset_t *base_cds = NULL;
    huffman_codeset_t *trace_cds = NULL;
    huffman_codeset_t *conf_cds = NULL;

    /* open all 3 filenames */
    if (NULL == (cp = strrchr(seq_file, '_')))
	return -1;

    sprintf(prb_file, "%.*s_prb.txt", (int)(cp-seq_file), seq_file);
    if (raw_mode)
	sprintf(sig_file, "../%.*s_int.txt", (int)(cp-seq_file), seq_file);
    else if (phased)
	sprintf(sig_file, "%.*s_sig2.txt", (int)(cp-seq_file), seq_file);
    else
	sprintf(sig_file, "%.*s_sig.txt", (int)(cp-seq_file), seq_file);

    if (NULL == (fp_seq = fopen(seq_file, "r"))) {
	return -1;
    }
    if (NULL == (fp_prb = fopen(prb_file, "r"))) {
	fclose(fp_seq);
	return -1;
    }
    if (NULL == (fp_sig = fopen(sig_file, "r"))) {
	fclose(fp_seq);
	fclose(fp_prb);
	return -1;
    }


    za = ArrayCreate(sizeof(ztr_t *), 0);
    la = ArrayCreate(sizeof(loc_t), 0);

    /* Fetch sequence */
    /*
     * Arguable we could claim that the first base is like a vector or adapter
     * sequence and should be included in the trace, but marked as clipped.
     * For now we simply omit this, but it's something to consider for the
     * next version.
     */

#ifdef USE_MODEL
    init_model();
#endif

    /* Cache all reads in memory */
    while ((seq = get_seq(fp_seq, skip, &l.lane, &l.tile, &l.x, &l.y))) {
	int   (*prb)[4] = get_prb(fp_prb, skip);
	float (*sig)[4] = get_sig(fp_sig, skip);
	Read *r;

	if (!prb || !sig) {
	    fprintf(stderr, "Couldn't load prb/sig for %s/%d\n",
		    seq_file, seq_num);
	    continue;
	}

	/* Create the ZTR file and output, if needed, the common header */
	if (NULL == (r = create_read(seq, prb, sig)))
	    continue;

	ARR(ztr_t *, za, nreads) = read2ztr(r);
	ARR(loc_t, la, nreads) = l;
	srf_compress_ztr(arr(ztr_t *, za, nreads), 0);
	read_deallocate(r);

	nreads++;
    }

#ifdef USE_MODEL
    average_model();
    subtract_models(za);
#endif

    /*
     * Compute averaged frequency metrics for HUFF encoding
     * For the 16-bit trace data we break this down into 1st and 2nd byte of
     * each value => trac1 and trac2 code sets.
     */
    base_cds  = ztr2codes(za, BASE_CODE, 1, ZTR_TYPE_BASE);

#ifdef SINGLE_HUFF
    conf_cds  = ztr2codes(za, CNF4_CODE, 1, ZTR_TYPE_CNF4);
    trace_cds = ztr2codes(za, SMP4_CODE, 1, ZTR_TYPE_SMP4);
#else
    conf_cds  = ztr2codes(za, CNF4_CODE, 4,         ZTR_TYPE_CNF4);
    trace_cds = ztr2codes(za, SMP4_CODE, N_TR_CODE, ZTR_TYPE_SMP4);
#endif

    /* Output traces */
    for (seq_num = 0; seq_num < nreads; seq_num++) {
	char name[1024], prefix[1024];
	int footer;
	mFILE *mf;
	srf_trace_body_t tb;
	ztr_t *z = arr(ztr_t *, za, seq_num);
	l = arr(loc_t, la, seq_num);

	ztr_add_hcode(z, base_cds, 0);
	ztr_add_hcode(z, trace_cds, 0);
	ztr_add_hcode(z, conf_cds, 0);

	srf_compress_ztr(z, 1);
	mf = encode_ztr(z, &footer);
	delete_ztr(z);

	format_name(prefix, prefix_fmt, l.lane, l.tile, l.x, l.y, seq_num);
	if (strcmp(prefix, last_prefix)) {
	    /* Prefix differs, so generate a new data block header */
	    srf_trace_hdr_t th;
	    strcpy(last_prefix, prefix);
	    srf_construct_trace_hdr(&th, prefix, mf->data, footer);
	    if (-1 == srf_write_trace_hdr(srf, &th))
		return -1;
	}

	/* Write out the variable element of a ZTR file */
	format_name(name, name_fmt, l.lane, l.tile, l.x, l.y, seq_num);
	srf_construct_trace_body(&tb, name, mf->data+footer, mf->size-footer);
	if (-1 == srf_write_trace_body(srf, &tb))
	    return -1;

	mfdestroy(mf);
    }

    huffman_codeset_destroy(base_cds);
    huffman_codeset_destroy(trace_cds);
    huffman_codeset_destroy(conf_cds);

    ArrayDestroy(za);
    ArrayDestroy(la);
    fclose(fp_seq);
    fclose(fp_prb);
    fclose(fp_sig);

    return 0;
}

int main(int argc, char **argv) {
    int i, ret = 0;
    FILE *outfp;
    int raw_mode = 0;
    char *outfn = "traces.srf";
    int skip = 0;
    int phased = 1;
    int quiet = 0;
    int dots = 0;
    char *name_fmt = "%x_%y";
    char *prefix_fmt = "s_%l_%t_";
    srf_t *srf;
    srf_cont_hdr_t *ch;

    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;
	} else if (!strcmp(argv[i], "-r")) {
	    raw_mode = 1;
	} else if (!strcmp(argv[i], "-p")) {
	    phased = 1;
	} else if (!strcmp(argv[i], "-u")) {
	    phased = 0;
	} else if (!strcmp(argv[i], "-q")) {
	    quiet = 1;
	} else if (!strcmp(argv[i], "-n")) {
	    name_fmt = argv[++i];
	} else if (!strcmp(argv[i], "-N")) {
	    prefix_fmt = argv[++i];
	} else if (!strcmp(argv[i], "-d")) {
	    dots = 1;
	} else if (!strcmp(argv[i], "-o")) {
	    if (i < argc) {
		outfn = argv[++i];
	    }
	} else if (!strcmp(argv[i], "-s")) {
	    if (i < argc) {
		skip = atoi(argv[++i]);
	    }
	} else {
	    fprintf(stderr, "Usage: solexa2ztr [-r] [-p] [-q] [-o file] "
		    "[-s n_bases] [-n name_format] [-N name_prefix_format] *_seq.txt ...\n");
	    fprintf(stderr, "Usage: solexa2ztr [options] *_seq.txt ...\n");
	    fprintf(stderr, "Options:\n");
	    fprintf(stderr, "    -r       Raw_mode (use ../*_int.txt)\n");
	    fprintf(stderr, "    -p       Phase corrected (use *_sig2.txt) - default\n");
	    fprintf(stderr, "    -u       Phase uncorrected (use *_sig.txt)\n");
	    fprintf(stderr, "    -d       Output 'dots' (default off)\n");
	    fprintf(stderr, "    -q       Quiet\n");
	    fprintf(stderr, "    -o file  Outputs to 'file' - default 'traces.srf'\n");
	    fprintf(stderr, "    -s sk    Skips 'sk' bases from start of each sequence - default 0\n");
	    fprintf(stderr, "    -n fmt   Format name suffix using the 'fmt' expansion rule as follows:\n");
	    fprintf(stderr, "             %%l = lane number\n");
	    fprintf(stderr, "             %%t = tile number\n");
	    fprintf(stderr, "             %%x = x coordinate\n");
	    fprintf(stderr, "             %%y = y coordinate\n");
	    fprintf(stderr, "             %%c = counter into file\n");
	    fprintf(stderr, "             %%L/%%T/%%X/%%Y/%%C = as above, but in hex.\n\n");
	    fprintf(stderr, "    -N fmt   Format name prefix; as above\n");
	    fprintf(stderr, "             eg \"-N Run10_%%l_%%t_ -n %%3X:%%3Y\"\n");
	    return 1;
	}
    }

    fds[FD_TRACE_RAW]=open("/tmp/srf_trace_raw",O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_TRACE_ORD]=open("/tmp/srf_trace_ord",O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_TRACE_CMP]=open("/tmp/srf_trace_cmp",O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_QUAL_RAW] =open("/tmp/srf_qual_raw", O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_QUAL_ORD] =open("/tmp/srf_qual_ord", O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_QUAL_RLE] =open("/tmp/srf_qual_rle", O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_QUAL_CMP] =open("/tmp/srf_qual_cmp", O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_BASE_RAW] =open("/tmp/srf_base_raw", O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_BASE_CMP] =open("/tmp/srf_base_cmp", O_CREAT|O_RDWR|O_TRUNC, 0666);

    if (NULL == (outfp = fopen(outfn, "wb"))) {
	perror(outfn);
    }

    /* SRF file header */
    srf = srf_create(outfp);
    ch = srf_construct_cont_hdr(NULL, "Bustard", "1.8.28");
    srf_write_cont_hdr(srf, ch);
    srf_destroy_cont_hdr(ch);

    /* Loop through remaining files appending to the archive */
    for (; i < argc; i++) {
	char c = '.';
	if (-1 == append(srf, argv[i], raw_mode, skip, phased, name_fmt,
			 prefix_fmt)) {
	    c = '!';
	    ret = 1;
	}
	if (!quiet && dots) {
	    putchar(c);
	    fflush(stdout);
	}
    }
    if (!quiet && dots)
	putchar('\n');

    srf_destroy(srf, 1);

    close(fds[FD_TRACE_RAW]);
    close(fds[FD_TRACE_ORD]);
    close(fds[FD_TRACE_CMP]);
    close(fds[FD_QUAL_RAW]);
    close(fds[FD_QUAL_ORD]);
    close(fds[FD_QUAL_RLE]);
    close(fds[FD_QUAL_CMP]);
    close(fds[FD_BASE_RAW]);
    close(fds[FD_BASE_CMP]);

    return ret;
}
