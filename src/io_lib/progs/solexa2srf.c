/* TO-DO:
 *
 * Replace the use of io_lib's Read struct in here with direct control of
 * ztr_t objects instead. This avoids needlessly calling read2ztr and it
 * should speed up the code somewhat.
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
 * Updated April 2007:   Added a format string for naming samples.
 * Updated June 2007:    Wrapped the output format in SRF (example).
 * Updated January 2008: Support for chastity filters via quahog files.
 * Updated April 2008:   Added -mf/-mr params from Camil Toma
 *                       Support for calibrated fastq (-qf/-qr).
 *
 * This code converts Solexa (aka Illumina) format text files into ZTR
 * format traces held within an SRF container. It reads (optionally) from
 * the *sig2.txt, *prb.txt, *seq.txt, *qhg.txt, *nse.txt, *int.txt and
 * calibrarted fastq solexa files.
 */

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <fcntl.h>
#include <zlib.h>
#include <math.h>

#include <io_lib/Read.h>
#include <io_lib/misc.h>
#include <io_lib/ztr.h>
#include <io_lib/array.h>
#include <io_lib/srf.h>

#define S2S_VERSION "1.6"

/* Move to autoconf */
#define HAVE_POPEN

/* #define SINGLE_HUFF */
/* #define DEBUG_OUT */

/* #define NO_ENTROPY_ENCODING */

/* A model saves approx 2% data size, vs 1.3% saving for delta */
/* Neither are used in the official SRF spec though. */
/* #define USE_MODEL */

#define N_TR_CODE 8

#define BASE_CODE (CODE_USER)
#define SIG4_CODE (CODE_USER+1)
#define CNF4_CODE (CODE_USER+2)
#define INT4_CODE (CODE_USER+3)
#define NSE4_CODE (CODE_USER+4)
#define CNF1_CODE (CODE_USER+5)

/*
 * io_lib is still based around the Read struct at its heart rather than
 * a more general ZTR style, so we have to shoehorn things in via
 * the private_data block.
 */
#define SIG_SIG 0
#define SIG_INT 1
#define SIG_NSE 2
typedef struct {
    float (*signal[3])[4]; /* raw intensities */
    char *prc;             /* Calibrated single confidence values */
    int baseline[3];
} read_pd;


/* --- Some wrappers around FILE * vs gzFile *, allowing for either --- */
/*
 * gzopen() works on both compressed and uncompressed data, but it has
 * a significant performance hit even for uncompressed data (tested as
 * 25s using FILE* to 46s via gzOpen and 66s via gzOpen when gzipped).
 *
 * Hence we use our own wrapper 'zfp' which is a FILE* when uncompressed
 * and gzFile* when compressed. This also means we could hide bzopen in
 * there too if desired.
 */

/*
 * Either a gzFile or a FILE.
 */
typedef struct {
    FILE   *fp;
    gzFile *gz;
} zfp;

/*
 * A wrapper for either fgets or gzgets depending on what has been
 * opened.
 */
char *zfgets(char *line, int size, zfp *zf) {
    if (zf->fp)
	return fgets(line, size, zf->fp);
    else
	return gzgets(zf->gz, line, size);
}

/* A replacement for either fopen or gzopen */
zfp *zfopen(const char *path, const char *mode) {
    char path2[1024];
    zfp *zf;

    if (!(zf = (zfp *)malloc(sizeof(*zf))))
	return NULL;
    zf->fp = NULL;
    zf->gz = NULL;

    /* Try normal fopen */
    if (zf->fp = fopen(path, mode)) {
	unsigned char magic[2];
	fread(magic, 1, 2, zf->fp);
	if (!(magic[0] == 0x1f &&
	      magic[1] == 0x8b)) {
	    fseek(zf->fp, 0, SEEK_SET);
	    return zf;
	}

	fclose(zf->fp);
	zf->fp = NULL;
    }

#ifdef HAVE_POPEN
    /*
     * I've no idea why, by gzgets is VERY slow, maybe because it handles
     * arbitrary seeks.
     * popen to gzip -cd is 3 times faster though.
     */
    if (access(path, R_OK) == 0) {
	sprintf(path2, "gzip -cd < %.*s", 1000, path);
	if (zf->fp = popen(path2, "r"))
	    return zf;
    }

    sprintf(path2, "gzip -cd < %.*s.gz", 1000, path);
    if (zf->fp = popen(path2, "r"))
	return zf;

    printf("Failed on %s\n", path);
#else
    /* Gzopen instead */
    if (zf->gz = gzopen(path, mode))
	return zf;

    sprintf(path2, "%.*s.gz", 1020, path);
    if (zf->gz = gzopen(path2, mode))
	return zf;
#endif

    perror(path);

    free(zf);
    return NULL;
}

int zfclose(zfp *zf) {
    int r = (zf->fp) ? fclose(zf->fp) : gzclose(zf->gz);
    free(zf);
    return r;
}

/* --- */

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
    char date[99];
    char machine[99];
    int run;
    int x, y;
    int tile;
    int lane;
} loc_t;

#define MAX_CYCLES 120

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

static char *ztr_encode_float_4(ztr_t *z, char *type,  int npoints,
				float (*sig)[4], int baseline,
				int *nbytes, char **mdata, int *mdbytes) {
    char *bytes;
    int i, j, k, t;

    if ((z->header.version_major > 1 ||
	z->header.version_minor >= 2) && baseline) {
	/* 1.2 onwards */
	char buf[256];
	int blen;
	blen = sprintf(buf, "%d", baseline);
	if (type) {
	    *mdata = (char *)malloc(6+blen+6+strlen(type));
	    *mdbytes = sprintf(*mdata, "TYPE%c%s%cOFFS%c%s",
			       0, type, 0,
			       0, buf) + 1;
	} else {
	    *mdata = (char *)malloc(6+blen);
	    *mdbytes = sprintf(*mdata, "OFFS%c%s", 0, buf) + 1;
	}
    } else {
	if (type) {
	    *mdata = (char *)malloc(6+strlen(type));
	    *mdbytes = sprintf(*mdata, "TYPE%c%s",
			       0, type) + 1;
	} else {
	    *mdata = NULL;
	    *mdbytes = 0;
	}
    }

    bytes = (char *)xmalloc(npoints * sizeof(TRACE)*4 + 2);
    for (k = 0, j = 2; k < 4; k++) {
	for (i = 0; i < npoints; i++) {
	    t = (int)sig[i][k];
	    bytes[j++] = (t >> 8) & 0xff;
	    bytes[j++] = (t >> 0) & 0xff;
	}
    }
    *nbytes = 4 * npoints * sizeof(TRACE) + 2;

    bytes[0] = ZTR_FORM_RAW;
    bytes[1] = 0;
    return bytes;
}

int add_smp4_chunk(ztr_t *z, char *type, int npoints, float (*sig)[4],
		   int baseline) {
    ztr_chunk_t *zc;
    char *data, *mdata;
    int dlen, mdlen;
    
    z->chunk = (ztr_chunk_t *)realloc(z->chunk,
				      ++z->nchunks * sizeof(ztr_chunk_t));

    zc = &z->chunk[z->nchunks-1];
    data = ztr_encode_float_4(z, type, npoints, sig, baseline,
			     &dlen, &mdata, &mdlen);
    zc->type = ZTR_TYPE_SMP4;
    zc->mdlength = mdlen;
    zc->mdata    = mdata;
    zc->dlength  = dlen;
    zc->data     = data;
    zc->ztr_owns = 1;

    return 0;
}

static char *ztr_encode_qual(ztr_t *z, int nbase, char *qual, char *scale,
			    int *nbytes, char **mdata, int *mdbytes) {
    char *bytes;
    int i;

    if (NULL == (bytes = xmalloc(1 + nbase)))
	return NULL;

    bytes[0] = ZTR_FORM_RAW;
    memcpy(&bytes[1], qual, nbase);

    /* Meta data, if appropriate */
    if (scale) {
	*mdbytes = strlen(scale) + 7;
	*mdata = (char *)malloc(*mdbytes);
	sprintf(*mdata, "SCALE%c%s", 0, scale);
    } else {
	*mdbytes = 0;
	*mdata = NULL;
    }

    *nbytes = 1+nbase;

    return bytes;
}

int add_cnf1_chunk(ztr_t *z, int nbase, char *qual, char *scale) {
    ztr_chunk_t *zc;
    char *data, *mdata;
    int dlen, mdlen;
    
    z->chunk = (ztr_chunk_t *)realloc(z->chunk,
				      ++z->nchunks * sizeof(ztr_chunk_t));

    zc = &z->chunk[z->nchunks-1];
    data = ztr_encode_qual(z, nbase, qual, scale, &dlen, &mdata, &mdlen);
    zc->type = ZTR_TYPE_CNF1;
    zc->mdlength = mdlen;
    zc->mdata    = mdata;
    zc->dlength  = dlen;
    zc->data     = data;
    zc->ztr_owns = 1;

    return 0;
}

int set_cnf4_metadata(ztr_t *z, char *scale) {
    int i;

    /* Find the chunk type */
    for (i = 0; i < z->nchunks; i++) {
	if (z->chunk[i].type == ZTR_TYPE_CNF4) {
	    break;
	}
    }
    if (i == z->nchunks)
	return -1;

    /* Reset the meta-data */
    if (z->chunk[i].mdata)
	free(z->chunk[i].mdata);
    z->chunk[i].mdlength = strlen(scale) + 7;
    z->chunk[i].mdata = (char *)malloc(z->chunk[i].mdlength);
    sprintf(z->chunk[i].mdata, "SCALE%c%s", 0, scale);    
}

/*
 * Adds a ZTR REGN chunk describing the paired-end structure. Ie which bit
 * is which. This is a simplified form of the more generic REGN contents.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int add_readpair_region(ztr_t *z, unsigned int cycle) {
    char *mdata = malloc(100);
    unsigned char *data = malloc(5);
    int mdlen;

    if (!data || !mdata)
	return -1;

    data[0] = 0;
    cycle--; /* we count from 0 */
    data[1] = (cycle >> 24) & 0xff;
    data[2] = (cycle >> 16) & 0xff;
    data[3] = (cycle >>  8) & 0xff;
    data[4] = (cycle >>  0) & 0xff;
    
    mdlen = sprintf(mdata, "NAME%cforward:P;reverse:P%c",0,0);
    return ztr_new_chunk(z, ZTR_TYPE_REGN, data, 5, mdata, mdlen) ? 0 : -1;
}

char *parse_4_int(char *str, int *val) {
    int minus = 0, count = 0;
    char c;
    enum state_t {BEFORE_NUM, IN_NUM} state = BEFORE_NUM;
    int ival = 0;

    val[0] = val[1] = val[2] = val[3] = 0;

    do {
	c = *str++;
	switch (state) {
	case BEFORE_NUM:
	    while (isspace(c))
	        c = *str++;
	    if (!c)
		return NULL;
	    state = IN_NUM;

	case IN_NUM:
	    if (c == '-' || c == '+') {
		minus = c == '-';
		c = *str++;
	    }
	    
	    while (c >= '0' && c <= '9') {
		ival = ival * 10 + c-'0';
		c = *str++;
	    }

	    switch(c) {
	    case '\n': case '\r': case '\t': case ' ': case '\0':
		if (minus)
		    ival = -ival;

		minus = 0;
		*val++ = ival;
		ival = 0;

		if (++count == 4) {
		    return str;
		} else {
		    state = BEFORE_NUM;
		}
		break;

	    default:
		goto error;
	    }
	    break;
	}

    } while (c);

    /*
     * If we get here we've run out of input before processing 4 values
     * (the while command terminated) or we've found malformed data and
     * jumped here directly.
     */
    
 error:
    fprintf(stderr, "Error: unexpected character '%c' during parsing\n", c);
    return NULL;
}

/*
 * Fast lookup table to convert integer fraction value to a floating point
 * value. This allows floats to be parsed to 2 decimal points with no
 * logs or divisions anywhere. Over-optimising? Maybe, but I have the bit
 * between my teeth and am running with it :-)
 */
static float f_lookup[100] = {
    0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,
    0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
    0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29,
    0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39,
    0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
    0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59,
    0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69,
    0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79,
    0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89,
    0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
};

char *parse_4_float(char *str, float *val, int *bin)
{
    int minus = 0, count = 0;
    char c;
    enum state_t {BEFORE_NUM, BEFORE_POINT, AFTER_POINT} state = BEFORE_NUM;
    double fval = 0;
    int ival1 = 0, ival2 = 0;

    val[0] = val[1] = val[2] = val[3] = 0;

    do {
	c = *str++;
	switch (state) {
	case BEFORE_NUM:
	    while (isspace(c))
	        c = *str++;
	    if (!c)
		return NULL;
	    state = BEFORE_POINT;

	case BEFORE_POINT:
	    if (c == '-' || c == '+') {
		minus = c == '-';
		c = *str++;
	    }

	    while (c >= '0' && c <= '9') {
		ival1 = ival1 * 10 + c-'0';
		c = *str++;
	    }

	    if (c == '.') {
		state = AFTER_POINT;
		break;
	    }

	    switch(c) {
	    case '\n': case '\r': case '\t': case ' ': case '\0':
		if (minus)
		    ival1 = -ival1;

		minus = 0;
		*val++ = ival1;

		if (bin) {
		    if (ival1 > 65535)
			bin[65535]++;
		    else if (ival1 < -65535)
			bin[-65535]++;
		    else
			bin[ival1]++;
		}

		ival1 = 0;

		if (++count == 4) {
		    return str;
		} else {
		    state = BEFORE_NUM;
		}
		break;

	    default:
		goto error;
	    }
	    break;
	
	case AFTER_POINT:
	    while (c >= '0' && c <= '9') {
		ival2 = ival2 * 10 + c-'0';
		c = *str++;
	    }

	    while (ival2 >= 100)
		ival2 /= 10;

	    fval = ival1 + f_lookup[ival2];
		
	    if (minus) {
		fval = -fval;
		ival1 = -ival1;
	    }

	    minus = 0;
	    *val++ = fval;

	    if (bin) {
		if (ival1 > 65535)
		    bin[65535]++;
		else if (ival1 < -65535)
		    bin[-65535]++;
		else
		    bin[ival1]++;
	    }

	    ival1 = ival2 = 0;

	    if (++count == 4) {
		return str;
	    } else {
		state = BEFORE_NUM;
	    }
	    break;
	}

    } while (c);

    /*
     * If we get here we've run out of input before processing 4 values
     * (the while command terminated) or we've found malformed data and
     * jumped here directly.
     */
    
 error:
    fprintf(stderr, "Error: unexpected character '%c' during parsing\n", c);
    return NULL;
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
char *get_seq(zfp *fp, int trim, int *lane, int *tile, int *x, int *y) {
    static char line[1024];
    char *cp, *end;
    int i4[4];

    if (NULL == zfgets(line, 1023, fp))
	return NULL;

    /* First 4 values */
    cp = parse_4_int(line, i4);
    *lane = i4[0];
    *tile = i4[1];
    *x    = i4[2];
    *y    = i4[3];

    /* Followed by the sequence */
    while (isspace(*cp))
	cp++;

    end = cp;
    while (*end && !isspace(*end))
	end++;
    *end = 0;

    return cp + trim;
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
int (*get_prb(zfp *fp, int trim))[4] {
    char line[MAX_CYCLES*20 +1];
    static int prb[MAX_CYCLES][4];
    char *cp;
    int c = 0; /* cycle */
    
    if (NULL == zfgets(line, MAX_CYCLES*20, fp))
	return NULL;

    cp = line;
    while (cp = parse_4_int(cp, prb[c])) {
	c++;
    }

    return &prb[trim];
}

/*
 * Generates data as per get_prb above but extracting it from a fastq
 * sequence and quality string already loaded to memory instead of a
 * _prb.txt file.  The scores returned are calibrated in the phred
 * scale with one value per base. For compatibility we still return 4
 * values per cycle (the remaining 3 being set to zero), hence the
 * requirement for the sequence to be known too.
 *
 * Returns point to Nx4 array of ints holding the scores.
 *         NULL on failure
 */
/* UNUSED */
int (*fastq_to_prb4(char *seq, char *qual, int trim))[4] {
    size_t nc = strlen(seq), i;
    static int prb[MAX_CYCLES][4];
    int qc, qo;

    for (i = 0; i < nc; i++) {
	switch(seq[i]) {
	case 'a': case 'A':
	    qc = prb[i][0] = qual[i] - 64;
	    /* FIXME: use a precomputed lookup table instead */
	    qo = -10 * log(3*(1+pow(10, qc/10.0))-1)/log(10) - 0.5;
	    prb[i][1] = prb[i][2] = prb[i][3] = qo;
	    break;

	case 'c': case 'C':
	    qc = prb[i][1] = qual[i] - 64;
	    qo = -10 * log(3*(1+pow(10, qc/10.0))-1)/log(10) - 0.5;
	    prb[i][0] = prb[i][2] = prb[i][3] = qo;
	    break;

	case 'g': case 'G':
	    qc = prb[i][2] = qual[i] - 64;
	    qo = -10 * log(3*(1+pow(10, qc/10.0))-1)/log(10) - 0.5;
	    prb[i][0] = prb[i][1] = prb[i][3] = qo;
	    break;

	case 't': case 'T':
	    qc = prb[i][3] = qual[i] - 64;
	    qo = -10 * log(3*(1+pow(10, qc/10.0))-1)/log(10) - 0.5;
	    prb[i][0] = prb[i][1] = prb[i][2] = qo;
	    break;

	default:
	    prb[i][0] = prb[i][1] = prb[i][2] = prb[i][3] = -5; /* -4.77 */
	    break;

	}
    }

    return &prb[trim];
}

/*
 * As for fastq_to_prb4, but instead using 1 quality value scaled in
 * the phred scale instead of log-odds scale.
 */
char *fastq_to_prb1(char *seq, unsigned char *qual, int trim) {
    size_t nc = strlen(seq), i;
    int j;
    signed char *prb;
    static int qlookup[256];
    static int qlookup_done = 0;

    if (!qlookup_done) {
	for (j = 0; j < 255; j++) {
	    qlookup[j] = (int)((10*log(1+pow(10, (j-64)/10.0))/log(10)+.499));
	}
	qlookup_done = 1;
    }

    prb = (signed char *)malloc(nc);
    for (i = 0; i < nc; i++) {
	prb[i] = qlookup[qual[i]];
    }

    return &prb[trim];
}

/*
 * Fetches the signal strength arrays as an Nx4 array from the *_sig.txt file.
 * The format of the lines is a series of 4-tuples (4 floating point numbers,
 * themselves separated by spaced) which are separated by tabs.
 *
 * The data returned is statically allocated, so it should not be freed and
 * the function is not rentrant.
 *
 * If 'bin' is true then we build a histogram of signal values truncated at
 * -65535 to +65535
 *
 * Returns point to Nx4 array of ints holding the signal strengths.
 *         NULL on failure
 */
float (*get_sig(zfp *fp, int trim, int *bin))[4] {
    char line[MAX_CYCLES*40 +1];
    int i4[4];
    char *cp;
    int c = 0;
    float (*sig)[4];

    if (NULL == (sig = malloc(MAX_CYCLES * sizeof(*sig))))
	return NULL;

    if (NULL == zfgets(line, MAX_CYCLES*40, fp))
	return NULL;

    /* Skip first 4 values */
    cp = parse_4_int(line, i4);
    while (cp = parse_4_float(cp, sig[c], bin)) {
	c++;
    }

    /* NOTE: Free using free(&sig[-trim]) */
    return &sig[trim];
}

float (*get_ipar_sig(float **ipar_data, int trim, int cluster,
		     int numberOfCycles, int *bin))[4] {
    float (*sig)[4] = malloc(MAX_CYCLES * sizeof(*sig));
    int cycle;

    if (NULL == sig)
    {
        return NULL;
    }
    for (cycle = 0; cycle < numberOfCycles; ++cycle)
    {
        int channel;
        for (channel = 0; channel < 4; ++channel)
        {
            sig[cycle][channel] = ipar_data[cycle][cluster * 4 + channel];
            if (bin)
            {
                int ival1 = trunc(sig[cycle][channel]);
                if (ival1 > 65535)
                    bin[65535]++;
                else if (ival1 < -65535)
                    bin[-65535]++;
                else
                    bin[ival1]++;
            }
        }
    }
    /* NOTE: Free using free(&sig[-trim]) */
    return &sig[trim];
}

/*
 * Extracts and purity and chastity values from a quahog output file.
 * The format is:
 * lane tile x y purity chastity similarity distance neighbours
 *
 * Returns 0 for success
 *        -1 for failure
 */
int get_chastity(zfp *fp,
		 float *purity, float *chastity,
		 float *similarity, float *distance) {
    char line[MAX_CYCLES*30 +1];
    int i4[4];
    char *cp;
    float scores[4];

    if (NULL == zfgets(line, MAX_CYCLES*30, fp))
	return -1;

    /* Skip first 4 values */
    if (NULL == (cp = parse_4_int(line, i4)))
	return -1;
    if (NULL == parse_4_float(cp, scores, NULL))
	return -1;

    if (purity)     *purity     = scores[0];
    if (chastity)   *chastity   = scores[1];
    if (similarity) *similarity = scores[2];
    if (distance)   *distance   = scores[3];

    return 0;
}

/*
 * Extracts the next sequence/quality from a fastq file.
 *
 * Returns 0 for success
 *        -1 for failure
 */
int get_fastq_entry(zfp *fp, char *name, char *seq, char *qual) {
    char name1[80], name2[80], *cp;

    if (NULL == zfgets(name1, 80, fp))
	return -1;

    if (cp = strrchr(name1, '/'))
	*cp = 0;
    if (cp = strrchr(name1, '\n'))
	*cp = 0;

    /* Copy name fixing "-0". Why minus zero as a coordinate? */
    cp = name1+1;
    while (*cp) {
	if (*cp == '-' && *(cp+1) == '0')
	    cp++;
	*name++ = *cp++;
    }
    *name++ = 0;

    if (NULL == zfgets(seq, MAX_CYCLES, fp))
	return -1;
    if (NULL == zfgets(name2, 80, fp))
	return -1;
    if (NULL == zfgets(qual, MAX_CYCLES, fp))
	return -1;

    return 0;
}


/*
 * Creates an io_lib Read object from sequence, probabilities and signal
 * strengths.
 * Note: we keep the signal values as floats in the private segment of the
 * read object as we'll adjust them later to generate the traceA/C/G/T
 * arrays.
 *
 * Returns: allocated Read on success
 *	    NULL on failure
 */
Read *create_read(char *seq, int (*prb)[4], read_pd *pd) {
    size_t nbases = strlen(seq), i;
    Read *r;
    int max = 0;
    
    if (NULL == (r = read_allocate(nbases, nbases)))
	return NULL;

    for (i = 0; i < nbases; i++) {
	/* Confidence values */
	r->prob_A[i] = prb[i][0];
	r->prob_C[i] = prb[i][1];
	r->prob_G[i] = prb[i][2];
	r->prob_T[i] = prb[i][3];

	/* Traces */
	r->private_data = pd;
	r->private_size = sizeof(*pd);

	/* Sequence & position */
	r->base[i] = seq[i];
	r->basePos[i] = i;
    }

    if (!pd->signal[SIG_SIG]) {
	if (r->traceA) {
	    free(r->traceA);
	    r->traceA = NULL;
	}
	if (r->traceC) {
	    free(r->traceC);
	    r->traceC = NULL;
	}
	if (r->traceG) {
	    free(r->traceG);
	    r->traceG = NULL;
	}
	if (r->traceT) {
	    free(r->traceT);
	    r->traceT = NULL;
	}
	r->NPoints = 0;
    }

    r->maxTraceVal = max;
    r->leftCutoff = 0;
    r->rightCutoff = nbases+1;

    /* Info: none for present */

    return r;
}

void format_name(char *name, char *fmt, loc_t *l, int count) {
    int n;

    for (; *fmt; fmt++) {
	switch(*fmt) {

	/* A format specifier */
	case '%': {
	    char *endp;
	    int l1 = 0;
	    char *fmt2;

	    /* Precision specifier */
	    l1 = strtol(++fmt, &endp, 10);
	    fmt = endp;
	    fmt2 = isupper(*fmt) ? "%0*x" : "%0*d";

	    switch(tolower(*fmt)) {
	    case 'd':
		n = sprintf(name, "%.*s", l1 ? l1 : 99, l->date);
		break;
	    case 'm':
		n = sprintf(name, "%.*s", l1 ? l1 : 99, l->machine);
		break;
	    case 'r':
		n = sprintf(name, fmt2, l1, l->run);
		break;
	    case 'l':
		n = sprintf(name, fmt2, l1, l->lane);
		break;
	    case 't':
		n = sprintf(name, fmt2, l1, l->tile);
		break;
	    case 'x':
		n = sprintf(name, fmt2, l1, l->x);
		break;
	    case 'y':
		n = sprintf(name, fmt2, l1, l->y);
		break;
	    case 'c':
		n = sprintf(name, fmt2, l1, count);
		break;
	    default:
		fprintf(stderr, "Invalid %% rule '%%%c'\n", *fmt);
		n = 0;
	    }
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
 * key/value are optional (may be NULL). When specified they further
 * restrict the chunks to process to those containing a specific meta-data
 * key/value pair.
 * A key/value pair with NULL value is taken as implying that the key must
 * not exist.
 *
 * Returns huffman_code_t on success
 *         NULL on failure
 */
huffman_codeset_t *ztr2codes(Array za, int code_set, int nc, int type,
			     char *key, char *value) {
    unsigned char *buf = NULL;
    size_t sz = 0, alloc = 0;
    int nreads = ArrayMax(za);
    int i, j;
    huffman_codeset_t *cds;
 
    /* Accumulated concatenated blocks of data from all ZTR files */
    for (i = 0; i < nreads; i++) {
	ztr_t *z = arr(ztr_t *, za, i);
	for (j = 0; j < z->nchunks; j++) {
	    int len;

	    if (z->chunk[j].type != type)
		continue;

	    /* Check meta-data if required */
	    if (key && value) {
		char *v = ztr_lookup_mdata_value(z, &z->chunk[j], key);
		if (!v || 0 != strcmp(v, value))
		    continue;
	    } else if (key) {
		if (ztr_lookup_mdata_value(z, &z->chunk[j], key))
		    continue;
	    }

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
 *     REGN chunk
 *     TEXT chunk (matrix, parameters, etc)
 *     SMP4 chunk header (type codes, but not data itself)
 * Footer:
 *     SMP4 data
 *     BASE chunk
 *     CNF4 chunk
 *     CNF1 chunk
 *     ...
 */
static void reorder_ztr(ztr_t *ztr) {
    int i, j, cnum;
    ztr_chunk_t tmp;
    int headers[] = {
	/* header */
	ZTR_TYPE_HUFF, ZTR_TYPE_BPOS, ZTR_TYPE_CLIP,
	ZTR_TYPE_REGN, ZTR_TYPE_TEXT,
	/* footer */
	ZTR_TYPE_SMP4, ZTR_TYPE_BASE, ZTR_TYPE_CNF4,
	ZTR_TYPE_CNF1
    };

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
	    sum[j] += be_int2(*tr);
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
	    double val = (be_int2(*tr));
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
	    sum[j] += be_int2(*tr);
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
	    double val = (be_int2(*tr));
	    val -= model[i][j] / r[j];
	    if (j == 0) /* experiment - better model only of called peak? */
		*tr = be_int2((int)(val));
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

/*
 * We have a histogram in 'bin' of values from -65536 to +65535.
 * From this we work out the optimal baseline value that trims the least
 * amount of data.
 *
 * Returns: baseline value.
 */
int compute_baseline(int quiet, int *bin) {
    int min_val = -65536, min_count = 0;
    int max_val = +65536, max_count = 0;

    /* Find used extents */
    while (!bin[++min_val]);
    while (!bin[--max_val]);

    /* Trim off both ends until extents fit in range */
    while (max_val - min_val > 65535) {
	if (max_count <= min_count) {
	    max_count += bin[max_val];
	    while (!bin[--max_val]);
	} else {
	    min_count += bin[min_val];
	    while (!bin[++min_val]);
	}
    }

    if (!quiet)
	printf("Trimmed %d min, %d max => dynamic range is %d..%d inclusive\n",
	       min_count, max_count, min_val, max_val);

    return -min_val;
}

/*
 * Applies a baseline correction to an array of signal values, trimming
 * anything remaining outside of 0 to 65535.
 */
void rescale_trace(Read *r, int sig_num, int baseline) {
    int i, mtv = 0;
    read_pd *pd = (read_pd *)r->private_data;
    float (*sig)[4] = NULL;

    if (!(sig = pd->signal[sig_num]))
	return;

    for (i = 0; i < r->NBases; i++) {
	sig[i][0] += baseline;
	sig[i][1] += baseline;
	sig[i][2] += baseline;
	sig[i][3] += baseline;

	if (sig[i][0] < 0) sig[i][0] = 0;
	if (sig[i][1] < 0) sig[i][1] = 0;
	if (sig[i][2] < 0) sig[i][2] = 0;
	if (sig[i][3] < 0) sig[i][3] = 0;

	if (sig[i][0] > 65535) sig[i][0] = 65535;
	if (sig[i][1] > 65535) sig[i][1] = 65535;
	if (sig[i][2] > 65535) sig[i][2] = 65535;
	if (sig[i][3] > 65535) sig[i][3] = 65535;

	if (sig_num != SIG_SIG)
	    continue;

	r->traceA[i] = (int)sig[i][0];
	r->traceC[i] = (int)sig[i][1];
	r->traceG[i] = (int)sig[i][2];
	r->traceT[i] = (int)sig[i][3];

	if (mtv < sig[i][0]) mtv = sig[i][0];
	if (mtv < sig[i][1]) mtv = sig[i][1];
	if (mtv < sig[i][2]) mtv = sig[i][2];
	if (mtv < sig[i][3]) mtv = sig[i][3];
    }

    pd->baseline[sig_num] = baseline;
    if (sig_num == SIG_SIG) {
	r->baseline = baseline;
	r->maxTraceVal = mtv;
    }
}

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

#ifdef DEBUG_OUT
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
#endif

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
	    char *key;
	case ZTR_TYPE_SMP4:
	    key = ztr_lookup_mdata_value(ztr, &ztr->chunk[i], "TYPE");
	    if (level == 0) {
#ifdef DEBUG_OUT
		if (key && 0 == strcmp(key, "SLXN"))
		write(fds[FD_TRACE_RAW],
		      ztr->chunk[i].data+2, ztr->chunk[i].dlength-2);
#endif
		if (!key || 0 != strcmp(key, "SLXN"))
		    if (-1 == compress_chunk(ztr, &ztr->chunk[i],
					     ZTR_FORM_TSHIFT, 0, 0))
			exit(1);
		/* delta_called(ztr, &ztr->chunk[i]); */
#ifdef USE_MODEL
		add_to_model(ztr, &ztr->chunk[i]);
#endif
	    }
	    if (level == 1) {
#ifdef DEBUG_OUT
		if (key && 0 == strcmp(key, "SLXN"))
		write(fds[FD_TRACE_ORD],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
#endif
#ifndef NO_ENTROPY_ENCODING
		if (key && 0 == strcmp(key, "SLXI"))
		    if (-1 == compress_chunk(ztr, &ztr->chunk[i],
					     ZTR_FORM_STHUFF, INT4_CODE, 0))
			exit(2);
		else if (key && 0 == strcmp(key, "SLXN"))
		    if (-1 == compress_chunk(ztr, &ztr->chunk[i],
					     ZTR_FORM_STHUFF, NSE4_CODE, 0))
			exit(3);
		else
		    if (-1 == compress_chunk(ztr, &ztr->chunk[i],
					     ZTR_FORM_STHUFF, SIG4_CODE, 0))
			exit(4);
#endif
#ifdef DEBUG_OUT
		if (key && 0 == strcmp(key, "SLXN"))
		write(fds[FD_TRACE_CMP],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
#endif
	    }
	    break;

	case ZTR_TYPE_CNF4:
	    /* uncalibrated 4x confidence values */
	    if (level == 0) {
#ifdef DEBUG_OUT
		write(fds[FD_QUAL_RAW],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
#endif
		if (-1 == compress_chunk(ztr, &ztr->chunk[i],
					 ZTR_FORM_QSHIFT, 0, 0))
		    exit(5);
#ifdef DEBUG_OUT
		write(fds[FD_QUAL_ORD],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
#endif
		if (-1 == compress_chunk(ztr, &ztr->chunk[i],
					 ZTR_FORM_XRLE2, 4, 0))
		    exit(6);
#ifdef DEBUG_OUT
		write(fds[FD_QUAL_RLE],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
#endif
	    } else {
#ifndef NO_ENTROPY_ENCODING
		if (-1 == compress_chunk(ztr, &ztr->chunk[i],
					 ZTR_FORM_STHUFF, CNF4_CODE, 0))
		    exit(7);
#endif
#ifdef DEBUG_OUT
		write(fds[FD_QUAL_CMP],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
#endif
	    }
	    break;

	case ZTR_TYPE_BPOS:
	    /*
	     * Skip compression of BPOS as it leads to memory leaks and
	     * inefficiencies due to being part of the block header.
	     *
	    if (level == 0) {
		compress_chunk(ztr, &ztr->chunk[i], ZTR_FORM_DELTA4, 1, 0);
		compress_chunk(ztr, &ztr->chunk[i], ZTR_FORM_32TO8,  0, 0);
		compress_chunk(ztr, &ztr->chunk[i], ZTR_FORM_XRLE2,   1, 0);
	    }
	    */
	    break;

	case ZTR_TYPE_BASE:
	    if (level == 1) {
#ifdef DEBUG_OUT
		write(fds[FD_BASE_RAW],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
#endif
#ifndef NO_ENTROPY_ENCODING
		if (-1 == compress_chunk(ztr, &ztr->chunk[i],
					 ZTR_FORM_STHUFF, BASE_CODE, 0))
		    exit(8);
#endif
#ifdef DEBUG_OUT
		write(fds[FD_BASE_CMP],
		      ztr->chunk[i].data, ztr->chunk[i].dlength);
#endif
	    }
	    break;

	case ZTR_TYPE_CNF1:
	    /* calibrated 1x confidence values */
	    if (level == 1) {
#ifndef NO_ENTROPY_ENCODING
		if (-1 == compress_chunk(ztr, &ztr->chunk[i],
					 ZTR_FORM_STHUFF, CNF1_CODE, 0))
		    exit(10);
#endif
	    }
	    break;
	}
    }
}

/*
 * Encodes a ztr file into an mFILE structure.
 * The 'footer' parameter, if supplied as non-NULL, is filled out with
 * the location several bytes into the SMP4 chunk or 0 if not found.
 *
 * If no_hcodes is true then we omit storing the interleaved deflate huffman
 * codes. In this case footer will be incorrect, but this is done for speed
 * purposes.
 *
 * returns NULL on failure
 */
mFILE *encode_ztr(ztr_t *ztr, int *footer, int no_hcodes) {
    mFILE *mf = mfcreate(NULL, 0);
    int pos, i;

    if (!no_hcodes)
	ztr_store_hcodes(ztr);
    reorder_ztr(ztr);
    ztr_mwrite_header(mf, &ztr->header);

    if (footer)
	*footer = 0;

    /* Write out chunks */
    for (i = 0; i < ztr->nchunks; i++) {
	pos = mftell(mf);
	ztr_mwrite_chunk(mf, &ztr->chunk[i]);
	if ((ztr->chunk[i].type == ZTR_TYPE_SMP4 ||
	     ztr->chunk[i].type == ZTR_TYPE_BASE) && footer && !*footer) {
	    /* allows meta-data up to 64k */
	    *footer = pos + 6;
	}
    }

    return mf;
}

/*
 * Slurps the entirety of a file into a malloced buffer and returns a pointer
 * to it.
 *
 * Returns: malloced buffer on success, *lenp equal to length
 *          NULL on failure
 */
static char *load(char *fn, int *lenp) {
    char *data = NULL;
    int dsize = 0;
    int dcurr = 0, len;
    int fd = 0;

    if (fn) {
	if (-1 == (fd = open(fn, O_RDONLY, 0))) {
	    perror(fn);
	    return NULL;
	}
    }

    do {
	if (dsize - dcurr < 8192) {
	    dsize = dsize ? dsize * 2 : 8192;
	    if (NULL == (data = realloc(data, dsize))) {
		if (fd)
		    close(fd);
		return NULL;
	    }
	}

	len = read(fd, data + dcurr, 8192);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
	if (fd)
	    close(fd);
	return NULL;
    }

    if (fd)
	close(fd);

    /* nul terminate; but not included in length */
    if (dsize - dcurr < 1) {
	dsize++;
	if (NULL == (data = realloc(data, dsize)))
	    return NULL;
    }
    data[dcurr] = 0;

    if (lenp)
	*lenp = dcurr;

    return data;
}

/*
 * Parameters:
 * - ipar_file (in): the name of the file, without the .gz extension
 *   (_int.txt.p or _nse.txt.p)
 * - bin (in): the bin used in parse_4_floats
 * - numberOfChannels (out): read from the IPAR file (supported only for
 *   4 channels)
 * - numberOfClusters (out): read from the IPAR file
 * - numberOfCycles   (out): inferred from reading the IPAR file
 *
 * Returns an array of cycle data. Each cycle data stores all the
 *            intensities for all the available clusters.
 * Returns NULL (and numberOfCycles == 0) on failure.
 */
#define BUFFER_SIZE (MAX_CYCLES*100)
float **load_ipar_data(char *ipar_file, int *bin, int *numberOfChannels,
		       int *numberOfClusters, int *numberOfCycles)
{
    zfp *file = NULL;
    float **result = NULL;
    char buffer[BUFFER_SIZE];
    char *getsResult = NULL;
    int resultCapacity   = 0;

    *numberOfChannels = 0;
    *numberOfClusters = 0;
    *numberOfCycles   = 0;

    if (NULL == (file = zfopen(ipar_file, "r")))
    {
        char errorMessage[1000];
        sprintf(errorMessage, "Failed to open %s.gz: errno: %d",
		ipar_file, errno);
        perror(errorMessage);
        return result;
    }

    /* get the layout of the data */
    if (NULL != (getsResult = zfgets(buffer, BUFFER_SIZE, file)) &&
        2 == sscanf(buffer, "#CH%d:OBJ%d",numberOfChannels,numberOfClusters) &&
        4 == *numberOfChannels)
    {
	int loadError;
        char errorMessage[4000];

        getsResult = zfgets(buffer, BUFFER_SIZE, file);

        /* skip comment lines and empty lines */
        while (NULL != getsResult &&
	       ('#' == buffer[0] || '\n' == buffer[0] || '\0' == buffer[0]))
        {
            getsResult = zfgets(buffer, BUFFER_SIZE, file);
        }

	/* we are expecting at least one cycle */
        loadError = (NULL == getsResult);
        errorMessage[0] = '\0';

        /* iterate over all the available cycles */
        while (!loadError && NULL != getsResult)
        {
	    float *currentCycle;
            unsigned int cluster;

            /* resize the result if necessary */
            if ( (*numberOfCycles) == resultCapacity)
            {
		float **tmp;

                resultCapacity += MAX_CYCLES;
		tmp = (float **)malloc(resultCapacity * sizeof(*result));

                if (NULL == tmp)
                {
                    loadError = 1;
                    strcpy(errorMessage, "Failed to allocate memory while "
			   "loading IPAR data");
                    getsResult = NULL;
                    break;
                }
                else
                {
                    unsigned int cycle = 0;
                    while (cycle < resultCapacity) 
                    {
                        tmp[cycle] = cycle < (*numberOfCycles)
			    ? result[cycle]
			    : NULL;
                        ++cycle;
                    }
                    free(result);
                    result = tmp;
                }
            }
            /* allocate the storage for the current cycle */
            currentCycle = malloc( (*numberOfChannels) *
				   (*numberOfClusters) *
				   sizeof(*currentCycle ));
            if (NULL == currentCycle)
            {
                loadError = 1;
                strcpy(errorMessage, "Failed to allocate memory while "
		       "loading IPAR data");
                getsResult = NULL;
                break;
            }

            /* store the data for all the clusters */
            for (cluster = 0; cluster < *numberOfClusters; ++cluster)
            {
		float *vals = currentCycle + ( (*numberOfChannels) * cluster);
                if (NULL == parse_4_float(getsResult, vals, NULL))
                {
                    free(currentCycle);
                    loadError = 1;
                    sprintf(errorMessage, "Failed to parse 4 floats from: %s",
			    ipar_file);
                    getsResult = NULL;
                    break;
                }
                getsResult = zfgets(buffer, BUFFER_SIZE, file);
            }

            /* skip the comment lines and empty lines */
            while (NULL != getsResult
		   && ('#'  == buffer[0] ||
		       '\n' == buffer[0] ||
		       '\0' == buffer[0]))
            {
                getsResult = zfgets(buffer, BUFFER_SIZE, file);
            }

            result[*numberOfCycles] = currentCycle;
            ++(*numberOfCycles);
        }

        if (loadError)
        {
            unsigned int cycle = 0;
            for (cycle = 0; cycle < (*numberOfCycles); ++cycle)
            {
                free (result[cycle]);
                result[cycle] = NULL;
            }
            free (result);
            result = NULL;
            resultCapacity = 0;
            perror(errorMessage);
            *numberOfCycles   = 0;
        }
    }
    else
    {
        char s[1000] = "load_ipar_data: unknown error";
        if (NULL == getsResult)
            sprintf(s, "couldn't read IPAR file %s", ipar_file);
        else if (0 == (*numberOfChannels) * (*numberOfClusters) )
            sprintf(s, "%s: invalid layout: %s", ipar_file, buffer);
        else if (4 != (*numberOfChannels))
            sprintf(s, "this application supports only 4 channels: %d",
		    *numberOfChannels);
        perror(s);
    }

    zfclose(file);
    return result;
}

/*
 * Looks at the current working directory in an attempt to guess the
 * date and run name from the run root.
 */
void parse_directory(loc_t *l) {
    char cwd[8192], *cp, *underscore;

    /* Defaults for when we fail */
    strcpy(l->date, "?");
    strcpy(l->machine, "?");

    if (NULL == getcwd(cwd, 8192))
	return;

    /* hop up three dirs from Bustard to root */
    if ((cp = strrchr(cwd, '/')))
	*cp = 0;
    if ((cp = strrchr(cwd, '/')))
	*cp = 0;
    if ((cp = strrchr(cwd, '/')))
	*cp = 0;
    if (!(cp = strrchr(cwd, '/')))
	return;

    /* cp is now the root dir. name */
    cp++;
    if ((underscore = strchr(cp, '_'))) {
	*underscore = 0;
	strncpy(l->date, cp, 99);
	cp = underscore+1;
	if ((underscore = strchr(cp, '_'))) {
	    *underscore = 0;
	    strncpy(l->machine, cp, 99);
	    cp = underscore+1;
	    l->run = atoi(cp);
	}
    }
}

/*
 * Takes a fastq name from gerald s_*_sequence.txt and extracts the lane
 * and tile number.
 *
 * Returns 0 (setting *lane, *tile) on success
 *        -1 on failure
 */
int get_tile_from_name(char *name, int *lane, int *tile) {
    char *cp;

    if (NULL == (cp = strchr(name, ':')))
	return -1;
    cp++;
    *lane = atoi(cp);

    if (NULL == (cp = strchr(cp, ':')))
	return -1;
    cp++;
    *tile = atoi(cp);

    return 0;
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
 * Returns: 0 written for success
 *	   -1 for failure
 */
int append(srf_t *srf, char *seq_file, char *fwd_fastq, char *rev_fastq,
	   int ipar, int raw, int proc, int skip,
	   int phased, float chastity, int quiet, int rev_cycle,
	   char *name_fmt, char *prefix_fmt,
	   char *matrix_f_name, char *matrix_r_name,
	   int *nr, int *nf) {
    char *cp, *matrix1 = NULL, *matrix2 = NULL, *params = NULL;
    char prb_file[1024], sig_file[1024], qhg_file[1024];
    char int_file[1024], nse_file[1024];
    char ipar_int_file[1024], ipar_nse_file[1024];
    zfp *fp_seq = NULL, *fp_prb = NULL;
    zfp *fp_sig = NULL, *fp_qhg = NULL;
    zfp *fp_int = NULL, *fp_nse = NULL;
    zfp *fp_qf  = NULL, *fp_qr  = NULL;
    char *seq, *slash;
    int seq_num = 0;
    char last_prefix[1024] = {'\0'};
    Array za = NULL, la = NULL, ra = NULL;
    loc_t l;
    int nreads = 0, filtered = 0;
    int err = -1;
    huffman_codeset_t *seq_cds = NULL;
    huffman_codeset_t *prb_cds = NULL;
    huffman_codeset_t *prc_cds = NULL;
    huffman_codeset_t *sig_cds = NULL;
    huffman_codeset_t *int_cds = NULL;
    huffman_codeset_t *nse_cds = NULL;
    int sig_bin_a[65536*2], *sig_bin = &sig_bin_a[65536];
    int int_bin_a[65536*2], *int_bin = &int_bin_a[65536];
    int nse_bin_a[65536*2], *nse_bin = &nse_bin_a[65536];
    int sig_baseline, int_baseline, nse_baseline;
    int last_lane = 0;
    int next_fastq = 1;
    char full_fmt[1024], fastq_name[1024];
    float **ipar_int_data = NULL;
    float **ipar_nse_data = NULL;
    int numberOfCycles; /* only for ipar */

    sprintf(full_fmt, "%s%s", prefix_fmt, name_fmt);

    memset(sig_bin_a, 0, 65536*2 * sizeof(*sig_bin_a));
    memset(int_bin_a, 0, 65536*2 * sizeof(*int_bin_a));
    memset(nse_bin_a, 0, 65536*2 * sizeof(*nse_bin_a));

    /* Handle full pathnames by cd to the directory first */
    if (slash = strrchr(seq_file, '/')) {
	*slash = 0;
	chdir(seq_file);
	seq_file = slash+1;
    }

    /* put run date, machine and run number in 'l' */
    parse_directory(&l);

    /* open the required files */
    if (NULL == (cp = strrchr(seq_file, '_')))
	return -1;

    sprintf(prb_file, "%.*s_prb.txt", (int)(cp-seq_file), seq_file);
    sprintf(qhg_file, "%.*s_qhg.txt", (int)(cp-seq_file), seq_file);
    sprintf(int_file, "../%.*s_int.txt", (int)(cp-seq_file), seq_file);
    sprintf(nse_file, "../%.*s_nse.txt", (int)(cp-seq_file), seq_file);
    sprintf(ipar_int_file, "../%.*s_int.txt.p", (int)(cp-seq_file), seq_file);
    sprintf(ipar_nse_file, "../%.*s_nse.txt.p", (int)(cp-seq_file), seq_file);

    if (phased) {
	sprintf(sig_file, "%.*s_sig2.txt", (int)(cp-seq_file), seq_file);
    } else {
	sprintf(sig_file, "%.*s_sig.txt", (int)(cp-seq_file), seq_file);
    }

    if (NULL == (fp_seq = zfopen(seq_file, "r")))
	goto error;
    if (NULL == (fp_prb = zfopen(prb_file, "r")))
	goto error;

    if (proc) {
	if (NULL == (fp_sig = zfopen(sig_file, "r")))
	    goto error;
    }

    if (raw) {
        if (ipar)
        {
            int numberOfChannels, numberOfClusters;
	    ipar_int_data = load_ipar_data(ipar_int_file,
					   int_bin,
					   &numberOfChannels,
					   &numberOfClusters,
					   &numberOfCycles);
	    ipar_nse_data = load_ipar_data(ipar_nse_file,
					   nse_bin,
					   &numberOfChannels,
					   &numberOfClusters,
					   &numberOfCycles);

            if (NULL == ipar_int_data || NULL == ipar_nse_data)
                goto error;
        }
        else
        {
            if (NULL == (fp_int = zfopen(int_file, "r")))
	        goto error;
            if (NULL == (fp_nse = zfopen(nse_file, "r")))
	        goto error;
        }
    }

    if (chastity > 0) {
	if (NULL == (fp_qhg = zfopen(qhg_file, "r")))
	    goto error;
    }

    if (fwd_fastq) {
	if (NULL == (fp_qf = zfopen(fwd_fastq, "r")))
	    goto error;
    }

    if (rev_fastq) {
	if (NULL == (fp_qr = zfopen(rev_fastq, "r")))
	    goto error;
    }

    za = ArrayCreate(sizeof(ztr_t *), 0);
    la = ArrayCreate(sizeof(loc_t), 0);
    ra = ArrayCreate(sizeof(Read *), 0);

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

    params = load("../.params", NULL);

    /* Cache all reads in memory */
    while ((seq = get_seq(fp_seq, skip, &l.lane, &l.tile, &l.x, &l.y))) {
	read_pd *pd = (read_pd *)calloc(1, sizeof(*pd));
	int   (*prb)[4];
	Read *r;
	float c;
	char fastq_seq[MAX_CYCLES+1];
	char fastq_qual[MAX_CYCLES+1];
	int eof = 0;

	/*
	 * If we're storing calibrated GERALD quality values instead of
	 * Bustard base-call output then we need to read the fwd and rev
	 * fastq files.
	 *
	 * There's probably a different number of entries in these though
	 * than there is in the *seq.txt (etc) files, so we need to skip
	 * entries. Logically speaking we'd be looping over the fastq
	 * data and then searching for the next sequence to match in
	 * the *seq.txt files, but the history of this code is such
	 * that the loop is the other way around. Hence for the odd
	 * use of a next_fastq variable.
	 */
	if (fp_qf || fp_qr) {
	    char name[1024];

	    format_name(name, full_fmt, &l, seq_num);

	    while (next_fastq) {
		int c = 0;
		int fq_lane, fq_tile;

		if (fp_qf) {
		    if (-1 == get_fastq_entry(fp_qf, fastq_name,
					      fastq_seq,
					      fastq_qual)) {
			eof=1;
			break;
		    }
		    c = strlen(fastq_seq);
		}
		if (fp_qr) {
		    if (-1 == get_fastq_entry(fp_qr, fastq_name,
					      &fastq_seq[c],
					      &fastq_qual[c])) {
			eof=1;
			break;
		    }
		}

		if (-1 == get_tile_from_name(fastq_name, &fq_lane, &fq_tile)) {
		    fprintf(stderr, "Couldn't parse name '%s'\n", fastq_name); 
		}
		if (l.lane == fq_lane && l.tile == fq_tile)
		    next_fastq = 0;
	    }

	    if (strcmp(fastq_name, name)) {
		char dummy[10240];
		next_fastq = 0;
		
		/* Skip our other files to keep all in sync */
		if (fp_qhg) zfgets(dummy, 10240, fp_qhg);
		if (fp_prb) zfgets(dummy, 10240, fp_prb);
		if (fp_sig) zfgets(dummy, 10240, fp_sig);
		if (fp_int) zfgets(dummy, 10240, fp_int);
		if (fp_nse) zfgets(dummy, 10240, fp_nse);

		filtered++;
		if (pd->signal[SIG_SIG]) free(&pd->signal[SIG_SIG][-skip]);
		if (pd->signal[SIG_INT]) free(&pd->signal[SIG_INT][-skip]);
		if (pd->signal[SIG_NSE]) free(&pd->signal[SIG_NSE][-skip]);
		if (pd->prc) free(&pd->prc[-skip]);
		free(pd);
		continue;

	    } else {

		next_fastq = 1;
	    }
	}

	if (eof)
	    break;

	if (fp_qhg) {
	    get_chastity(fp_qhg, NULL, &c, NULL, NULL);
	    if (c < chastity) {
		char dummy[10240];

		/*
		 * To keep in sync we still need to read from the other files,
		 * but we don't need to parse them. So we duplicate the work
		 * a bit here for efficiencies sake.
		 */
		if (fp_prb) zfgets(dummy, 10240, fp_prb);
		if (fp_sig) zfgets(dummy, 10240, fp_sig);
		if (fp_int) zfgets(dummy, 10240, fp_int);
		if (fp_nse) zfgets(dummy, 10240, fp_nse);

		filtered++;
		if (pd->signal[SIG_SIG]) free(&pd->signal[SIG_SIG][-skip]);
		if (pd->signal[SIG_INT]) free(&pd->signal[SIG_INT][-skip]);
		if (pd->signal[SIG_NSE]) free(&pd->signal[SIG_NSE][-skip]);
		if (pd->prc) free(&pd->prc[-skip]);
		free(pd);

		continue;
	    }
	}

	/* Load the approprate matrix file whenever changing lanes */
	if (l.lane != last_lane) {
	    char fn[100];

	    if (matrix1)
		free(matrix1);
	    sprintf(fn, matrix_f_name, l.lane);
	    last_lane = l.lane;
	    matrix1 = load(fn, NULL);

	    if (matrix2)
		free(matrix2);
	    if (rev_cycle) {
		sprintf(fn, matrix_r_name,
			l.lane, rev_cycle+1);
		last_lane = l.lane;
		matrix2 = load(fn, NULL);
	    }
	}

	/*
	 * We *either* read calibrated quality from the fastq file *or*
	 * from the prb.txt files.
	 */
	if (fp_qf || fp_qr) {
	    if (!(pd->prc = fastq_to_prb1(fastq_seq, fastq_qual, skip))) {
		fprintf(stderr, "Failed to generate prb from fastq\n");
		continue;
	    }
	} else {
	    pd->prc = NULL;
	}

	if (!(prb = get_prb(fp_prb, skip))) {
	    fprintf(stderr, "Couldn't load prb for %s/%d\n",
		    seq_file, seq_num);
	    continue;
	}


	if (fp_sig) {
	    if (!(pd->signal[SIG_SIG] = get_sig(fp_sig, skip, sig_bin))) {
		fprintf(stderr, "Couldn't load sig for %s/%d\n",
			seq_file, seq_num);
		continue;
	    }
	}

	if (fp_int) {
	    if (!(pd->signal[SIG_INT] = get_sig(fp_int, skip, int_bin))) {
		fprintf(stderr, "Couldn't load int for %s/%d\n",
			seq_file, seq_num);
		continue;
	    }
	}

        if (ipar_int_data) {
	    if (!(pd->signal[SIG_INT] = get_ipar_sig(ipar_int_data, skip,
						     nreads, numberOfCycles,
						     int_bin))) {
		fprintf(stderr, "Couldn't allocate memory for %s/%d\n",
			ipar_int_file, nreads);
		continue;
	    }            
        }

	if (fp_nse) {
	    if (!(pd->signal[SIG_NSE] = get_sig(fp_nse, skip, nse_bin))) {
		fprintf(stderr, "Couldn't load nse for %s/%d\n",
			seq_file, seq_num);
		continue;
	    }
	}


        if (ipar_nse_data) {
	    if (!(pd->signal[SIG_NSE] = get_ipar_sig(ipar_nse_data, skip,
						     nreads, numberOfCycles,
						     nse_bin))) {
		fprintf(stderr, "Couldn't allocate memory for %s/%d\n",
			ipar_nse_file, nreads);
		continue;
	    }            
        }

	/* Create the ZTR file and output, if needed, the common header */
	if (NULL == (r = create_read(seq, prb, pd)))
	    continue;

	ARR(Read *, ra, nreads) = r;
	ARR(loc_t, la, nreads) = l;

	nreads++;
    }

    if (nf)
	*nf = filtered;
    if (nr)
	*nr = nreads;

    if (!nreads)
	goto skip;

    if (fp_sig)
	sig_baseline = compute_baseline(quiet, sig_bin);
    if (fp_int || ipar_int_data)
	int_baseline = compute_baseline(quiet, int_bin);
    if (fp_nse || ipar_nse_data)
	nse_baseline = compute_baseline(quiet, nse_bin);
    
    /*
     * Now we know the min_sig and max_sig we can turn our float arrays in
     * the Read struct into integers and convert to ZTR.
     */
    for (seq_num = 0; seq_num < nreads; seq_num++) {
	ztr_t *z;
	ztr_chunk_t *tc;
	Read *r = arr(Read *, ra, seq_num);
	read_pd *pd = (read_pd *)r->private_data;

	if (fp_sig)
	    rescale_trace(r, SIG_SIG, sig_baseline);
	if (fp_int || ipar_int_data)
	    rescale_trace(r, SIG_INT, int_baseline);
	if (fp_nse || ipar_nse_data)
	    rescale_trace(r, SIG_NSE, nse_baseline);

	/* Standard conversion - seq, qual, processed trace */
	z = ARR(ztr_t *, za, seq_num) = read2ztr(r);
	set_cnf4_metadata(z, "LO");

	/* REGN chunk if needed */
	if (rev_cycle)
	    add_readpair_region(z, rev_cycle);

	/* Plus raw traces if desired */
	if (pd->signal[SIG_INT])
	    add_smp4_chunk(z, "SLXI", r->NBases, pd->signal[SIG_INT],
			   pd->baseline[SIG_INT]);
	if (pd->signal[SIG_NSE])
	    add_smp4_chunk(z, "SLXN", r->NBases, pd->signal[SIG_NSE],
			   pd->baseline[SIG_NSE]);

	if (pd->prc)
	    add_cnf1_chunk(z, r->NBases, pd->prc, "PH");

	/* Calibrated confidence values too */


	/* Miscellaneous text fields - only stored once per header */
	tc = ztr_add_text(z, NULL, "PROGRAM_ID", "solexa2srf v" S2S_VERSION);
	if (NULL == tc) {
	    fprintf(stderr, "Failed to add to TEXT chunk\n");
	    return -1;
	}

	if (matrix1) {
	    if (NULL == ztr_add_text(z, tc, "SOLEXA_MATRIX_FWD", matrix1)) {
		fprintf(stderr, "Failed to add to TEXT chunk\n");
		return -1;
	    }
	}

	if (matrix2) {
	    if (NULL == ztr_add_text(z, tc, "SOLEXA_MATRIX_REV", matrix2)) {
		fprintf(stderr, "Failed to add to TEXT chunk\n");
		return -1;
	    }
	}

	if (chastity > 0) {
	    char chas[100];
	    sprintf(chas, "%f", chastity);
	    if (NULL == ztr_add_text(z, tc, "SOLEXA_CHASTITY", chas)) {
		fprintf(stderr, "Failed to add to TEXT chunk\n");
		return -1;
	    }
	}

	if (params) {
	    if (NULL == ztr_add_text(z, tc, "SOLEXA_PARAMS", params)) {
		fprintf(stderr, "Failed to add to TEXT chunk\n");
		return -1;
	    }
	}

	srf_compress_ztr(arr(ztr_t *, za, seq_num), 0);

	if (pd->signal[SIG_SIG]) free(&pd->signal[SIG_SIG][-skip]);
	if (pd->signal[SIG_INT]) free(&pd->signal[SIG_INT][-skip]);
	if (pd->signal[SIG_NSE]) free(&pd->signal[SIG_NSE][-skip]);
	if (pd->prc) free(&pd->prc[-skip]);
	read_deallocate(r);
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
#ifndef NO_ENTROPY_ENCODING
    seq_cds = ztr2codes(za, BASE_CODE, 1, ZTR_TYPE_BASE, NULL,   NULL);
    prc_cds = ztr2codes(za, CNF1_CODE, 1, ZTR_TYPE_CNF1, NULL,   NULL);

#ifdef SINGLE_HUFF
    prb_cds = ztr2codes(za, CNF4_CODE, 1, ZTR_TYPE_CNF4, NULL,   NULL);
    if (fp_sig)
	sig_cds = ztr2codes(za, SIG4_CODE, 1, ZTR_TYPE_SMP4, "TYPE", NULL);
    if (fp_int || ipar_int_data)
	int_cds = ztr2codes(za, INT4_CODE, 1, ZTR_TYPE_SMP4, "TYPE", "SLXI");
    if (fp_nse || ipar_nse_data)
	nse_cds = ztr2codes(za, NSE4_CODE, 1, ZTR_TYPE_SMP4, "TYPE", "SLXN");
#else
    prb_cds = ztr2codes(za, CNF4_CODE, 4,         ZTR_TYPE_CNF4, NULL, NULL);
    if (fp_sig)
	sig_cds = ztr2codes(za, SIG4_CODE, N_TR_CODE,
			    ZTR_TYPE_SMP4, "TYPE", NULL);
    if (fp_int || ipar_int_data)
	int_cds = ztr2codes(za, INT4_CODE, N_TR_CODE,
			    ZTR_TYPE_SMP4, "TYPE", "SLXI");
    if (fp_nse || ipar_nse_data)
	nse_cds = ztr2codes(za, NSE4_CODE, 2,
			    ZTR_TYPE_SMP4, "TYPE", "SLXN");
#endif
#endif

    /* Output traces */
    for (seq_num = 0; seq_num < nreads; seq_num++) {
	char name[1024], prefix[1024];
	int footer;
	mFILE *mf;
	srf_trace_body_t tb;
	ztr_t *z = arr(ztr_t *, za, seq_num);
	l = arr(loc_t, la, seq_num);

	/*
	 * FIXME: copy chunks from one ztr file to the next, rather than
	 * re-encoding the HUFF chunks over and over again
	 */

#ifndef NO_ENTROPY_ENCODING
	ztr_add_hcode(z, seq_cds, 0);
	if (prb_cds)
	    ztr_add_hcode(z, prb_cds, 0);
	if (prc_cds)
	    ztr_add_hcode(z, prc_cds, 0);
	if (sig_cds)
	    ztr_add_hcode(z, sig_cds, 0);
	if (int_cds)
	    ztr_add_hcode(z, int_cds, 0);
	if (nse_cds)
	    ztr_add_hcode(z, nse_cds, 0);
#endif

	srf_compress_ztr(z, 1);

	format_name(prefix, prefix_fmt, &l, seq_num);
	if (strcmp(prefix, last_prefix)) {
	    /* Prefix differs, so generate a new data block header */
	    srf_trace_hdr_t th;
	    strcpy(last_prefix, prefix);

	    mf = encode_ztr(z, &footer, 0);
	    srf_construct_trace_hdr(&th, prefix, (unsigned char *)mf->data,
				    footer);
	    if (-1 == srf_write_trace_hdr(srf, &th))
		return -1;
	} else {
	    mf = encode_ztr(z, &footer, 1);
	}

	delete_ztr(z);

	/* Write out the variable element of a ZTR file */
	format_name(name, name_fmt, &l, seq_num);
	srf_construct_trace_body(&tb, name, (unsigned char *)mf->data+footer,
				 mf->size-footer);
	if (-1 == srf_write_trace_body(srf, &tb))
	    return -1;

	mfdestroy(mf);
    }

    if (seq_cds) huffman_codeset_destroy(seq_cds);
    if (sig_cds) huffman_codeset_destroy(sig_cds);
    if (prb_cds) huffman_codeset_destroy(prb_cds);
    if (prc_cds) huffman_codeset_destroy(prc_cds);
    if (int_cds) huffman_codeset_destroy(int_cds);
    if (nse_cds) huffman_codeset_destroy(nse_cds);

 skip:
    err = 0;
 error:

    if (ra)
	ArrayDestroy(ra);
    if (za)
	ArrayDestroy(za);
    if (la)
	ArrayDestroy(la);

    if (fp_seq)
	zfclose(fp_seq);
    if (fp_prb)
	zfclose(fp_prb);
    if (fp_sig)
	zfclose(fp_sig);
    if (fp_int)
	zfclose(fp_int);
    if (fp_nse)
	zfclose(fp_nse);
    if (fp_qhg)
	zfclose(fp_qhg);
    if (fp_qf)
	zfclose(fp_qf);
    if (fp_qr)
	zfclose(fp_qr);
    
    if (matrix1)
	free(matrix1);
    if (matrix2)
	free(matrix2);
    if (params)
	free(params);

    if (ipar_int_data)
    {
        unsigned int cycle = 0;
        for (cycle = 0; cycle < numberOfCycles; ++cycle) 
        {
            free(ipar_int_data[cycle]);
            ipar_int_data[cycle] = NULL;
        }
        free(ipar_int_data);
        ipar_int_data = NULL;
    }

    if (ipar_nse_data)
    {
        unsigned int cycle = 0;
        for (cycle = 0; cycle < numberOfCycles; ++cycle)
        {
            free(ipar_nse_data[cycle]);
            ipar_nse_data[cycle] = NULL;
        }
        free(ipar_nse_data);
        ipar_nse_data = NULL;
    }

    return err;
}

void usage(int code) {
    fprintf(stderr, "Solexa2srf v" S2S_VERSION "\n\n");
    fprintf(stderr, "Usage: solexa2srf [options] *_seq.txt ...\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -I       IPAR format for raw data\n");
    fprintf(stderr, "    -r       Add raw data:        ../_*{int,nse}.txt\n");
    fprintf(stderr, "    -R       Skip raw data:       ../_*{int,nse}.txt - default\n");
    fprintf(stderr, "    -p       Add processed data:  *sig2.txt - default\n");
    fprintf(stderr, "    -P       Skip processed data: *sig2.txt\n");
    fprintf(stderr, "    -u       Use phase uncorrected data (sig vs sig2)\n");

    fprintf(stderr, "    -c float Chastity filter (>= float)\n");
    fprintf(stderr, "    -d       Output 'dots' (default off)\n");
    fprintf(stderr, "    -q       Quiet\n");
    fprintf(stderr, "    -o file  Outputs to 'file' - default 'traces.srf'\n");
    fprintf(stderr, "    -s sk    Skips 'sk' bases from start of each sequence - default 0\n");
    fprintf(stderr, "    -2 cyc   The cycle number for the 2nd pair\n");
    fprintf(stderr, "    -qf fn   Get forward read calibrated quality from fastq file 'fn'\n");
    fprintf(stderr, "             (By default *_prb.txt is used instead)\n");
    fprintf(stderr, "    -qr fn   Get reverse read calibrated quality from fastq file 'fn'\n");
    fprintf(stderr, "             (By default *_prb.txt is used instead)\n");
    fprintf(stderr, "    -n fmt   Format name suffix using the 'fmt' expansion rule as follows:\n");
    fprintf(stderr, "             %%d = date (from run pwd)\n");
    fprintf(stderr, "             %%m = machine name (from run pwd)\n");
    fprintf(stderr, "             %%r = run number (from run pwd)\n");
    fprintf(stderr, "             %%l = lane number\n");
    fprintf(stderr, "             %%t = tile number\n");
    fprintf(stderr, "             %%x = x coordinate\n");
    fprintf(stderr, "             %%y = y coordinate\n");
    fprintf(stderr, "             %%c = counter into file\n");
    fprintf(stderr, "             %%L/%%T/%%X/%%Y/%%C = as above, but in hex.\n");
    fprintf(stderr, "    -N fmt   Format name prefix; as above\n");
    fprintf(stderr, "             eg \"-N Run10_%%l_%%t_ -n %%3X:%%3Y\"\n");
    fprintf(stderr, "    -mf name Name of the 1st pair matrix file.  By default\n"
	            "             ../Matrix/s_%%d_02_matrix.txt, where %%d is the lane.\n");
    fprintf(stderr, "    -mr name Name of the 2ndt pair matrix file.  By default\n"
	            "             ../Matrix/s_%%d_%%02d_matrix.txt, where %%d is the lane\n"
	            "             and %%02d is the first cycle in the 2nd pair.\n");
    fprintf(stderr, "\n");

    exit(code);
}

int main(int argc, char **argv) {
    int i, ret = 0;
    FILE *outfp;
    int ipar_mode = 0, raw_mode = 0, proc_mode = 1;
    char *outfn = "traces.srf";
    int skip = 0;
    int phased = 1;
    int quiet = 0;
    int dots = 0;
    char *name_fmt = "%x:%y";
    char *prefix_fmt = "%m_%r:%l:%t:";
    srf_t *srf;
    srf_cont_hdr_t *ch;
    float chastity = 0;
    int nreads = 0, nfiltered = 0;
    int rev_cycle = 0;
    char *fwd_fastq = NULL, *rev_fastq = NULL;
    char *matrix_f_name = "../Matrix/s_%d_02_matrix.txt";
    char *matrix_r_name = "../Matrix/s_%d_%02d_matrix.txt";

    /* Parse args */
    for (i = 1; i < argc && argv[i][0] == '-'; i++) {
	if (!strcmp(argv[i], "-")) {
	    break;
	} else if (!strcmp(argv[i], "-I")) {
	    ipar_mode = 1;
	} else if (!strcmp(argv[i], "-r")) {
	    raw_mode = 1;
	} else if (!strcmp(argv[i], "-R")) {
	    raw_mode = 0;
	} else if (!strcmp(argv[i], "-p")) {
	    proc_mode = 1;
	} else if (!strcmp(argv[i], "-P")) {
	    proc_mode = 0;
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
	} else if (!strcmp(argv[i], "-2")) {
	    if (i < argc) {
		rev_cycle = atoi(argv[++i]);
	    }
	} else if (!strcmp(argv[i], "-c")) {
	    if (i < argc) {
		chastity = atof(argv[++i]);
	    }	
	} else if (!strcmp(argv[i], "-qf")) {
	    if (i < argc) {
		fwd_fastq = argv[++i];
	    }
	} else if (!strcmp(argv[i], "-qr")) {
	    if (i < argc) {
		rev_fastq = argv[++i];
	    }
	} else if (!strcmp(argv[i], "-mf")) {
	    if (i < argc) {
		matrix_f_name = argv[++i];
	    }
	} else if (!strcmp(argv[i], "-mr")) {
	    if (i < argc) {
		matrix_r_name = argv[++i];
	    }
	} else if (!strcmp(argv[i], "-h")) {
	    usage(0);
	} else {
	    usage(1);
	}
    }

    if (i == argc)
	usage(0);

#ifdef DEBUG_OUT
    fds[FD_TRACE_RAW]=open("/tmp/srf_trace_raw",O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_TRACE_ORD]=open("/tmp/srf_trace_ord",O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_TRACE_CMP]=open("/tmp/srf_trace_cmp",O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_QUAL_RAW] =open("/tmp/srf_qual_raw", O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_QUAL_ORD] =open("/tmp/srf_qual_ord", O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_QUAL_RLE] =open("/tmp/srf_qual_rle", O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_QUAL_CMP] =open("/tmp/srf_qual_cmp", O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_BASE_RAW] =open("/tmp/srf_base_raw", O_CREAT|O_RDWR|O_TRUNC, 0666);
    fds[FD_BASE_CMP] =open("/tmp/srf_base_cmp", O_CREAT|O_RDWR|O_TRUNC, 0666);
#endif

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
	int nr, nf;
	if (!quiet) {
	    printf("Processing tile %s\n", argv[i]);
	}
	if (-1 == append(srf, argv[i], fwd_fastq, rev_fastq, ipar_mode,
			 raw_mode, proc_mode, skip, phased, chastity,
			 quiet, rev_cycle, name_fmt, prefix_fmt, 
			 matrix_f_name, matrix_f_name, &nr, &nf)) {
	    c = '!';
	    ret = 1;
	} else {
	    nreads += nr;
	    nfiltered += nf;
	}
	if (!quiet && dots) {
	    putchar(c);
	    fflush(stdout);
	}
    }
    if (!quiet && dots)
	putchar('\n');

    srf_destroy(srf, 1);

#ifdef DEBUG_OUT
    close(fds[FD_TRACE_RAW]);
    close(fds[FD_TRACE_ORD]);
    close(fds[FD_TRACE_CMP]);
    close(fds[FD_QUAL_RAW]);
    close(fds[FD_QUAL_ORD]);
    close(fds[FD_QUAL_RLE]);
    close(fds[FD_QUAL_CMP]);
    close(fds[FD_BASE_RAW]);
    close(fds[FD_BASE_CMP]);
#endif

    if (!quiet) {
	printf("Filtered %8d traces\n", nfiltered);
	printf("Wrote    %8d traces\n", nreads);
	printf("=>       %6.2f%% passed\n", 100.0*nreads/(nfiltered+nreads));
    }

    return ret;
}
