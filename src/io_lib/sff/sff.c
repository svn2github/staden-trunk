/*
 * Portions of this code have been derived from 454 Life Sciences Corporation's
 * getsff.c code (specifically the WriteSFFFile function.
 * It bears the following copyright notice:
 *
 * ------------------------------------------------------------
 * Copyright (c)[2001-2005] 454 Life Sciences Corporation. All Rights Reserved.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *
 * IN NO EVENT SHALL LICENSOR BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE.
 * ------------------------------------------------------------
 *
 * The remainder is Copyright Genome Research Limited (GRL) and is also
 * provided "AS IS" without warranty.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Read.h"
#include "xalloc.h"
#include "sff.h"

#define SFF_MAGIC_NUMBER (((unsigned int) '.') << 24) + \
                         (((unsigned int) 's') << 16) + \
                         (((unsigned int) 'f') << 8) + \
                         ((unsigned int) 'f')
#define SFF_VERSION "1.00"

/*
 * getuint2
 *
 * A function to convert a 2-byte TVF/SFF value into an integer.
 */
unsigned int getuint2(unsigned char *b)
{
    return ((unsigned int) b[0]) * 256 + ((unsigned int) b[1]);
}

/*
 * getuint4
 *
 * A function to convert a 4-byte TVF/SFF value into an integer.
 */
static unsigned int getuint4(unsigned char *b)
{
    return
        ((unsigned int) b[0]) * 256 * 256 * 256 +
        ((unsigned int) b[1]) * 256 * 256 +
        ((unsigned int) b[2]) * 256 +
        ((unsigned int) b[3]);
}


Read *sff2read(unsigned char *buf, size_t buflen) {
    int i;
    unsigned int signal, *flowSignals, *flowIndexes;
    unsigned int pos, magicNumber, numFlows, numBases, trimStart, trimEnd;
    unsigned int signalFormatCode, flowgramFormatCode;
    float *flowgram;
    char *bases, *flowChars, sffVersion[10];
    unsigned char *scores;
    Read *r = read_allocate(0,0);

    /*
     * Parse the SFF header bytes.
     */
    magicNumber = getuint4(buf);
    if (magicNumber != SFF_MAGIC_NUMBER) {
        fprintf(stderr, "Error:  Invalid SFF file format.\n");
	return NULL;
    }

    strncpy(sffVersion, (char *)buf+4, 4);
    sffVersion[4] = '\0';
    if (strcmp(sffVersion, SFF_VERSION) != 0) {
        fprintf(stderr, "Error: Unsupported SFF version:  %s\n", sffVersion);
	return NULL;
    }

    numFlows = getuint4(buf+8);

    signalFormatCode = getuint4(buf+12);
    flowgramFormatCode = getuint4(buf+16);

    flowChars = malloc(numFlows + 1);
    strncpy(flowChars, (char *)buf+20, numFlows);
    flowChars[numFlows] = '\0';

    /*
     * Parse out the SFF body values.
     */
    pos = 20 + numFlows;

    numBases = getuint4(buf + pos);
    pos += 4;
    trimStart = getuint2(buf + pos);
    pos += 2;
    trimEnd = getuint2(buf + pos);
    pos += 2;

    if (pos + numFlows * 4 + numBases * 4 != buflen) {
	fprintf(stderr, "Error:  Invalid SFF file format.\n");
	return NULL;
    }

    flowSignals = malloc(numFlows * sizeof(unsigned int));
    for (i=0; i < numFlows; i++) {
        flowSignals[i] = getuint2(buf + pos);
        pos += 2;
    }

    flowgram = malloc(numFlows * sizeof(float));
    for (i=0; i < numFlows; i++) {
        signal = getuint2(buf + pos);
        pos += 2;
        flowgram[i] = signal * 0.01;
    }

    flowIndexes = malloc(numBases * sizeof(unsigned int));
    for (i=0; i < numBases; i++) {
        flowIndexes[i] = getuint2(buf + pos);
        pos += 2;
    }

    bases = malloc(numBases + 1);
    scores = malloc(numBases);

    strncpy(bases, (char *)buf + pos, numBases);
    bases[numBases] = '\0';
    pos += numBases;

    strncpy((char *)scores, (char *)buf + pos, numBases);
    pos += numBases;

    /*
     * Output the SFF information.
     */
    r->nflows = numFlows;
    r->flow_order = flowChars;
    r->flow_raw = (unsigned int *)malloc(numFlows * sizeof(int));
    r->flow = (float *)malloc(numFlows * sizeof(float));
    
    r->leftCutoff = trimStart;
    r->rightCutoff = trimEnd;

    if (r->basePos) free(r->basePos);
    if (r->base)    free(r->base);
    if (r->prob_A)  free(r->prob_A);
    if (r->prob_C)  free(r->prob_C);
    if (r->prob_G)  free(r->prob_G);
    if (r->prob_T)  free(r->prob_T);

    r->NBases  = numBases;
    r->basePos = (uint_2 *)calloc(numBases, 2);
    r->base    = (char *)calloc(numBases, 1);
    r->prob_A  = (char *)calloc(numBases, 1);
    r->prob_C  = (char *)calloc(numBases, 1);
    r->prob_G  = (char *)calloc(numBases, 1);
    r->prob_T  = (char *)calloc(numBases, 1);

    for (i=0; i < numFlows; i++) {
	r->flow_raw[i] = flowSignals[i];
	r->flow[i] = flowgram[i];
    }

    for (i=0; i < numBases; i++) {
	r->base[i] = bases[i];
	r->prob_A[i] = 0;
	r->prob_C[i] = 0;
	r->prob_G[i] = 0;
	r->prob_T[i] = 0;
	switch (r->base[i]) {
	case 'A':
	case 'a':
	    r->prob_A[i] = scores[i];
	    break;
	case 'C':
	case 'c':
	    r->prob_C[i] = scores[i];
	    break;
	case 'G':
	case 'g':
	    r->prob_G[i] = scores[i];
	    break;
	case 'T':
	case 't':
	    r->prob_T[i] = scores[i];
	    break;
	}

	r->basePos[i] = flowIndexes[i];
    }

    free(flowSignals);
    free(flowgram);
    free(flowIndexes);
    free(bases);
    free(scores);

    return r;
}

Read *mfread_sff(mFILE *mf) {
    return sff2read((unsigned char *)mf->data, mf->size);
}
