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
#include "mach-io.h"
#include "misc.h"

void free_sff_common_header(sff_common_header *h) {
    if (!h)
	return;
    if (h->flow)
	xfree(h->flow);
    if (h->key)
	xfree(h->key);
    xfree(h);
}

sff_common_header *read_sff_common_header(mFILE *mf) {
    sff_common_header *h;

    if (NULL == (h = (sff_common_header *)xcalloc(1, sizeof(*h))))
	return NULL;

    if (!be_read_int_4(mf, &h->magic))
	return free_sff_common_header(h), NULL;
    if (4 != mfread(h->version, 1, 4, mf))
	return free_sff_common_header(h), NULL;
    if (!be_read_int_8(mf, &h->index_offset))
	return free_sff_common_header(h), NULL;
    if (!be_read_int_4(mf, &h->index_len))
	return free_sff_common_header(h), NULL;
    if (!be_read_int_4(mf, &h->nreads))
	return free_sff_common_header(h), NULL;
    if (!be_read_int_2(mf, &h->header_len))
	return free_sff_common_header(h), NULL;
    if (!be_read_int_2(mf, &h->key_len))
	return free_sff_common_header(h), NULL;
    if (!be_read_int_2(mf, &h->flow_len))
	return free_sff_common_header(h), NULL;
    if (!be_read_int_1(mf, &h->flowgram_format))
	return free_sff_common_header(h), NULL;

    if (h->magic != SFF_MAGIC) {
	xfree(h);
	return NULL;
    }
    
    if (NULL == (h->flow = (char *)xmalloc(h->flow_len)))
	return free_sff_common_header(h), NULL;
    if (NULL == (h->key  = (char *)xmalloc(h->key_len)))
	return free_sff_common_header(h), NULL;
    if (h->flow_len != mfread(h->flow, 1, h->flow_len, mf))
	return free_sff_common_header(h), NULL;
    if (h->key_len != mfread(h->key , 1, h->key_len,  mf))
	return free_sff_common_header(h), NULL;

    /* Pad to 8 chars */
    mfseek(mf, (mftell(mf) + 7)& ~7, SEEK_SET);

    return h;
}

void free_sff_read_header(sff_read_header *h) {
    if (!h)
	return;
    if (h->name)
	xfree(h->name);
    free(h);
}

sff_read_header *read_sff_read_header(mFILE *mf) {
    sff_read_header *h;

    if (NULL == (h = (sff_read_header *)xcalloc(1, sizeof(*h))))
	return NULL;

    if (!be_read_int_2(mf, &h->header_len))
	return free_sff_read_header(h), NULL;
    if (!be_read_int_2(mf, &h->name_len))
	return free_sff_read_header(h), NULL;
    if (!be_read_int_4(mf, &h->nbases))
	return free_sff_read_header(h), NULL;
    if (!be_read_int_2(mf, &h->clip_qual_left))
	return free_sff_read_header(h), NULL;
    if (!be_read_int_2(mf, &h->clip_qual_right))
	return free_sff_read_header(h), NULL;
    if (!be_read_int_2(mf, &h->clip_adapter_left))
	return free_sff_read_header(h), NULL;
    if (!be_read_int_2(mf, &h->clip_adapter_right))
	return free_sff_read_header(h), NULL;
    if (NULL == (h->name  = (char *)xmalloc(h->name_len)))
	return free_sff_read_header(h), NULL;
    if (h->name_len != mfread(h->name, 1, h->name_len, mf))
	return free_sff_read_header(h), NULL;
    
    /* Pad to 8 chars */
    mfseek(mf, (mftell(mf) + 7)& ~7, SEEK_SET);

    return h;
}

void free_sff_read_data(sff_read_data *d) {
    if (!d)
	return;
    if (d->flowgram)
	xfree(d->flowgram);
    if (d->flow_index)
	xfree(d->flow_index);
    if (d->bases)
	xfree(d->bases);
    if (d->quality)
	xfree(d->quality);
    xfree(d);
}

sff_read_data *read_sff_read_data(mFILE *mf, int nflows, int nbases) {
    sff_read_data *d;
    int i;

    if (NULL == (d = (sff_read_data *)xcalloc(1, sizeof(*d))))
	return NULL;

    if (NULL == (d->flowgram = (uint16_t *)xcalloc(nflows, 2)))
	return free_sff_read_data(d), NULL;
    if (nflows != mfread(d->flowgram, 2, nflows, mf))
	return free_sff_read_data(d), NULL;
    for (i = 0; i < nflows; i++)
	d->flowgram[i] = be_int2(d->flowgram[i]);

    if (NULL == (d->flow_index = (uint8_t*)xmalloc(nbases)))
	return free_sff_read_data(d), NULL;
    if (nbases != mfread(d->flow_index, 1, nbases, mf))
	return free_sff_read_data(d), NULL;

    if (NULL == (d->bases = (uint8_t*)xmalloc(nbases)))
	return free_sff_read_data(d), NULL;
    if (nbases != mfread(d->bases, 1, nbases, mf))
	return free_sff_read_data(d), NULL;

    if (NULL == (d->quality = (uint8_t*)xmalloc(nbases)))
	return free_sff_read_data(d), NULL;
    if (nbases != mfread(d->quality, 1, nbases, mf))
	return free_sff_read_data(d), NULL;

    return d;
}


Read *mfread_sff(mFILE *mf) {
    int i, bpos;
    Read *r;
    sff_common_header *ch;
    sff_read_header *rh;
    sff_read_data *rd;

    /* Load the SFF contents */
    if (NULL == (ch = read_sff_common_header(mf)))
	return NULL;
    if (NULL == (rh = read_sff_read_header(mf))) {
	free_sff_common_header(ch);
	return NULL;
    }
    if (NULL == (rd = read_sff_read_data(mf, ch->flow_len, rh->nbases))) {
	free_sff_common_header(ch);
	free_sff_read_header(rh);
	return NULL;
    }

    /* Convert to Read struct */
    r = read_allocate(0,0);
    if (r->basePos) free(r->basePos);
    if (r->base)    free(r->base);
    if (r->prob_A)  free(r->prob_A);
    if (r->prob_C)  free(r->prob_C);
    if (r->prob_G)  free(r->prob_G);
    if (r->prob_T)  free(r->prob_T);

    r->nflows = ch->flow_len;
    r->flow_order = ch->flow; ch->flow = NULL;
    r->flow_raw = NULL;
    r->flow = (float *)malloc(r->nflows * sizeof(float));
    for (i = 0; i < r->nflows; i++) {
	r->flow[i] = rd->flowgram[i] / 100.0;
    }

    r->NBases = rh->nbases;
    r->basePos = (uint_2 *)calloc(r->NBases, 2);
    r->base    = rd->bases; rd->bases = NULL;
    r->prob_A  = (char *)calloc(r->NBases, 1);
    r->prob_C  = (char *)calloc(r->NBases, 1);
    r->prob_G  = (char *)calloc(r->NBases, 1);
    r->prob_T  = (char *)calloc(r->NBases, 1);

    bpos = 0;
    for (i=0; i < r->NBases; i++) {
	r->prob_A[i] = 0;
	r->prob_C[i] = 0;
	r->prob_G[i] = 0;
	r->prob_T[i] = 0;
	switch (r->base[i]) {
	case 'A':
	case 'a':
	    r->prob_A[i] = rd->quality[i];
	    break;
	case 'C':
	case 'c':
	    r->prob_C[i] = rd->quality[i];
	    break;
	case 'G':
	case 'g':
	    r->prob_G[i] = rd->quality[i];
	    break;
	case 'T':
	case 't':
	    r->prob_T[i] = rd->quality[i];
	    break;
	}

	bpos += rd->flow_index[i];
	r->basePos[i] = bpos;
    }

    r->leftCutoff = MAX(rh->clip_qual_left, rh->clip_adapter_left);
    r->rightCutoff = MIN(rh->clip_qual_right
			 ? rh->clip_qual_right
			 : r->NBases+1,
			 rh->clip_adapter_right
			 ? rh->clip_adapter_right
			 : r->NBases+1);

    free_sff_common_header(ch);
    free_sff_read_header(rh);
    free_sff_read_data(rd);

    return r;
}
