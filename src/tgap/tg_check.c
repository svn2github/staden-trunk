/*
 * Checks the validity of data in the database.
 *
 * These can be split into actual corruptions and credibility warnings.
 * Eg a sequence length not matching the start/end range_t is a fatal error,
 * but a sequence being clipped to zero bases is valid, just rather
 * unbelievable.
 */
#include <staden_config.h>
#include <string.h>
#include <stdlib.h>

#include "tg_gio.h"
#include "tg_check.h"

/* Enable debugging, which may help track down the location of some errors */
//#define DEBUG_CHECK

#define NORM(x) (f_a * (x) + f_b)
#define NMIN(x,y) (MIN(NORM((x)),NORM((y))))
#define NMAX(x,y) (MAX(NORM((x)),NORM((y))))

typedef struct {
    /* Contig start/end */
    int cstart;
    int cend;
    
    /* Item counts */
    int nseq;
    int nanno;
    int nref;
} bin_stats;

/*
 * Checks a sequence struct is consistent, both internally and with the
 * range that points to it.
 */
int check_seq(GapIO *io, bin_index_t *bin, range_t *r) {
    int err = 0;
    int i, len;
    seq_t *s = cache_search(io, GT_Seq, r->rec);

    static char isbase[256], init_done = 0;

    if (!init_done) {
	memset(isbase, 0, 256);
	isbase['A'] = isbase['a'] = 1;
	isbase['C'] = isbase['c'] = 1;
	isbase['G'] = isbase['g'] = 1;
	isbase['T'] = isbase['t'] = 1;
	isbase['N'] = isbase['*'] = 1;
	isbase['0'] = isbase['1'] = 1;
	isbase['2'] = isbase['3'] = 1;
	init_done = 1;
    }

    if (!s)
	return 1;

    if (r->end - r->start + 1 != ABS(s->len)) {
	fprintf(stderr, "Seq %"PRIrec": length does not match bin-range\n",
		s->rec);
	err++;
    }

    if (s->bin != bin->rec) {
	fprintf(stderr, "Seq %"PRIrec": bin does not match observed bin\n",
		s->rec);
	err++;
    }

    if (!bin->rng ||
	s->bin_index != r - ArrayBase(range_t, bin->rng)) {
	fprintf(stderr, "Seq %"PRIrec": bin_index does not match range "
		"index\n", s->rec);
	err++;
    }

    if (s->left < 1 || s->right > ABS(s->len)) {
	fprintf(stderr, "Seq %"PRIrec": left/right clips outside of sequence "
		"bounds.\n", s->rec);
	err++;
    }

    if (s->right < s->left) {
	fprintf(stderr, "Seq %"PRIrec": right clip starts before left clip\n",
		s->rec);
	err++;
    }

    if (s->mapping_qual != r->mqual) {
	fprintf(stderr, "Seq %"PRIrec": mapping_qual disagrees with range\n",
		s->rec);
	err++;
    }

    /* TODO: Check r->pair_rec? */

    /* TODO: Check r->flags match s->flags */

    if (s->seq_tech > STECH_454) {
	fprintf(stderr, "Seq %"PRIrec": Unknown seq_tech value\n",
		s->rec);
	err++;
    }

    if (s->format > 3) {
	fprintf(stderr, "Seq %"PRIrec": Unknown format value\n",
		s->rec);
	err++;
    }

    if (s->format != SEQ_FORMAT_CNF1) {
	fprintf(stderr, "Seq %"PRIrec": Warning - outdated seq format\n",
		s->rec);
    }

    /* Parent rec / type */
    if (s->parent_rec && !cache_exists(io, s->parent_type, s->parent_rec)) {
	fprintf(stderr, "Seq %"PRIrec": parent_rec/type do not match\n",
		s->rec);
	err++;
    }


    /* Check valid sequence characters */
    len = ABS(s->len);
    for (i = 0; i < len; i++) {
	if (!isbase[(unsigned char)s->seq[i]]) {
	fprintf(stderr, "Seq %"PRIrec": seq has unexpected char '%c' (%d)\n",
		s->rec, isprint(s->seq[i]) ? s->seq[i] : '?', s->seq[i]);
	    err++;
	    break;
	}
    }

    if ((s->name && s->name_len != strlen(s->name)) ||
	(!s->name && s->name_len != 0)) {
	fprintf(stderr, "Seq %"PRIrec": name and name_len do not match\n",
		s->rec);
	err++;
    }

    return err;
}

/*
 * Walks a contig bin tree, executing callbacks per bin.
 *
 * Returns 0 on success
 *         number of errors on failure
 */
int bin_walk(GapIO *io, int fix, tg_rec rec, int offset, int complement,
	     bin_stats *bs) {
    bin_index_t *bin;
    int i, f_a, f_b, err = 0;
    bin_stats child_stats;
    int start, end, cstart, cend, nthis_seq = 0;

    if (!rec)
	return 1;

    if (NULL == (bin = cache_search(io, GT_Bin, rec)))
	return 1;

    /* Set f_a & f_b for NMIN/NMAX macros & handle complementing */
    if (bin->flags & BIN_COMPLEMENTED) {
	complement ^= 1;
    }
    if (complement) {
	f_a = -1;
	f_b = offset + bin->size-1;
    } else {
	f_a = +1;
	f_b = offset;
    }

#ifdef DEBUG_CHECK
    printf("Bin %"PRIrec"  %d..%d  nseq=%d\n",
	   bin->rec,
	   NMIN(bin->pos, bin->pos + bin->size-1),
	   NMAX(bin->pos, bin->pos + bin->size-1),
	   bin->nseqs);
#endif

    cache_incr(io, bin);


    /* Recurse */
    for (i = 0; i < 2; i++) {
	bin_index_t *ch;
	if (!bin->child[i])
	    continue;
	ch = get_bin(io, bin->child[i]);

	child_stats.cstart = INT_MAX;
	child_stats.cend   = INT_MIN;
	child_stats.nseq   = 0;
	child_stats.nanno  = 0;
	child_stats.nref   = 0;

	err += bin_walk(io, fix, bin->child[i],
			NMIN(ch->pos, ch->pos + ch->size-1) /* offset */,
			complement, &child_stats);

	bs->nseq  += child_stats.nseq;
	bs->nanno += child_stats.nanno;
	bs->nref  += child_stats.nref;

	if (child_stats.nseq) {
	    if (bs->cstart > child_stats.cstart)
		bs->cstart = child_stats.cstart;
	    if (bs->cend   < child_stats.cend)
		bs->cend   = child_stats.cend;
	}
    }


    /* Check range items */
    start  = INT_MAX;
    end    = INT_MIN;
    cstart = INT_MAX;
    cend   = INT_MIN;
    if (bin->rng) {
	for (i = 0; i < ArrayMax(bin->rng); i++) {
	    range_t *r = arrp(range_t, bin->rng, i);

	    if (r->flags & GRANGE_FLAG_UNUSED)
		continue;

#ifdef DEBUG_CHECK
	    printf("Range item %d: %d..%d (abs %d..%d)\n",
		   i, r->start, r->end,
		   NMIN(r->start, r->end), NMAX(r->start, r->end));
#endif
	    if (start > r->start)
		start = r->start;
	    if (end   < r->end)
		end   = r->end;

	    switch (r->flags & GRANGE_FLAG_ISMASK) {
	    case GRANGE_FLAG_ISSEQ: {
		bs->nseq++;
		nthis_seq++;
		err += check_seq(io, bin, r);

		if (cstart > r->start)
		    cstart = r->start;
		if (cend   < r->end)
		    cend   = r->end;
		break;
	    }

	    case GRANGE_FLAG_ISANNO:
		bs->nanno++;
		break;

	    case GRANGE_FLAG_ISREFPOS:
		bs->nref++;
		break;

	    case GRANGE_FLAG_ISUMSEQ:
	    case GRANGE_FLAG_ISCONS:
	    case GRANGE_FLAG_ISREF:
		break;

	    default:
		/* Unknown ISMASK flag */
		fprintf(stderr, "bin %"PRIrec" range item %d: "
			"Unknown GRANGE_FLAG_IS? flag: %d\n",
			bin->rec, i, r->flags & GRANGE_FLAG_ISMASK);
		err++;
	    }
	}
    }


    /* Check count validity to ensure this + children are correct */
    if (bin->nseqs != bs->nseq) {
	fprintf(stderr, "bin %"PRIrec": nseqs does not match observed "
		"counts\n", bin->rec);
	if (fix)
	    bin->nseqs = bs->nseq;
	err++;
    }
    if (io->db->version > 1 && bin->nanno != bs->nanno) {
	fprintf(stderr, "bin %"PRIrec": nanno does not match observed "
		"counts\n", bin->rec);
	if (fix)
	    bin->nanno = bs->nanno;
	err++;
    }
    if (io->db->version > 1 && bin->nrefpos != bs->nref) {
	fprintf(stderr, "bin %"PRIrec": nrefpos does not match observed "
		"counts\n", bin->rec);
	if (fix)
	    bin->nrefpos = bs->nref;
	err++;
    }


    /*
     * Check used start/end range, and accumulate absolute positions so we
     * can check the contig later.
     */
    if (bin->rng) {
	if (start != bin->start_used ||
	    end   != bin->end_used) {
	    fprintf(stderr, "bin %"PRIrec": used start/end range are "
		    "incorrect\n", bin->rec);
	    err++;
	}

	if (nthis_seq) {
	    if (bs->cstart > NMIN(cstart, cend))
		bs->cstart = NMIN(cstart, cend);
	    if (bs->cend   < NMAX(cstart, cend))
		bs->cend   = NMAX(cstart, cend);
	}
    } else {
	if (bin->start_used != 0 || bin->end_used != 0) {
	    fprintf(stderr, "bin %"PRIrec": used start/end are non-zero "
		    "in an empty bin\n", bin->rec);
	    err++;
	}
    }

    cache_decr(io, bin);

    return err;
}

int check_contig(GapIO *io, tg_rec crec) {
    contig_t *c;
    bin_stats bs;
    int err;

    c = cache_search(io, GT_Contig, crec);
    cache_incr(io, c);

    bs.cstart = INT_MAX;
    bs.cend   = INT_MIN;
    bs.nseq   = 0;
    bs.nanno  = 0;
    bs.nref   = 0;

    err = bin_walk(io, 0 /* fix */, c->bin, contig_offset(io, &c), 0, &bs);
    if (bs.cstart != c->start ||
	bs.cend   != c->end) {
	fprintf(stderr, "Contig %"PRIrec": used start/end range are "
		"incorrect\n", crec);
	err++;
    }
    
    cache_decr(io, c);

    return err;
}
